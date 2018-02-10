// This file is part of LeanSDR (c) <pabr@pabr.org>.
// See the toplevel README for more information.

// leansdrserv interfaces leansdr command pipelines with network sockets.

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <netinet/in.h>

void fatal(const char *s) { perror(s); exit(1); }

struct config {
  bool verbose;
  int data1_httpd;     // Port for data output HTTP server, or -1
  int info3_httpd;     // Port for info output HTTP server, or -1
  int control4_httpd;  // Port for control input HTTP server, or -1
  char *const *command;
  config()
    : verbose(false), data1_httpd(-1), info3_httpd(-1), control4_httpd(-1)
  { }
};

struct infobuffer {
  struct accumulator {
    const char *tag;
    int nlines;
    int nused;
    char **lines;
    accumulator *next;
    accumulator(const char *_tag, int _nlines) :
      tag(strdup(_tag)), nlines(_nlines), nused(0),
      lines(new char*[nlines]) { }
    void put(const char *line) {
      if ( nused == nlines ) {
	free(lines[0]);
	for ( int i=0; i<nused-1; ++i ) lines[i] = lines[i+1];
	--nused;
      }
      lines[nused++] = strdup(line);
    }
    void dump(FILE *f) {
      fprintf(f, "\"%s\":[", tag);
      for ( int i=0; i<nused; ++i ) fprintf(f, "%s%s", i?",":"", lines[i]);
      fprintf(f, "]");
      //nused = 0;
    }
  } *accumulators;
  infobuffer() : accumulators(NULL) { }
  accumulator *add_accumulator(const char *tag, int nlines) {
    accumulator *pa = new accumulator(tag, nlines);
    pa->next = accumulators;
    accumulators = pa;
    return pa;
  }
  void put(const char *tag, const char *line) {
    accumulator *pa;
    for ( pa=accumulators; pa; pa=pa->next )
      if ( ! strcmp(pa->tag, tag) ) break;
    if ( ! pa ) pa = add_accumulator(tag, 1);
    pa->put(line);
  }
  void dump(FILE *f) {
    fprintf(f, "{\n");
    for ( accumulator *pa=accumulators; pa; pa=pa->next ) {
      pa->dump(f);
      if ( pa->next ) fprintf(f, ",");
      fprintf(f, "\n");
    }
    fprintf(f, "}");
  }
};

int opt_listener(int port) {
  if ( port < 0 ) return -1;
  int fd = socket(AF_INET, SOCK_STREAM, 0);
  if ( fd < 0 ) fatal("socket");
  struct sockaddr_in addr;
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_port = ntohs(port);
  addr.sin_addr.s_addr = INADDR_ANY;
  int reuse = 1;
  setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));
  if ( bind(fd,(sockaddr*)&addr,sizeof(addr)) < 0 ) fatal("bind");
  if ( listen(fd, 2) < 0 ) fatal("listen");
  return fd;
}

int run(config &cfg) {
  infobuffer infobuf;
  bool intercept1=false, intercept3=false, control4=false;
  if ( cfg.data1_httpd >= 0 ) intercept1 = true;
  if ( cfg.info3_httpd >= 0 ) intercept3 = true;
  if ( cfg.control4_httpd >= 0 ) control4 = true;
  
  int fd1[2];  // pipe fd 1 from child process
  if ( intercept1 ) if ( pipe(fd1) ) fatal("pipe");
  int fd3[2];  // pipe fd 3 from child process
  if ( intercept3 ) if ( pipe(fd3) ) fatal("pipe");
  int fd4[2];  // pipe fd 4 to child process
  if ( control4 ) if ( pipe(fd4) ) fatal("pipe");

  pid_t child_pid = fork();
  
  if ( ! child_pid ) {
    // Execute child process
    if ( intercept1 ) {
      close(fd1[0]);
      dup2(fd1[1], 1);
      close(fd1[1]);
    } else fprintf(stderr, "Not intercepting stdout\n");
    if ( intercept3 ) {
      close(fd3[0]);
      dup2(fd3[1], 3);
      close(fd3[1]);
    } else fprintf(stderr, "Not intercepting fd 3\n");
    if ( control4 ) {
      close(fd4[1]);
      dup2(fd4[0], 4);
      close(fd4[0]);
    } else fprintf(stderr, "Not controlling fd 4\n");
    execvp(cfg.command[0], cfg.command);
    fatal("execvp"); 
  }

  if ( intercept1 ) close(fd1[1]);
  if ( intercept3 ) close(fd3[1]);
  if ( control4   ) close(fd4[0]);

  int fd_data1_httpd = opt_listener(cfg.data1_httpd);
  int fd_info3_httpd = opt_listener(cfg.info3_httpd);
  int fd_control4_httpd = opt_listener(cfg.control4_httpd);

  int fd_out1 = 1;  // Default: Forward to stdout

  char buf3[65536];
  int nbuf3 = 0;
  
  int fdmax = 0;
  if ( intercept1 && fd1[0]>fdmax ) fdmax = fd1[0];
  if ( intercept3 && fd3[0]>fdmax ) fdmax = fd3[0];
  if ( fd_data1_httpd > fdmax ) fdmax = fd_data1_httpd;
  if ( fd_info3_httpd > fdmax ) fdmax = fd_info3_httpd;
  if ( fd_control4_httpd > fdmax ) fdmax = fd_control4_httpd;
  
  while ( true ) {
    fd_set fds;
    FD_ZERO(&fds);
    if ( intercept1 ) FD_SET(fd1[0], &fds);
    if ( intercept3 ) FD_SET(fd3[0], &fds);
    if ( fd_data1_httpd >= 0 ) FD_SET(fd_data1_httpd, &fds);
    if ( fd_info3_httpd >= 0 ) FD_SET(fd_info3_httpd, &fds);
    if ( fd_control4_httpd >= 0 ) FD_SET(fd_control4_httpd, &fds);
    int ns = select(fdmax+1, &fds, NULL, NULL, NULL);
    if ( ns < 0 ) fatal("select");

    if ( intercept1 && FD_ISSET(fd1[0], &fds) ) {
      // Process stdout from child
      char buf[65536];
      int nr = read(fd1[0], buf, sizeof(buf));
      if ( nr < 0 ) fatal("read(stdout)");
      if ( ! nr ) exit(0);
      if ( fd_out1 >= 0 ) {
	int nw = write(fd_out1,buf,nr);
	if ( nw != nr ) {
	  if ( fd_out1 == 1 )
	    fatal("error or partial write on stdout");
	  else {
	    perror("error or partial write on client socket");
	    close(fd_out1);
	    fd_out1 = -1;
	  }
	}
      }
    }
    
    if ( intercept3 && FD_ISSET(fd3[0], &fds) ) {
      // Process fd3 from child
      int nr = read(fd3[0], buf3+nbuf3, sizeof(buf3)-nbuf3);
      if ( nr < 0 ) fatal("read(fd3)");
      if ( ! nr ) exit(0);
      nbuf3 += nr;
      char *tag=buf3, *eol;
      while ( tag<buf3+nbuf3 && (eol=strchr(tag,'\n')) ) {
	*eol = 0;
	char *space = strchr(tag,' ');
	if ( ! space ) fatal("Invalid format on fd 3");
	*space = 0;
	char *line = space+1;
	infobuf.put(tag, line);
	tag = eol+1;
      }
      nbuf3 = (buf3+nbuf3) - tag;
      memmove(buf3, tag, nbuf3);
    }
    
    if ( fd_data1_httpd>=0 && FD_ISSET(fd_data1_httpd,&fds) ) {
      int fd = accept(fd_data1_httpd, NULL, NULL);
      if ( fd < 0 ) fatal("accept(1)");
      FILE *f = fdopen(fd, "w");
      if ( ! f ) fatal("fdopen(1)");
      fprintf(f, "HTTP/1.0 200 OK\r\n");
      fprintf(f, "Content-Type: application/json\r\n");
      fprintf(f, "Access-Control-Allow-Origin: *\r\n");
      fprintf(f, "\r\n");
      fflush(f);
      // From now on, forward stdout from child to this socket
      if ( fd_out1 >= 0 ) close(fd_out1);
      fd_out1 = fd;
    }

    if ( fd_info3_httpd>=0 && FD_ISSET(fd_info3_httpd,&fds) ) {
      int fd = accept(fd_info3_httpd, NULL, NULL);
      if ( fd < 0 ) fatal("accept(3)");
      FILE *f = fdopen(fd, "w");
      if ( ! f ) fatal("fdopen(3)");
      fprintf(f, "HTTP/1.0 200 OK\r\n");
      fprintf(f, "Content-Type: application/json\r\n");
      fprintf(f, "Access-Control-Allow-Origin: *\r\n");
      fprintf(f, "\r\n");
      infobuf.dump(f);
      fflush(f);
      shutdown(fd, SHUT_RDWR);
      fclose(f);
    }

    if ( fd_control4_httpd>=0 && FD_ISSET(fd_control4_httpd,&fds) ) {
      int fd = accept(fd_control4_httpd, NULL, NULL);
      if ( fd < 0 ) fatal("accept(4)");
      FILE *f = fdopen(fd, "r+");
      if ( ! f ) fatal("fdopen(4)");
      char req[256];
      if ( ! fgets(req, sizeof(req), f) ) fatal("fgets(4)");
      if ( cfg.verbose ) fprintf(stderr, "Control request: %s\n", req);
      int nw = write(fd4[1], req, strlen(req));
      if ( nw != strlen(req) ) fatal("write(4)");
      char h[4096];
      while ( fgets(h,sizeof(h),f) && h[0] && h[0]!='\r' && h[0]!='\n' ) ;
      fprintf(f, "HTTP/1.0 200 OK\r\n");
      fprintf(f, "Content-Type: text/plain\r\n");
      fprintf(f, "Access-Control-Allow-Origin: *\r\n");
      fprintf(f, "\r\n");
      fprintf(f, "ECHO: %s", req);
      fflush(f);
      shutdown(fd, SHUT_RDWR);
      fclose(f);
    }

  }

}
  
// CLI

void usage(const char *name, FILE *f, int c) {
  fprintf(f, "Usage: %s [options] <command> ...\n", name);
  fprintf(f, "Run <program>, redirecting file descriptors as specified.\n");
  fprintf
    (f,
     "\nOptions:\n"
     "  -h                     Print this message\n"
     "  -v                     Verbose\n"
     "  --raw1-httpd PORT      Forward raw stream from FD 1 to HTTP\n"
     "  --info3-httpd PORT     Forward tagged data from FD 3 to HTTP\n"
     "  --control4-httpd PORT  Forward control from HTTP to FD 4\n"
     );
  exit(c);
}

int main(int argc, char *argv[]) {
  config cfg;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "-v") )
      cfg.verbose = true;
    else if ( ! strcmp(argv[i], "--data1-httpd") && i+1<argc )
      cfg.data1_httpd = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--info3-httpd") && i+1<argc )
      cfg.info3_httpd = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--control4-httpd") && i+1<argc )
      cfg.control4_httpd = atoi(argv[++i]);
    else if ( argv[i][0] == '-' )
      usage(argv[0], stderr, 1);
    else {
      cfg.command = &argv[i];
      break;
    }
  }

  return run(cfg);
}
