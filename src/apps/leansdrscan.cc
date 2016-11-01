#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/wait.h>

void fatal(const char *s) { perror(s); exit(1); }

struct field {
  int nvalues;
  char **values;
  int current;
  struct field *next;
  field(char *s) {
    int nsep = 0;
    for ( unsigned int i=0; i<strlen(s); ++i ) if ( s[i] == ',' ) ++nsep;
    nvalues = nsep+1;
    values = new char*[nvalues];
    values[0] = strtok(s, ",");
    for ( int i=1; i<nvalues; ++i )
      values[i] = strtok(NULL, ",");
    current = 0;
  }
  bool iterate() {
    ++current;
    if ( current == nvalues ) {
      current = 0;
      if ( next ) return next->iterate();
      return false;
    }
    return true;
  }
};

struct config {
  bool verbose;
  size_t maxsend;
  float timeout;
  bool rewind;
  field *fields;
  int nfields;
  config() :
    verbose(false), maxsend(16<<20), timeout(1.0),
    rewind(false), fields(NULL), nfields(0) { }
};

int do_write(int fd, char *buf, int count) {
  while ( count ) {
    int nw = write(fd, buf, count);
    if ( nw < 0 ) return nw;
    if ( ! nw ) fatal("eof");
    buf += nw;
    count -= nw;
  }
  return 0;
}

int run_program(config &cfg, char *const argv[]) {
  int fd0[2], fd1[2];
  if ( pipe(fd0) ) fatal("pipe");
  if ( pipe(fd1) ) fatal("pipe");
  pid_t child = fork();
  if ( ! child ) {
    // Child
    close(fd0[1]);
    close(fd1[0]);
    dup2(fd0[0], 0);
    dup2(fd1[1], 1);
    execvp(argv[0], argv);
    perror("execvp"); 
    exit(errno);
  }

  // Parent
  close(fd0[0]);
  close(fd1[1]);
  int nreceived = 0;
  size_t chunk = 65536;
  struct timeval latest;
  if ( gettimeofday(&latest, NULL) ) fatal("gettimeofday");
  
  size_t maxsend = cfg.maxsend;

  while ( true ) {
    fd_set fds;
    FD_ZERO(&fds);
    if ( !cfg.rewind || maxsend )
      FD_SET(0, &fds);
    FD_SET(fd1[0], &fds);
    struct timeval tv = { (int)cfg.timeout, (int)(cfg.timeout*1e6)%1000000 };
    int ns = select(fd1[0]+1, &fds, NULL, NULL, &tv);
    if ( ns < 0 ) fatal("select");

    // Timeout
    struct timeval now;
    if ( gettimeofday(&now, NULL) ) fatal("gettimeofday");
    float time_silent =
      (now.tv_sec -latest.tv_sec ) + (now.tv_usec-latest.tv_usec)*1e-6;
    if ( time_silent >= cfg.timeout ) {
      if ( cfg.verbose ) fprintf(stderr, "No output from child\n");
      break;
    }

    // Input data from our stdin

    if ( FD_ISSET(0, &fds) ) {
      char buf[chunk];
      if ( ! cfg.rewind ) {
	// Reading from live stream
	size_t nr = read(0, buf, sizeof(buf));
	if ( nr < 0 ) fatal("read");
	if ( ! nr ) {
	  if ( cfg.verbose ) fprintf(stderr, "End of stream, exiting\n");
	  exit(0);
	}
	if ( do_write(fd0[1], buf, nr) )
	  // Broken pipe, child has exited
	  break;
      } else {
	// Reading from file
	size_t maxread = maxsend;
	if ( maxread > sizeof(buf) ) maxread = sizeof(buf);
	ssize_t nr = read(0, buf, maxread);
	if ( nr < 0 ) fatal("read");
	if ( ! nr ) {
	  if ( cfg.verbose ) fprintf(stderr, "Sending EOF\n");
	  close(fd0[1]);
	  maxsend = 0;
	}
	if ( do_write(fd0[1], buf, nr) )
	  // Broken pipe, child has exited
	  break;
	maxsend -= nr;
      }
    }

    // Output data from stdout of child
    
    if ( FD_ISSET(fd1[0], &fds) ) {
      char buf[chunk];
      ssize_t nr = read(fd1[0], buf, sizeof(buf));
      if ( ! nr ) break;
      if ( nr < 0 ) fatal("read(child)");
      if ( ! cfg.rewind )
	// Live streaming
	if ( do_write(1, buf, nr) ) fatal("write");
      nreceived += nr;
      latest = now;
    }
  }

  close(fd0[1]);
  close(fd1[0]);
  kill(child, SIGKILL);
  int status;
  waitpid(child, &status, 0);
  return nreceived;
}


void print_command(char *argv[]) {
  for ( ; *argv; ++argv ) fprintf(stderr, " %s", *argv);
  fprintf(stderr, "\n");
}

int run(config &cfg) {
  // Don't die when child processes terminate
  signal(SIGPIPE, SIG_IGN);

  do {
    do {
      // Try the current combination of settings
      char *argv[cfg.nfields+1];
      int i = 0;
      for ( field *f=cfg.fields; f; ++i,f=f->next )
	argv[i] = f->values[f->current];
      argv[i] = NULL;
      if ( cfg.verbose ) {
	fprintf(stderr, "Trying command:");
	print_command(argv);
      }
      int nreceived = run_program(cfg, argv);
      // Seek to beginning of input file if not in live streaming mode
      if ( cfg.rewind )
	if ( lseek(0, 0, SEEK_SET) ) fatal("lseek");
      if ( nreceived ) {
	if ( cfg.verbose ) {
	  fprintf(stderr, "Got %d with command:", nreceived);
	  print_command(argv);
	}
	if ( cfg.rewind ) {
	  if ( cfg.verbose ) fprintf(stderr, "Now processing whole file.\n");
	  execvp(argv[0], argv);
	  exit(1);
	}
      }
      // Next combination of setting
    } while ( cfg.fields->iterate() );
    // Loop if in live streaming mode
  } while ( ! cfg.rewind );
  return 0;
}

// CLI

void usage(const char *name, FILE *f, int c) {
  fprintf(f, "Usage: %s [options] <program> [program settings]\n", name);
  fprintf(f, "Run <program>, cycling through combinations of settings.\n");
  fprintf(f, "Example: '%s -v cat -n,-e' will feed stdin through"
	  " 'cat -n' and 'cat -e' alternatively.\n", name);

  fprintf(f,
	  "\nOptions:\n"
	  "  -h              Print this message\n"
	  "  -v              Verbose\n"
	  "  --timeout N     Next settings if no output within N seconds\n"
	  "  --rewind        Rewind input (stdin must be a file)\n"
	  "  --probesize N   Forward only N bytes (with --rewind)\n"
	  );
  exit(c);
}

int main(int argc, const char *argv[]) {
  config cfg;

  int i;
  for ( i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "-v") )
      cfg.verbose = true;
    else if ( ! strcmp(argv[i], "--timeout") && i+1<argc )
      cfg.timeout = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--probesize") && i+1<argc )
      cfg.maxsend = atoll(argv[++i]);
    else if ( ! strcmp(argv[i], "--rewind") )
      cfg.rewind = true;
    else if ( argv[i][0] == '-' )
      usage(argv[0], stderr, 1);
    else
      break;
  }

  field **plast = &cfg.fields;
  for ( ; i<argc; ++i ) {
    field *f = new field(strdup(argv[i]));
    f->next = NULL;
    *plast = f;
    plast = &f->next;
    ++cfg.nfields;
  }

  if ( ! cfg.fields ) usage(argv[0], stderr, 1);
  
  if ( cfg.verbose ) {
    fprintf(stderr, "Fields:");
    for ( field *f=cfg.fields; f; f=f->next ) {
      fprintf(stderr, " ");
      if ( f->nvalues > 1 ) fprintf(stderr, "{");
      fprintf(stderr, "%s", f->values[0]);
      for ( int i=1; i<f->nvalues; ++i )
	fprintf(stderr, "|%s", f->values[i]);
      if ( f->nvalues > 1 ) fprintf(stderr, "}");
    }
    fprintf(stderr, "\n");
  }

  return run(cfg);
}
