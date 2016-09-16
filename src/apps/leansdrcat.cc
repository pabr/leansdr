#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>
#include <fcntl.h>
#include <string.h>

void fatal(const char *s) { perror(s); exit(1); }

void usage(const char *name, FILE *f, int c) {
  fprintf(f, "Usage: %s [options]\n", name);
  fprintf(f,
	  "Forward from stdin to stdout at constant rate.\n"
	  "\nOptions:\n"
	  "  --block      Pause when stdout is busy (default: '#' on stderr)\n"
	  "  --nonblock   Silently ignore when stdout is busy\n"
	  "  --cbr R      Set rate in bits per second\n"
	  "  --cbr8 R     Set rate in bytes per second\n"
	  "  --cbr16 R    Set rate in 16-bit words per second\n"
	  "  --cbr32 R    Set rate in 32-bit words per second\n"
	  "  --cbr64 R    Set rate in 64-bit words per second\n"
	  "  -h           Display this help message and exit\n"
	  );
  exit(c);
}

int main(int argc, const char *argv[]) {
  bool block=false, nonblock=false;
  int bytespersec = 2400000 * 2;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "--block") )
      block = true;
    else if ( ! strcmp(argv[i], "--nonblock") )
      nonblock = true;
    else if ( ! strcmp(argv[i], "--cbr") && i+1<argc )
      bytespersec = atoll(argv[++i]) / 8;
    else if ( ! strcmp(argv[i], "--cbr8") && i+1<argc )
      bytespersec = atoll(argv[++i]);
    else if ( ! strcmp(argv[i], "--cbr16") && i+1<argc )
      bytespersec = atoll(argv[++i]) * 2;
    else if ( ! strcmp(argv[i], "--cbr32") && i+1<argc )
      bytespersec = atoll(argv[++i]) * 4;
    else if ( ! strcmp(argv[i], "--cbr64") && i+1<argc )
      bytespersec = atoll(argv[++i]) * 8;
    else
      usage(argv[0], stderr, 1);
  }
  
  size_t blocksize = 4096;

  if ( bytespersec < blocksize )
    blocksize = bytespersec;

  if ( ! block ) {
    long flags = fcntl(1, F_GETFL);
    flags |= O_NONBLOCK;
    if ( fcntl(1, F_SETFL, flags) ) fatal("fcntl(F_SETFL)");
  }

  struct timeval tv0;
  if ( gettimeofday(&tv0, NULL) ) fatal("gettimeofday");
  
  unsigned long long current = 0;

  while ( 1 ) {
    struct timeval tv;
    if ( gettimeofday(&tv, NULL) ) fatal("gettimeofday");
    unsigned long long reltime =
      (tv.tv_sec -tv0.tv_sec )*1000000LL +
      (tv.tv_usec-tv0.tv_usec);
    unsigned long long target = reltime * bytespersec / 1000000;
    unsigned long long want = target - current;
    if ( want < blocksize ) {
      long long us = blocksize * 1000000LL / bytespersec;
      if ( us > 1000000 ) us = 1000000;
      usleep(us);
    } else {
      want = blocksize;
      unsigned char buf[want];
      ssize_t nr = read(0, buf, want);
      if ( nr < 0 ) fatal("read");
      if ( ! nr ) return 0;
      current += nr;
      for ( unsigned char *p=buf; nr; ) {
	ssize_t nw = write(1, p, nr);
	if ( nw < 0 ) {
	  if ( errno == EWOULDBLOCK ) {
	    if ( ! nonblock ) fprintf(stderr, "#");
	    nr = 0;
	  } else
	    fatal("write");
	} else if ( ! nw ) fatal("write: EOF");
	else {
	  // Try again.  If stdout is really busy we will get
	  // EWOULDBLOCK eventually, unless --block was set.
	  p += nw;
	  nr -= nw;
	}
      }
    }
  }    

}
