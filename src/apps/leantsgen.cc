#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#ifndef VERSION
#define VERSION "undefined"
#endif

void usage(const char *name, FILE *f, int c, const char *info=NULL) {
      fprintf(stderr, "Usage: %s [options]\n", name);
      fprintf(stderr, "Output numbered MPEG TS packets on stdout.\n");
      fprintf
	(f, "\nOptions:\n"
	 "  --version         Display version and exit\n"
	 "  -c INT            Number of packets to output\n"
	 );
      if ( info ) fprintf(f, "Error while processing '%s'\n", info);
      exit(c);
}

int main(int argc, char *argv[]) {
  const int SIZE = 188;
  int count = -1;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i],"-c") && i+1<argc )
      count = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--version") ) {
      printf("%s\n", VERSION);
      exit(0);
    }
    else
      usage(argv[0], stderr, 1, argv[i]);
  }
  
  for ( unsigned long t=0; count<0 || count--; ++t ) {
    unsigned char msg[SIZE];
    for ( int i=0; i+3<SIZE; i+=4 ) {
      // Format: 8-bit byte number, 24-bit packet number.
      msg[i] = i;
      msg[i+1] = t>>16;
      msg[i+2] = t>>8;
      msg[i+3] = t;
    }
    msg[0] = 0x47;
    if ( write(1,msg,SIZE) != SIZE ) { perror("write"); return 1; }
  }

  return 0;
}
