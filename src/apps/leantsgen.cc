#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char *argv[]) {
  const int SIZE = 188;
  int count = -1;

  for ( int i=1; i<argc; ++i ) {
    if ( !strcmp(argv[i],"-c") && i+1<argc )
      count = atoi(argv[++i]);
    else {
      fprintf(stderr, "Usage: %s [-c PACKETCOUNT]\n", argv[0]);
      fprintf(stderr, "Output numbered MPEG TS packets on stdout.\n");
      exit(1);
    }
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
