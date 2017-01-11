#ifndef LEANSDR_GENERIC_H
#define LEANSDR_GENERIC_H

#include <sys/types.h>
#include <unistd.h>

#include "leansdr/math.h"

namespace leansdr {

//////////////////////////////////////////////////////////////////////
// Simple blocks
//////////////////////////////////////////////////////////////////////

// [file_reader] reads raw data from a file descriptor into a [pipebuf].
// If the file descriptor is seekable, data can be looped.

template<typename T>
struct file_reader : runnable {
  file_reader(scheduler *sch, int _fdin, pipebuf<T> &_out)
    : runnable(sch, _out.name),
      loop(false),
      fdin(_fdin), out(_out),
      pos(0) {
  }
  void run() {
    size_t size = out.writable() * sizeof(T);
    if ( ! size ) return;
  again:
    ssize_t nr = read(fdin, out.wr(), size);
    if ( nr < 0 ) fatal("read");
    if ( !nr && !loop ) return;
    if ( ! nr ) {
      if ( sch->debug ) fprintf(stderr, "%s looping\n", name);
      off_t res = lseek(fdin, 0, SEEK_SET);
      if ( res == (off_t)-1 ) fatal("lseek");
      goto again;
    }
    if ( nr % sizeof(T) ) fatal("partial read");
    out.written(nr / sizeof(T));
  }
  bool loop;
private:
  int fdin;
  pipewriter<T> out;
  off_t pos;
};

// [file_writer] writes raw data from a [pipebuf] to a file descriptor.

template<typename T>
struct file_writer : runnable {
  file_writer(scheduler *sch, pipebuf<T> &_in, int _fdout) :
    runnable(sch, _in.name),
    in(_in), fdout(_fdout) {
  }
  void run() {
    int size = in.readable() * sizeof(T);
    if ( ! size ) return;
    int nw = write(fdout, in.rd(), size);
    if ( ! nw ) fatal("pipe");
    if ( nw < 0 ) fatal("write");
    if ( nw % sizeof(T) ) fatal("partial write");
    in.read(nw/sizeof(T));
  }
private:
  pipereader<T> in;
  int fdout;
};

// [file_printer] writes data from a [pipebuf] to a file descriptor,
// with printf-style formatting and optional scaling.

template<typename T>
struct file_printer : runnable {
  file_printer(scheduler *sch, const char *_format,
	      pipebuf<T> &_in, int _fdout) :
    runnable(sch, _in.name),
    scale(1), in(_in), format(_format), fdout(_fdout) {
  }
  void run() {
    int n = in.readable();
    T *pin=in.rd(), *pend=pin+n;
    for ( ; pin<pend; ++pin ) {
      char buf[256];
      int len = snprintf(buf, sizeof(buf), format, (*pin)*scale);
      if ( len < 0 ) fatal("obsolete glibc");
      int nw = write(fdout, buf, len);
      if ( nw != len ) fatal("partial write");
    }
    in.read(n);
  }
  T scale;
private:
  pipereader<T> in;
  const char *format;
  int fdout;
};

// [file_listprinter] writes all data available from a [pipebuf]
// to a file descriptor on a single line.
// Special case for complex.

template<typename T>
struct file_carrayprinter : runnable {
  file_carrayprinter(scheduler *sch,
		     const char *_head,
		     const char *_format,
		     const char *_tail,
		     pipebuf< complex<T> > &_in, int _fdout) :
    runnable(sch, _in.name),
    scale(1), in(_in),
    head(_head), format(_format), tail(_tail),
    fout(fdopen(_fdout,"w")) {
  }
  void run() {
    int n = in.readable();
    if ( n && fout ) {
      fprintf(fout, head, n);
      complex<T> *pin=in.rd(), *pend=pin+n;
      for ( ; pin<pend; ++pin )
	fprintf(fout, format, pin->re*scale, pin->im*scale);
      fprintf(fout, tail);
      fflush(fout);
    }
    in.read(n);
  }
  T scale;
private:
  pipereader< complex<T> > in;
  const char *head, *format, *tail;
  FILE *fout;
};

// [itemcounter] writes the number of input items to the output [pipebuf].
// [Tout] must be a numeric type.

template<typename Tin, typename Tout>
struct itemcounter : runnable {
  itemcounter(scheduler *sch, pipebuf<Tin> &_in, pipebuf<Tout> &_out)
    : runnable(sch, "itemcounter"),
      in(_in), out(_out) {
  }
  void run() {
    if ( out.writable() < 1 ) return;
    unsigned long count = in.readable();
    if ( ! count ) return;
    *out.wr() = count;
    in.read(count);
    out.written(1);
  }
private:
  pipereader<Tin> in;
  pipewriter<Tout> out;
};

// [decimator] forwards 1 in N sample.

template<typename T>
struct decimator : runnable {
  unsigned int d;

  decimator(scheduler *sch, int _d, pipebuf<T> &_in, pipebuf<T> &_out)
    : runnable(sch, "decimator"),
      d(_d),
      in(_in), out(_out) {
  }
  void run() {
    unsigned long count = min(in.readable()/d, out.writable());
    T *pin=in.rd(), *pend=pin+count*d, *pout=out.wr();
    for ( ; pin<pend; pin+=d, ++pout )
      *pout = *pin;
    in.read(count*d);
    out.written(count);
  }
private:
  pipereader<T> in;
  pipewriter<T> out;
};

  // [rate_estimator] accumulates counts of two quantities
  // and periodically outputs their ratio.

  template<typename T>
  struct rate_estimator : runnable {
    int sample_size;

    rate_estimator(scheduler *sch,
		   pipebuf<int> &_num, pipebuf<int> &_den,
		   pipebuf<float> &_rate)
      : runnable(sch, "rate_estimator"),
	sample_size(10000),
	num(_num), den(_den), rate(_rate),	
	acc_num(0), acc_den(0) {
    }
    
    void run() {
      if ( rate.writable() < 1 ) return;
      int count = min(num.readable(), den.readable());
      int *pnum=num.rd(), *pden=den.rd();
      for ( int n=count; n--; ++pnum,++pden ) {
	acc_num += *pnum;
	acc_den += *pden;
      }
      num.read(count);
      den.read(count);
      if ( acc_den >= sample_size ) {
	*rate.wr() = (float)acc_num / acc_den;
	rate.written(1);
	acc_num = acc_den = 0;
      }
    }
    
  private:
    pipereader<int> num, den;
    pipewriter<float> rate;
    T acc_num, acc_den;
  };

}  // namespace

#endif  // LEANSDR_GENERIC_H
