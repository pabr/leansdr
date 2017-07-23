#include <stdio.h>
#include <stdlib.h>

#include "leansdr/framework.h"
#include "leansdr/generic.h"
#include "leansdr/dsp.h"

using namespace leansdr;

typedef float f32;
typedef complex<f32> cf32;
typedef unsigned char u8;
typedef complex<u8> cu8;

template<typename T>
struct drifter : runnable {
  drifter(scheduler *sch,
	  pipebuf< complex<T> > &_in, pipebuf< complex<T> > &_out)
    : runnable(sch, "drifter"),
      in(_in), out(_out),
      t(0) {
    memset(drifts, 0, sizeof(drifts));
    for ( int i=0; i<65536; ++i ) {
      float a = 2*M_PI * i / 65536;
      lut_trig[i].re = cosf(a);
      lut_trig[i].im = sinf(a);
    }
  }
  
  static const int NCOMPONENTS = 3;

  struct component {
    float amp;       // Amplitude of fluctuation (Hz)
    float freq;      // Rate of fluctuation (Hz)
    signed long a;   // Phase at runtime (2pi/2^32)
  } drifts[NCOMPONENTS];
  
  void run() {
    unsigned long count = min(in.readable(), out.writable());
    complex<T> *pin = in.rd(), *pend = pin+count;
    complex<T> *pout = out.wr();
    signed short phase = 0;
    for ( ; pin<pend; ++pin,++pout,t+=1 ) {
      float f = 0;
      for ( int i=0; i<NCOMPONENTS; ++i ) {
	complex<float> *r = &lut_trig[(unsigned short)(drifts[i].a>>16)];
	f += drifts[i].amp * r->im;
	drifts[i].a += drifts[i].freq * 4294967296.0;
      }
      phase += f * 65536;
      complex<float> *r = &lut_trig[(unsigned short)phase];
      pout->re = pin->re*r->re - pin->im*r->im;
      pout->im = pin->re*r->im + pin->im*r->re;
    }
    in.read(count);
    out.written(count);
  }

private:
  complex<float> lut_trig[65536];
  pipereader< complex<T> > in;
  pipewriter< complex<T> > out;
  unsigned long t;
};

struct config {
  float loop_input;
  enum { IO_F32, IO_U8 } input_format, output_format;
  float scale;          // Input gain factor
  float awgn;           // White gaussian noise standard deviation
  bool deterministic;   // Pseudorandom noise
  float Fs;             // Sampling rate
  float Flo;            // Local oscillator freq
  float ppm;            // Local oscillator accuracy
  float drift_period, drift_rate;
  float drift2_amp, drift2_freq;
  // struct drifter::component drifts[5];
  
  config() :
    loop_input(false),
    input_format(IO_F32),
    output_format(IO_F32),
    scale(1),
    awgn(0),
    deterministic(false),
    Fs(0),
    Flo(0),
    ppm(-1), drift_period(0), drift_rate(0),
    drift2_amp(0), drift2_freq(0)
  {
  }
};

typedef complex<float> cf32;

int run(config &cfg) {
  scheduler sch;
  unsigned long BUF_BASEBAND = 4096;

  pipebuf<cf32> *pipe = NULL;

  switch ( cfg.input_format) {
  case config::config::IO_F32: {
    pipebuf<cf32> *p_stdin =
      new pipebuf<cf32>(&sch, "stdin", BUF_BASEBAND);
    file_reader<cf32> *r_stdin = new file_reader<cf32>(&sch, 0, *p_stdin);
    r_stdin->loop = cfg.loop_input;
    pipe = p_stdin;
    break;
  }
  case config::IO_U8: {
    pipebuf<cu8> *p_stdin = new pipebuf<cu8>(&sch, "stdin", BUF_BASEBAND);
    file_reader<cu8> *r_stdin = new file_reader<cu8>(&sch, 0, *p_stdin);
    r_stdin->loop = cfg.loop_input;
    pipebuf<cf32> *p_stdinf =
      new pipebuf<cf32>(&sch, "stdinf", BUF_BASEBAND);
    new cconverter<u8,128, f32,0, 1,1>(&sch, *p_stdin, *p_stdinf);
    pipe = p_stdinf;
    break;
  }
  }
  
  pipebuf<cf32> p_scaled(&sch, "scaled", BUF_BASEBAND);
  scaler<float,cf32,cf32> r_scale(&sch, cfg.scale, *pipe, p_scaled);
  pipe = &p_scaled;

  if ( ! cfg.deterministic )
    srand48(getpid());
  pipebuf<cf32> p_noise(&sch, "noise", BUF_BASEBAND);
  wgn_c<f32> r_noise(&sch, p_noise);
  r_noise.stddev = cfg.awgn;
  pipebuf<cf32> p_noisy(&sch, "noisy", BUF_BASEBAND);
  adder<cf32> r_addnoise(&sch, *pipe, p_noise, p_noisy);
  pipe = &p_noisy;

  pipebuf<cf32> p_drift(&sch, "drift", BUF_BASEBAND);
  drifter<float> r_drift(&sch, *pipe, p_drift);
  float maxoffs = cfg.Flo * cfg.ppm * 1e-6;
  r_drift.drifts[0].amp = maxoffs / cfg.Fs;
  if ( cfg.drift_period && cfg.drift_rate )
    fail("Specify only one of --drift-rate and --drift-period");
  if ( cfg.drift_period )
    r_drift.drifts[0].freq = (1.0/cfg.drift_period) / cfg.Fs;
  if ( cfg.drift_rate ) {
    if ( ! cfg.ppm ) fail("Need --ppm with --drift-rate");
    r_drift.drifts[0].freq = (cfg.drift_rate/(2*M_PI*cfg.ppm)) / cfg.Fs;
  }
  if ( cfg.drift2_amp && cfg.drift2_freq ) {
    r_drift.drifts[1].amp = cfg.drift2_amp / cfg.Fs;
    r_drift.drifts[1].freq = cfg.drift2_freq / cfg.Fs;
  }
  pipe = &p_drift;
  
  switch ( cfg.output_format ) {
  case config::IO_U8: {
    pipebuf<cu8> *p_stdout = new pipebuf<cu8>(&sch, "stdou", BUF_BASEBAND);
    new cconverter<f32,0, u8,128, 1,1>(&sch, *pipe, *p_stdout);
    new file_writer<cu8>(&sch, *p_stdout, 1);
    break;
  }
  case config::IO_F32: {
    new file_writer<cf32>(&sch, *pipe, 1);
    break;
  }
  }

  sch.run();

  return 0;
}

void usage(const char *name, FILE *f, int c) {
  fprintf(f, "Usage: %s [options]  < IQ.in  > IQ.out\n", name);
  fprintf(f, "Simulate an imperfect communication channel.\n");
  fprintf(f,
	  "\nInput options:\n"
	  "  --iu8              Interpret stdin as complex unsigned char\n"
	  "  --if32             Interpret stdin as complex float\n"
	  "  -f Hz              Specify sample rate\n"
	  "  --loop             Repeat (stdin must be a file)\n"
	  );
  fprintf(f,
	  "\nGain options:\n"
	  "  --scale FACTOR     Multiply by constant\n"
	  );
  fprintf(f,
	  "\nDrift options:\n"
	  "  --lo HZ            Specify nominal LO frequency\n"
	  "  --ppm PPM          Specify LO accuracy\n"
	  "  --drift-period S   Drift +-ppm every S seconds\n"
	  "  --drift-rate R     Drift with maximum rate R (Hz/s)\n"
	  "  --drift2-amp HZ    Add secondary drift (range in Hz)\n"
	  "  --drift2-freq HZ   Add secondary drift (rate in Hz)\n"
	  );
  fprintf(f,
	  "\nNoise options:\n"
	  "  --awgn STDDEV      Add white gaussian noise (dB)\n"
	  );
  fprintf(f,
	  "\nOutput options:\n"
	  "  --ou8              Output as complex unsigned char\n"
	  "  --of32             Output as complex float\n"
	  );
  exit(c);
}

int main(int argc, char *argv[]) {
  config cfg;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "--iu8") )
      cfg.input_format = config::IO_U8;
    else if ( ! strcmp(argv[i], "--if32") )
      cfg.input_format = config::IO_F32;
    else if ( ! strcmp(argv[i], "--loop") )
      cfg.loop_input = true;
    else if ( ! strcmp(argv[i], "--ou8") )
      cfg.output_format = config::IO_U8;
    else if ( ! strcmp(argv[i], "--of32") )
      cfg.output_format = config::IO_F32;
    else if ( ! strcmp(argv[i], "-f") && i+1<argc )
      cfg.Fs = atof(argv[++i]);
    // Scale
    else if ( ! strcmp(argv[i], "--scale") && i+1<argc )
      cfg.scale = atof(argv[++i]);
    // Noise
    else if ( ! strcmp(argv[i], "--awgn") && i+1<argc )
      cfg.awgn = expf(logf(10)*atof(argv[++i])/20);
    else if ( ! strcmp(argv[i], "--deterministic") )
      cfg.deterministic = true;
    // Drift
    else if ( ! strcmp(argv[i], "--lo") && i+1<argc )
      cfg.Flo = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--ppm") && i+1<argc )
      cfg.ppm = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--drift-period") && i+1<argc )
      cfg.drift_period = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--drift-rate") && i+1<argc )
      cfg.drift_rate = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--drift2-amp") && i+1<argc )
      cfg.drift2_amp = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--drift2-freq") && i+1<argc )
      cfg.drift2_freq = atof(argv[++i]);
    else
      usage(argv[0], stderr, 1);
  }

  return run(cfg);
}
