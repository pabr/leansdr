// This file is part of LeanSDR Copyright (C) 2016-2018 <pabr@pabr.org>.
// See the toplevel README for more information.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "leansdr/framework.h"
#include "leansdr/generic.h"
#include "leansdr/dsp.h"
#include "leansdr/sdr.h"
#include "leansdr/dvb.h"
#include "leansdr/rs.h"
#include "leansdr/filtergen.h"

using namespace leansdr;

template<typename T>
struct interpolator : runnable {
  unsigned int d;

  interpolator(scheduler *sch, int _d, pipebuf<T> &_in, pipebuf<T> &_out)
    : runnable(sch, "interpolator"),
      d(_d),
      in(_in), out(_out, d) { }
  void run() {
    unsigned long count = min(in.readable(), out.writable()/d);
    T *pin=in.rd(), *pend=pin+count, *pout=out.wr();
    for ( ; pin<pend; ++pin ) {
      *pout++ = *pin;
      for ( int skip=d-1; skip--; ) *pout++ = 0;
    }
    in.read(count);
    out.written(count*d);
  }
private:
  pipereader<T> in;
  pipewriter<T> out;
};

struct config {
  cstln_lut<256>::predef constellation;
  code_rate fec;
  float amp;       // Desired RMS constellation amplitude
  bool agc;
  int interp;
  int decim;
  float rolloff;
  float rrc_rej;
  enum { OUTPUT_F32, OUTPUT_S16 } output_format;
  bool fill;
  bool verbose, debug;
  config()
    : constellation(cstln_lut<256>::QPSK), fec(FEC12),
      amp(1.0), agc(false),
      interp(2), decim(1), rolloff(0.35), rrc_rej(10),
      output_format(OUTPUT_F32),
      fill(false),
      verbose(false), debug(false)
  { }
};
    
void run(config &cfg) {
  scheduler sch;
  sch.verbose = cfg.verbose;
  sch.debug = cfg.debug;

  int buf_factor = 2;
  unsigned long BUF_PACKETS = 12*buf_factor;  // TBD Reduce copying
  unsigned long BUF_BYTES = SIZE_RSPACKET*BUF_PACKETS;
  unsigned long BUF_SYMBOLS = BUF_BYTES*8 * 2;  // Worst case BPSK 1/2
  unsigned long BUF_BASEBAND = 4096;

  // TS PACKETS ON STDIN

  pipebuf<tspacket> p_tspackets(&sch, "TS packets", BUF_PACKETS);
  file_reader<tspacket> r_stdin(&sch, 0, p_tspackets);

  // RANDOMIZER

  pipebuf<tspacket> p_rtspackets(&sch, "rand TS packets", BUF_PACKETS);
  randomizer r_rand(&sch, p_tspackets, p_rtspackets);

  // RS-ENCODER

  pipebuf< rspacket<u8> > p_rspackets(&sch, "RS-enc packets", BUF_PACKETS);
  rs_encoder r_rsenc(&sch, p_rtspackets, p_rspackets);

  // INTERLEAVER

  pipebuf<u8> p_mpegbytes(&sch, "mpegbytes", BUF_BYTES);
  interleaver r_inter(&sch, p_rspackets, p_mpegbytes);

  // CONVOLUTIONAL CODER

  cstln_lut<256> *cstln = make_dvbs2_constellation(cfg.constellation, cfg.fec);
  int bits_per_symbol = log2i(cstln->nsymbols);

  if ( cfg.fec==FEC23 && (cstln->nsymbols==4 ||
			  cstln->nsymbols==64) ) {
    if ( cfg.verbose ) fprintf(stderr, "Handling rate 2/3 as 4/6\n");
    cfg.fec = FEC46;
  }
  pipebuf<u8> p_symbols(&sch, "symbols", BUF_SYMBOLS);
  dvb_convol r_convol(&sch, p_mpegbytes, p_symbols, cfg.fec, bits_per_symbol);

  // IQ MAPPER

  pipebuf<cf32> p_iqsymbols(&sch, "IQ symbols", BUF_SYMBOLS);
  cstln_transmitter<f32,0> r_mod(&sch, p_symbols, p_iqsymbols);
  r_mod.cstln = cstln;

  // RESAMPLER

  pipebuf<cf32> p_interp(&sch, "interpolated", BUF_BASEBAND);
  float Fm = 1.0 / cfg.interp;
  int order = cfg.interp * cfg.rrc_rej;
  float *coeffs;
  int ncoeffs = filtergen::root_raised_cosine
    (order, Fm, cfg.rolloff, &coeffs);
  // This yields about the desired power level even without AGC.
  filtergen::normalize_power(ncoeffs, coeffs, cfg.amp/cstln_amp);

  if ( sch.verbose )
    fprintf(stderr, "Interpolation: ratio %d/%d, rolloff %f, %d coeffs\n",
	    cfg.interp, cfg.decim, cfg.rolloff, ncoeffs);
  if ( sch.debug )
    filtergen::dump_filter("rrc", ncoeffs, coeffs);

  fir_resampler<cf32,float>
    r_resampler(&sch, ncoeffs, coeffs,
		p_iqsymbols, p_interp, cfg.interp, 1);

  // TBD Combine interp and decim
  pipebuf<cf32> p_resampled(&sch, "resampled", BUF_BASEBAND);
  decimator<cf32> r_decim(&sch, cfg.decim, p_interp, p_resampled);

  pipebuf<cf32> *tail = &p_resampled;

  // AGC

  if ( cfg.agc ) {
    pipebuf<cf32> *p_agc =
      new pipebuf<cf32>(&sch, "AGC", BUF_BASEBAND);
    simple_agc<f32> *r_agc =
      new simple_agc<f32>(&sch, *tail, *p_agc);
    r_agc->out_rms = cfg.amp / sqrtf((float)cfg.interp/cfg.decim);
    // Adjust bandwidth for large interpolation ratios.
    r_agc->bw = 0.001 * cfg.decim / cfg.interp;
    tail = p_agc;
  }

  // IQ ON STDOUT

  switch ( cfg.output_format ) {
  case config::OUTPUT_F32:
    (void)new file_writer<cf32>(&sch, *tail, 1);
    break;
  case config::OUTPUT_S16: {
    pipebuf<cs16> *p_stdout =
      new pipebuf<cs16>(&sch, "stdout", BUF_BASEBAND);
    (void)new cconverter<f32,0, int16_t,0, 32768,1>(&sch, *tail, *p_stdout);
    (void)new file_writer<cs16>(&sch, *p_stdout, 1);
    break;
  }
  default:
    fail("Output format not implemented");
  }

  if ( cfg.fill ) {
    if ( cfg.verbose ) fprintf(stderr, "Realtime mode\n");
    tspacket blank;
    memset(blank.data, 0, 188);
    blank.data[0] = 0x47;
    r_stdin.set_realtime(blank);
  }

  sch.run();
  sch.shutdown();
  if ( sch.verbose ) sch.dump();
}

// Command-line

void usage(const char *name, FILE *f, int c, const char *info=NULL) {
  fprintf(f, "Usage: %s [options]  < TS  > IQ\n", name);
  fprintf(f, "Modulate MPEG packets into a DVB-S baseband signal\n");
  fprintf(f, "Output float complex samples\n");
  fprintf
    (f, "\nOptions:\n"
     "  --const STRING           QPSK (default),\n"
     "                           BPSK .. 32APSK (DVB-S2),\n"
     "                           64APSKe (DVB-S2X),\n"
     "                           16QAM .. 256QAM (experimental)\n"
     "  --cr STRING              1/2, 2/3, 3/4, 5/6, 7/8\n"
     "  -f INTERP[/DECIM]        Samples per symbols (default: 2)\n"
     "  --roll-off FLOAT         RRC roll-off (default: 0.35)\n"
     "  --rrc-rej FLOAT          RRC filter rejection (defaut: 10)\n"
     "  --power FLOAT            Output power (dB, default: 0)\n"
     "  --agc                    Better regulation of output power\n"
     "  --f32                    Output 32-bit floats, range +-1.0 (default)\n"
     "  --s16                    Output 16-bit ints\n"
     "  --fill                   Insert blank packets\n"
     "  -v                       Output debugging info at startup and exit\n"
     "  -d                       Output debugging info during operation\n"
     "  --version                Display version and exit\n"
     );
  if ( info ) fprintf(f, "** Error while processing '%s'\n", info);
  exit(c);
}

int main(int argc, char *argv[]) {
  config cfg;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "-v") )
      cfg.verbose = true;
    else if ( ! strcmp(argv[i], "-d") )
      cfg.debug = true;
    else if ( ! strcmp(argv[i], "--version") ) {
      printf("%s\n", VERSION);
      exit(0);
    }
    else if ( ! strcmp(argv[i], "--cr") && i+1<argc ) {
      ++i;
      // DVB-S
      if      ( ! strcmp(argv[i], "1/2" ) ) cfg.fec = FEC12;
      else if ( ! strcmp(argv[i], "2/3" ) ) cfg.fec = FEC23;
      else if ( ! strcmp(argv[i], "3/4" ) ) cfg.fec = FEC34;
      else if ( ! strcmp(argv[i], "5/6" ) ) cfg.fec = FEC56;
      else if ( ! strcmp(argv[i], "7/8" ) ) cfg.fec = FEC78;
      // DVB-S2
      else if ( ! strcmp(argv[i], "4/5"  ) ) cfg.fec = FEC45;
      else if ( ! strcmp(argv[i], "8/9"  ) ) cfg.fec = FEC89;
      else if ( ! strcmp(argv[i], "9/10" ) ) cfg.fec = FEC910;
      else usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "--const") && i+1<argc ) {
      ++i;
      if      ( ! strcmp(argv[i], "BPSK" ) )
	cfg.constellation = cstln_lut<256>::BPSK;
      else if ( ! strcmp(argv[i], "QPSK" ) )
	cfg.constellation = cstln_lut<256>::QPSK;
      else if ( ! strcmp(argv[i], "8PSK" ) )
	cfg.constellation = cstln_lut<256>::PSK8;
      else if ( ! strcmp(argv[i], "16APSK" ) )
	cfg.constellation = cstln_lut<256>::APSK16;
      else if ( ! strcmp(argv[i], "32APSK" ) )
	cfg.constellation = cstln_lut<256>::APSK32;
      else if ( ! strcmp(argv[i], "64APSKe" ) )
	cfg.constellation = cstln_lut<256>::APSK64E;
      else if ( ! strcmp(argv[i], "16QAM" ) )
	cfg.constellation = cstln_lut<256>::QAM16;
      else if ( ! strcmp(argv[i], "64QAM" ) )
	cfg.constellation = cstln_lut<256>::QAM64;
      else if ( ! strcmp(argv[i], "256QAM" ) )
	cfg.constellation = cstln_lut<256>::QAM256;
      else usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "-f") && i+1<argc ) {
      ++i;
      cfg.decim = 1;
      if ( sscanf(argv[i], "%d/%d", &cfg.interp, &cfg.decim) < 1 )
	usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "--roll-off") && i+1<argc )
      cfg.rolloff = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--rrc-rej") && i+1<argc )
      cfg.rrc_rej = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--power") && i+1<argc )
      cfg.amp = expf(logf(10)*atof(argv[++i])/20);
    else if ( ! strcmp(argv[i], "--agc") )
      cfg.agc = true;
    else if ( ! strcmp(argv[i], "--f32") )
      cfg.output_format = config::OUTPUT_F32;
    else if ( ! strcmp(argv[i], "--s16") )
      cfg.output_format = config::OUTPUT_S16;
    else if ( ! strcmp(argv[i], "--fill") )
      cfg.fill = true;
    else 
      usage(argv[0], stderr, 1, argv[i]);
  }

  run(cfg);

  return 0;
}
