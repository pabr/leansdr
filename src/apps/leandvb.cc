// leandvb.cc copyright (c) 2016-2017 pabr@pabr.org
// http://www.pabr.org/radio/leandvb

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>

#include "leansdr/framework.h"
#include "leansdr/generic.h"
#include "leansdr/dsp.h"
#include "leansdr/sdr.h"
#include "leansdr/dvb.h"
#include "leansdr/rs.h"
#include "leansdr/gui.h"
#include "leansdr/filtergen.h"

using namespace leansdr;

// Main loop

struct config {
  bool verbose, debug;
  enum { INPUT_U8, INPUT_F32 } input_format;
  float float_scale; // Scaling factor for float data.
  bool loop_input;
  float Fs;          // Sampling frequency (Hz) 
  float Fderot;      // Shift the signal (Hz). Note: Ftune is faster
  int anf;           // Number of auto notch filters
  bool cnr;            // Measure CNR
  unsigned int decim;  // Decimation, 0=auto
  int fd_pp;         // FD for preprocessed data, or -1
  float awgn;        // Standard deviation of noise

  float Fm;          // QPSK symbol rate (Hz) 
  enum dvb_version { DVB_S, DVB_S2 } standard;
  cstln_lut<256>::predef constellation;
  code_rate fec;
  float Ftune;       // Bias frequency for the QPSK demodulator (Hz)
  bool allow_drift;
  bool fastlock;
  bool viterbi;
  bool hard_metric;
  bool resample;
  float resample_rej;  // Approx. filter rejection in dB
  bool rrc;            // Apply root raised-cosine filter
  float rolloff;       // Roll-off 0..1

  bool gui;          // Plot stuff
  float duration;    // Horizontal span of timeline GUI (s)
  bool linger;       // Keep GUI running after EOF
  int fd_info;       // FD for status information in text format, or -1
  int fd_const;      // FD for constellation and symbols, or -1

  config()
    : verbose(false),
      debug(false),
      input_format(INPUT_U8),
      float_scale(1.0),
      loop_input(false),
      Fs(2.4e6),
      Fderot(0),
      anf(1),
      cnr(false),
      decim(0),
      fd_pp(-1),
      awgn(0),
      Fm(2e6),
      standard(DVB_S),
      constellation(cstln_lut<256>::QPSK),
      fec(FEC12),
      Ftune(0),
      allow_drift(false),
      fastlock(false),
      viterbi(false),
      hard_metric(false),
      resample(false),
      resample_rej(10),
      rrc(false),
      rolloff(0.35),

      gui(false),
      duration(60),
      linger(false),
      fd_info(-1),
      fd_const(-1) {
  }
};

int run(config &cfg) {

  int w_timeline = 512, h_timeline = 256;
  int w_fft = 1024, h_fft = 256;
  int wh_const = 256;
  
  scheduler sch;
  sch.verbose = cfg.verbose;
  sch.debug = cfg.debug;

  int x0 = 100, y0 = 40;
  
  window_placement window_hints[] = {
    { "rawiq (iq)", x0, y0, wh_const,wh_const },
    { "rawiq (spectrum)", x0+300, y0, w_fft, h_fft },
    { "preprocessed (iq)", x0, y0+300, wh_const, wh_const },
    { "preprocessed (spectrum)", x0+300, y0+300, w_fft, h_fft },
    { "PSK symbols", x0, y0+600, wh_const, wh_const },
    { "timeline", x0+300, y0+600, w_timeline, h_timeline },
    { NULL, }
  };
  sch.windows = window_hints;

  int BUF_OVERSIZE = 4;
  // Min buffer size for baseband data
  //   scopes: 1024
  //   ss_estimator: 1024
  //   anf: 4096
  //   cstln_receiver: reads in chunks of 128+1
  unsigned long BUF_BASEBAND = 4096 * BUF_OVERSIZE;
  // Min buffer size for IQ symbols
  //   cstln_receiver: writes in chunks of 128/omega symbols (margin 128)
  //   deconv_sync: reads at least 64+32
  // A larger buffer improves performance significantly.
  unsigned long BUF_SYMBOLS = 1024 * BUF_OVERSIZE;
  // Min buffer size for unsynchronized bytes
  //   deconv_sync: writes 32 bytes
  //   mpeg_sync: reads up to 204*scan_syncs = 1632 bytes
  unsigned long BUF_BYTES = 2048 * BUF_OVERSIZE;
  // Min buffer size for synchronized (but interleaved) bytes
  //   mpeg_sync: writes 1 rspacket
  //   deinterleaver: reads 17*11*12+204 = 2448 bytes
  unsigned long BUF_MPEGBYTES = 2448 * BUF_OVERSIZE;
  // Min buffer size for packets: 1
  unsigned long BUF_PACKETS = BUF_OVERSIZE;
  // Min buffer size for misc measurements: 1
  unsigned long BUF_SLOW = BUF_OVERSIZE;

  // INPUT

  pipebuf<cf32> p_rawiq(&sch, "rawiq", BUF_BASEBAND);

  if ( cfg.input_format == config::INPUT_U8 ) {
    pipebuf<cu8> *p_stdin =
      new pipebuf<cu8>(&sch, "stdin", BUF_BASEBAND);
    file_reader<cu8> *r_stdin =
      new file_reader<cu8>(&sch, 0, *p_stdin);
    r_stdin->loop = cfg.loop_input;
    cconverter<u8,128, f32,0, 1,1> *r_convert =
      new cconverter<u8,128, f32,0, 1,1>(&sch, *p_stdin, p_rawiq);
  }
  if ( cfg.input_format == config::INPUT_F32 ) {
    pipebuf<cf32> *p_stdin =
      new pipebuf<cf32>(&sch, "stdin", BUF_BASEBAND);
    file_reader<cf32> *r_stdin =
      new file_reader<cf32>(&sch, 0, *p_stdin);
    r_stdin->loop = cfg.loop_input;
    scaler<float,cf32,cf32> *r_scale =
      new scaler<float,cf32,cf32>(&sch, cfg.float_scale, *p_stdin, p_rawiq);
  }

#ifdef GUI
  float amp = 128;

  if ( cfg.gui ) {
    cscope<f32> *r_cscope_raw =
      new cscope<f32>(&sch, p_rawiq, -amp, amp, "rawiq (iq)");
    spectrumscope<f32> *r_fft_raw =
      new spectrumscope<f32>(&sch, p_rawiq, amp, "rawiq (spectrum)");
    r_fft_raw->amax *= 0.25;
  }
#endif

  pipebuf<cf32> *p_preprocessed = &p_rawiq;

  // NOISE

  if ( cfg.awgn ) {
    if ( cfg.verbose ) 
      fprintf(stderr, "Adding noise with stddev %f\n", cfg.awgn);
    pipebuf<cf32> *p_noise =
      new pipebuf<cf32>(&sch, "noise", BUF_BASEBAND);
    wgn_c<f32> *r_noise =
      new wgn_c<f32>(&sch, *p_noise);
    r_noise->stddev = cfg.awgn;
    pipebuf<cf32> *p_noisy =
      new pipebuf<cf32>(&sch, "noisy", BUF_BASEBAND);
    adder<cf32> *r_addnoise =
      new adder<cf32>(&sch, *p_preprocessed, *p_noise, *p_noisy);
    p_preprocessed = p_noisy;
  }

  // NOTCH FILTER

  if ( cfg.anf ) {
    pipebuf<cf32> *p_autonotched =
      new pipebuf<cf32>(&sch, "autonotched", BUF_BASEBAND);
    auto_notch<f32> *r_auto_notch =
      new auto_notch<f32>(&sch, *p_preprocessed, *p_autonotched,
			  cfg.anf, 0);
    p_preprocessed = p_autonotched;
  } else {
    if ( cfg.verbose )
      fprintf(stderr, "ANF is disabled (requires a clean signal).\n");
  }

  // FREQUENCY CORRECTION

  if ( cfg.Fderot ) {
    if ( cfg.verbose )
      fprintf(stderr, "Derotating from %.3f kHz\n", cfg.Fderot/1e3);
    pipebuf<cf32> *p_derot =
      new pipebuf<cf32>(&sch, "derotated", BUF_BASEBAND);
    rotator<f32> *r_derot =
      new rotator<f32>(&sch, *p_preprocessed, *p_derot, -cfg.Fderot/cfg.Fs);
    p_preprocessed = p_derot;
  }

  // CNR ESTIMATION

  pipebuf<f32> p_cnr(&sch, "cnr", BUF_SLOW);
  cnr_fft<f32> *r_cnr = NULL;
  if ( cfg.cnr ) {
    if ( cfg.verbose )
      fprintf(stderr, "Measuring CNR\n");
    r_cnr = new cnr_fft<f32>(&sch, *p_preprocessed, p_cnr, cfg.Fm/cfg.Fs);
    r_cnr->decimation = 128*1024;  // Same as demod.decimation
  }

  // FILTERING

  if ( cfg.verbose ) fprintf(stderr, "Roll-off %g\n", cfg.rolloff);

  fir_filter<cf32,float> *r_resample = NULL;
  fir_filter<cf32,float> *r_rrc = NULL;
      
  int decim = 1;

  if ( cfg.resample ) {
    // Lowpass-filter and decimate.
    if ( cfg.decim )
      decim = cfg.decim;
    else {
      // Decimate to just above 4 samples per symbol
      float target_Fs = cfg.Fm * 4;
      decim = cfg.Fs / target_Fs;
      if ( decim < 1 ) decim = 1;
    }
    float transition = (cfg.Fm/2) * cfg.rolloff;
    int order = cfg.resample_rej * cfg.Fs / (22*transition);
    order = ((order+1)/2) * 2;  // Make even
    if ( cfg.verbose )
      fprintf(stderr, "Inserting filter: order %d, decimation %d.\n",
	      order, decim);
    pipebuf<cf32> *p_resampled =
      new pipebuf<cf32>(&sch, "resampled", BUF_BASEBAND);
    float *coeffs;
    float Fcut = (cfg.Fm/2) * (1+cfg.rolloff/2) / cfg.Fs;
    int ncoeffs = filtergen::lowpass(order, Fcut, &coeffs);
    if ( cfg.debug ) {
      for ( int i=0; i<order; ++i ) fprintf(stderr, "%f, ", coeffs[i]);
      fprintf(stderr, "\n");
    }
    r_resample = new fir_filter<cf32,float>
      (&sch, ncoeffs, coeffs, *p_preprocessed, *p_resampled, decim);
    p_preprocessed = p_resampled;
    cfg.Fs /= decim;
  }

  if ( cfg.rrc ) {
#if 1
    fprintf(stderr, "RRC not implemented (ignoring).\n");
#else
    float *coeffs;
    int ncoeffs = filtergen::root_raised_cosine
      (cfg.Fm/cfg.Fs, rolloff, &coeffs);
    if ( cfg.verbose )
      fprintf(stderr, "Inserting RRC filter, %d coeffs.\n", ncoeffs);
    pipebuf<cf32> *p_rrc =
      new pipebuf<cf32>(&sch, "rrc", BUF_BASEBAND);
    r_rrc = new fir_filter<cf32,float>
      (&sch, ncoeffs, coeffs, *p_preprocessed, *p_rrc);
    p_preprocessed = p_rrc;
#endif
  }

  // DECIMATION
  // (Unless already done in resampler)

  if ( !cfg.resample && cfg.decim>1 ) {
    decim = cfg.decim;
    if ( cfg.verbose )
      fprintf(stderr, "Inserting decimator 1/%u\n", decim);
    pipebuf<cf32> *p_decimated =
      new pipebuf<cf32>(&sch, "decimated", BUF_BASEBAND);
    decimator<cf32> *p_decim =
      new decimator<cf32>(&sch, decim, *p_preprocessed, *p_decimated);
    p_preprocessed = p_decimated;
    cfg.Fs /= decim;
  }

  if ( cfg.verbose )
    fprintf(stderr, "Fs after resampling/decimation: %f Hz\n", cfg.Fs);

#ifdef GUI
  if ( cfg.gui ) {
    cscope<f32> *r_cscope_pp =
      new cscope<f32>(&sch, *p_preprocessed, -amp, amp, "preprocessed (iq)");
    spectrumscope<f32> *r_fft_pp =
      new spectrumscope<f32>(&sch, *p_preprocessed, amp,
			     "preprocessed (spectrum)");
    r_fft_pp->amax *= 0.25;
    r_fft_pp->decimation /= decim;
  }
#endif

  // OUTPUT PREPROCESSED DATA

  if ( cfg.fd_pp >= 0 ) {
    if ( cfg.verbose )
      fprintf(stderr, "Writing preprocessed data to FD %d\n", cfg.fd_pp);
    file_writer<cf32> *r_ppout =
      new file_writer<cf32>(&sch, *p_preprocessed, cfg.fd_pp);
  }

  // QPSK

  pipebuf<softsymbol> p_symbols(&sch, "PSK soft-symbols", BUF_SYMBOLS);
  pipebuf<f32> p_freq(&sch, "freq", BUF_SLOW);
  pipebuf<f32> p_ss(&sch, "SS", BUF_SLOW);
  pipebuf<f32> p_mer(&sch, "MER", BUF_SLOW);
  pipebuf<cf32> p_sampled(&sch, "PSK symbols", BUF_BASEBAND);
  // TBD retype preprocess as unsigned char
  cstln_receiver<f32> demod(&sch, *p_preprocessed, p_symbols,
			    &p_freq, &p_ss, &p_mer, &p_sampled);
  cstln_lut<256> qpsk(cstln_lut<256>::QPSK);
  if ( cfg.hard_metric ) {
    if ( cfg.verbose )
      fprintf(stderr, "Using hard metric.\n");
    qpsk.harden();
  }
  demod.cstln = &qpsk;
  if ( cfg.standard == config::DVB_S2 ) {
    // For DVB-S2 testing only.
    // Constellation should be determined from PL signalling.
    fprintf(stderr, "DVB-S2: Testing symbol sampler only.\n");
    demod.cstln = make_dvbs2_constellation(cfg.constellation, cfg.fec);
  }
  demod.set_omega(cfg.Fs/cfg.Fm);
  if ( cfg.Ftune ) {
    if ( cfg.verbose )
      fprintf(stderr, "Biasing receiver to %.3f kHz\n", cfg.Ftune/1e3);
    demod.set_freq(cfg.Ftune/cfg.Fs);
  }
  if ( cfg.allow_drift ) {
    if ( cfg.verbose )
      fprintf(stderr, "Allowing unlimited drift.\n");
    demod.set_allow_drift(true);
  } else {
    if ( cfg.verbose )
      fprintf(stderr, "Frequency offset limits: %+.3f..%+.3f kHz.\n",
	      demod.min_freqw*cfg.Fs/65536/1000,
	      demod.max_freqw*cfg.Fs/65536/1000);
  }
  demod.meas_decimation = 128*1024;
  demod.meas_decimation /= decim;

  // TRACKING FILTERS

  if ( r_resample ) {
    r_resample->freq_tap = &demod.freq_tap;
    r_resample->tap_multiplier = 1.0 / decim;
    r_resample->freq_tol = cfg.Fm/(cfg.Fs*decim) * 0.1;
  }

  if ( r_rrc ) {
    r_rrc->freq_tap = &demod.freq_tap;
    r_rrc->freq_tol = cfg.Fm/cfg.Fs * 0.1;
  }

  if ( r_cnr ) {
    r_cnr->freq_tap = &demod.freq_tap;
    r_cnr->tap_multiplier = 1.0 / decim;
  }

#ifdef GUI
  if ( cfg.gui ) {
    cscope<f32> *r_scope_symbols =
      new cscope<f32>(&sch, p_sampled, -amp,amp);
    r_scope_symbols->decimation = 1;
    r_scope_symbols->cstln = &demod.cstln;
  }
#endif

  // NOT VITERBI (deconvolution only)
  // SYNCHRONIZATION

  //  pipebuf<u8> p_bits(&sch, "bits", BUF_DEINTERLEAVE*8);
  // EN 300 421, section 4.4.3, table 2 Punctured code, G1=0171, G2=0133
  //  deconvol r_deconv(&sch, p_symbols, p_bits, 0171, 0133, FEC78);
  //  deconvol_sync r_deconv(&sch, p_symbols, p_bits, FEC12);

  pipebuf<u8> p_bytes(&sch, "bytes", BUF_BYTES);
  pipebuf<int> p_lock(&sch, "lock", BUF_SLOW);

  deconvol_sync_simple *r_deconv = NULL;

  if ( cfg.viterbi ) {
    if ( cfg.fec != FEC12 )
      fail("Viterbi only for code rate 1/2");
    viterbi_sync *r_viterbi = new viterbi_sync(&sch, p_symbols, p_bytes);
    if ( cfg.fastlock ) r_viterbi->sync_decimation = 1;
  } else {
    r_deconv = make_deconvol_sync_simple(&sch, p_symbols, p_bytes, cfg.fec);
    r_deconv->fastlock = cfg.fastlock;
  }

  pipebuf<u8> p_mpegbytes(&sch, "mpegbytes", BUF_MPEGBYTES);
  mpeg_sync<u8,0> r_sync(&sch, p_bytes, p_mpegbytes, r_deconv, &p_lock);
  r_sync.fastlock = cfg.fastlock;

  // DEINTERLEAVING

  pipebuf< rspacket<u8> > p_rspackets(&sch, "RS-enc packets", BUF_PACKETS);
  deinterleaver<u8> r_deinter(&sch, p_mpegbytes, p_rspackets);

  // REED-SOLOMON

  pipebuf<tspacket> p_rtspackets(&sch, "rand TS packets", BUF_PACKETS);
  rs_decoder<u8,0> r_rsdec(&sch, p_rspackets, p_rtspackets);

  // DERANDOMIZATION

  pipebuf<tspacket> p_tspackets(&sch, "TS packets", BUF_PACKETS);
  derandomizer r_derand(&sch, p_rtspackets, p_tspackets);

  // OUTPUT

  file_writer<tspacket> r_stdout(&sch, p_tspackets, 1);

  // AUX OUTPUT

  if ( cfg.fd_info >= 0 ) {
    file_printer<f32> *r_printfreq =
      new file_printer<f32>(&sch, "FREQ %.0f\n", p_freq, cfg.fd_info);
    r_printfreq->scale = cfg.Fs;
    new file_printer<f32>(&sch, "SS %f\n", p_ss, cfg.fd_info);
    new file_printer<f32>(&sch, "MER %.1f\n", p_mer, cfg.fd_info);
    new file_printer<int>(&sch, "LOCK %d\n", p_lock, cfg.fd_info);
    new file_printer<f32>(&sch, "CNR %.1f\n", p_cnr, cfg.fd_info);
    // Output constants immediately
    FILE *f = fdopen(cfg.fd_info, "w");
    static const char *fec_names[] = { "1/2", "2/3", "3/4", "5/6", "7/8" };
    fprintf(f, "CR %s\n", fec_names[cfg.fec]);
    fprintf(f, "SR %f\n", cfg.Fm);
    fflush(f);
  }
  if ( cfg.fd_const >= 0 ) {
    cstln_lut<256> *c = demod.cstln;
    if ( c ) {
      // Output constellation immediately
      FILE *f = fdopen(cfg.fd_const, "w");
      fprintf(f, "CONST %d", c->nsymbols);
      for ( int i=0; i<c->nsymbols; ++i )
	fprintf(f, " %d,%d", c->symbols[i].re, c->symbols[i].im);
      fprintf(f, "\n");
      fflush(f);
    }
    new file_carrayprinter<f32>(&sch, "SYMBOLS %d", " %.0f,%.0f", "\n",
				p_sampled, cfg.fd_const);
  }

  // TIMELINE SCOPE

#ifdef GUI
  pipebuf<float> p_tscount(&sch, "packet counter", BUF_PACKETS*100);
  itemcounter<tspacket,float> r_tscounter(&sch, p_tspackets, p_tscount);
  float max_packet_rate = cfg.Fm / 8 / 204;
  float pixel_rate = cfg.Fs / demod.meas_decimation;
  float max_packets_per_pixel = max_packet_rate / pixel_rate;

  slowmultiscope<f32>::chanspec chans[] = {
    { &p_freq, "estimated frequency", "%3.3f kHz", {0,255,255},
      cfg.Fs*1e-3f,
      (cfg.Ftune-cfg.Fm/2)*1e-3f, (cfg.Ftune+cfg.Fm/2)*1e-3f,
      slowmultiscope<f32>::chanspec::WRAP },
    { &p_ss, "signal strength", "%3.0f", {255,0,0},
      1, 0,128, 
      slowmultiscope<f32>::chanspec::DEFAULT },
    { &p_mer, "MER", "%5.1f dB", {255,0,255},
      1, -10,20, 
      slowmultiscope<f32>::chanspec::DEFAULT },
    { &p_cnr, "CNR", "%5.1f dB", {255,255,0},
      1, -10,20, 
      (r_cnr?
       slowmultiscope<f32>::chanspec::DEFAULT:
       slowmultiscope<f32>::chanspec::DISABLED) },
    { &p_tscount, "TS recovery", "%3.0f %%", {255,255,0},
      110/max_packets_per_pixel, 0, 101,
      (slowmultiscope<f32>::chanspec::flag)
      (slowmultiscope<f32>::chanspec::ASYNC |
       slowmultiscope<f32>::chanspec::SUM) },
  };

  if ( cfg.gui ) {
    slowmultiscope<f32> *r_scope_timeline =
      new slowmultiscope<f32>(&sch, chans, sizeof(chans)/sizeof(chans[0]),
			      "timeline");
    r_scope_timeline->sample_freq = cfg.Fs / demod.meas_decimation;
    unsigned long nsamples = cfg.duration * cfg.Fs / demod.meas_decimation;
    r_scope_timeline->samples_per_pixel = (nsamples+w_timeline)/w_timeline;
  }
#endif  // GUI

  if ( cfg.debug )
    fprintf(stderr,
	    "Output:\n"
	    "  '_': packet received without errors\n"
	    "  '.': error-corrected packet\n"
	    "  '!': packet with remaining errors\n");
    
  sch.run();

  sch.shutdown();

  if ( cfg.verbose ) sch.dump();
  
  if ( cfg.gui && cfg.linger ) while ( 1 ) { sch.run(); usleep(10000); }

  return 0;
}

// Command-line
 
void usage(const char *name, FILE *f, int c) {
  fprintf(f, "Usage: %s [options]  < IQ  > TS\n", name);
  fprintf(f, "Demodulate DVB-S I/Q on stdin, output MPEG packets on stdout\n");
  fprintf(f,
	  "\nInput options:\n"
	  "  --u8           Input format is 8-bit unsigned (rtl_sdr, default)\n"
	  "  --f32          Input format is 32-bit float (gqrx)\n"
	  "  -f HZ          Input sample rate (default: 2.4e6)\n"
	  "  --loop         Repeat (stdin must be a file)\n");
  fprintf(f,
	  "\nPreprocessing options:\n"
	  "  --anf N        Number of birdies to remove (default: 1)\n"
	  "  --derotate HZ  For use with --fd-pp, otherwise use --tune\n"
	  "  --cnr          Measure CNR (CPU-intensive)\n"
	  "  --fd-pp NUM    Dump preprocessed IQ data to file descriptor\n"
	  );
  fprintf(f,
	  "\nDVB-S options:\n"
	  "  --sr HZ        Symbol rate (default: 2e6)\n"
	  "  --tune HZ      Bias frequency for demodulation\n"
	  "  --drift        Track frequency drift beyond safe limits\n"
	  "  --standard S   DVB-S (default), DVB-S2 (not implemented)\n"
	  "  --const C      QPSK (default), 8PSK .. 32APSK (DVB-S2 only)\n"
	  "  --cr N/D       Code rate 1/2 (default) .. 7/8 .. 9/10\n"
	  "  --fastlock     Synchronize more aggressively (CPU-intensive)\n"
	  "  --viterbi      Use Viterbi (CPU-intensive)\n"
	  "  --resample     Resample baseband (CPU-intensive)\n"
	  "  --resample-rej K  Aliasing rejection\n"
	  "  --decim N      Decimate baseband (with aliasing)\n"
	  "  --rrc          Apply root raised cosine filter (CPU-intensive)\n"
	  "  --hq           Enable all CPU-intensive features\n"
	  );
  fprintf(f,
	  "\nUI options:\n"
	  "  -h             Display this help message and exit\n"
	  "  -v             Output debugging info at startup and exit\n"
	  "  -d             Output debugging info during operation\n"
	  "  --fd-info NUM  Output demodulator status to file descriptor\n"
	  "  --fd-const NUM Output constellation and symbols to file descr\n"
	  );
#ifdef GUI
  fprintf(f,
	  "  --gui          Show constellation and spectrum\n"
	  "  --duration S   Width of timeline plot (default: 60)\n"
	  "  --linger       Keep GUI running after EOF\n"
	  );
#endif
  fprintf(f, "\nTesting options:\n"
	  "  --awgn STDDEV  Add white gaussian noise (slow)\n"
	  );
  exit(c);
}

int main(int argc, const char *argv[]) {
  config cfg;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "-v") )
      cfg.verbose = true;
    else if ( ! strcmp(argv[i], "-d") )
      cfg.debug = true;
    else if ( ! strcmp(argv[i], "-f") && i+1<argc )
      cfg.Fs = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--sr") && i+1<argc )
      cfg.Fm = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--standard") && i+1<argc ) {
      ++i;
      if      ( ! strcmp(argv[i], "DVB-S" ) )
	cfg.standard = config::DVB_S;
      else if ( ! strcmp(argv[i], "DVB-S2" ) )
	cfg.standard = config::DVB_S2;
      else usage(argv[0], stderr, 1);
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
      else usage(argv[0], stderr, 1);
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
      else usage(argv[0], stderr, 1);
    }
    else if ( ! strcmp(argv[i], "--fastlock") )
      cfg.fastlock = true;
    else if ( ! strcmp(argv[i], "--viterbi") )
      cfg.viterbi = true;
    else if ( ! strcmp(argv[i], "--hard-metric") )
      cfg.hard_metric = true;
    else if ( ! strcmp(argv[i], "--filter") ) {
      fprintf(stderr, "--filter is obsolete; use --resample.\n");
      cfg.resample = true;
    }
    else if ( ! strcmp(argv[i], "--resample") )
      cfg.resample = true;
    else if ( ! strcmp(argv[i], "--resample-rej") && i+1<argc )
      cfg.resample_rej = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--decim") && i+1<argc )
      cfg.decim = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--rrc") )
      cfg.rrc = true;
    else if ( ! strcmp(argv[i], "--roll-off") && i+1<argc )
      cfg.rolloff = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--hq") ) {
      cfg.fastlock = true;
      cfg.viterbi = true;
      cfg.resample = true;
      cfg.rrc = true;
      cfg.cnr = true;
    }
    else if ( ! strcmp(argv[i], "--anf") && i+1<argc )
      cfg.anf = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--cnr") )
      cfg.cnr = true;
    else if ( ! strcmp(argv[i], "--tune") && i+1<argc )
      cfg.Ftune = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--drift") )
      cfg.allow_drift = true;
    else if  ( ! strcmp(argv[i], "--gui") )
      cfg.gui = true;
    else if ( ! strcmp(argv[i], "--duration") && i+1<argc )
      cfg.duration = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--linger") )
      cfg.linger = true;
    else if ( ! strcmp(argv[i], "--f32") )
      cfg.input_format = config::INPUT_F32;
    else if ( ! strcmp(argv[i], "--u8") )
      cfg.input_format = config::INPUT_U8;
    else if ( ! strcmp(argv[i], "--float-scale") && i+1<argc )
      cfg.float_scale = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--loop") )
      cfg.loop_input = true;
    else if ( ! strcmp(argv[i], "--derotate") && i+1<argc )
      cfg.Fderot = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-pp") && i+1<argc )
      cfg.fd_pp = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--awgn") && i+1<argc )
      cfg.awgn = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-info") && i+1<argc )
      cfg.fd_info = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-const") && i+1<argc )
      cfg.fd_const = atoi(argv[++i]);
    else
      usage(argv[0], stderr, 1);
  }

  return run(cfg);
}
