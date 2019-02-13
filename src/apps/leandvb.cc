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
#include "leansdr/dvbs2.h"
#include "leansdr/rs.h"
#include "leansdr/gui.h"
#include "leansdr/filtergen.h"

#include "leansdr/hdlc.h"
#include "leansdr/iess.h"

using namespace leansdr;

// Main loop

struct config {
  bool verbose, debug, debug2;
  bool highspeed;      // Demodulate raw u8 I/Q without preprocessing
  enum {
    INPUT_U8, INPUT_S8,
    INPUT_S12, INPUT_S16,
    INPUT_F32
  } input_format;
  float float_scale;   // Scaling factor for float data.
  bool loop_input;
  int input_buffer;    // Extra input buffer size
  int buf_factor;      // Buffer sizing
  float Fs;            // Sampling frequency (Hz) 
  float Fderot;        // Shift the signal (Hz). Note: Ftune is faster
  int anf;             // Number of auto notch filters
  bool cnr;            // Measure CNR
  unsigned int decim;  // Decimation, 0=auto
  int fd_pp;           // FD for preprocessed data, or -1
  int fd_iqsymbols;    // FD for sampled symbols, or -1
  float awgn;          // Standard deviation of noise

  float Fm;            // QPSK symbol rate (Hz) 
  enum dvb_version { DVB_S, DVB_S2 } standard;
  cstln_predef constellation;
  bool strongpls;      // For S2 APSK, expect PLS symbols at maximum amplitude
  code_rate fec;
  float Ftune;         // Bias frequency for the QPSK demodulator (Hz)
  bool allow_drift;
  bool fastlock;
  bool viterbi;
  bool hard_metric;
  int ldpc_bf;
  const char *ldpc_helper;
  bool resample;
  float resample_rej;  // Approx. filter rejection in dB
  enum { SAMP_NEAREST, SAMP_LINEAR, SAMP_RRC } sampler;
  int rrc_steps;       // Discrete steps between symbols, 0=auto
  float rrc_rej;       // Approx. RRC filter rejection in dB
  float rolloff;       // Roll-off 0..1

  bool hdlc;           // Expect HDLC frames instead of MPEG packets
  bool packetized;     // Output frames with 16-bit BE length

  bool gui;            // Plot stuff
  float duration;      // Horizontal span of timeline GUI (s)
  bool linger;         // Keep GUI running after EOF
  int fd_info;         // FD for status information in text format, or -1
  float Finfo;         // Desired refresh rate on fd_info (Hz)
  int fd_const;        // FD for constellation and symbols, or -1
  int fd_spectrum;     // FD for spectrum data, or -1
  bool json;           // Use JSON syntax

  config()
    : verbose(false),
      debug(false),
      debug2(false),
      highspeed(false),
      input_format(INPUT_U8),
      float_scale(1.0),
      loop_input(false),
      input_buffer(0),
      buf_factor(4),
      Fs(2.4e6),
      Fderot(0),
      anf(1),
      cnr(false),
      decim(0),
      fd_pp(-1),
      fd_iqsymbols(-1),
      awgn(0),
      Fm(2e6),
      standard(DVB_S),
      constellation(QPSK),
      strongpls(false),
      fec(FEC12),
      Ftune(0),
      allow_drift(false),
      fastlock(false),
      viterbi(false),
      hard_metric(false),
      ldpc_bf(0),
      ldpc_helper(NULL),
      resample(false),
      resample_rej(10),
      sampler(config::SAMP_LINEAR),
      rrc_steps(0),
      rrc_rej(10),
      rolloff(0.35),

      hdlc(false),
      packetized(false),

      gui(false),
      duration(60),
      linger(false),
      fd_info(-1),
      Finfo(5),
      fd_const(-1),
      fd_spectrum(-1),
      json(false) {
  }
};

int decimation(float Fin, float Fout) {
  int d = Fin / Fout;
  return max(d, 1);
}

void output_initial_info(FILE *f, const config &cfg) {
  const char *quote = cfg.json ? "\"" : "";
  static const char *standard_names[] = {
    [config::DVB_S]="DVB-S",
    [config::DVB_S2]="DVB-S2",
  };
  fprintf(f, "STANDARD %s%s%s\n", quote, standard_names[cfg.standard], quote);
  if ( cfg.standard == config::DVB_S ) {
    fprintf(f, "CONSTELLATION %s%s%s\n",
	    quote, cstln_names[cfg.constellation], quote);
    fec_spec *fs = &fec_specs[cfg.fec];
    fprintf(f, "CR %s%d/%d%s\n", quote, fs->bits_in, fs->bits_out, quote);
  }
  fprintf(f, "SR %f\n", cfg.Fm);
}

struct runtime_common {
  scheduler *sch;
  pipebuf<cf32> *p_rawiq;
  fir_filter<cf32,float> *r_resample;
  pipebuf<cf32> *p_preprocessed;
  float Fspp;     // cfg.Fs adjusted after resampling/decimation
  int decimpp;    // cfg.decim or auto
  int rrc_steps;  // cfg.rrc_steps or auto-selected
  pipebuf<f32> *p_freq;
  pipebuf<f32> *p_ss;
  pipebuf<f32> *p_cnr;
  cnr_fft<f32> *r_cnr;
  pipebuf<f32> *p_mer;
  pipebuf<cf32> *p_cstln;
  pipebuf<cf32> *p_cstln_pls;
  sampler_interface<f32> *sampler;
  pipebuf<cf32> *p_iqsymbols;  // DIVERSITY
  pipebuf<int> *p_framelock;
  pipebuf<float> *p_vber;
  pipebuf<int> *p_lock;
  pipebuf<unsigned long> *p_locktime;
#ifdef GUI
  cscope<f32> *r_scope_cstln;
  cscope<f32> *r_scope_cstln_pls;
  pipebuf<float> *p_tscount;
#endif
  pipebuf<int> *p_vbitcount;
  pipebuf<int> *p_verrcount;

  pipebuf<tspacket> *p_tspackets;
  
  runtime_common(const config &cfg) {
    sch = new scheduler();
    sch->verbose = cfg.verbose;
    sch->debug = cfg.debug;

    int w_timeline = 512, h_timeline = 256;
    int w_fft = 1024, h_fft = 256;
    int wh_const = 256;

    int x0 = 100, y0 = 40;

    static window_placement window_hints[] = {
      { "rawiq (iq)", x0, y0, wh_const,wh_const },
      { "rawiq (spectrum)", x0+300, y0, w_fft, h_fft },
      { "preprocessed (iq)", x0, y0+300, wh_const, wh_const },
      { "preprocessed (spectrum)", x0+300, y0+300, w_fft, h_fft },
      { "PSK symbols", x0, y0+600, wh_const, wh_const },  // TBD obsolete
      { "cstln", x0, y0+600, wh_const, wh_const },
      { "PLS cstln", x0, y0+900, wh_const, wh_const },
      { "timeline", x0+300, y0+600, w_timeline, h_timeline },
      { NULL, }
    };
    sch->windows = window_hints;

    unsigned long S2_MAX_SYMBOLS = (90*(1+360)+36*((360-1)/16));
    // Min buffer size for baseband data
    //   scopes: 1024
    //   ss_estimator: 1024
    //   anf: 4096
    unsigned long BUF_BASEBAND;
    switch ( cfg.standard ) {
    case config::DVB_S:
      // cstln_receiver reads in chunks of 128+1
      BUF_BASEBAND = 4096 * cfg.buf_factor;
      break;
    case config::DVB_S2:
      // s2_frame_receiver need enough samples for one frame
      BUF_BASEBAND = S2_MAX_SYMBOLS * 2 * (cfg.Fs/cfg.Fm) * cfg.buf_factor;
      break;
    default:
      fail("not implemented");
    }

    // Min buffer size for TS packets: Up to 39 per BBFRAME
    unsigned long BUF_S2PACKETS = 39 * cfg.buf_factor;
    // Min buffer size for misc measurements: 1
    unsigned long BUF_SLOW = cfg.buf_factor;

    // INPUT

    p_rawiq = new pipebuf<cf32>(sch, "rawiq", BUF_BASEBAND);

    float amp = 1.0;

    switch ( cfg.input_format ) {
    case config::INPUT_U8: {
      pipebuf<cu8> *p_stdin =
	new pipebuf<cu8>(sch, "stdin", BUF_BASEBAND+cfg.input_buffer);
      file_reader<cu8> *r_stdin =
	new file_reader<cu8>(sch, 0, *p_stdin);
      r_stdin->loop = cfg.loop_input;
      cconverter<u8,128, f32,0, 1,1> *r_convert =
	new cconverter<u8,128, f32,0, 1,1>(sch, *p_stdin, *p_rawiq);
      amp = 128;
      break;
    }
    case config::INPUT_S8: {
      pipebuf<cs8> *p_stdin =
	new pipebuf<cs8>(sch, "stdin", BUF_BASEBAND+cfg.input_buffer);
      file_reader<cs8> *r_stdin =
	new file_reader<cs8>(sch, 0, *p_stdin);
      r_stdin->loop = cfg.loop_input;
      cconverter<s8,0, f32,0, 1,1> *r_convert =
	new cconverter<s8,0, f32,0, 1,1>(sch, *p_stdin, *p_rawiq);
      amp = 128;
      break;
    }
    case config::INPUT_S12:
    case config::INPUT_S16: {
      pipebuf<cs16> *p_stdin =
	new pipebuf<cs16>(sch, "stdin", BUF_BASEBAND+cfg.input_buffer);
      file_reader<cs16> *r_stdin =
	new file_reader<cs16>(sch, 0, *p_stdin);
      r_stdin->loop = cfg.loop_input;
      cconverter<s16,0, f32,0, 1,1> *r_convert =
	new cconverter<s16,0, f32,0, 1,1>(sch, *p_stdin, *p_rawiq);
      amp = (cfg.input_format==config::INPUT_S12 ? 2048 : 32768);
      break;
    }
    case config::INPUT_F32: {
      pipebuf<cf32> *p_stdin =
	new pipebuf<cf32>(sch, "stdin", BUF_BASEBAND+cfg.input_buffer);
      file_reader<cf32> *r_stdin =
	new file_reader<cf32>(sch, 0, *p_stdin);
      r_stdin->loop = cfg.loop_input;
      scaler<float,cf32,cf32> *r_scale =
	new scaler<float,cf32,cf32>(sch, cfg.float_scale, *p_stdin, *p_rawiq);
      amp = 2.0;
      break;
    }
    default:
      fail("Input format not implemented");
    }

    p_preprocessed = p_rawiq;

#ifdef GUI
    if ( cfg.gui ) {
      cscope<f32> *r_cscope_raw =
	new cscope<f32>(sch, *p_rawiq, -amp, amp, "rawiq (iq)");
      spectrumscope<f32> *r_fft_raw =
	new spectrumscope<f32>(sch, *p_rawiq, amp, "rawiq (spectrum)");
      r_fft_raw->amax *= 0.25;
    }
#else
    (void)amp;  // GCC warning
#endif

    // NOISE

    if ( cfg.awgn ) {
      if ( cfg.verbose )
	fprintf(stderr, "Adding noise with stddev %f\n", cfg.awgn);
      pipebuf<cf32> *p_noise =
	new pipebuf<cf32>(sch, "noise", BUF_BASEBAND);
      wgn_c<f32> *r_noise =
	new wgn_c<f32>(sch, *p_noise);
      r_noise->stddev = cfg.awgn;
      pipebuf<cf32> *p_noisy =
	new pipebuf<cf32>(sch, "noisy", BUF_BASEBAND);
      adder<cf32> *r_addnoise =
	new adder<cf32>(sch, *p_preprocessed, *p_noise, *p_noisy);
      p_preprocessed = p_noisy;
    }

    // NOTCH FILTER

    if ( cfg.anf ) {
      pipebuf<cf32> *p_autonotched =
	new pipebuf<cf32>(sch, "autonotched", BUF_BASEBAND);
      auto_notch<f32> *r_auto_notch =
	new auto_notch<f32>(sch, *p_preprocessed, *p_autonotched,
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
	new pipebuf<cf32>(sch, "derotated", BUF_BASEBAND);
      rotator<f32> *r_derot =
	new rotator<f32>(sch, *p_preprocessed, *p_derot, -cfg.Fderot/cfg.Fs);
      p_preprocessed = p_derot;
    }

    // CNR ESTIMATION

    p_cnr = new pipebuf<f32>(sch, "cnr", BUF_SLOW);
    if ( cfg.cnr ) {
      if ( cfg.verbose )
	fprintf(stderr, "Measuring CNR\n");
      r_cnr = new cnr_fft<f32>(sch, *p_preprocessed, *p_cnr, cfg.Fm/cfg.Fs);
      r_cnr->decimation = decimation(cfg.Fs, 1);  // 1 Hz
    } else {
      r_cnr = NULL;
    }

    // SPECTRUM

    pipebuf<f32[1024]> *p_spectrum = NULL;

    if ( cfg.fd_spectrum ) {
      if ( cfg.verbose )
	fprintf(stderr, "Measuring SPD\n");
      p_spectrum = new pipebuf<float[1024]>(sch, "spectrum", BUF_SLOW);
      spectrum<f32,1024> *r_spectrum =
	new spectrum<f32,1024>(sch, *p_preprocessed, *p_spectrum);
      r_spectrum->decimation = decimation(cfg.Fs, 1);  // 1 Hz
      r_spectrum->kavg = 0.5;
    }

    // FILTERING

    if ( cfg.verbose ) fprintf(stderr, "Roll-off %g\n", cfg.rolloff);

    Fspp = cfg.Fs;

    decimpp = 1;

    r_resample = NULL;
    if ( cfg.resample ) {
      // Lowpass-filter and decimate.
      if ( cfg.decim )
	decimpp = cfg.decim;
      else {
	// Decimate to just above 4 samples per symbol
	float target_Fs = cfg.Fm * 4;
	decimpp = cfg.Fs / target_Fs;
	if ( decimpp < 1 ) decimpp = 1;
      }
      float transition = (cfg.Fm/2) * cfg.rolloff;
      int order = cfg.resample_rej * cfg.Fs / (22*transition);
      order = ((order+1)/2) * 2;  // Make even
      if ( cfg.verbose )
	fprintf(stderr, "Inserting filter: order %d, decimation %d.\n",
		order, decimpp);
      pipebuf<cf32> *p_resampled =
	new pipebuf<cf32>(sch, "resampled", BUF_BASEBAND);
      float *coeffs;
#if 1  // Cut in middle of roll-off region
      float Fcut = (cfg.Fm/2) * (1+cfg.rolloff/2) / cfg.Fs;
#else  // Cut at beginning of roll-off region
      float Fcut = (cfg.Fm/2) / cfg.Fs;
#endif
      int ncoeffs = filtergen::lowpass(order, Fcut, &coeffs);
      filtergen::normalize_dcgain(ncoeffs, coeffs, 1);
      if ( cfg.debug ) filtergen::dump_filter("lowpass", ncoeffs, coeffs);
      r_resample = new fir_filter<cf32,float>
	(sch, ncoeffs, coeffs, *p_preprocessed, *p_resampled, decimpp);
      (void)r_resample;  // gcc warning
      p_preprocessed = p_resampled;
      Fspp /= decimpp;
    }

    // DECIMATION
    // (Unless already done in resampler)

    if ( !cfg.resample && cfg.decim>1 ) {
      decimpp = cfg.decim;
      if ( cfg.verbose )
	fprintf(stderr, "Inserting decimator 1/%u\n", decimpp);
      pipebuf<cf32> *p_decimated =
	new pipebuf<cf32>(sch, "decimated", BUF_BASEBAND);
      decimator<cf32> *p_decim =
	new decimator<cf32>(sch, decimpp, *p_preprocessed, *p_decimated);
      p_preprocessed = p_decimated;
      Fspp /= decimpp;
    }

    if ( cfg.verbose )
      fprintf(stderr, "Fs after resampling/decimation: %f Hz\n", Fspp);

#ifdef GUI
    if ( cfg.gui ) {
      cscope<f32> *r_cscope_pp =
	new cscope<f32>(sch, *p_preprocessed, -amp, amp, "preprocessed (iq)");
      spectrumscope<f32> *r_fft_pp =
	new spectrumscope<f32>(sch, *p_preprocessed, amp,
			       "preprocessed (spectrum)");
      r_fft_pp->amax *= 0.25;
      r_fft_pp->decimation /= decimpp;
    }
#endif

    // OUTPUT PREPROCESSED DATA

    if ( cfg.fd_pp >= 0 ) {
      if ( cfg.verbose )
	fprintf(stderr, "Writing preprocessed data to FD %d\n", cfg.fd_pp);
      file_writer<cf32> *r_ppout =
      new file_writer<cf32>(sch, *p_preprocessed, cfg.fd_pp);
    }

    // RECEIVER

    p_freq = new pipebuf<f32>(sch, "freq", BUF_SLOW);
    p_ss = new pipebuf<f32>(sch, "SS", BUF_SLOW);
    p_mer = new pipebuf<f32>(sch, "MER", BUF_SLOW);
    p_cstln = new pipebuf<cf32>(sch, "cstln", BUF_BASEBAND);
    p_cstln_pls = new pipebuf<cf32>(sch, "PLS cstln", BUF_BASEBAND);
    p_framelock = new pipebuf<int>(sch, "frame lock", BUF_SLOW);
    p_vber = new pipebuf<float>(sch, "VBER", BUF_SLOW);
    p_lock = new pipebuf<int>(sch, "lock", BUF_SLOW);
    p_locktime = new pipebuf<unsigned long>(sch, "locktime", BUF_S2PACKETS);

    rrc_steps = cfg.rrc_steps;

    switch ( cfg.sampler ) {
    case config::SAMP_NEAREST:
      sampler = new nearest_sampler<float>();
      break;
    case config::SAMP_LINEAR:
      sampler = new linear_sampler<float>();
      break;
    case config::SAMP_RRC: {
      float *coeffs;
      if ( rrc_steps == 0 ) {
	// At least 64 discrete sampling points between symbols
	rrc_steps = max(1, (int)(64*cfg.Fm / Fspp));
      }
      if ( cfg.verbose )
	fprintf(stderr, "RRC interpolator: %d steps\n", rrc_steps);
      float Frrc = Fspp * rrc_steps;  // Sample freq of the RRC filter
      float transition = (cfg.Fm/2) * cfg.rolloff;
      int order = cfg.rrc_rej * Frrc / (22*transition);
      int ncoeffs = filtergen::root_raised_cosine
	(order, cfg.Fm/Frrc, cfg.rolloff, &coeffs);
      // Total gain: rrc_steps * stride through the coeffs = 1
      filtergen::normalize_dcgain(ncoeffs, coeffs, rrc_steps);
      if ( cfg.verbose )
	fprintf(stderr, "RRC filter: %d coeffs.\n", ncoeffs);
      if ( cfg.debug2 ) filtergen::dump_filter("rrc", ncoeffs, coeffs);
      sampler = new fir_sampler<float,float>
	(ncoeffs, coeffs, rrc_steps);
      break;
    }
    default:
      fatal("Interpolator not implemented");
    }

    if ( cfg.fd_iqsymbols >= 0 ) {
      p_iqsymbols = new pipebuf<cf32> (sch, "cstln", BUF_BASEBAND);
      (void)new file_writer<cf32>(sch, *p_iqsymbols, cfg.fd_iqsymbols);
    } else {
      p_iqsymbols = NULL;
    }

#ifdef GUI
    float s_amp = 128;
    if ( cfg.gui ) {
      r_scope_cstln = new cscope<f32>(sch, *p_cstln, -s_amp,s_amp);
      r_scope_cstln->decimation = 1;
      if ( cfg.standard == config::DVB_S2 ) {
	// Dedicated scope for PLS symbols
	r_scope_cstln_pls = new cscope<f32>(sch, *p_cstln_pls, -s_amp,s_amp);
	r_scope_cstln_pls->decimation = 1;
      }
    } else {
      r_scope_cstln = NULL;
      r_scope_cstln_pls = NULL;
    }
#endif

    p_tspackets = new pipebuf<tspacket>(sch, "TS packets", BUF_S2PACKETS);

    p_vbitcount= new pipebuf<int>(sch, "Bits processed", BUF_S2PACKETS);
    p_verrcount = new pipebuf<int>(sch, "Bits corrected", BUF_S2PACKETS);

    // Standard-specific code would go here,
    // outputting into p_tspackets and into the measurements channels.

    new file_writer<tspacket>(sch, *p_tspackets, 1);

    // BER ESTIMATION

    p_vber = new pipebuf<float>(sch, "VBER", BUF_SLOW);
    rate_estimator<float> *r_vber =
      new rate_estimator<float>(sch, *p_verrcount, *p_vbitcount, *p_vber);
    // About twice per second, and slow enough for 2e-5 resolution
    //TBD r_vber->sample_size = cfg.Fm/2;  // About twice per second, depending on CR
    r_vber->sample_size = 1;
    if ( r_vber->sample_size < 50000 ) r_vber->sample_size = 50000;

    // AUX OUTPUT

    if ( cfg.fd_info >= 0 ) {
      file_printer<f32> *r_printfreq =
	new file_printer<f32>(sch, "FREQ %.0f\n", *p_freq, cfg.fd_info);
      r_printfreq->scale = cfg.Fs;
      new file_printer<f32>(sch, "SS %f\n", *p_ss, cfg.fd_info);
      new file_printer<f32>(sch, "MER %.1f\n", *p_mer, cfg.fd_info);
      new file_printer<int>(sch, "FRAMELOCK %d\n", *p_framelock, cfg.fd_info);
      new file_printer<int>(sch, "LOCK %d\n", *p_lock, cfg.fd_info);
      new file_printer<unsigned long>
	(sch, "LOCKTIME %lu\n", *p_locktime, cfg.fd_info,
	 decimation(cfg.Fm/8/204, cfg.Finfo));  // TBD CR
      new file_printer<f32>(sch, "CNR %.1f\n", *p_cnr, cfg.fd_info);
      new file_printer<float>(sch, "VBER %.6f\n", *p_vber, cfg.fd_info);
      // Output constants immediately
      FILE *f = fdopen(cfg.fd_info, "w");
      if ( ! f ) fatal("fdopen(fd_info)");
      output_initial_info(f, cfg);
      fflush(f);  // Do not close the FILE.
    }

#ifdef GUI
    p_tscount = new pipebuf<float>(sch, "packet counter", BUF_S2PACKETS);
    new itemcounter<tspacket,float>(sch, *p_tspackets, *p_tscount);

    float max_packet_rate = cfg.Fm / 8 / 204;
    float pixel_rate = cfg.Finfo;
    float max_packets_per_pixel = max_packet_rate / pixel_rate;

    slowmultiscope<f32>::chanspec chans[] = {
      { p_freq, "estimated frequency", "Offset %3.3f kHz", {0,255,255},
	cfg.Fm*1e-3f,  // TBD S2 specific
	(cfg.Ftune-cfg.Fm/2)*1e-3f, (cfg.Ftune+cfg.Fm/2)*1e-3f,
	slowmultiscope<f32>::chanspec::WRAP },
      { p_ss, "signal strength", "SS %3.3f", {255,0,0},
	1, 0,128,
	slowmultiscope<f32>::chanspec::DEFAULT },
      { p_mer, "MER", "MER %5.1f dB", {255,0,255},
	1, -10,20,
	slowmultiscope<f32>::chanspec::DEFAULT },
      { p_cnr, "CNR", "CNR %5.1f dB", {255,255,0},
	1, -10,20,
	(r_cnr?
	 slowmultiscope<f32>::chanspec::ASYNC:
	 slowmultiscope<f32>::chanspec::DISABLED) },
      { p_tscount, "TS recovery", "%3.0f %%", {255,255,0},
	110/max_packets_per_pixel, 0, 101,
	(slowmultiscope<f32>::chanspec::flag)
	(slowmultiscope<f32>::chanspec::ASYNC |
	 slowmultiscope<f32>::chanspec::SUM) },
    };

    if ( cfg.gui ) {
      slowmultiscope<f32> *r_scope_timeline =
	new slowmultiscope<f32>(sch, chans, sizeof(chans)/sizeof(chans[0]),
				"timeline");
      r_scope_timeline->sample_freq = cfg.Fs / cfg.Fm * cfg.Finfo;
      unsigned long nsamples = cfg.duration * cfg.Fs / cfg.Fm * cfg.Finfo;
      r_scope_timeline->samples_per_pixel = (nsamples+w_timeline)/w_timeline;
  }
#endif  // GUI

  }  // constructor

};  // runtime_common

#ifndef LEANSDR_EXTENSIONS
int run_dvbs2(config &cfg) {
  fail("DVB-S2 support not enabled at compile-time.");
  return 1;
}
#else
int run_dvbs2(config &cfg) {
  runtime_common run(cfg);

  // Min buffer size for slots: 4 for deinterleaver
  unsigned long BUF_SLOTS = modcod_info::MAX_SLOTS_PER_FRAME * cfg.buf_factor;
  pipebuf< plslot<llr_ss> > p_slots(run.sch, "PL slots", BUF_SLOTS);
  s2_frame_receiver<f32,llr_ss> demod(run.sch, run.sampler,
				      *run.p_preprocessed, p_slots,
				      run.p_freq, run.p_ss, run.p_mer,
				      run.p_cstln, run.p_cstln_pls,
				      run.p_iqsymbols, run.p_framelock);
  demod.omega = cfg.Fs/cfg.Fm;
  if ( cfg.Ftune ) {
    if ( cfg.verbose )
      fprintf(stderr, "Biasing frame receiver to %.3f kHz\n", cfg.Ftune/1e3);
    demod.Ftune = cfg.Ftune / cfg.Fm;  // Per symbol
  }
  demod.Fm = cfg.Fm;
  demod.meas_decimation = decimation(cfg.Fm, cfg.Finfo);
  demod.strongpls = cfg.strongpls;

#ifdef GUI
  if ( run.r_scope_cstln )
    run.r_scope_cstln->cstln = (cstln_base**)&demod.cstln;  // TBD variance
  if ( run.r_scope_cstln_pls )
    run.r_scope_cstln_pls->cstln = (cstln_base**)&demod.qpsk;  // TBD variance
#endif

  unsigned long BUF_FRAMES = cfg.buf_factor;
  pipebuf<bbframe> p_bbframes(run.sch, "BB frames", BUF_FRAMES);

  if ( ! cfg.ldpc_helper ) {
    // Bit-flipping mode.
    // Deinterleave into hard bits.
    pipebuf< fecframe<hard_sb> > *p_fecframes
      = new pipebuf< fecframe<hard_sb> >(run.sch, "FEC frames", BUF_FRAMES);
    new s2_deinterleaver<llr_ss,hard_sb>(run.sch, p_slots, *p_fecframes);
    // Decode FEC-protected frames into plain BB frames.
    s2_fecdec<bool,hard_sb> *r_fecdec =
      new s2_fecdec<bool,hard_sb>(run.sch, *p_fecframes, p_bbframes,
				  run.p_vbitcount, run.p_verrcount);
    r_fecdec->bitflips = cfg.ldpc_bf;
    if ( ! cfg.ldpc_bf )
      fprintf(stderr, "Warning: No LDPC error correction selected.\n");
  } else {
    // External LDPC decoder mode.
    // Deinterleave into soft bits.
    // TBD Latency
    pipebuf< fecframe<llr_sb> > *p_fecframes =
      new pipebuf< fecframe<llr_sb> >(run.sch, "FEC frames", BUF_FRAMES);
    new s2_deinterleaver<llr_ss,llr_sb>(run.sch, p_slots, *p_fecframes);
    // Decode FEC-protected frames into plain BB frames.
    new s2_fecdec_helper<llr_t,llr_sb>(run.sch, *p_fecframes, p_bbframes,
				       cfg.ldpc_helper,
				       run.p_vbitcount, run.p_verrcount);
  }

  // Deframe BB frames to TS packets
  s2_deframer deframer(run.sch, p_bbframes, *run.p_tspackets,
		       run.p_lock, run.p_locktime);

  run.sch->run();
  fprintf(stderr, "sch stopped\n");
  run.sch->shutdown();
  if ( cfg.debug ) run.sch->dump();

  return 0;
}  // run_dvbs2
#endif  // LEANSDR_EXTENSIONS

int run_dvbs(config &cfg) {
  runtime_common run(cfg);

  typedef eucl_ss softsymb;
  unsigned long BUF_SYMBOLS = 4096 * cfg.buf_factor;
  pipebuf<softsymb> p_symbols(run.sch, "PSK soft-symbols", BUF_SYMBOLS);
  cstln_receiver<f32,softsymb> demod(run.sch, run.sampler,
				     *run.p_preprocessed, p_symbols,
				     run.p_freq, run.p_ss, run.p_mer,
				     run.p_cstln);
  if ( cfg.constellation != QPSK &&
       cfg.constellation != BPSK )
    fprintf(stderr, "Warning: non-standard constellation for DVB-S\n");
  demod.cstln = make_dvbs2_constellation<softsymb>(cfg.constellation, cfg.fec);
#if 0  // Dump LUT as greymap
  demod.cstln->dump(stdout, 0);
  exit(0);
#endif
  if ( cfg.hard_metric ) {
    if ( cfg.verbose )
      fprintf(stderr, "Using hard metric.\n");
    demod.cstln->harden();
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
  if ( cfg.viterbi ) {
    if ( cfg.verbose ) fprintf(stderr, "PLL parameters for low SNR\n");
    demod.pll_adjustment /= 6;
  }
  demod.meas_decimation = decimation(cfg.Fs, cfg.Finfo);

#ifdef GUI
  if ( run.r_scope_cstln )
    run.r_scope_cstln->cstln = (cstln_base**)&demod.cstln;  // TBD variance
#endif

  // TRACKING FILTERS

  if ( run.r_resample ) {
    run.r_resample->freq_tap = &demod.freq_tap;
    run.r_resample->tap_multiplier = 1.0 / run.decimpp;
    run.r_resample->freq_tol = cfg.Fm/(run.Fspp*run.decimpp) * 0.1;
  }

  if ( run.r_cnr ) {
    run.r_cnr->freq_tap = &demod.freq_tap;
    run.r_cnr->tap_multiplier = 1.0 / run.decimpp;
  }

  // DECONVOLUTION AND SYNCHRONIZATION

  // Min buffer size for unsynchronized bytes
  //   deconv_sync: writes 32 bytes
  //   mpeg_sync: reads up to 204*scan_syncs = 1632 bytes
  unsigned long BUF_BYTES = 2048 * cfg.buf_factor;
  pipebuf<u8> p_bytes(run.sch, "bytes", BUF_BYTES);

  deconvol_sync_simple *r_deconv = NULL;

  if ( cfg.viterbi ) {
    if ( cfg.fec==FEC23 && (demod.cstln->nsymbols==4 ||
			    demod.cstln->nsymbols==64) ) {
      if ( cfg.verbose ) fprintf(stderr, "Handling rate 2/3 as 4/6\n");
      cfg.fec = FEC46;
    }
    viterbi_sync *r = new viterbi_sync(run.sch, p_symbols, p_bytes,
				       demod.cstln, cfg.fec);
    if ( cfg.fastlock ) r->resync_period = 1;
  } else {
    r_deconv = make_deconvol_sync_simple(run.sch, p_symbols, p_bytes, cfg.fec);
    r_deconv->fastlock = cfg.fastlock;
  }

  // Min buffer size for synchronized (but interleaved) bytes
  //   mpeg_sync: writes 1 rspacket
  //   deinterleaver: reads 17*11*12+204 = 2448 bytes
  unsigned long BUF_MPEGBYTES = 2448 * cfg.buf_factor;

  if ( cfg.hdlc ) {
    pipebuf<u8> *p_descrambled =
      new pipebuf<u8>(run.sch, "descrambled", BUF_MPEGBYTES);
    new etr192_descrambler(run.sch, p_bytes, *p_descrambled);
    pipebuf<u8> *p_frames =
      new pipebuf<u8>(run.sch, "frames", BUF_MPEGBYTES);
    hdlc_sync *r_sync = new hdlc_sync(run.sch, *p_descrambled, *p_frames, 2, 278);
    if ( cfg.fastlock ) r_sync->resync_period = 1;
    if ( cfg.packetized ) r_sync->header16 = true;
    new file_writer<u8>(run.sch, *p_frames, 1);
  }

  pipebuf<u8> p_mpegbytes(run.sch, "mpegbytes", BUF_MPEGBYTES);
  if ( ! cfg.hdlc ) {
    mpeg_sync<u8,0> *r_sync =
      new mpeg_sync<u8,0>(run.sch, p_bytes, p_mpegbytes, r_deconv,
			  run.p_lock, run.p_locktime);
    r_sync->fastlock = cfg.fastlock;
  }

  // DEINTERLEAVING

  // Min buffer size for TS packets: 1
  unsigned long BUF_PACKETS = cfg.buf_factor;
  pipebuf< rspacket<u8> > p_rspackets(run.sch, "RS-enc packets", BUF_PACKETS);
  deinterleaver<u8> r_deinter(run.sch, p_mpegbytes, p_rspackets);

  // REED-SOLOMON

  pipebuf<tspacket> p_rtspackets(run.sch, "rand TS packets", BUF_PACKETS);
  rs_decoder<u8,0> r_rsdec(run.sch, p_rspackets, p_rtspackets,
			   run.p_vbitcount, run.p_verrcount);

  // DERANDOMIZATION

  derandomizer r_derand(run.sch, p_rtspackets, *run.p_tspackets);

  if ( cfg.fd_const >= 0 ) {
    cstln_lut<eucl_ss,256> *c = demod.cstln;
    if ( c ) {
      // Output constellation immediately
      FILE *f = fdopen(cfg.fd_const, "w");
      if ( ! f ) fatal("fdopen(fd_const)");
      if ( cfg.json ) {
	fprintf(f, "CONST [");
	for ( int i=0; i<c->nsymbols; ++i )
	  fprintf(f, "%s[%d,%d]", i?",":"",
		  c->symbols[i].re, c->symbols[i].im);
	fprintf(f, "]\n");
      } else {
	fprintf(f, "CONST %d", c->nsymbols);
	for ( int i=0; i<c->nsymbols; ++i )
	  fprintf(f, " %d,%d", c->symbols[i].re, c->symbols[i].im);
	fprintf(f, "\n");
      }
      fflush(f);
    }
    file_carrayprinter<f32> *symbol_printer;
    if ( cfg.json )
      symbol_printer = new file_carrayprinter<f32>
	(run.sch, "SYMBOLS [", "[%.0f,%.0f]", ",", "]\n",
	 *run.p_cstln, cfg.fd_const);
    else
      symbol_printer =  new file_carrayprinter<f32>
	(run.sch, "SYMBOLS %d", " %.0f,%.0f", "", "\n",
	 *run.p_cstln, cfg.fd_const);
    symbol_printer->fixed_size = 128;
  }

  if ( cfg.fd_spectrum >= 0 ) {
  // TIMELINE SCOPE

#if 0
    float max_packet_rate = cfg.Fm / 8 / 204;
    float pixel_rate = cfg.Fs / demod.meas_decimation;
    float max_packets_per_pixel = max_packet_rate / pixel_rate;

    slowmultiscope<f32>::chanspec chans[] = {
      { &p_freq, "estimated frequency", "Offset %3.3f kHz", {0,255,255},
	cfg.Fs*1e-3f,
	(cfg.Ftune-cfg.Fm/2)*1e-3f, (cfg.Ftune+cfg.Fm/2)*1e-3f,
	slowmultiscope<f32>::chanspec::WRAP },
      { &p_ss, "signal strength", "SS %3.3f", {255,0,0},
	1, 0,128,
	slowmultiscope<f32>::chanspec::DEFAULT },
      { &p_mer, "MER", "MER %5.1f dB", {255,0,255},
	1, -10,20,
	slowmultiscope<f32>::chanspec::DEFAULT },
      { &p_cnr, "CNR", "CNR %5.1f dB", {255,255,0},
	1, -10,20,
	(r_cnr?
	 slowmultiscope<f32>::chanspec::ASYNC:
	 slowmultiscope<f32>::chanspec::DISABLED) },
      { &p_tscount, "TS recovery", "%3.0f %%", {255,255,0},
	110/max_packets_per_pixel, 0, 101,
	(slowmultiscope<f32>::chanspec::flag)
	(slowmultiscope<f32>::chanspec::ASYNC |
	 slowmultiscope<f32>::chanspec::SUM) },
    };

    if ( cfg.gui ) {
      slowmultiscope<f32> *r_scope_timeline =
	new slowmultiscope<f32>(run.sch, chans, sizeof(chans)/sizeof(chans[0]),
				"timeline");
      r_scope_timeline->sample_freq = cfg.Fs / demod.meas_decimation;
      unsigned long nsamples = cfg.duration * cfg.Fs / demod.meas_decimation;
      r_scope_timeline->samples_per_pixel = (nsamples+w_timeline)/w_timeline;
    }
#endif  // GUI
  }
  if ( cfg.debug ) {
    if ( ! cfg.hdlc )
      fprintf(stderr,
	      "Output:\n"
	      "  '_': packet received without errors\n"
	      "  '.': error-corrected packet\n"
	      "  '!': packet with remaining errors\n");
    else
      fprintf(stderr,
	      "Output:\n"
	      "  '_': HDLC frame with correct checksum\n"
	      "  '!': HDLC frame with invalid checksum\n"
	      "  '^': HDLC framing error\n");
  }

  run.sch->run();

  run.sch->shutdown();

  if ( cfg.debug ) run.sch->dump();

  if ( cfg.gui && cfg.linger ) while ( 1 ) { run.sch->run(); usleep(10000); }

  return 0;
}  // run_dvbs


#if 1  // TBD
int run_highspeed_s2(config &cfg) {
  fail("--hs is broken.");
  return 1;
}
#else
int run_highspeed_s2(config &cfg) {

  int w_timeline = 512, h_timeline = 256;
  int w_fft = 1024, h_fft = 256;
  int wh_const = 256;

  scheduler sch;
  sch.verbose = cfg.verbose;
  sch.debug = cfg.debug;

  int x0 = 100, y0 = 40;

  window_placement window_hints[] = {
    { "rawiq (iq)", x0, y0, wh_const,wh_const },
    { "PSK symbols", x0, y0+600, wh_const, wh_const },
    { "timeline", x0+300, y0+600, w_timeline, h_timeline },
    { NULL, }
  };
  sch.windows = window_hints;

  unsigned long MAX_SYMBOLS = (90*(1+360)+36*((360-1)/16)) * cfg.buf_factor;
  unsigned long BUF_BASEBAND = MAX_SYMBOLS * 2 * cfg.buf_factor;
  // Min buffer size for IQ symbols
  //   cstln_receiver: writes in chunks of 128/omega symbols (margin 128)
  //   deconv_sync: reads at least 64+32
  // A larger buffer improves performance significantly.
  unsigned long BUF_SYMBOLS = MAX_SYMBOLS * cfg.buf_factor;
  // Min buffer size for unsynchronized bytes
  //   deconv_sync: writes 32 bytes
  //   mpeg_sync: reads up to 204*scan_syncs = 1632 bytes
  unsigned long BUF_BYTES = 2048 * cfg.buf_factor;
  // Min buffer size for synchronized (but interleaved) bytes
  //   mpeg_sync: writes 1 rspacket
  //   deinterleaver: reads 17*11*12+204 = 2448 bytes
  unsigned long BUF_MPEGBYTES = 2448 * cfg.buf_factor;
  // Min buffer size for packets: 1
  unsigned long BUF_PACKETS = cfg.buf_factor;
 // Min buffer size for TS packets: Up to 39 per BBFRAME
  unsigned long BUF_S2PACKETS = 39 * cfg.buf_factor;
  // Min buffer size for misc measurements: 1
  unsigned long BUF_SLOW = cfg.buf_factor;

  // HIGHSPEED S2: INPUT

  if ( cfg.input_format != config::INPUT_S16 )
    fail("--hs requires --s16");

  pipebuf<cs16> p_rawiq(&sch, "rawiq", BUF_BASEBAND+cfg.input_buffer);
  file_reader<cs16> r_stdin(&sch, 0, p_rawiq);
  r_stdin.loop = cfg.loop_input;

#ifdef GUI
  float amp = 2048;
  if ( cfg.gui ) {
    cscope<s16> *r_cscope_raw =
      new cscope<s16>(&sch, p_rawiq, -amp, amp, "rawiq (iq)");
  }
#endif

  // HIGHSPEED S2: DEMOD

  pipebuf<cf32> p_sampled(&sch, "PSK symbols", BUF_SYMBOLS);
  pipebuf<plslot> p_slots(&sch, "slots", BUF_SYMBOLS/90);
  s2_frame_synchronizer frame_sync(&sch, p_rawiq, p_slots, &p_sampled);
  frame_sync.Fm = cfg.Fm;
  frame_sync.Ftune = cfg.Ftune / cfg.Fm;
  frame_sync.omega = cfg.Fs / cfg.Fm;

#ifdef GUI
  if ( cfg.gui ) {
    cscope<float> *r_scope_symbols =
      new cscope<float>(&sch, p_sampled, -amp, amp);
    r_scope_symbols->decimation = 1;
  }
#endif

  // HIGHSPEED S2: DECODING

  pipebuf<fecframe> p_fecframes(&sch, "FEC frames", BUF_PACKETS);
  s2_deinterleaver deint(&sch, p_slots, p_fecframes);
  pipebuf<bbframe> p_bbframes(&sch, "BB frames", BUF_PACKETS);
  s2_fecdec fecdec(&sch, p_fecframes, p_bbframes);
  pipebuf<tspacket> p_tspackets(&sch, "TS packets", BUF_S2PACKETS);
  s2_deframer deframer(&sch, p_bbframes, p_tspackets);
  file_writer<tspacket> r_stdout(&sch, p_tspackets, 1);
  fprintf(stderr, "Running S2 --hs\n");

  sch.run();

  sch.shutdown();

  if ( cfg.debug ) sch.dump();
  
  if ( cfg.gui && cfg.linger ) while ( 1 ) { sch.run(); usleep(10000); }

  return 0;
}
#endif

#if 1  // TBD
int run_highspeed(config &cfg) {
  fail("--hs is broken.");
  return 1;
}
#else
int run_highspeed(config &cfg) {

  int w_timeline = 512, h_timeline = 256;
  int w_fft = 1024, h_fft = 256;
  int wh_const = 256;
  
  scheduler sch;
  sch.verbose = cfg.verbose;
  sch.debug = cfg.debug;

  int x0 = 100, y0 = 40;
  
  window_placement window_hints[] = {
    { "rawiq (iq)", x0, y0, wh_const,wh_const },
    { "PSK symbols", x0, y0+600, wh_const, wh_const },
    { "timeline", x0+300, y0+600, w_timeline, h_timeline },
    { NULL, }
  };
  sch.windows = window_hints;

  // Min buffer size for baseband data
  //   scopes: 1024
  //   ss_estimator: 1024
  //   anf: 4096
  //   cstln_receiver: reads in chunks of 128+1
  unsigned long BUF_BASEBAND = 4096 * cfg.buf_factor;
  // Min buffer size for IQ symbols
  //   cstln_receiver: writes in chunks of 128/omega symbols (margin 128)
  //   deconv_sync: reads at least 64+32
  // A larger buffer improves performance significantly.
  unsigned long BUF_SYMBOLS = 1024 * cfg.buf_factor;
  // Min buffer size for unsynchronized bytes
  //   deconv_sync: writes 32 bytes
  //   mpeg_sync: reads up to 204*scan_syncs = 1632 bytes
  unsigned long BUF_BYTES = 2048 * cfg.buf_factor;
  // Min buffer size for synchronized (but interleaved) bytes
  //   mpeg_sync: writes 1 rspacket
  //   deinterleaver: reads 17*11*12+204 = 2448 bytes
  unsigned long BUF_MPEGBYTES = 2448 * cfg.buf_factor;
  // Min buffer size for packets: 1
  unsigned long BUF_PACKETS = cfg.buf_factor;
  // Min buffer size for misc measurements: 1
  unsigned long BUF_SLOW = cfg.buf_factor;

  // HIGHSPEED: INPUT
  
  if ( cfg.input_format != config::INPUT_U8 )
    fail("--hs requires --u8");

  pipebuf<cu8> p_rawiq(&sch, "rawiq", BUF_BASEBAND+cfg.input_buffer);
  file_reader<cu8> r_stdin(&sch, 0, p_rawiq);
  r_stdin.loop = cfg.loop_input;

#ifdef GUI
  float amp = 128;

  if ( cfg.gui ) {
    cscope<u8> *r_cscope_raw =
      new cscope<u8>(&sch, p_rawiq, 0, 2*amp, "rawiq (iq)");
  }
#endif

  // HIGHSPEED: QPSK

#define ALGEBRAIC_COMPAT 0  // Use legacy constellation receiver ?

  pipebuf<f32> p_freq(&sch, "freq", BUF_SLOW);
  pipebuf<cu8> p_sampled(&sch, "PSK symbols", BUF_BASEBAND);
#if ALGEBRAIC_COMPAT
  fprintf(stderr, "--hs: Using legacy receiver (slower)\n");
  pipebuf<cf32> p_rawiqf(&sch, "rawiq float", BUF_BASEBAND);
  cconverter<u8,128, f32,0, 1,1> r_convert_in(&sch, p_rawiq, p_rawiqf);

  pipebuf<softsymbol> p_symbols(&sch, "PSK soft-symbols", BUF_SYMBOLS);
  pipebuf<cf32> p_sampledf(&sch, "PSK fsymbols", BUF_BASEBAND);
  // TBD retype preprocess as unsigned char
  cstln_receiver<f32> demod(&sch, p_rawiqf, p_symbols,
			    &p_freq, NULL, NULL, &p_sampledf);
  cstln_lut<256> qpsk(QPSK);
  demod.cstln = &qpsk;
  // Convert the sampled symbols to cu8 for GUI
  cconverter<f32,0, u8,128, 1,1>
    r_convert_samp(&sch, p_sampledf, p_sampled);
#else
  // Use the new fixed-point receiver
  pipebuf<u8> p_symbols(&sch, "PSK hard symbols", BUF_SYMBOLS);
  fast_qpsk_receiver<u8> demod(&sch, p_rawiq, p_symbols,
			       &p_freq, &p_sampled);
#endif
  demod.set_omega(cfg.Fs/cfg.Fm);
  if ( cfg.Ftune ) {
    if ( cfg.verbose )
      fprintf(stderr, "Biasing receiver to %.3f kHz\n", cfg.Ftune/1e3);
    demod.set_freq(cfg.Ftune/cfg.Fs);
  }
  if ( cfg.allow_drift ) {
    if ( cfg.verbose )
      fprintf(stderr, "Allowing unlimited drift.\n");
    demod.allow_drift = true;
  } else {
    if ( cfg.verbose )
      fprintf(stderr, "Frequency offset limits: %+.3f..%+.3f kHz.\n",
	      demod.min_freqw*cfg.Fs/65536/1000,
	      demod.max_freqw*cfg.Fs/65536/1000);
  }
  demod.meas_decimation = decimation(cfg.Fs, cfg.Finfo);

#ifdef GUI
  if ( cfg.gui ) {
    cscope<u8> *r_scope_symbols =
      new cscope<u8>(&sch, p_sampled, 0,2*amp);
    r_scope_symbols->decimation = 1;
  }
#endif

  // HIGHSPEED: DECONVOLUTION

  if ( cfg.fec != FEC12 )
    fail("--hs currently supports code rate 1/2 only");

  pipebuf<u8> p_bytes(&sch, "bytes", BUF_BYTES);
#if ALGEBRAIC_COMPAT
  dvb_deconvol_sync_soft r_deconv(&sch, p_symbols, p_bytes);
#else
  dvb_deconvol_sync_hard r_deconv(&sch, p_symbols, p_bytes);
#endif
  r_deconv.resync_period = cfg.fastlock ? 1 : 32;

  // HIGHSPEED: SYNCHRONIZATION

  pipebuf<u8> p_mpegbytes(&sch, "mpegbytes", BUF_MPEGBYTES);
  pipebuf<int> p_lock(&sch, "lock", BUF_SLOW);
  pipebuf<u32> p_locktime(&sch, "locktime", BUF_PACKETS);
  mpeg_sync<u8,0> r_sync(&sch, p_bytes, p_mpegbytes, NULL,
			 &p_lock, &p_locktime);
  r_sync.fastlock = true;
  r_sync.resync_period = cfg.fastlock ? 1 : 32;
    
  // HIGHSPEED: DEINTERLEAVING

  pipebuf< rspacket<u8> > p_rspackets(&sch, "RS-enc packets", BUF_PACKETS);
  deinterleaver<u8> r_deinter(&sch, p_mpegbytes, p_rspackets);

  // HIGHSPEED: REED-SOLOMON

  pipebuf<int> p_vbitcount(&sch, "Bits processed", BUF_PACKETS);
  pipebuf<int> p_verrcount(&sch, "Bits corrected", BUF_PACKETS);
  pipebuf<tspacket> p_rtspackets(&sch, "rand TS packets", BUF_PACKETS);
  rs_decoder<u8,0> r_rsdec(&sch, p_rspackets, p_rtspackets,
			   &p_vbitcount, &p_verrcount);

  // HIGHSPEED: BER ESTIMATION

  pipebuf<float> p_vber(&sch, "VBER", BUF_SLOW);
  rate_estimator<float> r_vber(&sch, p_verrcount, p_vbitcount, p_vber);
  r_vber.sample_size = cfg.Fm/2;  // About twice per second, depending on CR
  // Require resolution better than 2E-5
  if ( r_vber.sample_size < 50000 ) r_vber.sample_size = 50000;

  // HIGHSPEED: DERANDOMIZATION

  pipebuf<tspacket> p_tspackets(&sch, "TS packets", BUF_PACKETS);
  derandomizer r_derand(&sch, p_rtspackets, p_tspackets);

  // HIGHSPEED: OUTPUT

  file_writer<tspacket> r_stdout(&sch, p_tspackets, 1);

  // HIGHSPEED: AUX OUTPUT

  if ( cfg.fd_info >= 0 ) {
    file_printer<f32> *r_printfreq =
      new file_printer<f32>(&sch, "FREQ %.0f\n", p_freq, cfg.fd_info);
    r_printfreq->scale = cfg.Fs;
    new file_printer<int>(&sch, "LOCK %d\n", p_lock, cfg.fd_info);
    new file_printer<u32>(&sch, "LOCKTIME %lu\n", p_locktime, cfg.fd_info,
			  decimation(cfg.Fm/8/204, cfg.Finfo));  // TBD CR
    new file_printer<float>(&sch, "VBER %.6f\n", p_vber, cfg.fd_info);
    // Output constants immediately
    FILE *f = fdopen(cfg.fd_info, "w");
    if ( ! f ) fatal("fdopen(fd_info)");
    output_initial_info(f, cfg);
    fflush(f);
  }
  if ( cfg.fd_const >= 0 ) {
    file_carrayprinter<u8> *symbol_printer;
    if ( cfg.json ) 
      symbol_printer = new file_carrayprinter<u8>
	(&sch, "SYMBOLS [", "[%.0f,%.0f]", ",", "]\n", p_sampled, cfg.fd_const);
    else
      symbol_printer =  new file_carrayprinter<u8>
	(&sch, "SYMBOLS %d", " %d,%d", "", "\n", p_sampled, cfg.fd_const);
    symbol_printer->fixed_size = 128;
  }

  // HIGHSPEED: TIMELINE SCOPE

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

  if ( cfg.debug ) sch.dump();
  
  if ( cfg.gui && cfg.linger ) while ( 1 ) { sch.run(); usleep(10000); }

  return 0;
}
#endif // TBD

// Command-line
 
void usage(const char *name, FILE *f, int c, const char *info=NULL) {
  fprintf(f, "Usage: %s [options]  < IQ  > TS\n", name);
  fprintf(f, "Demodulate DVB-S I/Q on stdin, output MPEG packets on stdout\n");
  fprintf
    (f,
     "\nInput options:\n"
     "  --u8           Input format is 8-bit unsigned (rtl_sdr, default)\n"
     "  --s12          Input format is 12/16-bit signed (PlutoSDR, LimeSDR)\n"
     "  --s16          Input format is 16-bit signed\n"
     "  --f32          Input format is 32-bit float (gqrx)\n"
     "  -f HZ          Input sample rate (Hz, default: 2.4e6)\n"
     "  --loop         Repeat (stdin must be a file)\n"
     "  --inbuf INT    Additional input buffering (samples)\n"
     );
  fprintf
    (f,
     "\nPreprocessing options:\n"
     "  --anf INT             Number of birdies to remove (default: 1)\n"
     "  --derotate HZ         For use with --fd-pp, otherwise use --tune\n"
     "  --resample            Resample baseband (CPU-intensive)\n"
     "  --resample-rej FLOAT  Aliasing rejection (default: 10)\n"
     "  --decim INT           Decimate baseband (causes aliasing)\n"
     "  --cnr                 Measure CNR (requires samprate>3*symbrate)\n"
     );
  fprintf
    (f,
     "\nDVB-S/S2 options:\n"
     "  --sr FLOAT        Symbol rate (Hz, default: 2e6)\n"
     "  --tune FLOAT      Bias frequency for demodulation (Hz)\n"
     "  --drift           Track frequency drift beyond safe limits\n"
     "  --standard S      DVB-S(default), DVB-S2\n"
     "  --const STRING    DVB-S constellation: QPSK(default), BPSK\n"
     "  --cr NUM/DEN      DVB-S code rate: 1/2(default) .. 7/8\n"
     "  --strongpls       DVB-S2: Expect PLS symbols at max amplitude\n"
     "  --fastlock        Synchronize more aggressively (CPU-intensive)\n"
     "  --sampler         Symbol estimation: nearest, linear, rrc\n"
     "  --rrc-steps INT   RRC interpolation factor\n"
     "  --rrc-rej FLOAT   RRC filter rejection (defaut: 10)\n"
     "  --roll-off FLOAT  RRC roll-off (default: 0.35)\n"
     "  --viterbi         DVB-S: Use Viterbi (CPU-intensive)\n"
     "  --hard-metric     Use Hamming distances with Viterbi\n"
     "  --ldpc-bf INT     Max number of LDPC bitflips (default: 0)\n"
     "  --ldpc-helper CMD Spawn external LDPC decoder:\n"
     "                    'CMD --standard DVB-S2 --modcod N [--shortframes]'\n"
     );
  fprintf
    (f,
     "\nCompatibility options:\n"
     "  --hdlc         Expect HDLC frames instead of MPEG packets\n"
     "  --packetized   Output 16-bit length prefix (default: as stream)\n"
     );
  fprintf
    (f,
     "\nGeneral options:\n"
     "  --buf-factor INT  Buffer size factor (default:4)\n"
     "  --hq              Maximize sensitivity\n"
     "                    (Enables all CPU-intensive features)\n"
     "  --hs              Maximize throughput (QPSK CR1/2 only)\n"
     "                    (Disables all preprocessing)\n"
     );
  fprintf
    (f,
     "\nUI options:\n"
     "  -h                   Display this help message and exit\n"
     "  -v                   Output debugging info at startup and exit\n"
     "  -d                   Output debugging info during operation\n"
     "  --version            Display version and exit\n"
     "  --fd-info FDNUM      Output demodulator status to file descriptor\n"
     "  --fd-const FDNUM     Output constellation and symbols to file descr\n"
     "  --fd-spectrum FDNUM  Output spectrum to file descr\n"
     "  --json               Use JSON syntax\n"
     );
#ifdef GUI
  fprintf
    (f,
     "  --gui             Show constellation and spectrum (X11)\n"
     "  --duration FLOAT  Width of timeline plot (s, default 60)\n"
     "  --linger          Keep GUI running after EOF\n"
     );
#endif
  fprintf
    (f, "\nTesting options:\n"
     "  --fd-pp FDNUM         Dump preprocessed IQ data to file descriptor\n"
     "  --fd-iqsymbols FDNUM  Dump sampled IQ symbols to file descriptor\n"
     "  --awgn FLOAT  Add white gaussian noise stddev (slow)\n"
     );
  if ( info ) fprintf(f, "** Error while processing '%s'\n", info);
  exit(c);
}

int main(int argc, const char *argv[]) {
  config cfg;

  for ( int i=1; i<argc; ++i ) {
    if      ( ! strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "-v") )
      cfg.verbose = true;
    else if ( ! strcmp(argv[i], "-d") ) {
      cfg.debug2 = cfg.debug;
      cfg.debug = true;
    }
    else if ( ! strcmp(argv[i], "--version") ) {
      printf("%s\n", VERSION);
      return 0;
    }
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
      else usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "--const") && i+1<argc ) {
      ++i;
      int c;
      for ( c=0; c<cstln_predef::COUNT; ++c )
	if ( ! strcmp(argv[i], cstln_names[c]) ) {
	  cfg.constellation = (cstln_predef)c;
	  break;
	}
      if ( c == cstln_predef::COUNT )
	usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "--strongpls") )
      cfg.strongpls = true;
    else if ( ! strcmp(argv[i], "--cr") && i+1<argc ) {
      ++i;
      int f;
      for ( f=0; f<FEC_COUNT; ++f )
	if ( ! strcmp(argv[i], fec_names[f]) ) {
	  cfg.fec = (code_rate)f;
	  break;
	}
      if ( f == FEC_COUNT )
	usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "--fastlock") )
      cfg.fastlock = true;
    else if ( ! strcmp(argv[i], "--viterbi") )
      cfg.viterbi = true;
    else if ( ! strcmp(argv[i], "--ldpc-bf") && i+1<argc )
      cfg.ldpc_bf = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--ldpc-helper") && i+1<argc )
      cfg.ldpc_helper = argv[++i];
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
    else if ( ! strcmp(argv[i], "--sampler") && i+1<argc ) {
      ++i;
      if      (!strcmp(argv[i],"nearest")) cfg.sampler = config::SAMP_NEAREST;
      else if (!strcmp(argv[i],"linear" )) cfg.sampler = config::SAMP_LINEAR;
      else if (!strcmp(argv[i],"rrc"    )) cfg.sampler = config::SAMP_RRC;
      else usage(argv[0], stderr, 1, argv[i]);
    }
    else if ( ! strcmp(argv[i], "--rrc-steps") && i+1<argc )
      cfg.rrc_steps = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--rrc-rej") && i+1<argc )
      cfg.rrc_rej = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--roll-off") && i+1<argc )
      cfg.rolloff = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--hq") ) {
      cfg.fastlock = true;
      cfg.viterbi = true;
      cfg.sampler = config::SAMP_RRC;
    }
    else if ( ! strcmp(argv[i], "--hs") )
      cfg.highspeed = true;
    else if ( ! strcmp(argv[i], "--anf") && i+1<argc )
      cfg.anf = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--cnr") )
      cfg.cnr = true;
    else if ( ! strcmp(argv[i], "--tune") && i+1<argc )
      cfg.Ftune = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--drift") )
      cfg.allow_drift = true;
    else if ( ! strcmp(argv[i], "--hdlc") )
      cfg.hdlc = true;
    else if ( ! strcmp(argv[i], "--packetized") )
      cfg.packetized = true;
#ifdef GUI
    else if  ( ! strcmp(argv[i], "--gui") )
      cfg.gui = true;
    else if ( ! strcmp(argv[i], "--duration") && i+1<argc )
      cfg.duration = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--linger") )
      cfg.linger = true;
#endif
    else if ( ! strcmp(argv[i], "--u8") )
      cfg.input_format = config::INPUT_U8;
    else if ( ! strcmp(argv[i], "--s8") )
      cfg.input_format = config::INPUT_S8;
    else if ( ! strcmp(argv[i], "--s12") )
      cfg.input_format = config::INPUT_S12;
    else if ( ! strcmp(argv[i], "--s16") )
      cfg.input_format = config::INPUT_S16;
    else if ( ! strcmp(argv[i], "--f32") )
      cfg.input_format = config::INPUT_F32;
    else if ( ! strcmp(argv[i], "--float-scale") && i+1<argc )
      cfg.float_scale = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--loop") )
      cfg.loop_input = true;
    else if ( ! strcmp(argv[i], "--inbuf")  && i+1<argc )
      cfg.input_buffer = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--buf-factor")  && i+1<argc )
      cfg.buf_factor = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--derotate") && i+1<argc )
      cfg.Fderot = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-pp") && i+1<argc )
      cfg.fd_pp = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--awgn") && i+1<argc )
      cfg.awgn = atof(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-info") && i+1<argc )
      cfg.fd_info = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-iqsymbols") && i+1<argc )
      cfg.fd_iqsymbols = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-const") && i+1<argc )
      cfg.fd_const = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--fd-spectrum") && i+1<argc )
      cfg.fd_spectrum = atoi(argv[++i]);
    else if ( ! strcmp(argv[i], "--json") )
      cfg.json = true;
    else
      usage(argv[0], stderr, 1, argv[i]);
  }

  if ( cfg.highspeed ) {
    fail("--hs mode is broken");
#if 0
    switch ( cfg.standard ) {
    case config::DVB_S:  return run_highspeed(cfg);
    case config::DVB_S2: return run_highspeed_s2(cfg);
    }
#endif
  }
  switch ( cfg.standard ) {
  case config::DVB_S:  return run_dvbs(cfg);
  case config::DVB_S2: return run_dvbs2(cfg);
  }
}
