// This file is part of LeanSDR Copyright (C) 2016-2020 <pabr@pabr.org>.
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

// TBD Interpolate frequency over whole frame, best-fit over pilots.
// TBD Interpolate gain between pilots / over whole frame.

#ifndef LEANSDR_DVBS2_H
#define LEANSDR_DVBS2_H

#include "leansdr/sdr.h"
#include "leansdr/softword.h"
#include "leansdr/dvb.h"
#include "leansdr/crc.h"
#include "leansdr/bch.h"
#include "leansdr/ldpc.h"

namespace leansdr {

#define DEBUG_CARRIER 0
#define DEBUG_LOOKAHEAD 0

  // S2 SOF
  // EN 302 307-1 section 5.5.2.1 SOF field

  template<typename T>
  struct s2_sof {
    static const uint32_t VALUE = 0x18d2e82;
    static const uint32_t MASK  = 0x3ffffff;
    static const int LENGTH = 26;
    complex<T> symbols[LENGTH];
    s2_sof() {
      for ( int s=0; s<LENGTH; ++s ) {
	int angle = ((VALUE>>(LENGTH-1-s))&1)*2 + (s&1);  // pi/2-BPSK
	symbols[s].re = cstln_amp * cosf(M_PI/4 + 2*M_PI*angle/4);
	symbols[s].im = cstln_amp * sinf(M_PI/4 + 2*M_PI*angle/4);
      }
    }
  };  // s2_sof


  // S2 PLS CODES
  // Precomputes the PLS code sequences.
  // EN 302 307-1 section 5.5.2.4 PLS code

  template<typename T>
  struct s2_plscodes {
    // PLS index format MODCOD[4:0]|SHORTFRAME|PILOTS
    static const int COUNT = 128;
    static const int LENGTH = 64;
    uint64_t codewords[COUNT];
    complex<T> symbols[COUNT][LENGTH];
    s2_plscodes() {
      uint32_t G[6] = { 0x55555555,
			0x33333333,
			0x0f0f0f0f,
			0x00ff00ff,
			0x0000ffff,
			0xffffffff };
      for ( int index=0; index<COUNT; ++index ) {
	uint32_t y = 0;
	for ( int row=0; row<6; ++row )
	  if ( (index>>(6-row)) & 1 ) y ^= G[row];
	uint64_t code = 0;
	for ( int bit=31; bit>=0; --bit ) {
	  int yi = (y>>bit) & 1;
	  if ( index & 1 ) code = (code<<2) | (yi<<1) | (yi^1);
	  else             code = (code<<2) | (yi<<1) | yi;
	}
	// Scrambling
	code ^= SCRAMBLING;
	// Store precomputed codeword.
	codewords[index] = code;
	// Also store as symbols.
	for ( int i=0; i<LENGTH; ++i ) {
	  int yi = (code>>(LENGTH-1-i)) & 1;
	  int nyi = yi ^ (i&1);
	  symbols[index][i].re = cstln_amp * (1-2*nyi) / sqrtf(2);
	  symbols[index][i].im = cstln_amp * (1-2* yi) / sqrtf(2);
	}
      }
    }
    static const uint64_t SCRAMBLING = 0x719d83c953422dfa;
  };  // s2_plscodes


  static const int PILOT_LENGTH = 36;

  // Date about pilots.
  // Mostly for consistency with s2_sof and s2_plscodes.

  template<typename T>
  struct s2_pilot {
    static const int LENGTH = PILOT_LENGTH;
    complex<T> symbol;
    s2_pilot() { symbol.re = symbol.im = cstln_amp*0.707107; }
  };  // s2_pilot

  // S2 SCRAMBLING
  // Precomputes the symbol rotations for PL scrambling.
  // EN 302 307-1 section 5.5.4 Physical layer scrambling

  struct s2_scrambling {
    uint8_t Rn[131072];  // 0..3  (* 2pi/4)
    s2_scrambling(int codenum=0) {
      uint32_t stx=0x00001, sty=0x3ffff;
      // x starts at codenum, wraps at index 2^18-1 by design
      for ( int i=0; i<codenum; ++i ) stx = lfsr_x(stx);
      // First half of sequence is LSB of scrambling angle
      for ( int i=0; i<131072; ++i ) {
	int zn = (stx ^ sty) & 1;
	Rn[i] = zn;
	stx = lfsr_x(stx);
	sty = lfsr_y(sty);
      }
      // Second half is MSB
      for ( int i=0; i<131072; ++i ) {
	int zn = (stx ^ sty) & 1;
	Rn[i] |= zn << 1;
	stx = lfsr_x(stx);
	sty = lfsr_y(sty);
      }
    }
    uint32_t lfsr_x(uint32_t X) {
      int bit = ((X>>7) ^ X) & 1;
      return ((bit<<18)|X) >> 1;
    }
    uint32_t lfsr_y(uint32_t Y) {
      int bit = ((Y>>10) ^ (Y>>7) ^ (Y>>5) ^ Y) & 1;
      return ((bit<<18)|Y) >> 1;
    }
  };  // s2_scrambling


  // S2 BBSCRAMBLING
  // Precomputes the xor pattern for baseband scrambling.
  // EN 302 307-1 section 5.2.2 BB scrambling

  struct s2_bbscrambling {
    s2_bbscrambling() {
      uint16_t st = 0x00a9;  // 000 0000 1010 1001  (Fig 5 reversed)
      for ( int i=0; i<sizeof(pattern); ++i ) {
	uint8_t out = 0;
	for ( int n=8; n--; ) {
	  int bit = ((st>>13) ^ (st>>14)) & 1;  // Taps
	  out = (out<<1) | bit;  // MSB first
	  st = (st<<1) | bit;    // Feedback
	}
	pattern[i] = out;
      }
    }
    void transform(const uint8_t *in, int bbsize, uint8_t *out) {
      for ( int i=0; i<bbsize; ++i )
	out[i] = in[i] ^ pattern[i];
    }
  private:
    uint8_t pattern[58192];  // Values 0..3
  };  // s2_bbscrambling


  // S2 PHYSICAL LAYER SIGNALLING

  struct s2_pls {
    int modcod;  // 0..31
    bool sf;
    bool pilots;
    int framebits() const { return sf ? 16200 : 64800; }
    bool is_dummy() { return (modcod==0); }
  };

  static const int PLSLOT_LENGTH = 90;

  template<typename SOFTSYMB>
  struct plslot {
    static const int LENGTH = PLSLOT_LENGTH;
    bool is_pls;
    union {
      s2_pls pls;
      SOFTSYMB symbols[LENGTH];
    };
  };


  // EN 302 307-1 section 5.5.2.2 MODCOD field
  // EN 302 307-1 section 6 Error performance

  const struct modcod_info {
    static const int MIN_SLOTS_PER_FRAME = 144;
    static const int MIN_SYMBOLS_PER_FRAME =
      (1+MIN_SLOTS_PER_FRAME) * PLSLOT_LENGTH;
    static const int MAX_SLOTS_PER_FRAME = 360;
    static const int MAX_SYMBOLS_PER_FRAME =
      (1+MAX_SLOTS_PER_FRAME) * PLSLOT_LENGTH +
      ((MAX_SLOTS_PER_FRAME-1)/16) * PILOT_LENGTH;
    int nslots_nf;     // Number of 90-symbol slots per normal frame
    int nsymbols;      // Symbols in the constellation
    cstln_base::predef c;
    code_rate rate;
    // Ideal Es/N0 for normal frames
    // EN 302 307 section 6 Error performance
    float esn0_nf;
    // Radii for APSK
    // EN 302 307, section 5.4.3, Table 9
    // EN 302 307, section 5.4.4, Table 10
    float g1, g2, g3;
  } modcod_infos[32] = {
    { 0, },
    // 1 - 11
    { 360,  4, cstln_base::QPSK,   FEC14,   -2.35 },
    { 360,  4, cstln_base::QPSK,   FEC13,   -1.24 },
    { 360,  4, cstln_base::QPSK,   FEC25,   -0.30 },
    { 360,  4, cstln_base::QPSK,   FEC12,    1.00 },
    { 360,  4, cstln_base::QPSK,   FEC35,    2.23 },
    { 360,  4, cstln_base::QPSK,   FEC23,    3.10 },
    { 360,  4, cstln_base::QPSK,   FEC34,    4.03 },
    { 360,  4, cstln_base::QPSK,   FEC45,    4.68 },
    { 360,  4, cstln_base::QPSK,   FEC56,    5.18 },
    { 360,  4, cstln_base::QPSK,   FEC89,    6.20 },
    { 360,  4, cstln_base::QPSK,   FEC910,   6.42 },
    // 12 - 17
    { 240,  8, cstln_base::PSK8,   FEC35,    5.50 },
    { 240,  8, cstln_base::PSK8,   FEC23,    6.62 },
    { 240,  8, cstln_base::PSK8,   FEC34,    7.91 },
    { 240,  8, cstln_base::PSK8,   FEC56,    9.35 },
    { 240,  8, cstln_base::PSK8,   FEC89,   10.69 },
    { 240,  8, cstln_base::PSK8,   FEC910,  10.98 },
    // 18 - 23
    { 180, 16, cstln_base::APSK16, FEC23,    8.97,   3.15 },
    { 180, 16, cstln_base::APSK16, FEC34,   10.21,   2.85 },
    { 180, 16, cstln_base::APSK16, FEC45,   11.03,   2.75 },
    { 180, 16, cstln_base::APSK16, FEC56,   11.61,   2.70 },
    { 180, 16, cstln_base::APSK16, FEC89,   12.89,   2.60 },
    { 180, 16, cstln_base::APSK16, FEC910,  13.13,   2.57 },
    // 24 - 28
    { 144, 32, cstln_base::APSK32, FEC34,   12.73,   2.84, 5.27 },
    { 144, 32, cstln_base::APSK32, FEC45,   13.64,   2.72, 4.87 },
    { 144, 32, cstln_base::APSK32, FEC56,   14.28,   2.64, 4.64 },
    { 144, 32, cstln_base::APSK32, FEC89,   15.69,   2.54, 4.33 },
    { 144, 32, cstln_base::APSK32, FEC910,  16.05,   2.53, 4.30 },
    // 29 - 31
    { 0, },
    { 0, },
    { 0, }
  };

  // Assert that a MODCOD number is valid
  const modcod_info *check_modcod(int m) {
    if ( m<0 || m>31 ) fail("Invalid MODCOD number");
    const modcod_info *r = &modcod_infos[m];
    if ( ! r->nslots_nf ) fail("Unsupported MODCOD");
    return r;
  }


  // S2 FRAME TRANSMITTER

  template<typename T>
  struct s2_frame_transmitter : runnable {
    s2_frame_transmitter(scheduler *sch,
			 pipebuf< plslot<hard_ss> > &_in,
			 pipebuf< complex<T> > &_out)
      : runnable(sch, "S2 frame transmitter"),
	in(_in), out(_out,modcod_info::MAX_SYMBOLS_PER_FRAME)
    {
      float amp = cstln_amp / sqrtf(2);
      qsymbols[0].re = +amp; qsymbols[0].im = +amp;
      qsymbols[1].re = +amp; qsymbols[1].im = -amp;
      qsymbols[2].re = -amp; qsymbols[2].im = +amp;
      qsymbols[3].re = -amp; qsymbols[3].im = -amp;
      // Clear the constellation cache.
      for ( int i=0; i<32; ++i ) pcsymbols[i] = NULL;
    }
    void run() {
      while ( in.readable() >= 1 ) {
	plslot<hard_ss> *pin = in.rd();
	if ( ! pin->is_pls ) fail("Expected PLS pseudo-slot");
	s2_pls *pls = &pin->pls;
	const modcod_info *mcinfo = check_modcod(pls->modcod);
	int nslots = (pls->sf ? mcinfo->nslots_nf/4 : mcinfo->nslots_nf);
	if ( in.readable() < 1+nslots ) break;
	// Require room for BBHEADER + slots + optional pilots.
	int nsymbols = ( (1+nslots) * pin->LENGTH +
			 (pls->pilots ? ((nslots-1)/16)*pilot.LENGTH : 0) );
	if ( out.writable() < nsymbols ) break;
	int nw = run_frame(pls, mcinfo, pin+1, nslots, out.wr());
	if ( nw != nsymbols ) fail("Bug: s2_frame_transmitter overflow");
	in.read(1+nslots);
	out.written(nsymbols);
      }
    }
    int run_frame(s2_pls *pls, const modcod_info *mcinfo,
		  const plslot<hard_ss> *pin, int nslots,
		  complex<T> *pout) {
      complex<T> *pout0 = pout;  // For sanity check
      // PLHEADER: SOF AND PLSCODE
      // EN 302 307-1 section 5.5.2 PL signalling
      memcpy(pout, sof.symbols, sof.LENGTH*sizeof(*pout));
      pout += sof.LENGTH;
      int pls_index = (pls->modcod<<2) | (pls->sf<<1) | pls->pilots;
      memcpy(pout, plscodes.symbols[pls_index], plscodes.LENGTH*sizeof(*pout));
      pout += plscodes.LENGTH;
      complex<T> *csymbols = get_csymbols(pls->modcod);
      // Slots and pilots
      int till_next_pilot = pls->pilots ? 16 : nslots;
      uint8_t *scr = &scrambling.Rn[0];
      for ( int S=0; S<nslots; ++S,++pin,--till_next_pilot ) {
	if ( till_next_pilot == 0 ) {
	  // Send pilot
	  for ( int s=0; s<pilot.LENGTH; ++s,++scr,++pout )
	    scramble(&qsymbols[0], *scr, pout);
	  till_next_pilot = 16;
	}
	// Send slot
	if ( pin->is_pls ) fail("s2_frame_transmitter: bad input sequence");
	const hard_ss *ps = pin->symbols;
	for ( int s=0; s<pin->LENGTH; ++s,++ps,++scr,++pout )
	  scramble(&csymbols[*ps], *scr, pout);
      }
      return pout - pout0;
    }
    inline void scramble(const complex<T> *src, uint8_t r, complex<T> *dst) {
      switch ( r ) {
      case 3: dst->re =  src->im; dst->im = -src->re; break;
      case 2: dst->re = -src->re; dst->im = -src->im; break;
      case 1: dst->re = -src->im; dst->im =  src->re; break;
      default: *dst = *src;
      }
    }
  private:
    pipereader< plslot<hard_ss> > in;
    pipewriter< complex<T> > out;
    complex<T> *pcsymbols[32];  // Constellations in use, indexed by modcod
    complex<T> *get_csymbols(int modcod) {
      if ( ! pcsymbols[modcod] ) {
	const modcod_info *mcinfo = check_modcod(modcod);
	if ( sch->debug )
	  fprintf(stderr, "Building constellation %s ratecode %d\n",
		  cstln_base::names[mcinfo->c], mcinfo->rate);
	// TBD Different Es/N0 for short frames ?
	cstln_lut<hard_ss,256> cstln(mcinfo->c, mcinfo->esn0_nf,
				     mcinfo->g1, mcinfo->g2, mcinfo->g3);
	pcsymbols[modcod] = new complex<T>[cstln.nsymbols];
	for ( int s=0; s<cstln.nsymbols; ++s ) {
	  pcsymbols[modcod][s].re = cstln.symbols[s].re;
	  pcsymbols[modcod][s].im = cstln.symbols[s].im;
	}
      }
      return pcsymbols[modcod];
    }
    complex<T> qsymbols[4];  // RMS cstln_amp
    s2_sof<T> sof;
    s2_pilot<T> pilot;
    s2_plscodes<T> plscodes;
    s2_scrambling scrambling;
  };  // s2_frame_transmitter


  // S2 FRAME RECEIVER

#define TEST_DIVERSITY 0

  template<typename T, typename SOFTSYMB>
  struct s2_frame_receiver : runnable {
    sampler_interface<T> *sampler;
    int meas_decimation;
    float Ftune;         // Tuning bias in cycles per symbol
    bool allow_drift;    // Unbounded carrier tracking
    float omega0;        // Samples per symbol
    bool strongpls;      // PL symbols at max amplitude (default: RMS)
    uint32_t modcods;    // Bitmask of desired modcods
    uint8_t framesizes;  // Bitmask of desired frame sizes
    bool fastlock;       // Synchronize more agressively
    bool fastdrift;      // Carrier drift faster than pilots
    float freq_tol;      // Tolerance on carrier frequency
    float sr_tol;        // Tolerance on symbol rate
    s2_frame_receiver(scheduler *sch,
		      sampler_interface<T> *_sampler,
		      pipebuf< complex<T> > &_in,
		      pipebuf< plslot<SOFTSYMB> > &_out,
		      pipebuf<float> *_freq_out=NULL,
		      pipebuf<float> *_ss_out=NULL,
		      pipebuf<float> *_mer_out=NULL,
		      pipebuf< complex<float> > *_cstln_out=NULL,
		      pipebuf< complex<float> > *_cstln_pls_out=NULL,
		      pipebuf< complex<float> > *_symbols_out=NULL,
		      pipebuf<int> *_state_out=NULL)
      : runnable(sch, "S2 frame receiver"),
	sampler(_sampler),
	meas_decimation(1048576),
	Ftune(0),
	allow_drift(false),
	strongpls(false),
	modcods(0xffffffff),
	framesizes(0x03),
	fastlock(false),
	fastdrift(false),
	freq_tol(0.25),
	sr_tol(100e-6),
	cstln(NULL),
	in(_in), out(_out,1+modcod_info::MAX_SLOTS_PER_FRAME),
	meas_count(0),
	freq_out(opt_writer(_freq_out)),
	ss_out(opt_writer(_ss_out)),
	mer_out(opt_writer(_mer_out)),
	cstln_out(opt_writer(_cstln_out,1024)),
	cstln_pls_out(opt_writer(_cstln_pls_out,1024)),
	symbols_out(opt_writer(_symbols_out,
			       modcod_info::MAX_SYMBOLS_PER_FRAME)),
	state_out(opt_writer(_state_out)),
	first_run(true),
	scrambling(0),
	pls_total_errors(0), pls_total_count(0)
    {
      float min_mer = -2.35;
      // Constellation for PLS
      bpsk = new cstln_lut<SOFTSYMB,256>(cstln_base::BPSK, min_mer);
      // Constellation for dummy frames
      qpsk = new cstln_lut<SOFTSYMB,256>(cstln_base::QPSK, min_mer);
      // Constellations for MODCODS.
      for ( int mc=0; mc<32; ++mc ) cstlns[mc] = NULL;

#if TEST_DIVERSITY
      fprintf(stderr, "** DEBUG: Diversity test mode (slower)\n");
#endif
    }

    enum {
      FRAME_DETECT,   // Looking for PLHEADER
      FRAME_PROBE,    // Aligned with PLHEADER, ready to recover carrier
      FRAME_LOCKED,   // Demodulating
    } state;

    // sampler_state holds the entire state of the PLL.
    // Useful for looking ahead (e.g. next pilots/SOF) then rewinding.

    struct sampler_state {
      complex<T> *p;  // Pointer to samples (update when entering run())
      float mu;       // Time of next symbol, counted from p
      float omega;    // Samples per symbol
      float gain;     // Scaling factor toward cstln_amp
      float ph16;     // Carrier phase at next symbol (cycles * 65536)
      float fw16;     // Carrier frequency (cycles per symbol * 65536)
      uint8_t *scr;   // Position in scrambling sequence for next symbol
      void normalize() {
	ph16 = fmodf(ph16, 65536.0f);  // Rounding direction irrelevant
      }
      void skip_symbols(int ns) {
	mu += omega * ns;
	p += (int)floorf(mu);
	mu -= floorf(mu);
	ph16 += fw16 * ns;
	normalize();
	scr += ns;
      }
      char *format() {
	static char buf[256];
	sprintf(buf, "%9.2lf %+6.0f ppm  %+3.0f °  %f",
		(double)((p-(complex<T>*)NULL)&262143)+mu,  // Arbitrary wrap
		fw16*1e6/65536, fmodfs(ph16,65536.0f)*360/65536, gain);
	return buf;
      }
    };

    float min_freqw16, max_freqw16;

    // State during FRAME_SEARCH and FRAME_LOCKED
    sampler_state ss_cache;

    void run() {
      if ( strongpls ) fail("--strongpls is broken.");
      // Require enough samples to detect one plheader,
      // TBD margin ?
      int min_samples = (1 + modcod_info::MAX_SYMBOLS_PER_FRAME +
			 sof.LENGTH+plscodes.LENGTH)*omega0 * 2;
      while ( in.readable() >= min_samples + sampler->readahead() &&
	      out.writable() >= 1+modcod_info::MAX_SLOTS_PER_FRAME &&
	      opt_writable(freq_out, 1) &&
	      opt_writable(ss_out, 1) &&
	      opt_writable(mer_out, 1) &&
	      opt_writable(symbols_out, modcod_info::MAX_SYMBOLS_PER_FRAME) &&
	      opt_writable(state_out, 1) ) {
	if ( first_run ) {
	  enter_frame_detect();
	  first_run = false;
	  continue;  // Must recheck writable().
	}
	switch ( state ) {
	case FRAME_DETECT: run_frame_detect(); break;
	case FRAME_PROBE: run_frame_probe(); break;
	case FRAME_LOCKED: run_frame_locked(); break;
	}
      }
    }

    // State transtion
    void enter_frame_detect() {
      state = FRAME_DETECT;
      // Setup sampler for interpolation only.
      ss_cache.fw16 = 65536 * Ftune;
      ss_cache.ph16 = 0;
      ss_cache.mu = 0;
      ss_cache.gain = 1;
      ss_cache.omega = omega0;
      // Set frequency tracking boundaries around tuning frequency.
      if ( allow_drift ) {
	// Track across the whole baseband.
	min_freqw16 = ss_cache.fw16 - omega0*65536;
	max_freqw16 = ss_cache.fw16 + omega0*65536;
      } else {
	min_freqw16 = ss_cache.fw16 - freq_tol*65536.0;
	max_freqw16 = ss_cache.fw16 + freq_tol*65536.0;
      }

      opt_write(state_out, 0);
      if ( sch->debug ) fprintf(stderr, "DETECT\n");

      if ( fastlock || first_run )
	discard = 0;
      else {
	// Discard some data so that CPU usage during PLHEADER detection
	// is at same level as during steady-state demodulation.
	// This has no effect if the first detection is successful.
	float duty_factor = 5;
	discard = ( modcod_info::MAX_SYMBOLS_PER_FRAME * omega0 *
		    (duty_factor+drand48()-0.5) );
      }
    }

    size_t discard;

    void run_frame_detect() {
      if ( discard ) {
	size_t d = min(discard, in.readable());
	in.read(d);
	discard -= d;
	return;
      }

      sampler->update_freq(ss_cache.fw16/omega0);

      const int search_range = modcod_info::MIN_SYMBOLS_PER_FRAME;
      ss_cache.p = in.rd();
      float q = find_plheader(&ss_cache, search_range);

#if DEBUG_CARRIER
      fprintf(stderr, "CARRIER diffcorr: %.0f%%  %s\n",
	      q*100, ss_cache.format());
#endif

      in.read(ss_cache.p-in.rd());

      enter_frame_probe();
    }

    void enter_frame_probe() {
      if ( sch->debug ) fprintf(stderr, "PROBE\n");
      state = FRAME_PROBE;
    }

    void run_frame_probe() {
      return run_frame_probe_locked();
    }

    void enter_frame_locked() {
      state = FRAME_LOCKED;
      if ( sch->debug ) fprintf(stderr, "LOCKED\n");
      opt_write(state_out, 1);
    }

    void run_frame_locked() {
      return run_frame_probe_locked();
    }

    // Process one frame.
    // Perform additional carrier estimation if state==FRAME_PROBE.

    void run_frame_probe_locked() {
      complex<float> *psampled;  // Data symbols (one per slot)
      if ( cstln_out && cstln_out->writable()>=1024 )
	psampled = cstln_out->wr();
      else
	psampled = NULL;
      complex<float> *psampled_pls;  // PLHEADER symbols
      if ( cstln_pls_out && cstln_pls_out->writable()>=1024 )
	psampled_pls = cstln_pls_out->wr();
      else
	psampled_pls = NULL;

#if TEST_DIVERSITY
      complex<float> *psymbols = symbols_out ? symbols_out->wr() : NULL;
      float scale_symbols = 1.0 / cstln_amp;
#endif

      sampler_state ss = ss_cache;
      ss.p = in.rd();
      ss.normalize();
      sampler->update_freq(ss.fw16/omega0);

#if DEBUG_CARRIER
      fprintf(stderr, "CARRIER frame: %s\n", ss.format());
#endif

      // Interpolate PLHEADER.

      const int PLH_LENGTH = sof.LENGTH + plscodes.LENGTH;
      complex<float> plh_symbols[PLH_LENGTH];
      llr_t plh_llrs[PLH_LENGTH];
      for ( int s=0; s<PLH_LENGTH; ++s ) {
	complex<float> p = interp_next(&ss) * ss.gain;
	plh_symbols[s] = p;
	if ( psampled_pls ) *psampled_pls++ = p;
	// Soft-decode as pi/2-BPSK.
	typename cstln_lut<SOFTSYMB,256>::result *cr =
	  (s&1) ? bpsk->lookup(p.im,-p.re) : bpsk->lookup(p.re,p.im);
	plh_llrs[s] = cr->ss.bits[0];   // TBD Don't assume SOFTSYMB==llr_ss.
      }

      // Decode SOF.

      int llr_sof = 0;
      {
	static const uint32_t highbit = (uint32_t)1 << (sof.LENGTH-1);
	for ( int i=0; i<sof.LENGTH; ++i ) {
	  llr_t l = plh_llrs[i];
	  if ( (sof.VALUE<<i) & highbit ) llr_sof -= l;  // Expected 1
	  else                            llr_sof += l;  // Expected 0
	}
      }

      // Decode PLSCODE.

      int plscode_index = 0;  // Avoid compiler warning
      int llr_plscode = -128 * plscodes.LENGTH;  // TBD Don't assume llr_ss.
      // Compute LLR for each candidate.
      for ( int c=0; c<plscodes.COUNT; ++c ) {
	uint64_t cw = plscodes.codewords[c];
	int llr = 0;
	static const uint64_t highbit = (uint64_t)1 << (plscodes.LENGTH-1);
	for ( int i=0; i<plscodes.LENGTH; ++i ) {
	  llr_t l = plh_llrs[sof.LENGTH+i];
	  if ( (cw<<i) & highbit ) llr -= l;  // Expected 1
	  else                     llr += l;  // Expected 0
	}
	if ( llr > llr_plscode ) { llr_plscode=llr; plscode_index=c; }
      }

#if DEBUG_CARRIER
      // Perform ML decoding of PLHEADER to confirm LLR-based result.
      {
	int plscode_index2 = 0;
	float best_m2 = -1;
	complex<T> expected[PLH_LENGTH];
	for ( int i=0; i<sof.LENGTH; ++i )
	  expected[i] = sof.symbols[i];
	for ( int c=0; c<plscodes.COUNT; ++c ) {
	  for ( int i=0; i<plscodes.LENGTH; ++i )
	    expected[sof.LENGTH+i] = plscodes.symbols[c][i];
	  sampler_state ssc = ss;
	  complex<T> symbs[PLH_LENGTH];  // Altered in-place
	  memcpy(symbs, plh_symbols, sizeof(symbs));
	  float m2 = match_ph_amp(expected, symbs, PLH_LENGTH, &ssc);
	  if ( m2 > best_m2 ) { best_m2=m2; plscode_index2=c; }
	}
	fprintf(stderr, "CARRIER MLPLH %d %d %4.1f dB %s\n",
		plscode_index, plscode_index2, 10*log10f(best_m2),
		(plscode_index==plscode_index2)?"OK":"NOK");
      }
#endif

      // Sanity-check on PLHEADER plausibility.
      // Possibly useful to detect a bad lock or end of carrier.
      if ( llr_sof+llr_plscode < 0 ) {
	if ( sch->debug2 )
	  fprintf(stderr, "Weak PLHEADER (%d+%d)\n", llr_sof, llr_plscode);
	in.read(ss.p-in.rd());
	enter_frame_detect();
	return;
      }

      // ss now points to first data slot.
      ss.scr = scrambling.Rn;

      complex<float> plh_expected[PLH_LENGTH];
      for ( int i=0; i<sof.LENGTH; ++i )
	plh_expected[i] = sof.symbols[i];
      for ( int i=0; i<plscodes.LENGTH; ++i )
	plh_expected[sof.LENGTH+i] = plscodes.symbols[plscode_index][i];

      if ( state == FRAME_PROBE ) {
	// Carrier frequency from differential detector is still not reliable.
	// Use known PLH symbols to improve.
	match_freq(plh_expected, plh_symbols, PLH_LENGTH, &ss);
#if DEBUG_CARRIER
	fprintf(stderr, "CARRIER freq: %s\n", ss.format());
#endif
      }

      // Use known PLH symbols to estimate carrier phase and amplitude.

      float mer2_plh = match_ph_amp(plh_expected, plh_symbols, PLH_LENGTH, &ss);
      float mer_plh_db = 10 * log10f(mer2_plh);
#if DEBUG_CARRIER
      fprintf(stderr, "CARRIER plheader: %s MER %.1f dB\n",
	      ss.format(), mer_plh_db);
#endif

      // Parse PLSCODE.

      s2_pls pls;
      pls.modcod = plscode_index >> 2;  // Guaranteed 0..31
      pls.sf = plscode_index & 2;
      pls.pilots = plscode_index & 1;

      plslot<SOFTSYMB> *pout=out.wr(), *pout0=pout;

      if  ( sch->debug2 )
	fprintf(stderr, "PLS: mc=%2d, sf=%d, pilots=%d llr=%3d %4.1f dB ",
		pls.modcod, pls.sf, pls.pilots,
		llr_sof+llr_plscode, mer_plh_db);

      // Determine contents of frame.

      int S;                            // Data slots in this frame
      cstln_lut<SOFTSYMB,256> *dcstln;  // Constellation for data slots

      if ( pls.is_dummy() ) {
	S = 36;
	dcstln = qpsk;
      } else {
	const modcod_info *mcinfo = &modcod_infos[pls.modcod];
	if ( ! mcinfo->nslots_nf ) {
	  fprintf(stderr, "Unsupported or corrupted MODCOD\n");
	  in.read(ss.p-in.rd());
	  enter_frame_detect();
	  return;
	}
#if 0
	if ( mer_plh_db < mcinfo->esn0_nf - 1.0f ) {
	  // TBD Try to decode anyway, or drop frame, or unlock ?
	  if ( sch->debug ) fprintf(stderr, "Insufficient MER\n");
	  in.read(ss.p-in.rd());
	  enter_frame_detect();
	  return;
	}
#endif
	if ( pls.sf && mcinfo->rate==FEC910 ) {  // TBD use fec_infos
	  fprintf(stderr, "Unsupported or corrupted FEC\n");
	  in.read(ss.p-in.rd());
	  enter_frame_detect();
	  return;
	}
	S = pls.sf ? mcinfo->nslots_nf/4 : mcinfo->nslots_nf;
	// Constellation for data slots.
	dcstln = get_cstln(pls.modcod);
	cstln = dcstln;  // Used by GUI
	// Output special slot with PLS information.
	pout->is_pls = true;
	pout->pls = pls;
	++pout;
      }

      // Now we know the frame structure.
      // Estimate carrier over known symbols.

#if DEBUG_LOOKAHEAD
      static int plh_counter = 0;  // For debugging only
      fprintf(stderr, "\nLOOKAHEAD %d PLH sr %+3.0f  %s  %.1f dB\n",
	      plh_counter, (1-ss.omega/omega0)*1e6,
	      ss.format(), 10*log10f(mer2_plh));
#endif
      // Next SOF
      sampler_state ssnext;
      float mer2_next;
      {
	ssnext = ss;
	int ns = ( S*PLSLOT_LENGTH +
		   (pls.pilots?((S-1)/16)*pilot.LENGTH:0) );
	// Find next SOP at expected position +- tolerance on symbol rate.
	int ns_tol = lrintf(ns*sr_tol);
	ssnext.omega = omega0;
	ssnext.skip_symbols(ns-ns_tol);
	float q = find_plheader(&ssnext, ns_tol*2);
	// Don't trust frequency from differential correlator.
	// Our current estimate should be better.
	ssnext.fw16 = ss.fw16;
	mer2_next = interp_match_sof(&ssnext);
#if DEBUG_CARRIER
	fprintf(stderr, "\nCARRIER next: %.0f%%  %s  %.1f dB\n",
		q*100, ssnext.format(), 10*log10f(mer2_next));
#else
	(void)mer2_next;
#endif
	// Estimate symbol rate (considered stable over whole frame)
	float dist = (ssnext.p-ss.p) + (ssnext.mu-ss.mu);
	// Set symbol period accordingly.
	ss.omega = dist / (ns+sof.LENGTH);
      }
      // Pilots
      int npilots = (S-1) / 16;
      sampler_state sspilots[npilots];
      float mer2_pilots = 0;
      int llr_pilots = 0;
      // Detect pilots
      if ( pls.pilots ) {
	sampler_state ssp = ss;
	for ( int i=0; i<npilots; ++i ) {
	  ssp.skip_symbols(16*PLSLOT_LENGTH);
	  float mer2 = interp_match_pilot(&ssp);
	  sspilots[i] = ssp;
	  mer2_pilots += mer2;
#if DEBUG_LOOKAHEAD
	  fprintf(stderr, "LOOKAHEAD %d PILOT%02d  sr %+3.0f  %s  %.1f dB\n",
		  plh_counter, i, (1-ssp.omega/omega0)*1e6,
		  ssp.format(), 10*log10f(mer2));
#endif
	}
      }

      if ( pls.pilots ) {
	// Measure average frequency, using pilots to unwrap.
	float totalph = 0;
	float prevph = ss.ph16;
	int span = 16*PLSLOT_LENGTH + pilot.LENGTH;
	// Wrap phase around current freq estimate, pilot by pilot.
	for ( int i=0; i<npilots; ++i ) {
	  float dph = sspilots[i].ph16 - (prevph+ss.fw16*span);
	  totalph += fmodfs(dph, 65536.0f);
	  prevph = sspilots[i].ph16;
	}
	// TBD Use data between last pilot and next SOF ?
	// Tricky when there is still an integer ambiguity.
	ss.fw16 += totalph / (span*npilots);
      } else {
	// Measure average frequency over whole frame.
	// Mostly useful for null frames.
	int span = S*plslot<SOFTSYMB>::LENGTH + sof.LENGTH;
	float dph = ssnext.ph16 - (ss.ph16+ss.fw16*span);
	ss.fw16 += fmodfs(dph,65536.0f) / span;
      }
#if DEBUG_LOOKAHEAD
      fprintf(stderr, "LOOKAHEAD %d NEXTSOF  sr %+3.0f  %s  %.1f dB\n",
	      plh_counter, (1-ssnext.omega/omega0)*1e6,
	      ssnext.format(), 10*log10f(mer2_pilots/npilots));
      ++plh_counter;
#endif

      if ( state == FRAME_PROBE ) {
	float fw0 = ss.fw16;
	match_frame(&ss, &pls, S, dcstln);
	// Apply retroactively from midpoint of pilots and next SOF.
	float fw_adj = ss.fw16 - fw0;
	if ( pls.pilots )
	  for ( int i=0; i<npilots; ++i )
	    sspilots[i].ph16 += fw_adj * PILOT_LENGTH / 2;
	ssnext.ph16 += fw_adj * sof.LENGTH / 2;
#if DEBUG_CARRIER
	fprintf(stderr, "CARRIER disambiguated:  %s\n", ss.format());
#endif
      }

      // TBD In FRAME_LOCKED, match_frame with sliprange +-1
      // to avoid broken locks ?

      // Per-frame statistics on carrier frequency.
      // Useful at very low symbol rates to detect situations where
      // the carrier fluctuates too fast for pilot-aided recovery.
      statistics<float> freq_stats;

      // Read slots and pilots

#if DEBUG_CARRIER
      fprintf(stderr, "CARRIER data: %s\n", ss.format());
#endif

      int pilot_errors = 0;

      for ( int slot=0; slot<S; ++slot,++pout ) {
	// Pilot before this slot ?
	if ( pls.pilots && !(slot&15) && slot ) {
	  ss.skip_symbols(pilot.LENGTH);
	  ss.ph16 = sspilots[slot/16-1].ph16;
	}

	// Time for pilot-aided carrier recovery ?
	if ( pls.pilots && !(slot&15) && slot+16<S ) {
	  // Sequence of data slots followed by pilots
	  sampler_state *ssp = &sspilots[slot/16];
	  int span = 16*pout->LENGTH + pilot.LENGTH;
	  float dph = ssp->ph16 - (ss.ph16+ss.fw16*span);
	  ss.fw16 += fmodfs(dph,65536.0f) / span;
	} else if ( ( pls.pilots && !(slot&15) && slot+16>=S ) ||
		    ( !pls.pilots && !slot ) ) {
	  // Sequence of data slots followed by SOF
	  int span = (S-slot)*pout->LENGTH + sof.LENGTH;
	  float dph = ssnext.ph16 - (ss.ph16+ss.fw16*span);
	  ss.fw16 += fmodfs(dph,65536.0f) / span;
	}

	// Read slot.
	freq_stats.add(ss.fw16);
	pout->is_pls = false;
	complex<float> p;  // Export last symbols for cstln_out
	for ( int s=0; s<pout->LENGTH; ++s ) {
	  p = interp_next(&ss) * ss.gain;
#if TEST_DIVERSITY
	  if ( psymbols )
	    *psymbols++ = p * scale_symbols;
#endif
	  if ( !pls.pilots || fastdrift )
	    (void)track_symbol(&ss, p, dcstln);  // SLOW
	  complex<float> d = descramble(&ss, p);
#if 1  // Slow
	  SOFTSYMB *symb = &dcstln->lookup(d.re, d.im)->ss;
#else  // Avoid scaling floats. May wrap at very low SNR.
	  SOFTSYMB *symb = &dcstln->lookup((int)d.re, (int)d.im)->ss;
#endif
	  pout->symbols[s] = *symb;
	}
	ss.normalize();
	if ( psampled ) *psampled++ = p;
      }  // slot

      if ( sch->debug2 )
	fprintf(stderr, "sr%+.0f fs=%.0f %4.1f dB\n",
		(1-ss.omega/omega0)*1e6,
		(freq_stats.max()-freq_stats.min())*1e6/65536.0f,
		pls.pilots?10*log10f(mer2_pilots/npilots):-99);

      if ( ss.fw16<min_freqw16 || ss.fw16>max_freqw16 ) {
	if ( sch->debug ) fprintf(stderr, "Carrier out of bounds\n");
	in.read(ss.p-in.rd());
	enter_frame_detect();
	return;
      }

      // Commit whole frame after final SOF.
      if ( ! pls.is_dummy() ) {
	if ( (modcods&(1<<pls.modcod)) && (framesizes&(1<<pls.sf)) )
	  out.written(pout-pout0);
	else {
	  if ( sch->debug2 )
	    fprintf(stderr, "modcod %d size %d rejected\n", pls.modcod, pls.sf);
	}
      }

      // Write back sampler progress
      meas_count += ss.p - in.rd();
      in.read(ss.p-in.rd());
      ss_cache = ss;

      // Measurements
      if ( psampled ) cstln_out->written(psampled-cstln_out->wr());
      if ( psampled_pls ) cstln_pls_out->written(psampled_pls-cstln_pls_out->wr());
#if TEST_DIVERSITY
      if ( psymbols ) symbols_out->written(psymbols-symbols_out->wr());
#endif
      if ( meas_count >= meas_decimation ) {
	opt_write(freq_out, ss_cache.fw16/65536/ss_cache.omega);
	opt_write(ss_out, cstln_amp / ss_cache.gain);
	opt_write(mer_out, mer2_plh ? 10*log10f(mer2_plh) : -99);
	meas_count -= meas_decimation;
      }

      if ( state == FRAME_PROBE ) {
	// First frame completed successfully.  Validate the lock.
	enter_frame_locked();
      }
    }

    // Find most likely PLHEADER between *pss and *pss+search_range
    // by differential correlation.
    // Align symbol timing.
    // Estimate carrier frequency (to about +-10kppm).
    // Estimate carrier phase (to about +-PI/4).
    // Adjust *ss accordingly.
    // Initialize AGC.
    // Return confidence (unbounded, 0=bad, 1=nominal).

    // TBD: When the buffer contains multiple equally-good candidates
    //      (this can happen with QPSK short frames or 8PSK or 16/32APSK),
    //      select the earliest one, in order to maximize data recovery.

    float find_plheader(sampler_state *pss, int search_range) {
      complex<T> best_corr = 0;
      int best_imu = 0;  // Avoid compiler warning.
      int best_pos = 0;  // Avoid compiler warning.

      // Symbol clock is not recovered yet, so we try fractional symbols.
      const int interp = 8;

      for ( int imu=0; imu<interp; ++imu ) {
	const int ndiffs = search_range + sof.LENGTH + plscodes.LENGTH;
	complex<T> diffs[ndiffs];
	sampler_state ss = *pss;
	ss.mu += imu * ss.omega / interp;
	// Compute rotation between consecutive symbols.
	complex<T> prev = 0;
	for ( int i=0; i<ndiffs; ++i ) {
	  complex<T> p = interp_next(&ss);
	  diffs[i] = conjprod(prev, p);
	  prev = p;
	}
	// Find best PLHEADER candidate position.
	const int ncorrs = search_range;
	for ( int i=0; i<ncorrs; ++i ) {
	  complex<T> c = correlate_plheader_diff(&diffs[i]);
	  //if ( cnorm2(c) > cnorm2(best_corr) ) {
	  // c.im>0 enforces frequency error +-Fm/4
	  if ( cnorm2(c)>cnorm2(best_corr) && c.im>0 ) {
	    best_corr = c;
	    best_imu = imu;
	    best_pos = i;
	  }
	}
      }  // imu

      // Setup sampler according to best match.
      pss->mu += best_imu * pss->omega / interp;
      pss->skip_symbols(best_pos);
      pss->normalize();

#if 0
      // Lowpass-filter the differential correlator.
      // (Does not help much.)
      static complex<float> acc = 0;
      static const float k = 0.05;
      acc = best_corr*k + acc*(1-k);
      best_corr = acc;
#endif

      // Get rough estimate of carrier frequency from differential correlator.
      // (best_corr is nominally +j).
      float freqw = atan2f(-best_corr.re, best_corr.im);
      pss->fw16 += freqw * 65536 / (2*M_PI);
      // Force refresh because correction may be large.
      sampler->update_freq(pss->fw16/omega0);

      // Interpolate SOF with initial frequency estimation.
      // Estimate phase by correlation.
      // Initialize AGC (naive).
      float q;
      {
	sampler_state ss = *pss;
	float power = 0;
	complex<T> symbs[sof.LENGTH];
	for ( int i=0; i<sof.LENGTH; ++i ) {
	  symbs[i] = interp_next(&ss);
	  power += cnorm2(symbs[i]);
	}
	complex<float> c = conjprod(sof.symbols, symbs, sof.LENGTH);
	c *= 1.0f / sof.LENGTH;
	align_phase(pss, c);
	float signal_amp = sqrtf(power/sof.LENGTH);
	q = sqrtf(cnorm2(c)) / (cstln_amp*signal_amp);
	pss->gain = cstln_amp / signal_amp;
      }
      return q;
    }

    // Correlate PLHEADER.

    complex<float> correlate_plheader_diff(complex<T> *diffs) {
      complex<T> csof = correlate_sof_diff(diffs);
      complex<T> cplsc = correlate_plscode_diff(&diffs[sof.LENGTH]);
      // Use csof+cplsc or csof-cplsc, whichever maximizes likelihood.
      complex<T> c0 = csof + cplsc;  // Best when b7==0 (pilots off)
      complex<T> c1 = csof - cplsc;  // Best when b7==1 (pilots on)
      complex<T> c = (cnorm2(c0)>cnorm2(c1)) ? c0 : c1;
      return c * (1.0f/(26-1+64/2));
    }

    // Correlate 25 differential transitions in SOF.

    complex<float> correlate_sof_diff(complex<T> *diffs) {
      complex<T> c = 0;
      const uint32_t dsof = sof.VALUE ^ (sof.VALUE>>1);
      for ( int i=0; i<sof.LENGTH; ++i ) {
	// Constant  odd bit => +PI/4
	// Constant even bit => -PI/4
	// Toggled   odd bit => -PI/4
	// Toggled  even bit => +PI/4
	if ( ((dsof>>(sof.LENGTH-1-i)) ^ i) & 1 ) c += diffs[i];
	else                                      c -= diffs[i];
      }
      return c;
    }

    // Correlate 32 data-independent transitions in PLSCODE.

    complex<float> correlate_plscode_diff(complex<T> *diffs) {
      complex<T> c = 0;
      uint64_t dscr = plscodes.SCRAMBLING ^ (plscodes.SCRAMBLING>>1);
      for ( int i=1; i<plscodes.LENGTH; i+=2 ) {
	if ( (dscr>>(plscodes.LENGTH-1-i)) & 1 ) c -= diffs[i];
	else                                     c += diffs[i];
      }
      return c;
    }

    // Adjust carrier frequency in *pss to match target phase after ns symbols.

    void adjust_freq(sampler_state *pss, int ns, const sampler_state *pnext) {
      // Note: Minimum deviation from current estimate.
      float adj = fmodfs(pnext->ph16-(pss->ph16+pss->fw16*ns),65536.0f) / ns;
      // 2sps vcm avec adj 10711K
      // 2sps vcm avec adj/2 10265k
      // 2sps vcm sans adj 8060K
      // 5sps 2MS1 avec adj 2262K
      // 5sps 2MS1 avec adj/2 3109K
      // 5sps 2MS1 avec adj/3 3254K
      // 5sps 2MS1 sans adj 3254K
      // 1.2sps oldbeacon avec adj 5774K
      // 1.2sps oldbeacon avec adj/2 15M
      // 1.2sps oldbeacon avec adj/3 15.6M max 17M
      // 1.2sps oldbeacon sans adj 14M malgré drift
      pss->fw16 += adj;
    }

    // Estimate frequency from known symbols by differential correlation.
    // Adjust *ss accordingly.
    // Retroactively frequency-shift the samples in-place.
    //
    // *ss must point after the sampled symbols.
    // expect[] must be PSK with amplitude cstln_amp.
    // Idempotent.
    // Not affected by phase error nor by gain mismatch.
    // Can handle +-25% error.
    // Result is typically such that residual phase error over a PLHEADER
    // spans less than 90°.

    void match_freq(complex<float> *expect, complex<float> *recv, int ns,
		    sampler_state *ss) {
      complex<float> diff = 0;
      for ( int i=0; i<ns-1; ++i ) {
	complex<float> de = conjprod(expect[i], expect[i+1]);
	complex<float> dr = conjprod(recv[i], recv[i+1]);
	diff += conjprod(de, dr);
      }
      float dfw16 = atan2f(diff.im,diff.re) * 65536 / (2*M_PI);
      // Derotate.
      for ( int i=0; i<ns; ++i )
	recv[i] *= trig.expi(-dfw16*i);
      ss->fw16 += dfw16;
      ss->ph16 += dfw16 * ns;  // Retroactively
    }

    // Interpolate a pilot block.
    // Perform ML match.

    float interp_match_pilot(sampler_state *pss) {
      complex<T> symbols[pilot.LENGTH];
      complex<T> expected[pilot.LENGTH];
      for ( int i=0; i<pilot.LENGTH; ++i ) {
	complex<float> p = interp_next(pss) * pss->gain;
	symbols[i] = descramble(pss, p);
	expected[i] = pilot.symbol;
	//fprintf(stderr, "%f %f\n", symbols[i].re, symbols[i].im);
      }
      return match_ph_amp(expected, symbols, pilot.LENGTH, pss);
    }

    // Interpolate a SOF.
    // Perform ML match.

    float interp_match_sof(sampler_state *pss) {
      complex<T> symbols[pilot.LENGTH];
      for ( int i=0; i<sof.LENGTH; ++i ) {
	symbols[i] = interp_next(pss) * pss->gain;
      }
      return match_ph_amp(sof.symbols, symbols, sof.LENGTH, pss);
    }

    // Estimate phase and amplitude from known symbols.
    // Adjust *ss accordingly.
    // Retroactively derotate and scale the samples in-place.
    // Return MER^2.
    //
    // *ss must point after the sampled symbols, with gain applied.
    // expect[] must be PSK (not APSK) with amplitude cstln_amp.
    // Idempotent.

    float match_ph_amp(complex<float> *expect, complex<float> *recv, int ns,
			  sampler_state *ss) {
      complex<float> rr = 0;
      for ( int i=0; i<ns; ++i )
	rr += conjprod(expect[i], recv[i]);
      rr *= 1.0f / (ns*cstln_amp);
      float dph16 = atan2f(rr.im,rr.re) * 65536 / (2*M_PI);
      ss->ph16 += dph16;
      rr *= trig.expi(-dph16);
      // rr.re is now the modulation amplitude.
      float dgain = cstln_amp / rr.re;
      ss->gain *= dgain;
      // Rotate and scale.  Compute error power.
      complex<float> adj = trig.expi(-dph16) * dgain;
      float ev2 = 0;
      for ( int i=0; i<ns; ++i ) {
	recv[i] *= adj;
	ev2 += cnorm2(recv[i]-expect[i]);
      }
      // Return MER^2.
      ev2 *= 1.0f / (ns*cstln_amp*cstln_amp);
      return 1.0f / ev2;
    }

    // Adjust frequency in *pss by an integer number of cycles per data block
    // so that phases align best on data symbols.
    //
    // *pss must point to the beginning of a frame.

    void match_frame(sampler_state *pss, const s2_pls *pls, int S,
		     cstln_lut<SOFTSYMB,256> *dcstln) {
      // With pilots: Use first block of data slots.
      // Without pilots: Use whole frame.
      int ns = pls->pilots ? 16*90 : S*90;
      // Pilots: steps of ~700 ppm (= 1 cycle between pilots)
      // No pilots, normal frames: steps of ~30(QPSK) - ~80(32APSK) ppm
      // No pilots, short frames: steps of ~120(QPSK) - 300(32APSK) ppm
      int nwrap = pls->pilots ? 16*90+pilot.LENGTH : S*90+sof.LENGTH;
      // Frequency search range.
      int sliprange = pls->pilots ? 10 : 50;  // TBD Customizable ?
      float besterr = 1e99;
      float bestslip = 0;  // Avoid compiler warning
      for ( int slip=-sliprange; slip<=sliprange; ++slip ) {
	sampler_state ssl = *pss;
	float dfw = slip * 65536.0f / nwrap;
	ssl.fw16 += dfw;
	// Apply retroactively from midpoint of preceeding PLHEADER,
	// where the phase from match_ph_amp is most reliable.
	ssl.ph16 += dfw * (plscodes.LENGTH+sof.LENGTH)/2;
	float err = 0;
	for ( int s=0; s<ns; ++s ) {
	  complex<float> p = interp_next(&ssl) * ssl.gain;
	  typename cstln_lut<SOFTSYMB,256>::result *cr =
	    dcstln->lookup(p.re, p.im);
	  complex<int8_t> &cp = dcstln->symbols[cr->symbol];
	  complex<float> ev(p.re-cp.re, p.im-cp.im);
	  err += cnorm2(ev);
	}
	err /= ns * cstln_amp * cstln_amp;
	if ( err < besterr ) { besterr=err; bestslip=slip; }
#if DEBUG_CARRIER
	fprintf(stderr, "slip %+3d  %6.0f ppm  err=%f\n",
		slip, ssl.fw16*1e6/65536, err);
#endif
      }
      pss->fw16 += bestslip * 65536.0f / nwrap;
    }  // match_frame

    void shutdown() {
      if ( sch->verbose )
	fprintf(stderr, "PL errors: %d/%d (%.0f ppm)\n",
		pls_total_errors, pls_total_count,
		1e6 * pls_total_errors / pls_total_count);
    }

    complex<float> descramble(sampler_state *ss, const complex<float> &p) {
      int r = *ss->scr++;
      complex<float> res;
      switch ( r ) {
      case 3: res.re = -p.im; res.im =  p.re; break;
      case 2: res.re = -p.re; res.im = -p.im; break;
      case 1: res.re =  p.im; res.im = -p.re; break;
      default: res = p;
      }
      return res;
    }

    // Interpolator

    inline complex<float> interp_next(sampler_state *ss) {
      // Skip to next sample
      while ( ss->mu >= 1 ) { ++ss->p; ss->mu-=1.0f; }
      // Interpolate
#if 0
      // Interpolate linearly then derotate.
      // This will fail with large carrier offsets (e.g. --tune).
      float cmu = 1.0f - ss->mu;
      complex<float> s(ss->p[0].re*cmu + ss->p[1].re*ss->mu,
		       ss->p[0].im*cmu + ss->p[1].im*ss->mu);
      ss->mu += ss->omega;
      // Derotate
      const complex<float> &rot = trig.expi(-ss->ph16);
      ss->ph16 += ss->fw16;
      return rot * s;
#else
      // Use generic interpolator
      complex<float> s = sampler->interp(ss->p, ss->mu, ss->ph16);
      ss->mu += ss->omega;
      ss->ph16 += ss->fw16;
      return s;
#endif
    }

    // Adjust phase in [ss] to cancel offset observed as [c].

    void align_phase(sampler_state *ss, const complex<float> &c) {
      float err = atan2f(c.im,c.re) * 65536 / (2*M_PI);
      ss->ph16 += err;
    }

    inline uint8_t track_symbol(sampler_state *ss, const complex<float> &p,
				cstln_lut<SOFTSYMB,256> *c) {
      static const float kph = 4e-2;
      static const float kfw = 1e-4;
      // Decision
      typename cstln_lut<SOFTSYMB,256>::result *cr = c->lookup(p.re, p.im);
      // Carrier tracking
      ss->ph16 += cr->phase_error * kph;
      ss->fw16 += cr->phase_error * kfw;
      if ( ss->fw16 < min_freqw16 ) ss->fw16 = min_freqw16;
      if ( ss->fw16 > max_freqw16 ) ss->fw16 = max_freqw16;
      return cr->symbol;
    }

  public:
    cstln_lut<SOFTSYMB,256> *bpsk;  // For PLHEADER
    cstln_lut<SOFTSYMB,256> *qpsk;  // For dummy frames
    s2_plscodes<T> plscodes;
    cstln_lut<SOFTSYMB,256> *cstlns[32];  // Constellations in use, by modcod
    cstln_lut<SOFTSYMB,256> *get_cstln(int modcod) {
      if ( ! cstlns[modcod] ) {
	const modcod_info *mcinfo = &modcod_infos[modcod];
	if ( sch->debug )
	  fprintf(stderr, "Creating LUT for %s ratecode %d\n",
		  cstln_base::names[mcinfo->c], mcinfo->rate);
	// Assume worst-case SNR for this MODCOD.
	// TBD Can this be detrimental at high SNR ?
	// TBD Needs a different value for short frames ?
	cstlns[modcod] = new cstln_lut<SOFTSYMB,256>
	  (mcinfo->c, mcinfo->esn0_nf, mcinfo->g1,mcinfo->g2,mcinfo->g3);
#if 0
	fprintf(stderr, "Dumping constellation LUT to stdout.\n");
	cstlns[modcod]->dump(stdout);
	exit(0);
#endif
      }
      return cstlns[modcod];
    }
    cstln_lut<SOFTSYMB,256> *cstln;  // Last seen, or NULL (legacy)

    trig16 trig;
    pipereader< complex<T> > in;
    pipewriter< plslot<SOFTSYMB> > out;
    int meas_count;
    pipewriter<float> *freq_out, *ss_out, *mer_out;
    pipewriter< complex<float> > *cstln_out;
    pipewriter< complex<float> > *cstln_pls_out;
    pipewriter< complex<float> > *symbols_out;
    pipewriter<int> *state_out;
    bool first_run;
    // S2 constants
    s2_scrambling scrambling;
    s2_sof<T> sof;
    s2_pilot<T> pilot;
    // Performance stats on PL signalling
    uint32_t pls_total_errors, pls_total_count;
  };  // s2_frame_receiver


  template<typename SOFTBYTE>
  struct fecframe {
    s2_pls pls;
    SOFTBYTE bytes[64800/8];  // Contains 16200/8 or 64800/8 bytes.
  };


  // S2 INTERLEAVER
  // EN 302 307-1 section 5.3.3 Bit Interleaver

  struct s2_interleaver : runnable {
    s2_interleaver(scheduler *sch,
		   pipebuf< fecframe<hard_sb> > &_in,
		   pipebuf< plslot<hard_ss> > &_out)
      : runnable(sch, "S2 interleaver"),
	in(_in), out(_out,1+360)
    {
    }
    void run() {
      while ( in.readable() >= 1 ) {
	const s2_pls *pls = &in.rd()->pls;
	const modcod_info *mcinfo = check_modcod(pls->modcod);
	int nslots = pls->sf ? mcinfo->nslots_nf/4 : mcinfo->nslots_nf;
	if ( out.writable() < 1+nslots ) return;
	const hard_sb *pbytes = in.rd()->bytes;
	// Output pseudo slot with PLS.
	plslot<hard_ss> *ppls = out.wr();
	ppls->is_pls = true;
	ppls->pls = *pls;
	out.written(1);
	// Interleave
	plslot<hard_ss> *pout = out.wr();
	if ( mcinfo->nsymbols == 4 )
	  serialize_qpsk(pbytes, nslots, pout);
	else {
	  int bps = log2(mcinfo->nsymbols);
	  int rows = pls->framebits() / bps;
	  if ( mcinfo->nsymbols==8 && mcinfo->rate==FEC35 )
	    interleave(bps, rows, pbytes, nslots, false, pout);
	  else
	    interleave(bps, rows, pbytes, nslots, true, pout);
	}
	in.read(1);
	out.written(nslots);
      }
    }
  private:

    // Fill slots with serialized QPSK symbols, MSB first.
    static void serialize_qpsk(const hard_sb *pin, int nslots,
			       plslot<hard_ss> *pout) {
#if 0  // For reference
      hard_sb acc;
      int nacc = 0;
      for ( ; nslots; --nslots,++pout ) {
	pout->is_pls = false;
	hard_ss *ps = pout->symbols;
	for ( int ns=pout->LENGTH; ns--; ++ps ) {
	  if ( nacc < 2 ) { acc=*pin++; nacc=8; }
	  *ps = acc>>6;
	  acc <<= 2;
	  nacc -= 2;
	}
      }
      if ( nacc ) fail("Bug: s2_interleaver");
#else
      if ( nslots % 2 ) fatal("Bug: Truncated byte");
      for ( ; nslots; nslots-=2 ) {
	hard_sb b;
	hard_ss *ps;
	// Slot 0 (mod 2)
	pout->is_pls = false;
	ps = pout->symbols;
	for ( int i=0; i<22; ++i ) {
	  b = *pin++;
	  *ps++ = (b>>6);
	  *ps++ = (b>>4) & 3;
	  *ps++ = (b>>2) & 3;
	  *ps++ = (b   ) & 3;
	}
	b = *pin++;
	*ps++ = (b>>6);
	*ps++ = (b>>4) & 3;
	// Slot 1 (mod 2)
	++pout;
	pout->is_pls = false;
	ps = pout->symbols;
	*ps++ = (b>>2) & 3;
	*ps++ = (b   ) & 3;
	for ( int i=0; i<22; ++i ) {
	  b = *pin++;
	  *ps++ = (b>>6);
	  *ps++ = (b>>4) & 3;
	  *ps++ = (b>>2) & 3;
	  *ps++ = (b   ) & 3;
	}
	++pout;
      }
#endif
    }

    // Fill slots with interleaved symbols.
    // EN 302 307-1 figures 7 and 8
#if 0  // For reference
    static void interleave(int bps, int rows,
			   const hard_sb *pin, int nslots,
			   bool msb_first, plslot<hard_ss> *pout) {
      if ( bps==4 && rows==4050 && msb_first )
	return interleave4050(pin, nslots, pout);
      if ( rows % 8 ) fatal("modcod/framesize combination not supported\n");
      int stride = rows/8;  // Offset to next column, in bytes
      hard_sb accs[bps];    // One accumulator per column
      int nacc = 0;         // Bits in each column accumulator
      for ( ; nslots; --nslots,++pout ) {
	pout->is_pls = false;
	hard_ss *ps = pout->symbols;
	for ( int ns=pout->LENGTH; ns--; ++ps ) {
	  if ( ! nacc ) {
	    const hard_sb *pi = pin;
	    for ( int b=0; b<bps; ++b,pi+=stride ) accs[b] = *pi;
	    ++pin;
	    nacc = 8;
	  }
	  hard_ss symb = 0;
	  if ( msb_first )
	    for ( int b=0; b<bps; ++b ) {
	      symb = (symb<<1) | (accs[b]>>7);
	      accs[b] <<= 1;
	    }
	  else
	    for ( int b=bps; b--; ) {
	      symb = (symb<<1) | (accs[b]>>7);
	      accs[b] <<= 1;
	    }
	  --nacc;
	  *ps = symb;
	}
      }
      if ( nacc ) fail("Bug: s2_interleaver");
    }
#else  // reference
    static void interleave(int bps, int rows,
			   const hard_sb *pin, int nslots,
			   bool msb_first, plslot<hard_ss> *pout) {
      void (*func)(int rows, const hard_sb *pin, int nslots,
		   plslot<hard_ss> *pout) = 0;
      if ( msb_first )
	switch ( bps ) {
	case 2: func = interleave<1,2>; break;
	case 3: func = interleave<1,3>; break;
	case 4: func = interleave<1,4>; break;
	case 5: func = interleave<1,5>; break;
	default: fail("Bad bps");
	}
      else
	switch ( bps ) {
	case 2: func = interleave<0,2>; break;
	case 3: func = interleave<0,3>; break;
	case 4: func = interleave<0,4>; break;
	case 5: func = interleave<0,5>; break;
	default: fail("Bad bps");
	}
      (*func)(rows, pin, nslots, pout);
    }
    template<int MSB_FIRST, int BPS>
    static void interleave(int rows, const hard_sb *pin, int nslots,
			   plslot<hard_ss> *pout) {
      if ( BPS==4 && rows==4050 && MSB_FIRST )
	return interleave4050(pin, nslots, pout);
      if ( rows % 8 ) fatal("modcod/framesize combination not supported\n");
      int stride = rows/8;  // Offset to next column, in bytes
      if ( nslots % 4 ) fatal("Bug: Truncated byte");
      // plslot::symbols[] are not packed across slots,
      // so we need tos split bytes at boundaries.
      for ( ; nslots; nslots-=4 ) {
	const hard_sb *pi;
	hard_sb accs[BPS];    // One accumulator per column
	hard_ss *ps;
	// Slot 0 (mod 4): 88+2
	pout->is_pls = false;
	ps = pout->symbols;
	for ( int i=0; i<11; ++i ) {
	  split_byte<BPS>(pin++, stride, accs);
	  pop_symbols<MSB_FIRST,BPS>(accs, &ps, 8);
	}
	split_byte<BPS>(pin++, stride, accs);
	pop_symbols<MSB_FIRST,BPS>(accs, &ps, 2);
	++pout;
	// Slot 1 (mod 4): 6+80+4
	pout->is_pls = false;
	ps = pout->symbols;
	pop_symbols<MSB_FIRST,BPS>(accs, &ps, 6);
	for ( int i=0; i<10; ++i ) {
	  split_byte<BPS>(pin++, stride, accs);
	  pop_symbols<MSB_FIRST,BPS>(accs, &ps, 8);
	}
	split_byte<BPS>(pin++, stride, accs);
	pop_symbols<MSB_FIRST,BPS>(accs, &ps, 4);
	++pout;
	// Slot 2 (mod 4): 4+80+6
	pout->is_pls = false;
	ps = pout->symbols;
	pop_symbols<MSB_FIRST,BPS>(accs, &ps, 4);
	for ( int i=0; i<10; ++i ) {
	  split_byte<BPS>(pin++, stride, accs);
	  pop_symbols<MSB_FIRST,BPS>(accs, &ps, 8);
	}
	split_byte<BPS>(pin++, stride, accs);
	pop_symbols<MSB_FIRST,BPS>(accs, &ps, 6);
	++pout;
	// Slot 3 (mod 4): 2+88
	pout->is_pls = false;
	ps = pout->symbols;
	pop_symbols<MSB_FIRST,BPS>(accs, &ps, 2);
	for ( int i=0; i<11; ++i ) {
	  split_byte<BPS>(pin++, stride, accs);
	  pop_symbols<MSB_FIRST,BPS>(accs, &ps, 8);
	}
	++pout;
      }
    }
    template<int BPS>
    static inline void split_byte(const hard_sb *pi, int stride,
				  hard_sb accs[BPS]) {
      // TBD Pass stride as template parameter.
      for ( int b=0; b<BPS; ++b,pi+=stride ) accs[b] = *pi;
    }
    template<int MSB_FIRST, int BPS>
    static void pop_symbols(hard_sb accs[BPS], hard_ss **ps, int ns) {
      for ( int i=0; i<ns; ++i ) {
	hard_ss symb = 0;
	// Check unrolling and constant propagation.
	for ( int b=0; b<BPS; ++b )
	  if ( MSB_FIRST ) symb = (symb<<1) | (accs[b]      >>7);
	  else             symb = (symb<<1) | (accs[BPS-1-b]>>7);
	for ( int b=0; b<BPS; ++b ) accs[b] <<= 1;
	*(*ps)++ = symb;
      }
    }
#endif  // reference

    // Special case for 16APSK short frames.
    // 4050 rows is not a multiple of 8.
    static void interleave4050(const hard_sb *pin, int nslots,
			       plslot<hard_ss> *pout) {
      const hard_sb *pin0 = pin;
      int rows = 4050;
      hard_sb accs[4];    // One accumulator per column
      int nacc = 0;         // Bits in each column accumulator
      for ( ; nslots; --nslots,++pout ) {
	pout->is_pls = false;
	hard_ss *ps = pout->symbols;
	for ( int ns=pout->LENGTH; ns--; ++ps ) {
	  if ( ! nacc ) {
	    if ( nslots==1 && ns==1 ) {
	      // Special case just to avoid reading beyond end of buffer
	      accs[0] =  pin[   0];
	      accs[1] = (pin[ 506]<<2) | (pin[ 507]>>6);
	      accs[2] = (pin[1012]<<4) | (pin[1013]>>4);
	      accs[3] = (pin[1518]<<6);
	    } else {
	      accs[0] =  pin[   0];
	      accs[1] = (pin[ 506]<<2) | (pin[ 507]>>6);
	      accs[2] = (pin[1012]<<4) | (pin[1013]>>4);
	      accs[3] = (pin[1518]<<6) | (pin[1519]>>2);
	    }
	    ++pin;
	    nacc = 8;
	  }
	  hard_ss symb = 0;
	  for ( int b=0; b<4; ++b ) {
	    symb = (symb<<1) | (accs[b]>>7);
	    accs[b] <<= 1;
	  }
	  --nacc;
	  *ps = symb;
	}
      }
    }

    pipereader< fecframe<hard_sb> > in;
    pipewriter< plslot<hard_ss> > out;
  };  // s2_interleaver


  // S2 DEINTERLEAVER
  // EN 302 307-1 section 5.3.3 Bit Interleaver

  template<typename SOFTSYMB, typename SOFTBYTE>
  struct s2_deinterleaver : runnable {
    s2_deinterleaver(scheduler *sch,
		     pipebuf< plslot<SOFTSYMB> > &_in,
		     pipebuf< fecframe<SOFTBYTE> > &_out)
      : runnable(sch, "S2 deinterleaver"),
	in(_in), out(_out)
    {
    }
    void run() {
      while ( in.readable()>=1 && out.writable()>=1 ) {
	plslot<SOFTSYMB> *pin = in.rd();
	if ( ! pin->is_pls ) fail("s2_deinterleaver: bad input sequence");
	s2_pls *pls = &pin->pls;
	const modcod_info *mcinfo = check_modcod(pls->modcod);
	int nslots = pls->sf ? mcinfo->nslots_nf/4 : mcinfo->nslots_nf;
	if ( in.readable() < 1+nslots ) return;
	fecframe<SOFTBYTE> *pout = out.wr();
	pout->pls = *pls;
	SOFTBYTE *pbytes = pout->bytes;
	if ( mcinfo->nsymbols == 4 )
	  deserialize_qpsk(pin+1, nslots, pbytes);
	else {
	  int bps = log2(mcinfo->nsymbols);
	  int rows = pls->framebits() / bps;
	  if ( mcinfo->nsymbols==8 && mcinfo->rate==FEC35 )
	    deinterleave(bps, rows, pin+1, nslots, false, pbytes);
	  else
	    deinterleave(bps, rows, pin+1, nslots, true, pbytes);
	}
	in.read(1+nslots);
	out.written(1);
      }
    }

  private:

    // Deserialize slots of QPSK symbols, MSB first.
    static void deserialize_qpsk(plslot<SOFTSYMB> *pin, int nslots,
				 SOFTBYTE *pout) {
      SOFTBYTE acc;
      softword_clear(&acc);  // gcc warning
      int nacc = 0;
      for ( ; nslots; --nslots,++pin ) {
	SOFTSYMB *ps = pin->symbols;
	for ( int ns=pin->LENGTH; ns--; ++ps ) {
	  pack_qpsk_symbol(*ps, &acc, nacc);
	  nacc += 2;
	  if ( nacc == 8 ) {  // TBD unroll
	    *pout++ = acc;
	    nacc = 0;
	  }
	}
      }
    }

    // Deinterleave slots of symbols.
    // EN 302 307-1 figures 7 and 8
#if 0  // For reference
    static void deinterleave(int bps, int rows,
			     const plslot<SOFTSYMB> *pin, int nslots,
			     bool msb_first, SOFTBYTE *pout) {
      if ( bps==4 && rows==4050 && msb_first )
	return deinterleave4050(pin, nslots, pout);
      if ( rows % 8 ) fatal("modcod/framesize combination not supported\n");
      int stride = rows/8;  // Offset to next column, in bytes
      SOFTBYTE accs[bps];
      for ( int b=0; b<bps; ++b ) softword_clear(&accs[b]);  // gcc warning
      int nacc = 0;
      for ( ; nslots; --nslots,++pin ) {
	const SOFTSYMB *ps = pin->symbols;
	for ( int ns=pin->LENGTH; ns--; ++ps ) {
	  split_symbol(*ps, bps, accs, nacc, msb_first);
	  ++nacc;
	  if ( nacc == 8 ) {
	    SOFTBYTE *po = pout;
	    for ( int b=0; b<bps; ++b,po+=stride ) *po = accs[b];
	    ++pout;
	    nacc = 0;
	  }
	}
      }
      if ( nacc ) fail("Bug: s2_deinterleaver");
    }
#else  // reference
    static void deinterleave(int bps, int rows,
			     const plslot<SOFTSYMB> *pin, int nslots,
			     bool msb_first, SOFTBYTE *pout) {
      void (*func)(int rows, const plslot<SOFTSYMB> *pin, int nslots,
		   SOFTBYTE *pout) = 0;
      if ( msb_first )
	switch ( bps ) {
	case 2: func = deinterleave<1,2>; break;
	case 3: func = deinterleave<1,3>; break;
	case 4: func = deinterleave<1,4>; break;
	case 5: func = deinterleave<1,5>; break;
	default: fail("Bad bps");
	}
      else
	switch ( bps ) {
	case 2: func = deinterleave<0,2>; break;
	case 3: func = deinterleave<0,3>; break;
	case 4: func = deinterleave<0,4>; break;
	case 5: func = deinterleave<0,5>; break;
	default: fail("Bad bps");
	}
      (*func)(rows, pin, nslots, pout);
    }

    template<int MSB_FIRST, int BPS>
    static void deinterleave(int rows, const plslot<SOFTSYMB> *pin, int nslots,
			     SOFTBYTE *pout) {
      if ( BPS==4 && rows==4050 && MSB_FIRST )
	return deinterleave4050(pin, nslots, pout);
      if ( rows % 8 ) fatal("modcod/framesize combination not supported\n");
      int stride = rows/8;  // Offset to next column, in bytes
      SOFTBYTE accs[BPS];
      for ( int b=0; b<BPS; ++b ) softword_clear(&accs[b]);  // gcc warning
      int nacc = 0;
      for ( ; nslots; --nslots,++pin ) {
	const SOFTSYMB *ps = pin->symbols;
	for ( int ns=pin->LENGTH; ns--; ++ps ) {
	  split_symbol(*ps, BPS, accs, nacc, MSB_FIRST);
	  ++nacc;
	  if ( nacc == 8 ) {  // TBD Unroll, same as interleave()
	    SOFTBYTE *po = pout;
	    // TBD Pass stride as template parameter.
	    for ( int b=0; b<BPS; ++b,po+=stride ) *po = accs[b];
	    ++pout;
	    nacc = 0;
	  }
	}
      }
      if ( nacc ) fail("Bug: s2_deinterleaver");
    }
#endif  // reference

    // Special case for 16APSK short frames.
    // 4050 rows is not a multiple of 8
    // so we process rows one at a time rather than in chunks of 8.
    static void deinterleave4050(const plslot<SOFTSYMB> *pin, int nslots,
				 SOFTBYTE *pout) {
      const int rows = 4050;
      SOFTBYTE accs[4];
      for ( int b=0; b<4; ++b ) softword_clear(&accs[b]);  // gcc warning
      int nacc = 0;
      for ( ; nslots; --nslots,++pin ) {
	const SOFTSYMB *ps = pin->symbols;
	for ( int ns=pin->LENGTH; ns--; ++ps ) {
	  split_symbol(*ps, 4, accs, nacc, true);
	  ++nacc;
	  if ( nacc == 8 ) {
	    for ( int b=0; b<8; ++b ) {
	      softwords_set(pout, rows*0+b, softword_get(accs[0],b));
	      softwords_set(pout, rows*1+b, softword_get(accs[1],b));
	      softwords_set(pout, rows*2+b, softword_get(accs[2],b));
	      softwords_set(pout, rows*3+b, softword_get(accs[3],b));
	    }
	    ++pout;
	    nacc = 0;
	  }
	}
      }
      if ( nacc != 2 ) fatal("Bug: Expected 2 leftover rows\n");
      // Pad with random symbol so that we can use accs[].
      for ( int b=nacc; b<8; ++b )
	split_symbol(pin->symbols[0], 4, accs, b, true);
      for ( int b=0; b<nacc; ++b ) {
	softwords_set(pout, rows*0+b, softword_get(accs[0],b));
	softwords_set(pout, rows*1+b, softword_get(accs[1],b));
	softwords_set(pout, rows*2+b, softword_get(accs[2],b));
	softwords_set(pout, rows*3+b, softword_get(accs[3],b));
      }
    }

    // Spread LLR symbol across hard columns.
    // Must call 8 times before using result because we use bit shifts.
    static inline void split_symbol(const llr_ss &ps, int bps,
				    hard_sb accs[/*bps*/], int nacc,
				    bool msb_first) {
      if ( msb_first ) {
	for ( int b=0; b<bps; ++b )
	  accs[b] = (accs[b]<<1) | llr_harden(ps.bits[bps-1-b]);
      } else {
	for ( int b=0; b<bps; ++b )
	  accs[b] = (accs[b]<<1) | llr_harden(ps.bits[b]);
      }
    }
    // Fast variant
    template<int MSB_FIRST, int BPS>
    static inline void split_symbol(const llr_ss &ps,
				    hard_sb accs[/*bps*/], int nacc) {
      if ( MSB_FIRST ) {
	for ( int b=0; b<BPS; ++b )
	  accs[b] = (accs[b]<<1) | llr_harden(ps.bits[BPS-1-b]);
      } else {
	for ( int b=0; b<BPS; ++b )
	  accs[b] = (accs[b]<<1) | llr_harden(ps.bits[b]);
      }
    }

    // Spread LLR symbol across LLR columns.
    static inline void split_symbol(const llr_ss &ps, int bps,
				    llr_sb accs[/*bps*/], int nacc,
				    bool msb_first) {
      if ( msb_first ) {
	for ( int b=0; b<bps; ++b ) accs[b].bits[nacc] = ps.bits[bps-1-b];
      } else {
	for ( int b=0; b<bps; ++b ) accs[b].bits[nacc] = ps.bits[b];
      }
    }
    // Fast variant
    template<int MSB_FIRST, int BPS>
    static inline void split_symbol(const llr_ss &ps,
				    llr_sb accs[/*bps*/], int nacc) {
      if ( MSB_FIRST ) {
	for ( int b=0; b<BPS; ++b ) accs[b].bits[nacc] = ps.bits[BPS-1-b];
      } else {
	for ( int b=0; b<BPS; ++b ) accs[b].bits[nacc] = ps.bits[b];
      }
    }

    // Merge QPSK LLR symbol into hard byte.
    static inline void pack_qpsk_symbol(const llr_ss &ps,
					hard_sb *acc, int nacc) {
      // TBD Must match LLR law, see softsymb_harden.
      uint8_t s = llr_harden(ps.bits[0]) | (llr_harden(ps.bits[1])<<1);
      *acc = (*acc<<2) | s;
    }

    // Merge QPSK LLR symbol into LLR byte.
    static inline void pack_qpsk_symbol(const llr_ss &ps,
					llr_sb *acc, int nacc) {
      acc->bits[nacc]   = ps.bits[1];
      acc->bits[nacc+1] = ps.bits[0];
    }


    pipereader< plslot<SOFTSYMB> > in;
    pipewriter< fecframe<SOFTBYTE> > out;
  };  // s2_deinterleaver


  typedef ldpc_table<uint16_t> s2_ldpc_table;
  typedef ldpc_engine<bool,hard_sb,8,uint16_t> s2_ldpc_engine;

  #include "dvbs2_data.h"

  static const struct fec_info {
    static const int KBCH_MAX = 58192;
    int Kbch;       // BCH message size (bits)
    int kldpc;      // LDPC message size (= BCH codeword size) (bits)
    int t;          // BCH error correction
    const s2_ldpc_table *ldpc;
  } fec_infos[2][FEC_COUNT] = {
    { // Normal frames
      [FEC12]  = { 32208, 32400, 12, &ldpc_nf_fec12  },
      [FEC23]  = { 43040, 43200, 10, &ldpc_nf_fec23  },
      [FEC46]  = {     0                             },
      [FEC34]  = { 48408, 48600, 12, &ldpc_nf_fec34  },
      [FEC56]  = { 53840, 54000, 10, &ldpc_nf_fec56  },
      [FEC78]  = {     0                             },
      [FEC45]  = { 51648, 51840, 12, &ldpc_nf_fec45  },
      [FEC89]  = { 57472, 57600,  8, &ldpc_nf_fec89  },
      [FEC910] = { 58192, 58320,  8, &ldpc_nf_fec910 },
      [FEC14]  = { 16008, 16200, 12, &ldpc_nf_fec14  },
      [FEC13]  = { 21408, 21600, 12, &ldpc_nf_fec13  },
      [FEC25]  = { 25728, 25920, 12, &ldpc_nf_fec25  },
      [FEC35]  = { 38688, 38880, 12, &ldpc_nf_fec35  },
    },
    { // Short frames
      [FEC12]  = {  7032,  7200, 12, &ldpc_sf_fec12 },
      [FEC23]  = { 10632, 10800, 12, &ldpc_sf_fec23 },
      [FEC46]  = {                                  },
      [FEC34]  = { 11712, 11880, 12, &ldpc_sf_fec34 },
      [FEC56]  = { 13152, 13320, 12, &ldpc_sf_fec56 },
      [FEC78]  = {                                  },
      [FEC45]  = { 12432, 12600, 12, &ldpc_sf_fec45 },
      [FEC89]  = { 14232, 14400, 12, &ldpc_sf_fec89 },
      [FEC910] = {                                  },
      [FEC14]  = {  3072,  3240, 12, &ldpc_sf_fec14 },
      [FEC13]  = {  5232,  5400, 12, &ldpc_sf_fec13 },
      [FEC25]  = {  6312,  6480, 12, &ldpc_sf_fec25 },
      [FEC35]  = {  9552,  9720, 12, &ldpc_sf_fec35 },
    },
  };

  struct bbframe {
    s2_pls pls;
    uint8_t bytes[58192/8];  // Kbch/8 max
  };

  // S2_LDPC_ENGINES
  // Initializes LDPC engines for all DVB-S2 FEC settings.

  template<typename SOFTBIT, typename SOFTBYTE>
  struct s2_ldpc_engines {
    typedef ldpc_engine<SOFTBIT,SOFTBYTE,8,uint16_t> s2_ldpc_engine;
    s2_ldpc_engine *ldpcs[2][FEC_COUNT];  // [shortframes][fec]
    s2_ldpc_engines() {
      memset(ldpcs, 0, sizeof(ldpcs));
      for ( int sf=0; sf<=1; ++sf ) {
	for ( int fec=0; fec<FEC_COUNT; ++fec ) {
	  const fec_info *fi = &fec_infos[sf][fec];
	  if ( ! fi->ldpc ) {
	    ldpcs[sf][fec] = NULL;
	  } else {
	    int n = (sf ? 64800/4 : 64800);
	    int k = fi->kldpc;
	    ldpcs[sf][fec] = new s2_ldpc_engine(fi->ldpc, k, n);
	  }
	}
      }
    }
    void print_node_stats() {
      for ( int sf=0; sf<=1; ++sf )
	for ( int fec=0; fec<FEC_COUNT; ++fec ) {
	  s2_ldpc_engine *ldpc = ldpcs[sf][fec];
	  if ( ldpc ) ldpc->print_node_stats();
	}
    }
  };  // s2_ldpc_engines


  // S2_BCH_ENGINES
  // Initializes BCH engines for all DVB-S2 FEC settings.

  struct s2_bch_engines {
    bch_interface *bchs[2][FEC_COUNT];
    // N=t*m
    // The generator of GF(2^m) is always g1.
    // Normal frames with 8, 10 or 12 polynomials.
    typedef bch_engine<uint32_t,192,17,16,uint16_t,0x002d> s2_bch_engine_nf12;
    typedef bch_engine<uint32_t,160,17,16,uint16_t,0x002d> s2_bch_engine_nf10;
    typedef bch_engine<uint32_t,128,17,16,uint16_t,0x002d> s2_bch_engine_nf8;
    // Short frames with 12 polynomials.
    typedef bch_engine<uint32_t,168,17,14,uint16_t,0x002b> s2_bch_engine_sf12;
    s2_bch_engines() {
      bitvect<uint32_t,17> bch_polys[2][12];  // [shortframes][polyindex]
      // EN 302 307-1 5.3.1 Table 6a (polynomials for normal frames)
      bch_polys[0][ 0] = bitvect<uint32_t,17>(0x1002d);  // g1
      bch_polys[0][ 1] = bitvect<uint32_t,17>(0x10173);  // g2
      bch_polys[0][ 2] = bitvect<uint32_t,17>(0x10fbd);  // g3
      bch_polys[0][ 3] = bitvect<uint32_t,17>(0x15a55);  // g4
      bch_polys[0][ 4] = bitvect<uint32_t,17>(0x11f2f);  // g5
      bch_polys[0][ 5] = bitvect<uint32_t,17>(0x1f7b5);  // g6
      bch_polys[0][ 6] = bitvect<uint32_t,17>(0x1af65);  // g7
      bch_polys[0][ 7] = bitvect<uint32_t,17>(0x17367);  // g8
      bch_polys[0][ 8] = bitvect<uint32_t,17>(0x10ea1);  // g9
      bch_polys[0][ 9] = bitvect<uint32_t,17>(0x175a7);  // g10
      bch_polys[0][10] = bitvect<uint32_t,17>(0x13a2d);  // g11
      bch_polys[0][11] = bitvect<uint32_t,17>(0x11ae3);  // g12
      // EN 302 307-1 5.3.1 Table 6b (polynomials for short frames)
      bch_polys[1][ 0] = bitvect<uint32_t,17>(0x402b);  // g1
      bch_polys[1][ 1] = bitvect<uint32_t,17>(0x4941);  // g2
      bch_polys[1][ 2] = bitvect<uint32_t,17>(0x4647);  // g3
      bch_polys[1][ 3] = bitvect<uint32_t,17>(0x5591);  // g4
      bch_polys[1][ 4] = bitvect<uint32_t,17>(0x6b55);  // g5
      bch_polys[1][ 5] = bitvect<uint32_t,17>(0x6389);  // g6
      bch_polys[1][ 6] = bitvect<uint32_t,17>(0x6ce5);  // g7
      bch_polys[1][ 7] = bitvect<uint32_t,17>(0x4f21);  // g8
      bch_polys[1][ 8] = bitvect<uint32_t,17>(0x460f);  // g9
      bch_polys[1][ 9] = bitvect<uint32_t,17>(0x5a49);  // g10
      bch_polys[1][10] = bitvect<uint32_t,17>(0x5811);  // g11
      bch_polys[1][11] = bitvect<uint32_t,17>(0x65ef);  // g12
      // Redundant with fec_infos[], but needs static template argument.
      memset(bchs, 0, sizeof(bchs));
      bchs[0][FEC12]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[0][FEC23]  = new s2_bch_engine_nf10(bch_polys[0],10);
      bchs[0][FEC34]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[0][FEC56]  = new s2_bch_engine_nf10(bch_polys[0],10);
      bchs[0][FEC45]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[0][FEC89]  = new s2_bch_engine_nf8 (bch_polys[0], 8);
      bchs[0][FEC910] = new s2_bch_engine_nf8 (bch_polys[0], 8);
      bchs[0][FEC14]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[0][FEC13]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[0][FEC25]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[0][FEC35]  = new s2_bch_engine_nf12(bch_polys[0],12);
      bchs[1][FEC12]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC23]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC34]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC56]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC45]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC89]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC14]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC13]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC25]  = new s2_bch_engine_sf12(bch_polys[1],12);
      bchs[1][FEC35]  = new s2_bch_engine_sf12(bch_polys[1],12);
    }
  };  // s2_bch_engines


  // S2 BASEBAND DESCRAMBLER AND FEC ENCODER
  // EN 302 307-1 section 5.2.2
  // EN 302 307-1 section 5.3

  struct s2_fecenc : runnable {
    typedef ldpc_engine<bool,hard_sb,8,uint16_t> s2_ldpc_engine;
    s2_fecenc(scheduler *sch,
	      pipebuf<bbframe> &_in, pipebuf< fecframe<hard_sb> > &_out)
      : runnable(sch, "S2 fecenc"),
	in(_in), out(_out)
    {
      if ( sch->debug ) s2ldpc.print_node_stats();
    }
    void run() {
      while ( in.readable()>=1 && out.writable()>=1 ) {
	bbframe *pin = in.rd();
	fecframe<hard_sb> *pout = out.wr();
	run_frame(in.rd(), out.wr());
	in.read(1);
	out.written(1);
      }
    }
  private:
    void run_frame(const bbframe *pin, fecframe<hard_sb> *pout) {
      const modcod_info *mcinfo = check_modcod(pin->pls.modcod);
      const fec_info *fi = &fec_infos[pin->pls.sf][mcinfo->rate];
      pout->pls = pin->pls;
      hard_sb *pbytes = pout->bytes;
      bbscrambling.transform(pin->bytes, fi->Kbch/8, pbytes);
      { // BCH
	size_t msgbytes = fi->Kbch / 8;
	size_t cwbytes = fi->kldpc / 8;
	bch_interface *bch = s2bch.bchs[pin->pls.sf][mcinfo->rate];
	bch->encode(pbytes, msgbytes, pbytes+msgbytes);
      }
      { // LDPC
	size_t msgbits = fi->kldpc;
	size_t cwbits = pin->pls.framebits();
	s2_ldpc_engine *ldpc = s2ldpc.ldpcs[pin->pls.sf][mcinfo->rate];
	ldpc->encode(fi->ldpc, pbytes, msgbits, cwbits, pbytes+msgbits/8);
      }
    }
    pipereader<bbframe> in;
    pipewriter< fecframe<hard_sb> > out;
    s2_bbscrambling bbscrambling;
    s2_bch_engines s2bch;
    s2_ldpc_engines<bool,hard_sb> s2ldpc;
  };  // s2_fecenc


  // S2 FEC DECODER AND BASEBAND DESCRAMBLER
  // EN 302 307-1 section 5.3
  // EN 302 307-1 section 5.2.2

  template<typename SOFTBIT, typename SOFTBYTE>
  struct s2_fecdec : runnable {
    int bitflips;
    s2_fecdec(scheduler *sch,
	      pipebuf< fecframe<SOFTBYTE> > &_in, pipebuf<bbframe> &_out,
	      pipebuf<int> *_bitcount=NULL,
	      pipebuf<int> *_errcount=NULL)
      : runnable(sch, "S2 fecdec"),
	bitflips(0),
	in(_in), out(_out),
	bitcount(opt_writer(_bitcount,1)),
	errcount(opt_writer(_errcount,1))
    {
      if ( sch->debug ) s2ldpc.print_node_stats();
    }
    void run() {
      while ( in.readable()>=1 && out.writable()>=1 &&
	      opt_writable(bitcount,1) && opt_writable(errcount,1) ) {
	fecframe<SOFTBYTE> *pin = in.rd();
	const modcod_info *mcinfo = check_modcod(pin->pls.modcod);
	const fec_info *fi = &fec_infos[pin->pls.sf][mcinfo->rate];
	bool corrupted = false;
	bool residual_errors;
	if ( true ) {
	  // LDPC decode
	  size_t cwbits = pin->pls.framebits();
	  size_t msgbits = fi->kldpc;
	  size_t chkbits = cwbits - msgbits;
	  s2_ldpc_engine *ldpc = s2ldpc.ldpcs[pin->pls.sf][mcinfo->rate];
	  int ncorr = ldpc->decode_bitflip
	    (fi->ldpc, pin->bytes, msgbits, cwbits, bitflips);
	  if ( sch->debug2 ) fprintf(stderr, "LDPCCORR = %d\n", ncorr);
	}
	uint8_t *hardbytes = softbytes_harden(pin->bytes, fi->kldpc/8, bch_buf);
	if ( true ) {
	  // BCH decode
	  size_t cwbytes = fi->kldpc / 8;
	  size_t msgbytes = fi->Kbch / 8;
	  size_t chkbytes = cwbytes - msgbytes;
	  // Decode with suitable BCH decoder for this MODCOD
	  bch_interface *bch = s2bch.bchs[pin->pls.sf][mcinfo->rate];
	  int ncorr = bch->decode(hardbytes, cwbytes);
	  if ( sch->debug2 ) fprintf(stderr, "BCHCORR = %d\n", ncorr);
	  corrupted = (ncorr < 0);
	  residual_errors = (ncorr != 0);
	  // Report VER
	  opt_write(bitcount, fi->Kbch);
	  opt_write(errcount, (ncorr>=0) ? ncorr : fi->Kbch);
	}
	int bbsize = fi->Kbch / 8;
	// TBD Some decoders want the bad packets.
#if 0
	if ( corrupted ) {
	  fprintf(stderr, "Passing bad frame\n");
	  corrupted = false;
	}
#endif
	if ( ! corrupted ) {
	  // Descramble and output
	  bbframe *pout = out.wr();
	  pout->pls = pin->pls;
	  bbscrambling.transform(hardbytes, bbsize, pout->bytes);
	  out.written(1);
	}
	if ( sch->debug )
	  fprintf(stderr, "%c", corrupted ? ':' : residual_errors ? '.' : '_');
	in.read(1);
      }
    }
  private:
    s2_ldpc_engines<SOFTBIT,SOFTBYTE> s2ldpc;
    uint8_t bch_buf[64800/8];  // Temp storage for hardening before BCH
    s2_bch_engines s2bch;
    s2_bbscrambling bbscrambling;
    pipereader< fecframe<SOFTBYTE> > in;
    pipewriter<bbframe> out;
    pipewriter<int> *bitcount, *errcount;
  };  // s2_fecdec

  // External LDPC decoder
  // Spawns a user-specified command, FEC frames on stdin/stdout.

  template<typename T, int _SIZE>
  struct simplequeue {
    static const int SIZE = _SIZE;
    simplequeue() { rd = wr = count = 0; }
    bool full() { return count == SIZE; }
    T *put() { T *res=&q[wr]; wr=(wr+1)%SIZE; ++count; return res; }
    bool empty() { return count == 0; }
    const T *peek() { return &q[rd]; }
    const T *get() { const T *res=&q[rd]; rd=(rd+1)%SIZE; --count; return res; }
    //  private:
    int rd, wr, count;
    T q[SIZE];
  };

  template<typename SOFTBIT, typename SOFTBYTE>
  struct s2_fecdec_helper : runnable {
    int batch_size;
    int nhelpers;
    bool must_buffer;
    s2_fecdec_helper(scheduler *sch,
		     pipebuf< fecframe<SOFTBYTE> > &_in,
		     pipebuf<bbframe> &_out,
		     const char *_command,
		     pipebuf<int> *_bitcount=NULL,
		     pipebuf<int> *_errcount=NULL)
      : runnable(sch, "S2 fecdec io"),
	batch_size(32),
	nhelpers(1),
	must_buffer(false),
	in(_in), out(_out),
	command(_command),
	bitcount(opt_writer(_bitcount,1)),
	errcount(opt_writer(_errcount,1))
    {
      for ( int mc=0; mc<32; ++mc )
	for ( int sf=0; sf<2; ++sf )
	  pools[mc][sf].procs = NULL;
    }
    void run() {
      bool work_done = false;
      // Send work until all helpers block.
      bool all_blocked = false;
      while ( in.readable() >= 1 && !jobs.full() ) {
	if ( ! send_frame(in.rd()) ) { all_blocked=true; break; }
	in.read(1);
	work_done = true;
      }
      // Risk blocking on read() only when we have nothing else to do
      // and we know a result is coming.
      while ( (all_blocked || !work_done || jobs.full()) &&
	      !jobs.empty() &&
	      jobs.peek()->h->b_out &&
	      out.writable()>=1 &&
	      opt_writable(bitcount,1) && opt_writable(errcount,1) ) {
	receive_frame(jobs.get());
      }
    }
  private:
    struct helper_instance {
      int fd_tx;       // To helper
      int fd_rx;       // From helper
      int batch_size;  // Latency
      int b_in;        // Jobs in input queue
      int b_out;       // Jobs in output queue
    };
    struct pool {
      helper_instance *procs;  // NULL or [nprocs]
      int nprocs;
    } pools[32][2];  // [modcod][sf]
    struct helper_job {
      s2_pls pls;
      helper_instance *h;
    };
    simplequeue<helper_job,1024> jobs;
    // Try to send a frame.  Return false all helpers were busy.
    bool send_frame(fecframe<SOFTBYTE> *pin) {
      pool *p = get_pool(&pin->pls);
      for ( int i=0; i<p->nprocs; ++i ) {
	helper_instance *h = &p->procs[i];
	size_t iosize = (pin->pls.framebits()/8) * sizeof(SOFTBYTE);
	// fprintf(stderr, "Writing %lu to fd %d\n", iosize, h->fd_tx);
	int nw = write(h->fd_tx, pin->bytes, iosize);
	if ( nw<0 && errno==EWOULDBLOCK ) continue;
	if ( nw < 0 ) fatal("write(LDPC helper");
	if ( nw != iosize ) fatal("partial write(LDPC helper)");
	helper_job *job = jobs.put();
	job->pls = pin->pls;
	job->h = h;
	++h->b_in;
	if ( h->b_in >= h->batch_size ) {
	  h->b_in -= h->batch_size;
	  h->b_out += h->batch_size;
	}
	return true;
      }
      return false;
    }
    // Return a pool of running helpers for a given modcod.
    pool *get_pool(const s2_pls *pls) {
      pool *p = &pools[pls->modcod][pls->sf];
      if ( ! p->procs ) {
	p->procs = new helper_instance[nhelpers];
	for ( int i=0; i<nhelpers; ++i ) spawn_helper(&p->procs[i], pls);
	p->nprocs = nhelpers;
      }
      return p;
    }
    // Spawn a helper process.
    void spawn_helper(helper_instance *h, const s2_pls *pls) {
      if ( sch->debug )
	fprintf(stderr, "Spawning LDPC helper: modcod=%d sf=%d\n",
		pls->modcod, pls->sf);
      int tx[2], rx[2];
      if ( pipe(tx) || pipe(rx) ) fatal("pipe");
      // Size the pipes so that the helper never runs out of work to do.
      int pipesize = 64800 * batch_size;
      if ( fcntl(tx[0], F_SETPIPE_SZ, pipesize) < 0 ||
	   fcntl(rx[0], F_SETPIPE_SZ, pipesize) < 0 ||
	   fcntl(tx[1], F_SETPIPE_SZ, pipesize) < 0 ||
	   fcntl(rx[1], F_SETPIPE_SZ, pipesize) < 0 ) {
	fprintf(stderr,
		"*** Failed to increase pipe size.\n"
		"*** Try echo %d > /proc/sys/fs/pipe-max-size\n", pipesize);
	if ( must_buffer )
	  fatal("F_SETPIPE_SZ");
	else
	  fprintf(stderr, "*** Throughput will be suboptimal.\n");
      }
      int child = vfork();
      if ( ! child ) {
	// Child process
	close(tx[1]); dup2(tx[0], 0);
	close(rx[0]); dup2(rx[1], 1);
	char mc_arg[16];
	sprintf(mc_arg, "%d", pls->modcod);
	const char *sf_arg = pls->sf ? "--shortframes" : NULL;
	const char *argv[] = { command, "--modcod", mc_arg, sf_arg, NULL };
	execve(command, (char *const*)argv, NULL);
	fatal(command);
      }
      h->fd_tx = tx[1];  close(tx[0]);
      h->fd_rx = rx[0];  close(rx[1]);
#if 1
      h->batch_size = 32;  // TBD
#else
      fprintf(stderr, "DEBUG: Assuming LDPC helper with no latency\n");
      h->batch_size = 1;  // TBD
#endif
      h->b_in = h->b_out = 0;
      int flags = fcntl(h->fd_tx, F_GETFL);
      if ( fcntl(h->fd_tx, F_SETFL, flags|O_NONBLOCK) ) fatal("fcntl(helper)");
    }

    // Receive a finished job.
    void receive_frame(const helper_job *job) {
      // Read corrected frame from helper
      const s2_pls *pls = &job->pls;
      size_t iosize = (pls->framebits()/8) * sizeof(ldpc_buf[0]);
      int nr = read(job->h->fd_rx, ldpc_buf, iosize);
      if ( nr < 0 ) fatal("read(LDPC helper)");
      if ( nr != iosize ) fatal("partial read(LDPC helper)");
      --job->h->b_out;
      // Decode BCH.
      const modcod_info *mcinfo = check_modcod(job->pls.modcod);
      const fec_info *fi = &fec_infos[job->pls.sf][mcinfo->rate];
      uint8_t *hardbytes = softbytes_harden(ldpc_buf, fi->kldpc/8, bch_buf);
      size_t cwbytes = fi->kldpc / 8;
      size_t msgbytes = fi->Kbch / 8;
      size_t chkbytes = cwbytes - msgbytes;
      bch_interface *bch = s2bch.bchs[job->pls.sf][mcinfo->rate];
      int ncorr = bch->decode(hardbytes, cwbytes);
      if ( sch->debug2 ) fprintf(stderr, "BCHCORR = %d\n", ncorr);
      bool corrupted = (ncorr < 0);
      // Report VBER
      opt_write(bitcount, fi->Kbch);
      opt_write(errcount, (ncorr>=0) ? ncorr : fi->Kbch);
#if 0
      // TBD Some decoders want the bad packets.
	if ( corrupted ) {
	  fprintf(stderr, "Passing bad frame\n");
	  corrupted = false;
	}
#endif
      if ( ! corrupted ) {
	// Descramble and output
	bbframe *pout = out.wr();
	pout->pls = job->pls;
	bbscrambling.transform(hardbytes, fi->Kbch/8, pout->bytes);
	out.written(1);
      }
      if ( sch->debug )
	fprintf(stderr, "%c", corrupted ? '!' : ncorr ? '.' : '_');
    }
    pipereader< fecframe<SOFTBYTE> > in;
    pipewriter<bbframe> out;
    const char *command;
    SOFTBYTE ldpc_buf[64800/8];
    uint8_t bch_buf[64800/8];  // Temp storage for hardening before BCH
    s2_bch_engines s2bch;
    s2_bbscrambling bbscrambling;
    pipewriter<int> *bitcount, *errcount;
  };  // s2_fecdec_helper

  // S2 FRAMER
  // EN 302 307-1 section 5.1 Mode adaptation

  struct s2_framer : runnable {
    uint8_t rolloff_code;  // 0=0.35, 1=0.25, 2=0.20, 3=reserved
    // User must provide pls_seq[n_pls_seq].
    // For ACM, user can change pls_seq[0] at runtime.
    // For VCM with a repeating pattern, use n_pls_seq>=2.
    s2_pls *pls_seq;
    int n_pls_seq;

    s2_framer(scheduler *sch, pipebuf<tspacket> &_in, pipebuf<bbframe> &_out)
      : runnable(sch, "S2 framer"),
	n_pls_seq(0),
	pls_index(0),
	in(_in), out(_out)
    {
      nremain = 0;
      remcrc = 0;  // CRC for nonexistent previous packet
    }
    void run() {
      while ( out.writable() >= 1 ) {
	if ( ! n_pls_seq ) fail("PLS not specified");
	s2_pls *pls = &pls_seq[pls_index];
	const modcod_info *mcinfo = check_modcod(pls->modcod);
	const fec_info *fi = &fec_infos[pls->sf][mcinfo->rate];
	int framebytes = fi->Kbch / 8;
	if ( ! framebytes ) fail("MODCOD/framesize combination not allowed");
	if ( 10+nremain+188*in.readable() < framebytes )
	  break;  // Not enough data to fill a frame
	bbframe *pout = out.wr();
	pout->pls = *pls;
	uint8_t *buf = pout->bytes;
	uint8_t *end = buf + framebytes;
	// EN 302 307-1 section 5.1.6 Base-Band Header insertion
	uint8_t *bbheader = buf;
	uint8_t matype1 = 0;
	matype1 |= 0xc0;  // TS
	matype1 |= 0x20;  // SIS
	matype1 |= (n_pls_seq==1) ? 0x10 : 0x00;  // CCM/ACM
	// TBD ISSY/NPD required for ACM ?
	matype1 |= rolloff_code;
	*buf++ = matype1;              // MATYPE-1
	*buf++ = 0;                    // MATYPE-2
	uint16_t upl = 188 * 8;
	*buf++ = upl >> 8;             // UPL MSB
	*buf++ = upl;                  // UPL LSB
	uint16_t dfl = (framebytes-10) * 8;
	*buf++ = dfl >> 8;             // DFL MSB
	*buf++ = dfl;                  // DFL LSB
	*buf++ = 0x47;                 // SYNC
	uint16_t syncd = nremain * 8;
	*buf++ = syncd >> 8;           // SYNCD MSB
	*buf++ = syncd;                // SYNCD LSB
	*buf++ = crc8.compute(bbheader, 9);
	// Data field
	memcpy(buf, rembuf, nremain);  // Leftover from previous runs
	buf += nremain;
	while ( buf < end ) {
	  tspacket *tsp = in.rd();
	  if ( tsp->data[0] != MPEG_SYNC ) fail("Invalid TS");
	  *buf++ = remcrc;  // Replace SYNC with CRC of previous.
	  remcrc = crc8.compute(tsp->data+1, tspacket::SIZE-1);
	  int nused = end - buf;
	  if ( nused > tspacket::SIZE-1 ) nused = tspacket::SIZE-1;
	  memcpy(buf, tsp->data+1, nused);
	  buf += nused;
	  if ( buf == end ) {
	    nremain = (tspacket::SIZE-1) - nused;
	    memcpy(rembuf, tsp->data+1+nused, nremain);
	  }
	  in.read(1);
	}
	if ( buf != end ) fail("Bug: s2_framer");
	out.written(1);
	++pls_index;
	if ( pls_index == n_pls_seq ) pls_index = 0;
      }
    }
  private:
    int pls_index;  // Next slot to use in pls_seq
    pipereader<tspacket> in;
    pipewriter<bbframe> out;
    crc8_engine crc8;
    int nremain;
    uint8_t rembuf[tspacket::SIZE];
    uint8_t remcrc;
  };  // s2_framer


  // S2 DEFRAMER
  // EN 302 307-1 section 5.1 Mode adaptation

  struct s2_deframer : runnable {
    int fd_bbf;  // FD for raw BB frames, or -1
    bool bbf_l3sync;  // Output 0xb8 headers to fd_bbf
    int fd_gse;  // FD for generic streams, or -1
    s2_deframer(scheduler *sch, pipebuf<bbframe> &_in, pipebuf<tspacket> &_out,
		pipebuf<int> *_state_out=NULL,
		pipebuf<unsigned long> *_locktime_out=NULL)
      : runnable(sch, "S2 deframer"),
      fd_bbf(-1), bbf_l3sync(false), fd_gse(-1),
	nleftover(-1),
	in(_in), out(_out,MAX_TS_PER_BBFRAME),
	current_state(false),
	state_out(opt_writer(_state_out,2)),
	report_state(true),
	locktime(0),
	locktime_out(opt_writer(_locktime_out,MAX_TS_PER_BBFRAME))
    { }
    void run() {
      while ( in.readable()>=1 && out.writable()>=MAX_TS_PER_BBFRAME &&
	      opt_writable(state_out,2) &&
	      opt_writable(locktime_out,MAX_TS_PER_BBFRAME) ) {
	if ( report_state ) {
	  // Report unlocked state on first invocation.
	  opt_write(state_out, 0);
	  report_state = false;
	}
	run_bbframe(in.rd());
	in.read(1);
      }
    }
  private:
    void run_bbframe(bbframe *pin) {
      uint8_t *bbh = pin->bytes;
      // EN 302 307 section 5.1.6 Base-Band Header Insertion
      uint8_t streamtype = bbh[0] >> 6;
      bool sis = bbh[0] & 32;
      bool ccm = bbh[0] & 16;
      bool issyi = bbh[0] & 8;
      bool npd = bbh[0] & 4;
      int ro_code = bbh[0] & 3;
      uint8_t isi = bbh[1];  // if !sis
      uint16_t upl = (bbh[2]<<8) | bbh[3];
      uint16_t dfl = (bbh[4]<<8) | bbh[5];
      uint8_t sync = bbh[6];
      uint16_t syncd = (bbh[7]<<8) | bbh[8];
      uint8_t crcexp = crc8.compute(bbh, 9);
      uint8_t crc = bbh[9];
      uint8_t *data = bbh + 10;
      if ( sch->debug2 ) {
	static const char *stnames[] = { "GP", "GC", "??", "TS" };
	static float ro_values[] = { 0.35, 0.25, 0.20, 0 };
	fprintf(stderr, "BBH: crc %02x/%02x(%s) %s %s(ISI=%d) %s%s%s ro=%.2f"
		" upl=%d dfl=%d sync=%02x syncd=%d\n",
		crc, crcexp, (crc==crcexp)?"OK":"KO",
		stnames[streamtype], (sis?"SIS":"MIS"), isi,
		(ccm?"CCM":"ACM"), (issyi?" ISSYI":""), (npd?" NPD":""),
		ro_values[ro_code],
		upl, dfl, sync, syncd);
      }
      if ( crc!=crcexp || dfl>fec_info::KBCH_MAX ) {
	// Note: Maybe accept syncd=65535
	fprintf(stderr, "Bad bbframe\n");
	nleftover = -1;
	info_unlocked();
	return;  // Max one state_out per loop
      }
      // TBD: Supporting byte-oriented payloads only.
      if ( (dfl&7) || (syncd&7) ) {
	fprintf(stderr, "Unsupported bbframe\n");
	nleftover = -1;
	info_unlocked();
	return;  // Max one state_out per loop
      }
      if ( fd_bbf >= 0 ) {
	if ( bbf_l3sync ) {
	  uint8_t l3sync = 0xb8;
	  ssize_t nw = write(fd_bbf, &l3sync, 1);
	  if ( nw != 1 ) fatal("write(l3sync)");
	}
	size_t n = 10 + dfl/8;
	ssize_t nw = write(fd_bbf, bbh, n);
	if ( nw < 0 ) fatal("write(bbf)");
	if ( nw != n ) fail("partial write(bbf)");
      }
      if ( streamtype==3 && upl==188*8 && sync==0x47 && syncd<=dfl)
	handle_ts(data, dfl, syncd, sync);
      else if ( streamtype == 1 ) {
	if ( fd_gse >= 0 ) {
	  ssize_t nw = write(fd_gse, data, dfl/8);
	  if ( nw < 0 ) fatal("write(gse)");
	  if ( nw != dfl/8 ) fail("partial write(gse");
	}
      } else {
	fprintf(stderr, "Unrecognized bbframe\n");
      }
    }
    void handle_ts(uint8_t *data, uint16_t dfl, uint16_t syncd, uint8_t sync) {
      int pos;  // Start of useful data in this bbframe
      if ( nleftover < 0 ) {
	// Not synced. Skip unusable data at beginning of bbframe.
	pos = syncd/8;
	if ( sch->debug ) fprintf(stderr, "Start TS at %d\n", pos);
	nleftover = 0;
      } else {
	// Sanity check
	if ( syncd/8 != 188-nleftover ) {
	  if ( sch->debug ) fprintf(stderr, "Lost a bbframe ?\n");
	  nleftover = -1;
	  info_unlocked();
	  return;  // Max one state_out per loop
	}
	pos = 0;
      }
      while ( pos+(188-nleftover)+1 <= dfl/8 ) {
	// Enough data available for one packet and its CRC.
	tspacket *pout = out.wr();
	memcpy(pout->data, leftover, nleftover);  // NOP most of the time
	memcpy(pout->data+nleftover, data+pos, 188-nleftover);
	pout->data[0] = sync;  // Replace CRC
	uint8_t crc = crc8.compute(pout->data+1, 188-1);
	if ( data[pos+(188-nleftover)] == crc ) {
	  info_good_packet();
	} else {
	  pout->data[1] |= 0x80;  // Set TEI bit
	  if ( sch->debug ) fprintf(stderr, "C");
	}
	out.written(1);
	pos += 188 - nleftover;
	nleftover = 0;
      }
      int remain = dfl/8 - pos;
      if ( nleftover+remain > sizeof(leftover) ) fail("Bug: TS deframer");
      memcpy(leftover+nleftover, data+pos, remain);
      nleftover += remain;
    }

    void info_unlocked() {
      info_is_locked(false);
      locktime = 0;
    }
    void info_good_packet() {
      info_is_locked(true);
      ++locktime;
      opt_write(locktime_out, locktime);
    }
    void info_is_locked(bool newstate) {
      if ( newstate != current_state ) {
	opt_write(state_out, newstate?1:0);
	current_state = newstate;
      }
    }

    crc8_engine crc8;
    int nleftover;  // Bytes in leftover[]:
                    // -1      = not synced.
                    // 0       = no leftover data
                    // 1 - 187 = incomplete packet
                    // 188     =  waiting for CRC
    uint8_t leftover[188];
    static const int MAX_TS_PER_BBFRAME = fec_info::KBCH_MAX/8/188 + 1;
    bool locked;
    pipereader<bbframe> in;
    pipewriter<tspacket> out;
    int current_state;
    pipewriter<int> *state_out;
    bool report_state;
    unsigned long locktime;
    pipewriter<unsigned long> *locktime_out;
  };  // s2_deframer

}  // namespace

#endif  // LEANSDR_DVBS2_H
