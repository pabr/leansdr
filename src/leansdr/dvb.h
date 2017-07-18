#ifndef LEANSDR_DVB_H
#define LEANSDR_DVB_H

#include <stdint.h>

#include "leansdr/viterbi.h"
#include "leansdr/convolutional.h"
#include "leansdr/sdr.h"
#include "leansdr/rs.h"

namespace leansdr {

  static const int SIZE_RSPACKET = 204;
  static const int MPEG_SYNC = 0x47;
  static const int MPEG_SYNC_INV = (MPEG_SYNC^0xff);
  static const int MPEG_SYNC_CORRUPTED = 0x55;

  // Generic deconvolution

  enum code_rate {
    FEC12, FEC23, FEC34, FEC56, FEC78,  // DVB-S
    FEC45, FEC89, FEC910,               // DVB-S2
  };

  // Customize APSK radii according to code rate

  cstln_lut<256> *make_dvbs2_constellation(cstln_lut<256>::predef c,
					   code_rate r) {
    float gamma1=1, gamma2=1;
    switch ( c ) {
    case cstln_lut<256>::APSK16:
      // EN 302 307, section 5.4.3, Table 9
      switch ( r ) {
      case FEC23:  gamma1 = 3.15; break;
      case FEC34:  gamma1 = 2.85; break;
      case FEC45:  gamma1 = 2.75; break;
      case FEC56:  gamma1 = 2.70; break;
      case FEC89:  gamma1 = 2.60; break;
      case FEC910: gamma1 = 2.57; break;
      default: fail("Code rate not supported with APSK16");
      }
      break;
    case cstln_lut<256>::APSK32:
      // EN 302 307, section 5.4.4, Table 10
      switch ( r ) {
      case FEC34:  gamma1 = 2.84; gamma2 = 5.27; break;
      case FEC45:  gamma1 = 2.72; gamma2 = 4.87; break;
      case FEC56:  gamma1 = 2.64; gamma2 = 4.64; break;
      case FEC89:  gamma1 = 2.54; gamma2 = 4.33; break;
      case FEC910: gamma1 = 2.53; gamma2 = 4.30; break;
      default: fail("Code rate not supported with APSK32");
      }
      break;
    default:
      break;
    }
    return new cstln_lut<256>(c, gamma1, gamma2);
  }

  // EN 300 421, section 4.4.3, table 2 Punctured code, G1=0171, G2=0133
  static const int DVBS_G1 = 0171;
  static const int DVBS_G2 = 0133;

//  G1 = 0b1111001 
//  G2 = 0b1011011
//
//  G1 = [ 1 1 1 1 0 0 1 ]
//  G2 = [ 1 0 1 1 0 1 1 ]
//
//  C = [ G2     ;
//        G1     ;
//        0 G2   ;
//        0 G1   ;
//        0 0 G2 ;
//        0 0 G1 ]
//
//  C = [ 1 0 1 1 0 1 1 0 0 0 0 0 0 ;
//        1 1 1 1 0 0 1 0 0 0 0 0 0 ;
//        0 1 0 1 1 0 1 1 0 0 0 0 0 ;
//        0 1 1 1 1 0 0 1 0 0 0 0 0 ;
//        0 0 1 0 1 1 0 1 1 0 0 0 0 ;
//        0 0 1 1 1 1 0 0 1 0 0 0 0 ;
//        0 0 0 1 0 1 1 0 1 1 0 0 0 ;
//        0 0 0 1 1 1 1 0 0 1 0 0 0 ;
//        0 0 0 0 1 0 1 1 0 1 1 0 0 ;
//        0 0 0 0 1 1 1 1 0 0 1 0 0 ;
//        0 0 0 0 0 1 0 1 1 0 1 1 0 ;
//        0 0 0 0 0 1 1 1 1 0 0 1 0 ;
//        0 0 0 0 0 0 1 0 1 1 0 1 1 ;
//        0 0 0 0 0 0 1 1 1 1 0 0 1 ]
//
//  IQ = [ Q1; I1; ... Q10; I10 ] = C * S
//
//  D * C == [ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ]
// 
//  D = [ 0 1 0 1 1 1 0 1 1 1 0 0 0 0]
//  D = 0x3ba 

  template<typename Tbyte, Tbyte BYTE_ERASED>
  struct deconvol_sync : runnable {
    deconvol_sync(scheduler *sch,
		  pipebuf<softsymbol> &_in,
		  pipebuf<Tbyte> &_out,
		  uint32_t gX, uint32_t gY,
		  uint32_t pX, uint32_t pY)
      : runnable(sch, "deconvol_sync"),
	fastlock(false),
	in(_in), out(_out,SIZE_RSPACKET),
	skip(0) {
      conv = new uint32_t[2];
      conv[0] = gX;
      conv[1] = gY;
      nG = 2;
      punct = new uint32_t[2];
      punct[0] = pX;
      punct[1] = pY;
      punctperiod = 0;
      punctweight = 0;
      for ( int i=0; i<2; ++i ) {
	int nbits = log2(punct[i]) + 1;
	if ( nbits > punctperiod ) punctperiod = nbits;
	punctweight += hamming_weight(punct[i]);
      }
      if ( sch->verbose ) 
	fprintf(stderr, "puncturing %d/%d\n", punctperiod, punctweight);
      deconv = new iq_t[punctperiod];
      deconv2 = new iq_t[punctperiod];
      inverse_convolution();
      init_syncs();
      locked = &syncs[0];
    }

    typedef uint64_t signal_t;
    typedef uint64_t iq_t;

    static int log2(uint64_t x) {
      int n = -1;
      for ( ; x; ++n,x>>=1 ) ;
      return n;
    }

    iq_t convolve(signal_t s) {
      int sbits = log2(s) + 1;
      iq_t iq = 0;
      unsigned char state = 0;
      for ( int b=sbits-1; b>=0; --b ) {  // Feed into convolver, MSB first
	unsigned char bit = (s>>b) & 1;
	state = (state>>1) | (bit<<6);  // Shift register
	for ( int j=0; j<nG; ++j ) {
	  unsigned char xy = parity(state&conv[j]);  // Taps
	  if ( punct[j] & (1<<(b%punctperiod)) )
	    iq = (iq<<1) | xy;
	}
      }
      return iq;
    }
    
    void run() {
      run_decoding();
    }
    
    void next_sync() {
      if ( fastlock ) fatal("Bug: next_sync() called with fastlock");
      ++locked;
      if ( locked == &syncs[NSYNCS] ) {
	locked = &syncs[0];
	// Try next symbol alignment (for FEC other than 1/2)
	skip = 1;
      }
    }

    bool fastlock;

  private:

    static const int maxsbits = 64;
    iq_t response[maxsbits];

    //static const int traceback = 48;  // For code rate 7/8
    static const int traceback = 64;  // For code rate 7/8 with fastlock

    void solve_rec(iq_t prefix, int nprefix, signal_t exp, iq_t *best) {
      if ( prefix > *best ) return;
      if ( nprefix > sizeof(prefix)*8 ) return;
      int solved = 1;
      for ( int b=0; b<maxsbits; ++b ) {
	if ( parity(prefix&response[b]) != ((exp>>b)&1) ) {
	  // Current candidate does not solve this column.
	  if ( (response[b]>>nprefix) == 0 )
	    // No more bits to trace back.
	    return;
	  solved = 0;
	}
      }
      if ( solved ) { *best = prefix; return; }
      solve_rec(prefix,                    nprefix+1, exp, best);
      solve_rec(prefix|((iq_t)1<<nprefix), nprefix+1, exp, best);
    }

    static const int LATENCY = 0;

    void inverse_convolution() {
      for ( int sbit=0; sbit<maxsbits; ++sbit ) {
	response[sbit] = convolve((iq_t)1<<sbit);
	//fprintf(stderr, "response %d = %x\n", sbit, response[sbit]);
      }
      for ( int b=0; b<punctperiod; ++b ) {
	deconv[b] = -(iq_t)1;
	solve_rec(0, 0, 1<<(LATENCY+b), &deconv[b]);
      }
      
      // Alternate polynomials for fastlock
      for ( int b=0; b<punctperiod; ++b ) {
	uint64_t d=deconv[b], d2=d;
	// 1/2
	if ( d == 0x00000000000003baLL ) d2 = 0x0000000000038ccaLL;
	// 2/3
	if ( d == 0x0000000000000f29LL ) d2 = 0x000000003c569329LL;
	if ( d == 0x000000000003c552LL ) d2 = 0x00000000001dee1cLL;
	if ( d == 0x0000000000007948LL ) d2 = 0x00000001e2b49948LL;
	if ( d == 0x00000000000001deLL ) d2 = 0x00000000001e2a90LL;
	// 3/4
	if ( d == 0x000000000000f247LL ) d2 = 0x000000000fd6383bLL;
	if ( d == 0x00000000000fd9eeLL ) d2 = 0x000000000fd91392LL;
	if ( d == 0x0000000000f248d8LL ) d2 = 0x00000000fd9eef18LL;
	// 5/6
	if ( d == 0x0000000000f5727fLL ) d2 = 0x000003d5c909758fLL;
	if ( d == 0x000000003d5c90aaLL ) d2 = 0x0f5727f0229c90aaLL;
	if ( d == 0x000000003daa371cLL ) d2 = 0x000003d5f45630ecLL;
	if ( d == 0x0000000f5727ff48LL ) d2 = 0x0000f57d28260348LL;
	if ( d == 0x0000000f57d28260LL ) d2 = 0xf5727ff48128260LL;
	// 7/8
	if ( d == 0x0000fbeac76c454fLL ) d2 = 0x00fb11d6ba045a8fLL;
	if ( d == 0x00000000fb11d6baLL ) d2 = 0xfbea3c7d930e16baLL;
	if ( d == 0x0000fb112d5038dcLL ) d2 = 0x00fb112d5038271cLL;
	if ( d == 0x000000fbea3c7d68LL ) d2 = 0x00fbeac7975462a8LL;
	if ( d == 0x00000000fb112d50LL ) d2 = 0x00fbea3c86793290LL;
	if ( d == 0x0000fb112dabd2e0LL ) d2 = 0x00fb112d50c3cd20LL;
	if ( d == 0x00000000fb11d640LL ) d2 = 0x00fbea3c8679c980LL;
	if ( d2 == d ) fail("Alt polynomial not provided");
	deconv2[b] = d2;
      }

      if ( sch->debug ) {
	for ( int b=0; b<punctperiod; ++b )
	  fprintf(stderr, "deconv[%d]=0x%016llx %d taps / %d bits\n",
		  b, (unsigned long long)deconv[b],
		  hamming_weight(deconv[b]), log2(deconv[b])+1);
      }

      // Sanity check
      for ( int b=0; b<punctperiod; ++b ) {
	for ( int i=0; i<maxsbits; ++i ) {
	  iq_t iq = convolve((iq_t)1<<(LATENCY+i));
	  int expect = (b==i) ? 1 : 0;
	  int d = parity(iq&deconv[b]);
	  if ( d != expect )
	    fail("Failed to inverse convolutional coding");
	  int d2 = parity(iq&deconv2[b]);
	  if ( d2 != expect )
	    fail("Failed to inverse convolutional coding (alt)");
	}
	if ( traceback > sizeof(iq_t)*8 )
	  fail("Bug: traceback exceeds register size");
	if ( log2(deconv[b])+1 > traceback )
	  fail("traceback insufficient for deconvolution");
	if ( log2(deconv2[b])+1 > traceback )
	  fail("traceback insufficient for deconvolution (alt)");
      }
    }

    static const int NSYNCS = 4;

    struct sync_t {
      u8 lut[2][2];  // lut[(re>0)?1:0][(im>0)?1:0] = 0b000000IQ
      iq_t in;
      int n_in;
      signal_t out;
      int n_out;
      // Auxiliary shift register for fastlock
      iq_t in2;
      int n_in2, n_out2;
    } syncs[NSYNCS];

    void init_syncs() {
      // EN 300 421, section 4.5, Figure 5 QPSK constellation
      // Four rotations * two conjugations.
      // 180° rotation is detected as polarity inversion in mpeg_sync.
      for ( int sync_id=0; sync_id<NSYNCS; ++sync_id ) {
	for ( int re_pos=0; re_pos<=1; ++re_pos )
	  for ( int im_pos=0; im_pos<=1; ++im_pos ) {
	    int re_neg = !re_pos, im_neg = !im_pos;
	    int I, Q;
	    switch ( sync_id ) {
	    case 0:  // Direct 0°
	      I = re_pos ? 0 : 1;
	      Q = im_pos ? 0 : 1;
	      break;
	    case 1:  // Direct 90°
	      I = im_pos ? 0 : 1;
	      Q = re_neg ? 0 : 1;
	      break;
	    case 2:  // Conj 0°
	      I = re_pos ? 0 : 1;
	      Q = im_pos ? 1 : 0;
	      break;
	    case 3:  // Conj 90°
	      I = im_pos ? 1 : 0;
	      Q = re_neg ? 0 : 1;
	      break;
#if 0
	    case 4:  // Direct 180°
	      I = re_neg ? 0 : 1;
	      Q = im_neg ? 0 : 1;
	      break;
	    case 5:  // Direct 270°
	      I = im_neg ? 0 : 1;
	      Q = re_pos ? 0 : 1;
	      break;
	    case 6:  // Conj 180°
	      I = re_neg ? 0 : 1;
	      Q = im_neg ? 1 : 0;
	      break;
	    case 7:  // Conj 270°
	      I = im_neg ? 1 : 0;
	      Q = re_pos ? 0 : 1;
	      break;
#endif
	    }
	    syncs[sync_id].lut[re_pos][im_pos] = (I<<1) | Q;
	  }
	syncs[sync_id].n_in = 0;
	syncs[sync_id].n_out = 0;
	syncs[sync_id].n_in2 = 0;
	syncs[sync_id].n_out2 = 0;
      }
    }

    // TODO: Unroll for each code rate setting.
    // 1/2: 8 symbols -> 1 byte
    // 2/3 12 symbols -> 2 bytes
    // 3/4 16 symbols -> 3 bytes
    // 5/6 24 symbols -> 5 bytes
    // 7/8 32 symbols -> 7 bytes

    inline Tbyte readbyte(sync_t *s, softsymbol *&p) {
      while ( s->n_out < 8 ) {
	iq_t iq = s->in;
	while ( s->n_in < traceback ) {
	  u8 iqbits = s->lut[(p->symbol&2)?1:0][p->symbol&1];
	  ++p;
	  iq = (iq<<2) | iqbits;
	  s->n_in += 2;
	}
	s->in = iq;
	for ( int b=punctperiod-1; b>=0; --b ) {
	  u8 bit = parity(iq&deconv[b]);
	  s->out = (s->out<<1) | bit;
	}
	s->n_out += punctperiod;
	s->n_in -= punctweight;
      }
      Tbyte res = (s->out >> (s->n_out-8)) & 255;
      s->n_out -= 8;
      return res;
    }

    inline unsigned long readerrors(sync_t *s, softsymbol *&p) {
      unsigned long res = 0;
      while ( s->n_out2 < 8 ) {
	iq_t iq = s->in2;
	while ( s->n_in2 < traceback ) {
	  u8 iqbits = s->lut[(p->symbol&2)?1:0][p->symbol&1];
	  ++p;
	  iq = (iq<<2) | iqbits;
	  s->n_in2 += 2;
	}
	s->in2 = iq;
	for ( int b=punctperiod-1; b>=0; --b ) {
	  u8 bit  = parity(iq&deconv[b]);
	  u8 bit2 = parity(iq&deconv2[b]);
	  if ( bit2 != bit ) ++res;
	}
	s->n_out2 += punctperiod;
	s->n_in2 -= punctweight;
      }
      s->n_out2 -= 8;
      return res;
    }

    void run_decoding() {
      in.read(skip);
      skip = 0;

      // 8 byte margin to fill the deconvolver
      if ( in.readable() < 64 ) return;
      int maxrd = (in.readable()-64) / (punctweight/2) * punctperiod / 8;
      int maxwr = out.writable();
      int n = (maxrd<maxwr) ? maxrd : maxwr;
      if ( ! n ) return;
      // Require enough symbols to discriminate in fastlock mode
      // (threshold must be less than size of rspacket)
      if ( n < 32 ) return;
      
      if ( fastlock ) {
	// Try all sync alignments
	unsigned long errors_best = 1 << 30;
	sync_t *best = &syncs[0];
	for ( sync_t *s=syncs; s<syncs+NSYNCS; ++s ) {
	  softsymbol *pin = in.rd();
	  unsigned long errors = 0;
	  for ( int c=n; c--; )
	    errors += readerrors(s, pin);
	  if ( errors < errors_best ) {
	    errors_best = errors;
	    best = s;
	  }
	}
	if ( best != locked ) {
	  // Another alignment produces fewer bit errors
	  if ( sch->debug )
	    fprintf(stderr, "{%d->%d}\n",
		    (int)(locked-syncs), (int)(best-syncs));
	  locked = best;
	}
	// If deconvolution bit error rate > 33%, try next sample alignment
	if ( errors_best > n*8/3 ) {
	  // fprintf(stderr, ">");
	  skip = 1;
	}
      }

      softsymbol *pin=in.rd(), *pin0=pin;
      Tbyte *pout=out.wr(), *pout0=pout;
      while ( n-- )
	*pout++ = readbyte(locked, pin);
      in.read(pin-pin0);
      out.written(pout-pout0);
    }    

    pipereader<softsymbol> in;
    pipewriter<Tbyte> out;
    // DECONVOL
    int nG;
    uint32_t *conv;  // [nG] Convolution polynomials; MSB is newest
    uint32_t *punct;  // [nG] Puncturing pattern
    int punctperiod, punctweight;
    iq_t *deconv;  // [punctperiod] Deconvolution polynomials
    iq_t *deconv2;  // [punctperiod] Alternate polynomials (for fastlock)
    sync_t *locked;
    int skip;

  };

  typedef deconvol_sync<u8,0> deconvol_sync_simple;

  deconvol_sync_simple *make_deconvol_sync_simple(scheduler *sch,
						  pipebuf<softsymbol> &_in,
						  pipebuf<u8> &_out,
						  enum code_rate rate) {
    // EN 300 421, section 4.4.3 Inner coding
    uint32_t pX, pY;
    switch ( rate ) {
    case FEC12:
      pX = 0x1;  // 1
      pY = 0x1;  // 1
      break;
    case FEC23:
      pX = 0xa;  // 1010  (Handle as FEC4/6, no half-symbols)
      pY = 0xf;  // 1111
      break;
    case FEC34:
      pX = 0x5;  // 101
      pY = 0x6;  // 110
      break;
    case FEC56:
      pX = 0x15;  // 10101
      pY = 0x1a;  // 11010
      break;
    case FEC78:
      pX = 0x45;  // 1000101
      pY = 0x7a;  // 1111010
      break;
    default:
      //fail("Code rate not implemented");
      // For testing DVB-S2 constellations.
      fprintf(stderr, "Code rate not implemented; proceeding anyway\n");
      pX = pY = 1;
    }
    return new deconvol_sync_simple(sch, _in, _out, DVBS_G1, DVBS_G2, pX, pY);
  }


  // CONVOLUTIONAL ENCODER

  struct dvb_convol : runnable {
    typedef u8 uncoded_byte;
    typedef u8 hardsymbol;
    dvb_convol(scheduler *sch,
	       pipebuf<uncoded_byte> &_in,
	       pipebuf<hardsymbol> &_out)
      : runnable(sch, "dvb_convol"),
	in(_in), out(_out,8)
    {
    }
    void run() {
      int count = min(in.readable(), out.writable()/8);
      static const u8 remap[] = { 0, 1, 2, 3 };
      convol.run(in.rd(), remap, out.wr(), count);
      in.read(count);
      out.written(count*8);
    }
  private:
    pipereader<uncoded_byte> in;
    pipewriter<hardsymbol> out;
    convol_poly2<uint32_t, 0171, 0133> convol;
  };  // dvb_convol


  // NEW ALGEBRAIC DECONVOLUTION

  // QPSK 1/2 only;
  // With DVB-S polynomials hardcoded.

  template<typename Tin>
  struct dvb_deconvol_sync : runnable {
    typedef u8 decoded_byte;

    int resync_period;

    static const int chunk_size = 64;  // At least 2*sizeof(Thist)/8

    dvb_deconvol_sync(scheduler *sch,
		      pipebuf<Tin> &_in,
		      pipebuf<decoded_byte> &_out)
      : runnable(sch, "deconvol_sync_multipoly"),
	resync_period(32),
	in(_in), out(_out,chunk_size),
	resync_phase(0)
    {
      init_syncs();
      locked = &syncs[0];
    }

    void run() {

      while ( in.readable() >= chunk_size*8 &&
	      out.writable() >= chunk_size ) {
	int errors_best = 1 << 30;
	sync_t *best = NULL;
	for ( sync_t *s=syncs; s<syncs+NSYNCS; ++s ) {
	  if ( resync_phase!=0 && s!=locked )
	    // Decode only the currently-locked alignment
	    continue;
	  Tin *pin = in.rd();
	  static decoded_byte dummy[chunk_size];
	  decoded_byte *pout = (s==locked) ? out.wr() : dummy;
	  int nerrors = s->deconv.run(pin, s->lut, pout, chunk_size);
	  if ( nerrors < errors_best ) { errors_best=nerrors; best=s; }
	}
	in.read(chunk_size*8);
	out.written(chunk_size);
	if ( best != locked ) {
	  if ( sch->debug ) fprintf(stderr, "%%%d", (int)(best-syncs));
	  locked = best;
	}	
	if ( ++resync_phase >= resync_period ) resync_phase = 0;
      } // Work to do

    }  // run()

  private:
    pipereader<Tin> in;
    pipewriter<decoded_byte> out;
    int resync_phase;
    
    static const int NSYNCS = 4;

    struct sync_t {
      deconvol_poly2<Tin, uint32_t, uint64_t, 0x3ba, 0x38f70> deconv;
      u8 lut[4];  // TBD Swap and flip bits in the polynomials instead.
    } syncs[NSYNCS];

    sync_t *locked;

    void init_syncs() {
      for ( int s=0; s<NSYNCS; ++s )
      // EN 300 421, section 4.5, Figure 5 QPSK constellation
      // Four rotations * two conjugations.
      // 180° rotation is detected as polarity inversion in mpeg_sync.
      //         2 | 0
      //         --+--
      //         3 | 1
      // 0°
      syncs[0].lut[0] = 0;
      syncs[0].lut[1] = 1;
      syncs[0].lut[2] = 2;
      syncs[0].lut[3] = 3;
      // 90°
      syncs[1].lut[0] = 2;
      syncs[1].lut[1] = 0;
      syncs[1].lut[2] = 3;
      syncs[1].lut[3] = 1;
      // 0° conjugated
      syncs[2].lut[0] = 1;
      syncs[2].lut[1] = 0;
      syncs[2].lut[2] = 3;
      syncs[2].lut[3] = 2;
      // 90° conjugated
      syncs[3].lut[0] = 0;
      syncs[3].lut[1] = 2;
      syncs[3].lut[2] = 1;
      syncs[3].lut[3] = 3;
    }

  };  // dvb_deconvol_sync


  typedef dvb_deconvol_sync<softsymbol> dvb_deconvol_sync_soft;
  typedef dvb_deconvol_sync<u8> dvb_deconvol_sync_hard;


  // BIT ALIGNMENT AND MPEG SYNC DETECTION

  template<typename Tbyte, Tbyte BYTE_ERASED>
  struct mpeg_sync : runnable {
    int scan_syncs, want_syncs;
    unsigned long lock_timeout;
    bool fastlock;
    int resync_period;

    mpeg_sync(scheduler *sch,
	      pipebuf<Tbyte> &_in,
	      pipebuf<Tbyte> &_out,
	      deconvol_sync<Tbyte,0> *_deconv,
	      pipebuf<int> *_state_out=NULL,
	      pipebuf<unsigned long> *_locktime_out=NULL)
      : runnable(sch, "sync_detect"),
	scan_syncs(8), want_syncs(4),
	lock_timeout(4),
	fastlock(false),
	resync_period(1),
	in(_in), out(_out, SIZE_RSPACKET*(scan_syncs+1)),
	deconv(_deconv),
	polarity(0),
	resync_phase(0),
	bitphase(0), synchronized(false),
	next_sync_count(0),
	report_state(true) {
      state_out = _state_out ? new pipewriter<int>(*_state_out) : NULL;
      locktime_out =
	_locktime_out ? new pipewriter<unsigned long>(*_locktime_out) : NULL;
    }
    
    void run() {
      if ( report_state && state_out && state_out->writable()>=1 ) {
	// Report unlocked state on first invocation.
	state_out->write(0);
	report_state = false;
      }
      if ( synchronized )
	run_decoding();
      else {
	if ( fastlock ) run_searching_fast(); else run_searching();
      }
    }

    void run_searching() {
      bool next_sync = false;
      int chunk = SIZE_RSPACKET * scan_syncs;
      while ( in.readable() >= chunk+1 &&  // Need 1 ahead for bit shifting
	      out.writable() >= chunk &&  // Use as temp buffer
	      ( !state_out || state_out->writable()>=1 ) ) {
	if ( search_sync() ) return;
	in.read(chunk);
	// Switch to next bit alignment
	++bitphase;
	if ( bitphase == 8 ) {
	  bitphase = 0;
	  next_sync = true;
	}
      }

      if ( next_sync ) {
	// No lock this time
	++next_sync_count;
	if ( next_sync_count >= 3 ) {
	  // After a few cycles without a lock, resync the deconvolver.
	  next_sync_count = 0;
	  if ( deconv ) deconv->next_sync();
	}
      }
    }

    void run_searching_fast() {
      int chunk = SIZE_RSPACKET * scan_syncs;
      while ( in.readable() >= chunk+1 &&  // Need 1 ahead for bit shifting
	      out.writable() >= chunk &&  // Use as temp buffer
	      ( !state_out || state_out->writable()>=1 ) ) {
	if ( resync_phase == 0 ) {
	  // Try all bit alighments
	  for ( bitphase=0; bitphase<=7; ++bitphase ) {
	    if ( search_sync() ) return;
	  }
	}
	in.read(SIZE_RSPACKET);
	if ( ++resync_phase >= resync_period ) resync_phase = 0;
      }
    }

    bool search_sync() {
      int chunk = SIZE_RSPACKET * scan_syncs;
      // Bit-shift [scan_sync] packets according to current [bitphase]
      Tbyte *pin = in.rd(), *pend = pin+chunk;
      Tbyte *pout = out.wr();
      unsigned short w = *pin++;
      for ( ; pin<=pend; ++pin,++pout ) {
	w = (w<<8) | *pin;
	*pout = w >> bitphase;
      }
      // Search for [want_sync] start codes at all 204 offsets
      for ( int i=0; i<SIZE_RSPACKET; ++i ) {
	int nsyncs_p=0, nsyncs_n=0;  // # start codes assuming pos/neg polarity
	int phase8_p=-1, phase8_n=-1;  // Position in sequence of 8 packets
	Tbyte *p = &out.wr()[i];
	for ( int j=0; j<scan_syncs; ++j,p+=SIZE_RSPACKET ) {
	  Tbyte b = *p;
	  if ( b==MPEG_SYNC )     { ++nsyncs_p; phase8_n=(8-j)&7; }
	  if ( b==MPEG_SYNC_INV ) { ++nsyncs_n; phase8_p=(8-j)&7; }
	}
	// Detect most likely polarity
	int nsyncs;
	if ( nsyncs_p > nsyncs_n)
	  { polarity=0;  nsyncs=nsyncs_p; phase8=phase8_p; }
	else
	  { polarity=-1; nsyncs=nsyncs_n; phase8=phase8_n; }
	if ( nsyncs>=want_syncs && phase8>=0 ) {
	  if ( sch->debug ) fprintf(stderr, "Locked\n");
	  if ( ! i ) {  // Avoid fixpoint detection in scheduler
	    i = SIZE_RSPACKET;
	    phase8 = (phase8+1) & 7;
	  }
	  in.read(i);  // Skip to first start code
	  synchronized = true;
	  lock_timeleft = lock_timeout;
	  locktime = 0;
	  if ( state_out )
	    state_out->write(1);
	  return true;
	}
      }
      return false;
    }

    void run_decoding() {
      while ( in.readable() >= SIZE_RSPACKET+1 &&  // +1 for bit shifting
	      out.writable() >= SIZE_RSPACKET &&
	      ( !state_out || state_out->writable()>=1 ) &&
	      ( !locktime_out || locktime_out->writable()>=1 ) ) {
	Tbyte *pin = in.rd(), *pend = pin+SIZE_RSPACKET;
	Tbyte *pout = out.wr();
	unsigned short w = *pin++;
	for ( ; pin<=pend; ++pin,++pout ) {
	  w = (w<<8) | *pin;
	  *pout = (w >> bitphase) ^ polarity;
	}
	in.read(SIZE_RSPACKET);
	Tbyte syncbyte = *out.wr();
	out.written(SIZE_RSPACKET);
	++locktime;
	if ( locktime_out )
	  locktime_out->write(locktime);
	// Reset timer if sync byte is correct
	Tbyte expected = phase8 ? MPEG_SYNC : MPEG_SYNC_INV;
	if ( syncbyte == expected ) lock_timeleft = lock_timeout;
	phase8 = (phase8+1) & 7;
	--lock_timeleft;
	if ( ! lock_timeleft ) {
	  if ( sch->debug ) fprintf(stderr, "Unlocked\n");
	  synchronized = false;
	  next_sync_count = 0;
	  if ( state_out )
	    state_out->write(0);
	  return;
	}
      }
    }

  private:
    pipereader<Tbyte> in;
    pipewriter<Tbyte> out;
    deconvol_sync<Tbyte,0> *deconv;
    unsigned char polarity;  // XOR mask, 0 or 0xff
    int resync_phase;
    int bitphase;
    bool synchronized;
    int next_sync_count;
    int phase8;  // Position in 8-packet cycle, -1 if not synchronized
    unsigned long lock_timeleft;
    unsigned long locktime;
    pipewriter<int> *state_out;
    pipewriter<unsigned long> *locktime_out;
    bool report_state;
  };


  template<typename Tbyte>
  struct rspacket { Tbyte data[SIZE_RSPACKET]; };


  // INTERLEAVER

  struct interleaver : runnable {
    interleaver(scheduler *sch, pipebuf< rspacket<u8> > &_in,
		pipebuf< u8 > &_out)
      : runnable(sch, "interleaver"),
	in(_in), out(_out,SIZE_RSPACKET) {
    }
    void run() {
      while ( in.readable() >= 12 &&
	      out.writable() >= SIZE_RSPACKET ) {
	rspacket<u8> *pin=in.rd();
	u8 *pout = out.wr();
	int delay = 0;
	for ( int i=0; i<SIZE_RSPACKET; ++i,++pout,delay=(delay+1)%12 )
	  *pout = pin[11-delay].data[i];
	in.read(1);
	out.written(SIZE_RSPACKET);
      }
    }
  private:
    pipereader< rspacket<u8> > in;
    pipewriter<u8> out;
  };  // interleaver


  // DEINTERLEAVER

  template<typename Tbyte>
  struct deinterleaver : runnable {
    deinterleaver(scheduler *sch, pipebuf<Tbyte> &_in,
		  pipebuf< rspacket<Tbyte> > &_out)
      : runnable(sch, "deinterleaver"),
	in(_in), out(_out) {
    }
    void run() {
      while ( in.readable() >= 17*11*12+SIZE_RSPACKET &&
	      out.writable() >= 1 ) {
	Tbyte *pin = in.rd()+17*11*12, *pend=pin+SIZE_RSPACKET;
	Tbyte *pout= out.wr()->data;
	for ( int delay=17*11; pin<pend;
	      ++pin,++pout,delay=(delay-17+17*12)%(17*12) )
	  *pout = pin[-delay*12];
	in.read(SIZE_RSPACKET);
	out.written(1);
      }
    }
  private:
    pipereader<Tbyte> in;
    pipewriter< rspacket<Tbyte> > out;
  };  // deinterleaver


  static const int SIZE_TSPACKET = 188;
  struct tspacket { u8 data[SIZE_TSPACKET]; };


  // RS ENCODER

  struct rs_encoder : runnable {
    rs_encoder(scheduler *sch,
	       pipebuf<tspacket> &_in, pipebuf< rspacket<u8> > &_out)
      : runnable(sch, "RS encoder"),
	in(_in), out(_out) { }

    void run() {
      while ( in.readable()>=1 && out.writable()>=1 ) {
	u8 *pin = in.rd()->data;
	u8 *pout = out.wr()->data;
	// The first 188 bytes are the uncoded message P(X)
	memcpy(pout, pin, SIZE_TSPACKET);
	// Append 16 RS parity bytes R(X) = - (P(X)*X^16 mod G(X))
	// so that G(X) divides the coded message S(X) = P(X)*X^16 - R(X).
	rs.encode(pout);
	in.read(1);
	out.written(1);
      }
    }
  private:
    rs_engine rs;
    pipereader<tspacket> in;
    pipewriter< rspacket<u8> > out;
  };  // rs_encoder


  // RS DECODER

  template<typename Tbyte, int BYTE_ERASED>
  struct rs_decoder : runnable {
    rs_engine rs;
    rs_decoder(scheduler *sch,
	       pipebuf< rspacket<Tbyte> > &_in,
	       pipebuf<tspacket> &_out,
	       pipebuf<int> *_bitcount=NULL,
	       pipebuf<int> *_errcount=NULL)
      : runnable(sch, "RS decoder"),
	in(_in), out(_out) {
      bitcount = _bitcount ? new pipewriter<int>(*_bitcount) : NULL;
      errcount = _errcount ? new pipewriter<int>(*_errcount) : NULL;
    }
    void run() {
      if ( bitcount && bitcount->writable()<1 ) return;
      if ( errcount && errcount->writable()<1 ) return;

      int nbits=0, nerrs=0;

      while ( in.readable()>=1 && out.writable()>=1 ) {
	Tbyte *pin = in.rd()->data;
	u8 *pout = out.wr()->data;

	nbits += SIZE_RSPACKET * 8;

	// The message is the first 188 bytes.
	if ( sizeof(Tbyte) == 1 )
	  memcpy(pout, pin, SIZE_TSPACKET);      
	else
	  fail("Erasures not implemented");

	u8 synd[16];
	bool corrupted = rs.syndromes(pin, synd);

#if 0
	if ( ! corrupted ) {
	  // Test BM
	  fprintf(stderr, "Simulating errors\n");
	  pin[203] ^= 42;
	  pin[202] ^= 99;
	  corrupted = rs.syndromes(pin, synd);
	}
#endif
	if ( ! corrupted ) {
	  if ( sch->debug ) 
	    fprintf(stderr, "_");  // Packet received without errors.
	} else {
	  corrupted = rs.correct(synd, pout, pin, &nerrs);
	  if ( sch->debug ) {
	    if ( ! corrupted )
	      fprintf(stderr, ".");  // Errors were corrected.
	    else
	      fprintf(stderr, "!");  // Packet still corrupted.
	  }
	}

	in.read(1);
      
	// Output corrupted packets (with a special mark)
	// otherwise the derandomizer will lose synchronization.
	if ( corrupted ) pout[0] ^= MPEG_SYNC_CORRUPTED;
	out.written(1);

      }
      if ( nbits ) {
	if ( bitcount ) bitcount->write(nbits);
	if ( errcount ) errcount->write(nerrs);
      }
    }
  private:
    pipereader< rspacket<Tbyte> > in;
    pipewriter<tspacket> out;
    pipewriter<int> *bitcount, *errcount;
  };  // rs_decoder


  // RANDOMIZER

  struct randomizer : runnable {
    randomizer(scheduler *sch,
	       pipebuf<tspacket> &_in, pipebuf<tspacket> &_out)
      : runnable(sch, "derandomizer"),
	in(_in), out(_out) {
      precompute_pattern();
      pos = pattern;
      pattern_end = pattern + sizeof(pattern)/sizeof(pattern[0]);
    }
    void precompute_pattern() {
      // EN 300 421, section 4.4.1 Transport multiplex adaptation
      pattern[0] = 0xff;  // Invert one in eight sync bytes
      unsigned short st = 000251;  // 0b 000 000 010 101 001 (Fig 2 reversed)
      for ( int i=1; i<188*8; ++i ) {
	u8 out = 0;
	for ( int n=8; n--; ) {
	  int bit = ((st>>13) ^ (st>>14)) & 1;  // Taps
	  out = (out<<1) | bit;  // MSB first
	  st = (st<<1) | bit;  // Feedback
	}
	pattern[i] = (i%188) ? out : 0;  // Inhibit on sync bytes
      }
    }
    void run() {
      while ( in.readable()>=1 && out.writable()>=1 ) {
	u8 *pin = in.rd()->data, *pend = pin+SIZE_TSPACKET;
	u8 *pout= out.wr()->data;
	if ( pin[0] != MPEG_SYNC ) 
	  fprintf(stderr, "randomizer: bad MPEG sync %02x\n", pin[0]);
	for ( ; pin<pend; ++pin,++pout,++pos ) *pout = *pin ^ *pos;
	if ( pos == pattern_end ) pos = pattern;
	in.read(1);
	out.written(1);
      }
    }
  private:
    u8 pattern[188*8], *pattern_end, *pos;
    pipereader<tspacket> in;
    pipewriter<tspacket> out;
  };  // randomizer


  // DERANDOMIZER

  struct derandomizer : runnable {
    derandomizer(scheduler *sch,
		 pipebuf<tspacket> &_in, pipebuf<tspacket> &_out)
      : runnable(sch, "derandomizer"),
	in(_in), out(_out) {
      precompute_pattern();
      pos = pattern;
      pattern_end = pattern + sizeof(pattern)/sizeof(pattern[0]);
    }
    void precompute_pattern() {
      // EN 300 421, section 4.4.1 Transport multiplex adaptation
      pattern[0] = 0xff;  // Restore the inverted sync byte
      unsigned short st = 000251;  // 0b 000 000 010 101 001 (Fig 2 reversed)
      for ( int i=1; i<188*8; ++i ) {
	u8 out = 0;
	for ( int n=8; n--; ) {
	  int bit = ((st>>13) ^ (st>>14)) & 1;  // Taps
	  out = (out<<1) | bit;  // MSB first
	  st = (st<<1) | bit;  // Feedback
	}
	pattern[i] = (i%188) ? out : 0;  // Inhibit on sync bytes
      }
    }
    void run() {
      while ( in.readable()>=1 && out.writable()>=1 ) {
	u8 *pin = in.rd()->data, *pend = pin+SIZE_TSPACKET;
	u8 *pout= out.wr()->data;
	if ( pin[0] == MPEG_SYNC_INV ||
	     pin[0] == (MPEG_SYNC_INV^MPEG_SYNC_CORRUPTED) ) {
	  if ( pos != pattern ) {
	    if ( sch->debug )
	      fprintf(stderr, "derandomizer: resynchronizing\n");
	    pos = pattern;
	  }
	}
	for ( ; pin<pend; ++pin,++pout,++pos ) *pout = *pin ^ *pos;
	if ( pos == pattern_end ) pos = pattern;
	in.read(1);

	u8 sync = out.wr()->data[0];
	if ( sync == MPEG_SYNC ) {
	  out.written(1);
	} else {
	  if ( sync != (MPEG_SYNC^MPEG_SYNC_CORRUPTED) )
	    if ( sch->debug ) fprintf(stderr, "(%02x)", sync);
	  out.wr()->data[1] |= 0x80;  // Set the Transport Error Indicator bit
	  // We could output corrupted packets here, in case the
	  // MPEG decoder can use them somehow.
	  //out.written(1);  
	}
      }
    }
  private:
    u8 pattern[188*8], *pattern_end, *pos;
    pipereader<tspacket> in;
    pipewriter<tspacket> out;
  };  // derandomizer


  // VITERBI DECODING
  // QPSK 1/2 only.

  struct viterbi_sync : runnable {
    static const int TRACEBACK = 32;  // Suitable for QPSK 1/2
    typedef uint8_t TS, TCS, TUS;
    typedef uint16_t TBM;
    typedef uint32_t TPM;
    typedef bitpath<uint32_t,TUS,1,TRACEBACK> dvb_path;
    typedef viterbi_dec<TS,64, TUS,2, TCS,4, TBM, TPM, dvb_path> dvb_dec;
    typedef trellis<TS,64, TUS,2, 4> dvb_trellis;

  private:
    pipereader<softsymbol> in;
    pipewriter<unsigned char> out;
    static const int NSYNCS = 4;
    dvb_dec *syncs[NSYNCS];
    int current_sync;
    static const int chunk_size = 128;
    int resync_phase;
  public:
    int resync_period;

    viterbi_sync(scheduler *sch,
		 pipebuf<softsymbol> &_in,
		 pipebuf<unsigned char> &_out)
      : runnable(sch, "viterbi_sync"),
	in(_in), out(_out, chunk_size),
	current_sync(0),
	resync_phase(0),
	resync_period(32)   // 1/32 = 9% synchronization overhead
    {
      dvb_trellis *trell = new dvb_trellis();
      uint64_t dvb_polynomials[] = { DVBS_G1, DVBS_G2 };
      trell->init_convolutional(dvb_polynomials);
      for ( int s=0; s<NSYNCS; ++s ) syncs[s] = new dvb_dec(trell);
    }

    inline TUS update_sync(int s, TBM m[4], TPM *discr) {
      // EN 300 421, section 4.5 Baseband shaping and modulation
      // EN 302 307, section 5.4.1
      //              
      //    IQ=10=(2) | IQ=00=(0)
      //    ----------+----------
      //    IQ=11=(3) | IQ=01=(1)
      //
      TBM vm[4];
      switch ( s ) {
      case 0:  // Mapping for 0°
	vm[0]=m[0]; vm[1]=m[1]; vm[2]=m[2]; vm[3]=m[3]; break;
      case 1:  // Mapping for 90°
	vm[0]=m[2]; vm[2]=m[3]; vm[3]=m[1]; vm[1]=m[0]; break;
      case 2:  // Mapping for 0° conjugated
	vm[0]=m[1]; vm[1]=m[0]; vm[2]=m[3]; vm[3]=m[2]; break;
      case 3:  // Mapping for 90° conjugated
	vm[0]=m[3]; vm[2]=m[2]; vm[3]=m[0]; vm[1]=m[1]; break;
      default:
	return 0;  // Avoid compiler warning
      }
      return syncs[s]->update(vm, discr);
    }

    void run() {
      while ( in.readable()>=8*chunk_size && out.writable()>=chunk_size ) {
	unsigned long totaldiscr[NSYNCS];
	for ( int s=0; s<NSYNCS; ++s ) totaldiscr[s] = 0;
	for ( int bytenum=0; bytenum<chunk_size; ++bytenum ) {
	  // Decode one byte
	  unsigned char byte = 0;
	  softsymbol *pin = in.rd();
	  if ( ! resync_phase ) {
	    // Every one in [resync_period] chunks, run all decoders
	    // and compute average quality metrics
	    for ( int b=0; b<8; ++b,++pin ) {
	      TUS bits[NSYNCS];
	      for ( int s=0; s<NSYNCS; ++s ) {
		TPM discr;
		bits[s] = update_sync(s, pin->metrics4, &discr);
		if ( bytenum*8 > TRACEBACK ) totaldiscr[s] += discr;
	      }
	      byte = (byte<<1) | bits[current_sync];
	    }
	  } else {
	    // Otherwise run only the selected decoder
	    for ( int b=0; b<8; ++b,++pin ) {
	      TUS bit = update_sync(current_sync, pin->metrics4, NULL);
	      byte = (byte<<1) | bit;
	    }
	  }
	  in.read(8);
	  out.write(byte);
	}  // chunk_size
	if ( ! resync_phase ) {
	  // Switch to another decoder ?
	  int best = current_sync;
	  for ( int s=0; s<NSYNCS; ++s )
	    if ( totaldiscr[s] > totaldiscr[best] ) best = s;
	  if ( best != current_sync ) {
	    if ( sch->debug ) fprintf(stderr, "{%d->%d}", current_sync, best);
	    current_sync = best;
	  }
	}
	if ( ++resync_phase >= resync_period ) resync_phase = 0;
      }
    }
    
  };  // viterbi_sync


  // VITERBI DECODING - BPSK
  // ETSI TR 101 198: The BPSK variant of DVB-S transmits R=I,Q,I,Q,...
  // Code rate 1/2 only.

  struct viterbi_sync_bpsk : runnable {
    static const int TRACEBACK = 32;  // Suitable for BPSK 1/2
    typedef uint8_t TS, TCS, TUS;
    typedef uint16_t TBM;
    typedef uint32_t TPM;
    typedef bitpath<uint32_t,TUS,1,TRACEBACK> dvb_path;
    typedef viterbi_dec<TS,64, TUS,2, TCS,4, TBM, TPM, dvb_path> dvb_dec;
    typedef trellis<TS,64, TUS,2, 4> dvb_trellis;

  private:
    pipereader<softsymbol> in;
    pipewriter<unsigned char> out;
    static const int NSYNCS = 2;  // Alignment of symbol pairs
    dvb_dec *syncs[2];
    int current_sync;
    static const int chunk_size = 128;
    int resync_phase;
  public:
    int resync_period;

    viterbi_sync_bpsk(scheduler *sch,
		      pipebuf<softsymbol> &_in,
		      pipebuf<unsigned char> &_out)
      : runnable(sch, "viterbi_sync_bpsk"),
	in(_in), out(_out, chunk_size),
	current_sync(0),
	resync_phase(0),
	resync_period(32)
    {
      dvb_trellis *trell = new dvb_trellis();
      uint64_t dvb_polynomials[] = { DVBS_G1, DVBS_G2 };
      trell->init_convolutional(dvb_polynomials);
      for ( int s=0; s<NSYNCS; ++s ) syncs[s] = new dvb_dec(trell);
    }

    inline TUS update_sync(int s, TBM mX[4], TBM mY[4], TPM *discr) {
      // Reconstruct QPSK metrics from pairs of BPSK symbols.
      // EN 300 421, section 4.5 Baseband shaping and modulation
      // EN 302 307, section 5.4.1
      //              
      //    IQ=10=(2) | IQ=00=(0)
      //    ----------+----------
      //    IQ=11=(3) | IQ=01=(1)
      //
      TBM vm[4];
      vm[0] = mX[0] + mY[0];
      vm[1] = mX[0] + mY[1];
      vm[2] = mX[1] + mY[0];
      vm[3] = mX[1] + mY[1];
      return syncs[s]->update(vm, discr);
    }

    void run() {
      while ( in.readable()>=2*8*chunk_size+1 &&  // +1 for pair alignment
	      out.writable()>=chunk_size ) {
	unsigned long totaldiscr[NSYNCS];
	for ( int s=0; s<NSYNCS; ++s ) totaldiscr[s] = 0;
	for ( int bytenum=0; bytenum<chunk_size; ++bytenum ) {
	  // Decode one byte
	  unsigned char byte = 0;
	  softsymbol *pin = in.rd();
	  if ( ! resync_phase ) {
	    // Every one in [resync_period] chunks, run all decoders
	    // and compute average quality metrics
	    for ( int b=0; b<8; ++b,pin+=2 ) {
	      TUS bits[NSYNCS];
	      for ( int s=0; s<NSYNCS; ++s ) {
		TPM discr;
		bits[s] = update_sync(s,
				      (pin+s)->metrics4, (pin+s+1)->metrics4,
				      &discr);
		if ( bytenum*8 > TRACEBACK ) totaldiscr[s] += discr;
	      }
	      byte = (byte<<1) | bits[current_sync];
	    }
	  } else {
	    // Otherwise run only the selected decoder
	    for ( int b=0; b<8; ++b,pin+=2 ) {
	      TUS bit = update_sync(current_sync,
				    (pin+current_sync)->metrics4,
				    (pin+current_sync+1)->metrics4, NULL);
	      byte = (byte<<1) | bit;
	    }
	  }
	  in.read(16);
	  out.write(byte);
	}  // chunk_size
 	if ( ! resync_phase ) {
	  // Switch to another decoder ?
	  int best = current_sync;
	  for ( int s=0; s<NSYNCS; ++s )
	    if ( totaldiscr[s] > totaldiscr[best] ) best = s;
	  if ( best != current_sync ) {
	    if ( sch->debug ) fprintf(stderr, "{%d->%d}", current_sync, best);
	    current_sync = best;
	  }
	}
	if ( ++resync_phase >= resync_period ) resync_phase = 0;
      }
    }

  };  // viterbi_sync_bpsk

}  // namespace

#endif  // LEANSDR_DVB_H
