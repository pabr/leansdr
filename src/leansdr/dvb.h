#ifndef LEANSDR_DVB_H
#define LEANSDR_DVB_H

namespace leansdr {

  static const int SIZE_RSPACKET = 204;
  static const int MPEG_SYNC = 0x47;
  static const int MPEG_SYNC_INV = (MPEG_SYNC^0xff);
  static const int MPEG_SYNC_CORRUPTED = 0x55;

  // Generic deconvolution

  enum code_rate { FEC12, FEC23, FEC34, FEC56, FEC78 };

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
		  unsigned long gX, unsigned long gY,
		  unsigned long pX, unsigned long pY)
      : runnable(sch, "deconvol_sync"),
	in(_in), out(_out,SIZE_RSPACKET),
	skip(0) {
      conv = new unsigned long[2];
      conv[0] = gX;
      conv[1] = gY;
      nG = 2;
      punct = new unsigned long[2];
      punct[0] = pX;
      punct[1] = pY;
      punctperiod = 0;
      punctweight = 0;
      for ( int i=0; i<2; ++i ) {
	int nbits = log2(punct[i]) + 1;
	if ( nbits > punctperiod ) punctperiod = nbits;
	punctweight += hamming(punct[i]);
      }
      if ( sch->verbose ) 
	fprintf(stderr, "puncturing %d/%d\n", punctperiod, punctweight);
      deconv = new iq_t[punctperiod];
      inverse_convolution();
      init_syncs();
      locked = &syncs[0];
    }

    unsigned char hamming(unsigned long x) {
      int h = 0;
      for ( ; x; x>>=1 ) h += x&1;
      return h;
    }
    unsigned char parity(unsigned long x) {
      unsigned char parity = 0;
      for ( ; x; x>>=1 ) parity ^= x&1;
      return parity;
    }
    unsigned char parity(unsigned long long x) {
      unsigned char parity = 0;
      for ( ; x; x>>=1 ) parity ^= x&1;
      return parity;
    }

    typedef unsigned long long signal_t;
    typedef unsigned long long iq_t;

    static int log2(unsigned long long x) {
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
      ++locked;
      if ( locked == &syncs[NSYNCS] ) {
	locked = &syncs[0];
	// Try next symbol alignment (for FEC other than 1/2)
	skip = 1;
      }
    }

  private:

    static const int maxsbits = 64;
    iq_t response[maxsbits];

    static const int traceback = 48;  // For code rate 7/8

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

    void inverse_convolution() {
      for ( int sbit=0; sbit<maxsbits; ++sbit ) {
	response[sbit] = convolve((iq_t)1<<sbit);
	//fprintf(stderr, "response %d = %x\n", sbit, response[sbit]);
      }
      for ( int b=0; b<punctperiod; ++b ) {
	deconv[b] = -(iq_t)1;
	solve_rec(0, 0, 1<<b, &deconv[b]);
      }

      // Sanity check
      for ( int b=0; b<punctperiod; ++b ) {
	for ( int i=0; i<maxsbits; ++i ) {
	  iq_t iq = convolve((iq_t)1<<i);
	  unsigned long d = parity(iq&deconv[b]);
	  unsigned long expect = (b==i) ? 1 : 0;
	  if ( d != expect )
	    fail("Failed to inverse convolutional coding");
	}
	if ( log2(deconv[b])+1 > traceback )
	  fail("traceback exceeds limit");
      }
    }

    static const int NSYNCS = 8;

    struct sync_t {
      u8 lut[2][2];  // lut[(re>0)?1:0][(im>0)?1:0] = 0b000000IQ
      iq_t in;
      int n_in;
      signal_t out;
      int n_out;
    } syncs[NSYNCS];

    void init_syncs() {
      // EN 300 421, section 4.5, Figure 5 QPSK constellation
      // Four rotations * two conjugations.
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
	    case 2:  // Direct 180°
	      I = re_neg ? 0 : 1;
	      Q = im_neg ? 0 : 1;
	      break;
	    case 3:  // Direct 270°
	      I = im_neg ? 0 : 1;
	      Q = re_pos ? 0 : 1;
	      break;
	    case 4:  // Conj 0°
	      I = re_pos ? 0 : 1;
	      Q = im_pos ? 1 : 0;
	      break;
	    case 5:  // Conj 90°
	      I = im_pos ? 1 : 0;
	      Q = re_neg ? 0 : 1;
	      break;
	    case 6:  // Conj 180°
	      I = re_neg ? 0 : 1;
	      Q = im_neg ? 1 : 0;
	      break;
	    case 7:  // Conj 270°
	      I = im_neg ? 1 : 0;
	      Q = re_pos ? 0 : 1;
	      break;
	    }
	    syncs[sync_id].lut[re_pos][im_pos] = (I<<1) | Q;
	  }
	syncs[sync_id].n_in = 0;
	syncs[sync_id].n_out = 0;
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
	while ( s->n_in < traceback ) {
	  u8 iq = s->lut[(p->symbol&2)?1:0][p->symbol&1];
	  ++p;
	  s->in = (s->in<<2) | iq;
	  s->n_in += 2;
	}
	iq_t iq = s->in >> (s->n_in-40);
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

    void run_decoding() {
      in.read(skip);
      skip = 0;

      // 8 byte margin to fill the deconvolver
      int maxrd = (in.readable()-64) / (punctweight/2) * punctperiod / 8;
      int maxwr = out.writable();
      int n = (maxrd<maxwr) ? maxrd : maxwr;
      if ( n <= 0 ) return;
      
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
    unsigned long *conv;  // [nG] Convolution polynomials; MSB is newest
    unsigned long *punct;  // [nG] Puncturing pattern
    int punctperiod, punctweight;
    iq_t *deconv;  // [punctperiod] Deconvolution polynomials
    sync_t *locked;
    int skip;

  };

  typedef deconvol_sync<u8,0> deconvol_sync_simple;

  deconvol_sync_simple make_deconvol_sync_simple(scheduler *sch,
					  pipebuf<softsymbol> &_in,
					  pipebuf<u8> &_out,
					  enum code_rate rate) {
    // EN 300 421, section 4.4.3 Inner coding
    unsigned long pX, pY;
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
      fail("Code rate not implemented");
    }
    return deconvol_sync_simple(sch, _in, _out, DVBS_G1, DVBS_G2, pX, pY);
  }

  template<typename Tbyte, Tbyte BYTE_ERASED>
  struct mpeg_sync : runnable {
    int scan_syncs, want_syncs;
    unsigned long lock_timeout;

    mpeg_sync(scheduler *sch,
	      pipebuf<Tbyte> &_in,
	      pipebuf<Tbyte> &_out,
	      deconvol_sync<Tbyte,0> *_deconv,
	      pipebuf<int> *_state_out=NULL)
      : runnable(sch, "sync_detect"),
	scan_syncs(4), want_syncs(2),
	lock_timeout(4),
	in(_in), out(_out, SIZE_RSPACKET*(scan_syncs+1)),
	deconv(_deconv),
	bitphase(0), synchronized(false),
	report_state(true) {
      state_out = _state_out ? new pipewriter<int>(*_state_out) : NULL;
    }
    
    void run() {
      if ( report_state && state_out && state_out->writable()>=1 ) {
	*state_out->wr() = 0;
	state_out->written(1);
	report_state = false;
      }
      if ( ! synchronized ) run_searching(); else run_decoding();
    }

    void run_searching() {
      int chunk = SIZE_RSPACKET * (scan_syncs+1);
      while ( in.readable() >= chunk+1 &&
	      out.writable() >= chunk &&
	      ( !state_out || state_out->writable()>=1 ) ) {
	Tbyte *pin = in.rd(), *pend = pin+chunk;
	Tbyte *pout = out.wr();
	for ( ; pin<pend; ++pin,++pout ) {
	  unsigned short w = ((unsigned short)pin[0]<<8) | pin[1];
	  *pout = w >> bitphase;
	}
	for ( int i=0; i<SIZE_RSPACKET; ++i ) {
	  int nsyncs = 0;
	  Tbyte *p = &out.wr()[i];
	  int phase8 = -1;
	  for ( int j=0; j<scan_syncs; ++j,p+=SIZE_RSPACKET ) {
	    Tbyte b = *p;
	    if ( b==MPEG_SYNC )
	      ++nsyncs;
	    if ( b==MPEG_SYNC_INV ) phase8 = (8-j)&7;
	  }
	  if ( nsyncs>=want_syncs && phase8>=0 ) {
	    if ( sch->debug ) fprintf(stderr, "Locked\n");
	    if ( ! i ) {  // Avoid fixpoint detection
	      i = SIZE_RSPACKET;
	      phase8 = (phase8+1) & 7;
	    }
	    in.read(i);  // Skip until beginning
	    synchronized = true;
	    lock_timeleft = lock_timeout;
	    if ( state_out ) {
	      *state_out->wr() = 1;
	      state_out->written(1);
	    }
	    return;
	  }
	}
	in.read(chunk);
	++bitphase;
	if ( bitphase == 8 ) {
	  bitphase = 0;
	  deconv->next_sync();
	}
      }
    }

    void run_decoding() {
      while ( in.readable() >= SIZE_RSPACKET+1 &&
	      out.writable() >= SIZE_RSPACKET &&
	       ( !state_out || state_out->writable()>=1 ) ) {
	Tbyte *pin = in.rd(), *pend = pin+SIZE_RSPACKET;
	Tbyte *pout = out.wr();
	for ( ; pin<pend; ++pin,++pout ) {
	  unsigned short w = ((unsigned short)pin[0]<<8) | pin[1];
	  *pout = w >> bitphase;
	}
	in.read(SIZE_RSPACKET);
	Tbyte syncbyte = *out.wr();
	out.written(SIZE_RSPACKET);
	// Reset timer if sync byte is correct
	Tbyte expected = phase8 ? MPEG_SYNC : MPEG_SYNC_INV;
	if ( syncbyte == expected ) lock_timeleft = lock_timeout;
	phase8 = (phase8+1) & 7;
	--lock_timeleft;
	if ( ! lock_timeleft ) {
	  if ( sch->debug ) fprintf(stderr, "Unlocked\n");
	  synchronized = false;
	  if ( state_out ) {
	    *state_out->wr() = 0;
	    state_out->written(1);
	  }
	  return;
	}
      }
    }

  private:
    pipereader<Tbyte> in;
    pipewriter<Tbyte> out;
    deconvol_sync<Tbyte,0> *deconv;
    int bitphase;
    bool synchronized;
    int phase8;
    unsigned long lock_timeleft;
    pipewriter<int> *state_out;
    bool report_state;
  };

  // DEINTERLEAVING

  template<typename Tbyte>
  struct rspacket { Tbyte data[SIZE_RSPACKET]; };

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
  };

      static const int SIZE_TSPACKET = 188;
  struct tspacket { u8 data[SIZE_TSPACKET]; };

  // DERANDOMIZATION

  struct derandomizer : runnable {
    derandomizer(scheduler *sch, pipebuf<tspacket> &_in, pipebuf<tspacket> &_out)
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
  };

}  // namespace

#endif  // LEANSDR_DVB_H
