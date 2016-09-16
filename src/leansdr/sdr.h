#ifndef LEANSDR_SDR_H
#define LEANSDR_SDR_H

namespace leansdr {

  //////////////////////////////////////////////////////////////////////
  // SDR blocks
  //////////////////////////////////////////////////////////////////////
  
  typedef unsigned char u8;
  typedef unsigned short u16;
  typedef signed char s8;
  typedef float f32;
  
  typedef complex<f32> cf32;
  typedef complex<f32> iqsymbol;
  typedef complex<u8> cu8;
  typedef complex<s8> cs8;
  
  template<typename T>
  struct auto_notch : runnable {
    int decimation;
    float k;
    auto_notch(scheduler *sch, pipebuf< complex<T> > &_in,
	       pipebuf< complex<T> > &_out, int _nslots,
	       T _agc_rms_setpoint)
      : runnable(sch, "auto_notch"),
	decimation(1024), k(0.002),   // k(0.01)
	fft(4096),
	in(_in), out(_out,fft.n),
	nslots(_nslots), slots(new slot[nslots]),
	phase(0), gain(1), agc_rms_setpoint(_agc_rms_setpoint) {
      for ( int s=0; s<nslots; ++s ) {
	slots[s].i = -1;
	slots[s].expj = new complex<float>[fft.n];
      }
    }
    void run() {
      while ( in.readable()>=fft.n && out.writable()>=fft.n ) {
	if ( ! phase ) detect();
	if ( ++phase >= decimation ) phase = 0;
	process();
	in.read(fft.n);
	out.written(fft.n);
      }
    }
    void detect() {
      complex<T> *pin = in.rd();
      complex<float> data[fft.n];
      float m0=0, m2=0;
      for ( int i=0; i<fft.n; ++i ) {
	data[i].re = pin[i].re;
	data[i].im = pin[i].im;
	m2 += (float)pin[i].re*pin[i].re + (float)pin[i].im*pin[i].im;
	if ( gen_abs(pin[i].re) > m0 ) m0 = gen_abs(pin[i].re);
	if ( gen_abs(pin[i].im) > m0 ) m0 = gen_abs(pin[i].im);
      }
      if ( agc_rms_setpoint && m2 ) {
	float rms = gen_sqrt(m2/fft.n);
	if ( sch->debug ) fprintf(stderr, "(pow %f max %f)", rms, m0);
	float new_gain = agc_rms_setpoint / rms;
	gain = gain*0.9 + new_gain*0.1;
      }
      fft.inplace(data, true);
      float amp[fft.n];
      for ( int i=0; i<fft.n; ++i ) amp[i] = hypotf(data[i].re, data[i].im);
      for ( slot *s=slots; s<slots+nslots; ++s ) {
	int iamax = 0;
	for ( int i=0; i<fft.n; ++i )
	  if ( amp[i] > amp[iamax] ) iamax=i;
	if ( iamax != s->i ) {
	  if ( sch->debug ) 
	    fprintf(stderr, "%s: slot %d new peak %d -> %d\n",
		    name, (int)(s-slots), s->i, iamax);
	  s->i = iamax;
	  s->estim.re = 0;
	  s->estim.im = 0;
	  s->estt = 0;
	  for ( int i=0; i<fft.n; ++i ) {
	    float a = 2 * M_PI * s->i * i / fft.n;
	    s->expj[i].re = cosf(a);
	    s->expj[i].im = sinf(a);
	  }
	}
	amp[iamax] = 0;
	if ( iamax-1 >= 0 ) amp[iamax-1] = 0;
	if ( iamax+1 < fft.n ) amp[iamax+1] = 0;
      }
    }
    void process() {
      complex<T> *pin=in.rd(), *pend=pin+fft.n, *pout=out.wr();
      for ( slot *s=slots; s<slots+nslots; ++s ) s->ej = s->expj;
      for ( ; pin<pend; ++pin,++pout ) {
	complex<float> out = *pin;
	// TODO Optimize for nslots==1 ?
	for ( slot *s=slots; s<slots+nslots; ++s->ej,++s ) {
	  complex<float> bb(pin->re*s->ej->re + pin->im*s->ej->im,
			    -pin->re*s->ej->im + pin->im*s->ej->re);
	  s->estim.re = bb.re*k + s->estim.re*(1-k);
	  s->estim.im = bb.im*k + s->estim.im*(1-k);
	  complex<float> sub(s->estim.re*s->ej->re - s->estim.im*s->ej->im,
			     s->estim.re*s->ej->im + s->estim.im*s->ej->re);
	  out.re -= sub.re;
	  out.im -= sub.im;
	}
	pout->re = gain * out.re;
	pout->im = gain * out.im;
      }
    }
    
  private:
    cfft_engine<float> fft;
    pipereader< complex<T> > in;
    pipewriter< complex<T> > out;
    int nslots;
    struct slot {
      int i;
      complex<float> estim;
      complex<float> *expj, *ej;
      int estt;
    } *slots;
    int phase;
    float gain;
    T agc_rms_setpoint;
  };
  
  template<typename T>
  struct ss_estimator : runnable {
    unsigned long window_size;  // Samples per estimation
    unsigned long decimation;  // Output rate
    ss_estimator(scheduler *sch, pipebuf< complex<T> > &_in, pipebuf<T> &_out)
      : runnable(sch, "SS estimator"),
	window_size(1024), decimation(1), 
	in(_in), out(_out),
	phase(0) {
    }
    void run() {
      while ( in.readable()>=window_size && out.writable()>=1 ) {
	if ( ! phase ) {
	  complex<T> *p=in.rd(), *pend=p+window_size;
	  float s = 0;
	  for ( ; p<pend; ++p )
	    s += (float)p->re*p->re + (float)p->im*p->im;
	  *out.wr() = sqrtf(s/window_size);
	  out.written(1);
	}
	in.read(window_size);
	if ( ++phase >= decimation ) phase = 0;
      }
    }
  private:  
    pipereader< complex<T> > in;
    pipewriter<T> out;
    unsigned long phase;
  };
  
  typedef unsigned short u_angle;  //  [0,2PI[ in 65536 steps
  typedef signed short s_angle;  // [-PI,PI[ in 65536 steps

  // GENERIC CONSTELLATION DECODING BY LOOK-UP TABLE.
  // R must be a power of 2.
  // Up to 256 symbols.
  
  struct softsymbol {
    unsigned char symbol; // 000000IQ for QPSK
    unsigned char metric;
  };

  // Target RMS amplitude for AGC
  //const float cstln_amp = 73;  // Best for 32APSK 9/10
  //const float cstln_amp = 90;  // Best for QPSK
  const float cstln_amp = 75;  // Trade-off

  template<int R>
  struct cstln_lut {
    complex<signed char> *symbols;
    int nsymbols;
    enum predef { BPSK, QPSK, PSK8, APSK16, APSK32 };
    cstln_lut(predef type, float gamma1=1, float gamma2=1) {
      switch ( type ) {
      case BPSK:
	nsymbols = 2;
	symbols = new complex<signed char>[nsymbols];
	symbols[0] = polar(1, 2, 0);
	symbols[1] = polar(1, 2, 1);
	make_lut_from_symbols();
	break;
      case QPSK:
	// EN 300 421, section 4.5 Baseband shaping and modulation
	// EN 302 307, section 5.4.1
	nsymbols = 4;
	symbols = new complex<signed char>[nsymbols];
	symbols[0] = polar(1, 4, 0.5);
	symbols[1] = polar(1, 4, 3.5);
	symbols[2] = polar(1, 4, 1.5);
	symbols[3] = polar(1, 4, 2.5);
	make_lut_from_symbols();
	break;
      case PSK8:
	// EN 302 307, section 5.4.2
	nsymbols = 8;
	symbols = new complex<signed char>[nsymbols];
	symbols[0] = polar(1, 8, 0);
	symbols[1] = polar(1, 8, 7);
	symbols[2] = polar(1, 8, 4);
	symbols[3] = polar(1, 8, 5);
	symbols[4] = polar(1, 8, 2);
	symbols[5] = polar(1, 8, 7);
	symbols[6] = polar(1, 8, 3);
	symbols[7] = polar(1, 8, 6);
	make_lut_from_symbols();
	break;
      case APSK16: {
	// EN 302 307, section 5.4.3
	float r1 = sqrtf(4 / (1+3*gamma1*gamma1));
	float r2 = gamma1 * r1;
	nsymbols = 16;
	symbols = new complex<signed char>[nsymbols];
	symbols[0]  = polar(r2, 12,  1.5);
	symbols[1]  = polar(r2, 12, 10.5);
	symbols[2]  = polar(r2, 12,  4.5);
	symbols[3]  = polar(r2, 12,  7.5);
	symbols[4]  = polar(r2, 12,  0.5);
	symbols[5]  = polar(r2, 12, 11.5);
	symbols[6]  = polar(r2, 12,  5.5);
	symbols[7]  = polar(r2, 12,  6.5);
	symbols[8]  = polar(r2, 12,  2.5);
	symbols[9]  = polar(r2, 12,  9.5);
	symbols[10] = polar(r2, 12,  3.5);
	symbols[11] = polar(r2, 12,  8.5);
	symbols[12] = polar(r1, 4,   0.5);
	symbols[13] = polar(r1, 4,   3.5);
	symbols[14] = polar(r1, 4,   1.5);
	symbols[15] = polar(r1, 4,   2.5);
	make_lut_from_symbols();
	break;
      }
      case APSK32: {
	// EN 302 307, section 5.4.3
	float r1 = sqrtf(8 / (1+3*gamma1*gamma1+4*gamma2*gamma2));
	float r2 = gamma1 * r1;
	float r3 = gamma2 * r1;
	nsymbols = 32;
	symbols = new complex<signed char>[nsymbols];
	symbols[0]  = polar(r2, 12,  1.5);
	symbols[1]  = polar(r2, 12,  2.5);
	symbols[2]  = polar(r2, 12, 10.5);
	symbols[3]  = polar(r2, 12,  9.5);
	symbols[4]  = polar(r2, 12,  4.5);
	symbols[5]  = polar(r2, 12,  3.5);
	symbols[6]  = polar(r2, 12,  7.5);
	symbols[7]  = polar(r2, 12,  8.5);
	symbols[8]  = polar(r3, 16,  1  );
	symbols[9]  = polar(r3, 16,  3  );
	symbols[10] = polar(r3, 16, 14  );
	symbols[11] = polar(r3, 16, 12  );
	symbols[12] = polar(r3, 16,  6  );
	symbols[13] = polar(r3, 16,  4  );
	symbols[14] = polar(r3, 16,  9  );
	symbols[15] = polar(r3, 16, 11  );
	symbols[16] = polar(r2, 12,  0.5);
	symbols[17] = polar(r1,  4,  0.5);
	symbols[18] = polar(r2, 12, 11.5);
	symbols[19] = polar(r1,  4,  3.5);
	symbols[20] = polar(r2, 12,  5.5);
	symbols[21] = polar(r1,  4,  1.5);
	symbols[22] = polar(r2, 12,  6.5);
	symbols[23] = polar(r1,  4,  2.5);
	symbols[24] = polar(r3, 16,  0  );
	symbols[25] = polar(r3, 16,  2  );
	symbols[26] = polar(r3, 16, 15  );
	symbols[27] = polar(r3, 16, 13  );
	symbols[28] = polar(r3, 16,  7  );
	symbols[29] = polar(r3, 16,  5  );
	symbols[30] = polar(r3, 16,  8  );
	symbols[31] = polar(r3, 16, 10  );
	make_lut_from_symbols();
	break;
      }
      default:
	fail("Constellation not implemented");
      }
    }
    struct result {
      struct softsymbol ss;
      s_angle phase_error;
    };
    inline result *lookup(int I, int Q) {
      return &lut[(unsigned char)I][(unsigned char)Q];
    }
  private:
    complex<signed char> polar(float r, int n, float i) {
      float a = i * 2*M_PI / n;
      return complex<signed char>(r*cosf(a)*cstln_amp, r*sinf(a)*cstln_amp);
    }
    result lut[R][R];
    void make_lut_from_symbols() {
      for ( int I=-R/2; I<R/2; ++I )
	for ( int Q=-R/2; Q<R/2; ++Q ) {
	  unsigned int dmin = R*2;
	  unsigned char smin = 0;
	  for ( int s=0; s<nsymbols; ++s ) {
	    unsigned int d = hypotf(I-symbols[s].re, Q-symbols[s].im);
	    if ( d < dmin ) { dmin=d; smin=s; }
	  }
	  float ph_symbol = atan2f(symbols[smin].im,symbols[smin].re);
	  float ph_err = atan2f(Q,I) - ph_symbol;
	  result *pr = &lut[I&(R-1)][Q&(R-1)];
	  if ( dmin > 255 ) fail("dmin overflow");
	  pr->ss.symbol = smin;
	  pr->ss.metric = dmin;
	  //pr->ss.metric = 255 * (int)dmin / dmin2;
	  pr->phase_error = ph_err * 65536 / (2*M_PI);
	}
    }
  };
  
  // CONSTELLATION RECEIVER
  
  template<typename T>
  struct cstln_receiver : runnable {
    cstln_lut<256> *cstln;
    unsigned long meas_decimation;      // Measurement rate
    float omega, min_omega, max_omega;  // Samples per symbol
    u_angle freqw, min_freqw, max_freqw;  // Freq offset in angle per sample
    static const unsigned int chunk_size = 128;
    float kest;
    
    cstln_receiver(scheduler *sch,
		   pipebuf< complex<T> > &_in,
		   pipebuf<softsymbol> &_out,
		   pipebuf<float> *_freq_out=NULL,
		   pipebuf<float> *_ss_out=NULL,
		   pipebuf<float> *_mer_out=NULL,
		   pipebuf<cf32> *_cstln_out=NULL)
      : runnable(sch, "Constellation receiver"),
	cstln(NULL),
	meas_decimation(1048576),
	kest(0.01),
	in(_in), out(_out, chunk_size),
	est_insp(cstln_amp*cstln_amp), agc_gain(1),
	mu(0), phase(0),
	est_sp(0), est_ep(0),
	meas_count(0)  {
      set_omega(1);
      set_freq(0);
      freq_out = _freq_out ? new pipewriter<float>(*_freq_out) : NULL;
      ss_out = _ss_out ? new pipewriter<float>(*_ss_out) : NULL;
      mer_out = _mer_out ? new pipewriter<float>(*_mer_out) : NULL;
      cstln_out = _cstln_out ? new pipewriter<cf32>(*_cstln_out) : NULL;
      memset(hist, 0, sizeof(hist));
      init_trig_tables();
    }
    
    void set_omega(float _omega, float tol=10e-6) {
      omega = _omega;
      min_omega = omega * (1-tol);
      max_omega = omega * (1+tol);
      update_freq_limits();
    }
    
    void set_freq(float freq) {
      freqw = freq * 65536;
      update_freq_limits();
    }
    void update_freq_limits() {
      // Keep PLL away from +-symbolrate/4
      // Note: Don't cast from negative float to unsigned on ARM.
      min_freqw = (s_angle)(freqw - 65536/max_omega/8);
      max_freqw = (s_angle)(freqw + 65536/max_omega/8);
    }
    
    void run() {
      if ( ! cstln ) fail("constellation not set");
      
      // Magic constants that work with the qa recordings.
      const signed long freq_alpha = 0.04 * 65536;
      const signed long freq_beta = 0.001 * 65536;
      float gain_mu = 0.02 / (cstln_amp*cstln_amp) * 2;
      
      int max_meas = chunk_size/meas_decimation + 1;
      // Largin margin on output_size because mu adjustments
      // can lead to more than chunk_size/min_omega symbols.
      while ( in.readable() >= chunk_size+1 &&
	      out.writable() >= chunk_size &&
	      ( !freq_out  || freq_out ->writable()>=max_meas ) &&
	      ( !ss_out    || ss_out   ->writable()>=max_meas ) &&
	      ( !mer_out   || mer_out  ->writable()>=max_meas ) &&
	      ( !cstln_out || cstln_out->writable()>=max_meas ) ) {
	
	complex<T> *pin=in.rd(), *pin0=pin, *pend=pin+chunk_size;
	softsymbol *pout=out.wr(), *pout0=pout;
	
	// These are scoped outside the loop for SS and MER estimation.
	complex<float> sg; // Symbol before AGC;
	complex<float> s;  // For MER estimation and constellation viewer
	complex<signed char> *cstln_point = NULL;
	
	while ( pin < pend ) {
	  // Here mu is the time of the next symbol counted from 0 at pin.
	  if ( mu < 1 ) {
	    // Here 0<=mu<1 is the fractional time of the next symbol
	    // between pin and pin+1.
	    
	    // Derotate pin[0] and pin[1]
	    float cosph, sinph;
	    cosph = fastcos(-phase);
	    sinph = fastsin(-phase);
	    complex<float> s0(pin[0].re*cosph - pin[0].im*sinph,
			      pin[0].re*sinph + pin[0].im*cosph);
	    cosph = fastcos(-(phase+freqw));
	    sinph = fastsin(-(phase+freqw));
	    complex<float> s1(pin[1].re*cosph - pin[1].im*sinph,
			      pin[1].re*sinph + pin[1].im*cosph);
	    
	    // Interpolate linearly
	    sg = s0*(1-mu) + s1*mu;
	    s = sg * agc_gain;
	    
	    // Constellation look-up
	    cstln_lut<256>::result *cr = cstln->lookup(s.re, s.im);
	    *pout = cr->ss;
	    ++pout;
	    
	    // PLL
#if 0
	    signed short c1 = (cr->phase_error * freq_alpha + 1) >> 16;
	    signed short c2 = (cr->phase_error * freq_alpha) / 65536;
	  //	  if ( c1 != c2 ) fprintf(stderr, "\n### %d %d %d\n", cr->phase_error, c1, c2);
	    phase += (cr->phase_error * freq_alpha) / 65536;
	    freqw += (cr->phase_error * freq_beta) / 65536;
#else
	    phase += (cr->phase_error * freq_alpha + 32768) >> 16;
	    freqw += (cr->phase_error * freq_beta + 32768) >> 16;
#endif  
	    
	    // Modified Mueller and Müller
	    // mu[k]=real(c[k]-c[k-2])*conj(p[k-1])-(p[k]-p[k-2])*conj(c[k-1]))
	    //      =dot(c[k]-c[k-2],p[k-1]) - dot(p[k]-p[k-2],c[k-1])
	    // p = received signals
	    // c = decisions (constellation points)
	    hist[2] = hist[1];
	    hist[1] = hist[0];
	    hist[0].p.re = s.re;
	    hist[0].p.im = s.im;
	    cstln_point = &cstln->symbols[cr->ss.symbol];
	    hist[0].c.re = cstln_point->re;
	    hist[0].c.im = cstln_point->im;
	    float muerr =
	      ( (hist[0].p.re-hist[2].p.re)*hist[1].c.re +
		(hist[0].p.im-hist[2].p.im)*hist[1].c.im ) -
	      ( (hist[0].c.re-hist[2].c.re)*hist[1].p.re +
		(hist[0].c.im-hist[2].c.im)*hist[1].p.im );
	    float mucorr = muerr * gain_mu;
	    const float max_mucorr = 0.1;
	    // TBD Optimize out statically
	    if ( mucorr < -max_mucorr ) mucorr = -max_mucorr;
	    if ( mucorr >  max_mucorr ) mucorr =  max_mucorr;
	    mu += mucorr;
	    mu += omega;  // Next symbol time;
	  } // mu<1
	  
	  // Next sample
	  ++pin;
	  --mu;
	  phase += freqw;
	}  // chunk_size
	
	in.read(pin-pin0);
	out.written(pout-pout0);
	
	// Output the last interpolated PSK symbol, once per chunk_size
	if ( cstln_out ) {
	  *cstln_out->wr() = s;
	  cstln_out->written(1);
	}
	
	// AGC
	// For APSK we must do AGC on the symbols, not the whole signal.
	float insp = sg.re*sg.re + sg.im*sg.im;
	est_insp = insp*kest + est_insp*(1-kest);
	if ( est_insp )
	  agc_gain = cstln_amp / gen_sqrt(est_insp);
	
	// SS and MER
	float sig_power = s.re*s.re+s.im*s.im;
	est_sp = sig_power*kest + est_sp*(1-kest);
	if ( ! cstln_point ) fail("No sample");
	complex<float> errvect(s.re-cstln_point->re, s.im-cstln_point->im);
	float errvect_power = errvect.re*errvect.re + errvect.im*errvect.im;
	est_ep = errvect_power*kest + est_ep*(1-kest);
	
	// This is best done periodically ouside the inner loop,
	// but will cause non-deterministic output.
	
	if ( (signed short)(freqw-min_freqw)<0 ||
	     (signed short)(max_freqw-freqw)<0 )
	  freqw = min_freqw + (unsigned short)(max_freqw-min_freqw) / 2;
      
	// Output measurements
	
	meas_count += pin-pin0;
	while ( meas_count >= meas_decimation ) {
	  meas_count -= meas_decimation;
	  if ( freq_out ) {
	    *freq_out->wr() = (float)(signed short)freqw / 65536;
	    freq_out->written(1);
	  }
	  if ( ss_out ) {
	    *ss_out->wr() = sqrtf(est_insp);
	    ss_out->written(1);
	  }
	  if ( mer_out ) {
	    float mer = est_ep ? 10*logf(est_sp/est_ep)/logf(10) : 0;
	    *mer_out->wr() = mer;
	    mer_out->written(1);
	  }
	  
	}
	
      }  // Work to do
    }
    
  private:
    struct {
      complex<float> p;  // Received symbol
      complex<float> c;  // Matched constellation point
    } hist[3];
    pipereader< complex<T> > in;
    pipewriter<softsymbol> out;
    float est_insp, agc_gain;
    float mu;  // PSK time expressed in clock ticks
    u_angle phase;
    float est_sp;  // Estimated RMS signal power
    float est_ep;  // Estimated RMS error vector power
    unsigned long meas_count;
    pipewriter<float> *freq_out, *ss_out, *mer_out;
    pipewriter<cf32> *cstln_out;

    float lut_cos[65536];
    float fastcos(u_angle a) { return lut_cos[a]; }
    float fastsin(u_angle a) { return lut_cos[(u_angle)(a-16384)]; }

    void init_trig_tables() {
      for ( int a=0; a<65536; ++a )
	lut_cos[a] = cosf(a*2*M_PI/65536);
    }
  };
  
  // FREQUENCY SHIFTER
  // Resolution is sample_freq/65536
  
  template<typename T>
  struct rotator : runnable {
    rotator(scheduler *sch, pipebuf< complex<T> > &_in,
	    pipebuf< complex<T> > &_out, float freq)
      : runnable(sch, "rotator"),
	in(_in), out(_out), index(0) {
      int ifreq = freq * 65536;
      if ( sch->debug )
	fprintf(stderr, "Rotate: req=%f real=%f\n", freq, ifreq/65536.0);
      for ( int i=0; i<65536; ++i ) {
	lut_cos[i] = cosf(2*M_PI * i * ifreq / 65536);
	lut_sin[i] = sinf(2*M_PI * i * ifreq / 65536);
      }
    }
    void run() {
      unsigned long count = min(in.readable(), out.writable());
      complex<T> *pin = in.rd(), *pend = pin+count;
      complex<T> *pout = out.wr();
      for ( ; pin<pend; ++pin,++pout,++index ) {
	float c = lut_cos[index];
	float s = lut_sin[index];
	pout->re = pin->re*c - pin->im*s;
	pout->im = pin->re*s + pin->im*c;
      }
      in.read(count);
      out.written(count);
    }
  private:
    pipereader< complex<T> > in;
    pipewriter< complex<T> > out;
    float lut_cos[65536];
    float lut_sin[65536];
    unsigned short index;  // Current phase
  };
  
}  // namespace

#endif  // LEANSDR_SDR_H
