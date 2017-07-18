#ifndef LEANSDR_SDR_H
#define LEANSDR_SDR_H

#include "leansdr/math.h"
#include "leansdr/dsp.h"

namespace leansdr {

  // Abbreviations for floating-point types

  typedef float f32;

  typedef complex<u8> cu8;
  typedef complex<s8> cs8;
  typedef complex<u16> cu16;
  typedef complex<s16> cs16;
  typedef complex<f32> cf32;


  //////////////////////////////////////////////////////////////////////
  // SDR blocks
  //////////////////////////////////////////////////////////////////////
  
  // AUTO-NOTCH FILTER

  // Periodically detects the [nslots] strongest peaks with a FFT,
  // removes them with a first-order filter.

  template<typename T>
  struct auto_notch : runnable {
    int decimation;
    float k;
    auto_notch(scheduler *sch, pipebuf< complex<T> > &_in,
	       pipebuf< complex<T> > &_out, int _nslots,
	       T _agc_rms_setpoint)
      : runnable(sch, "auto_notch"),
	decimation(1024*4096), k(0.002),   // k(0.01)
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
	phase += fft.n;
	if ( phase >= decimation ) {
	  phase -= decimation;
	  detect();
	}
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


  // SIGNAL STRENGTH ESTIMATOR

  // Outputs RMS values.

  template<typename T>
  struct ss_estimator : runnable {
    unsigned long window_size;  // Samples per estimation
    unsigned long decimation;  // Output rate
    ss_estimator(scheduler *sch, pipebuf< complex<T> > &_in, pipebuf<T> &_out)
      : runnable(sch, "SS estimator"),
	window_size(1024), decimation(1024),
	in(_in), out(_out),
	phase(0) {
    }
    void run() {
      while ( in.readable()>=window_size && out.writable()>=1 ) {
	phase += window_size;
	if ( phase >= decimation ) {
	  phase -= decimation;
	  complex<T> *p=in.rd(), *pend=p+window_size;
	  float s = 0;
	  for ( ; p<pend; ++p )
	    s += (float)p->re*p->re + (float)p->im*p->im;
	  out.write(sqrtf(s/window_size));
	}
	in.read(window_size);
      }
    }
  private:  
    pipereader< complex<T> > in;
    pipewriter<T> out;
    unsigned long phase;
  };

  template<typename T>
  struct ss_amp_estimator : runnable {
    unsigned long window_size;  // Samples per estimation
    unsigned long decimation;  // Output rate
    ss_amp_estimator(scheduler *sch, pipebuf< complex<T> > &_in,
		     pipebuf<T> &_out_ss,
		     pipebuf<T> &_out_ampmin, pipebuf<T> &_out_ampmax)
      : runnable(sch, "SS estimator"),
	window_size(1024), decimation(1024),
	in(_in), out_ss(_out_ss),
	out_ampmin(_out_ampmin), out_ampmax(_out_ampmax),
	phase(0) {
    }
    void run() {
      while ( in.readable() >= window_size &&
	      out_ss.writable() >= 1 &&
	      out_ampmin.writable() >= 1 &&
	      out_ampmax.writable() >= 1 ) {
	phase += window_size;
	if ( phase >= decimation ) {
	  phase -= decimation;
	  complex<T> *p=in.rd(), *pend=p+window_size;
	  float s2 = 0;
	  float amin=1e38, amax=0;
	  for ( ; p<pend; ++p ) {
	    float mag2 = (float)p->re*p->re + (float)p->im*p->im;
	    s2 += mag2;
	    float mag = sqrtf(mag2);
	    if ( mag < amin ) amin = mag;
	    if ( mag > amax ) amax = mag;
	  }
	  out_ss.write(sqrtf(s2/window_size));
	  out_ampmin.write(amin);
	  out_ampmax.write(amax);
	}
	in.read(window_size);
      }
    }
  private:  
    pipereader< complex<T> > in;
    pipewriter<T> out_ss, out_ampmin, out_ampmax;
    unsigned long phase;
  };
  
  // AGC

  template<typename T>
  struct simple_agc : runnable {
    float out_rms;    // Desired RMS output power
    float bw;         // Bandwidth
    float estimated;  // Input power
    simple_agc(scheduler *sch,
	       pipebuf< complex<T> > &_in,
	       pipebuf< complex<T> > &_out)
      : runnable(sch, "AGC"),
	out_rms(1), bw(0.001), estimated(0),
	in(_in), out(_out) {
    }
  private:
    pipereader< complex<T> > in;
    pipewriter< complex<T> > out;
    static const int chunk_size = 128;
    void run() {
      while ( in.readable() >= chunk_size &&
	      out.writable() >= chunk_size ) {
	complex<T> *pin=in.rd(), *pend=pin+chunk_size;
	float amp2 = 0;
	for ( ; pin<pend; ++pin ) amp2 += pin->re*pin->re + pin->im*pin->im;
	amp2 /= chunk_size;
	if ( ! estimated ) estimated = amp2;
	estimated = estimated*(1-bw) + amp2*bw;
	float gain = estimated ? out_rms / sqrtf(estimated) : 0;
	pin = in.rd();
	complex<T> *pout = out.wr();
	float bwcomp = 1 - bw;
	for ( ; pin<pend; ++pin,++pout ) {
	  pout->re = pin->re * gain;
	  pout->im = pin->im * gain;
	}
	in.read(chunk_size);
	out.written(chunk_size);
      }
    }
  };  // simple_agc


  typedef unsigned short u_angle;  //  [0,2PI[ in 65536 steps
  typedef signed short s_angle;  // [-PI,PI[ in 65536 steps


  // GENERIC CONSTELLATION DECODING BY LOOK-UP TABLE.

  // Metrics and phase errors are pre-computed on a RxR grid.
  // R must be a power of 2.
  // Up to 256 symbols.
  
  struct softsymbol {
    unsigned short metrics4[4];  // For Viterbi QPSK
    unsigned char symbol; // 000000IQ for QPSK
  };

  // Target RMS amplitude for AGC
  //const float cstln_amp = 73;  // Best for 32APSK 9/10
  //const float cstln_amp = 90;  // Best for QPSK
  //const float cstln_amp = 64;  // Best for BPSK
  //const float cstln_amp = 75;  // Best for BPSK at 45°
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
#if 0  // BPSK at 0°
	symbols[0] = polar(1, 2, 0);
	symbols[1] = polar(1, 2, 1);
#else  // BPSK at 45°
	symbols[0] = polar(1, 8, 1);
	symbols[1] = polar(1, 8, 5);
#endif
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
    inline result *lookup(float I, float Q) {
      // Handling of overflows beyond the lookup table:
      // - For BPSK/QPSK/8PSK we only care about the phase,
      //   so the following is fine.
      // - For amplitude modulations this is not appropriate.
      //   However, if there is enough noise to cause overflow,
      //   demodulation would probably fail anyway.
      //
      // Comment-out for better throughput at high SNR.
#if 1
      while ( I<-128 || I>127 || Q<-128 || Q>127 ) {
	I *= 0.5;
	Q *= 0.5;
      }
#endif
      return &lut[(u8)(s8)I][(u8)(s8)Q];
    }
    inline result *lookup(int I, int Q) {
      // Ignore wrapping modulo 256
      return &lut[(u8)I][(u8)Q];
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
	  result *pr = &lut[I&(R-1)][Q&(R-1)];
	  unsigned int dmin = R*2;
	  unsigned char smin = 0;
	  for ( int s=0; s<nsymbols; ++s ) {
	    unsigned int d2 =
	      (I-symbols[s].re)*(I-symbols[s].re) +
	      (Q-symbols[s].im)*(Q-symbols[s].im);
	    if ( d2 > 65535 ) fail("Unexpected constellation");
	    unsigned int d = sqrtf(d2);
	    if ( d < dmin ) { dmin=d; smin=s; }
	    if ( nsymbols <= 4 ) {
	      pr->ss.metrics4[s] = d2;
	    }
	  }
	  float ph_symbol = atan2f(symbols[smin].im,symbols[smin].re);
	  float ph_err = atan2f(Q,I) - ph_symbol;
	  if ( dmin > 255 ) fail("dmin overflow");
	  pr->ss.symbol = smin;
	  pr->phase_error = (s32)(ph_err * 65536 / (2*M_PI));  // Mod 65536
	}
    }

  public:
    // Convert soft metric to Hamming distance
    void harden() {
      for ( int i=0; i<R; ++i )
	for ( int q=0; q<R; ++q ) {
	  unsigned short *m = lut[i][q].ss.metrics4;
	  int best;
	  if ( nsymbols == 2 ) {
	    best = (m[0]<m[1]) ? 0 : 1;
	    m[0] = (best==0) ? 0 : 1;
	    m[1] = (best==1) ? 0 : 1;
	  } else {
	    if ( m[0]<=m[1] && m[0]<=m[2] && m[0]<=m[3] ) best = 0;
	    else if ( m[1]<=m[2] && m[1]<=m[3] ) best = 1;
	    else if ( m[2]<=m[3] ) best = 2;
	    else best = 3;
	    for ( int s=0; s<4; ++s )
	      m[s] = hamming_weight((uint8_t)(s^best));
	  }
	}
    }

  };  // cstln_lut


  // SAMPLER INTERFACE FOR CSTLN_RECEIVER
  
  template<typename T>
  struct sampler_interface {
    virtual complex<T> interp(const complex<T> *pin, float mu, float phase) = 0;
    virtual void update_freq(float freqw) { }  // 65536 = 1 Hz
    virtual int readahead() { return 0; }
  };


  // NEAREST-SAMPLE SAMPLER FOR CSTLN_RECEIVER
  // Suitable for bandpass-filtered, oversampled signals only

  template<typename T>
  struct nearest_sampler : sampler_interface<T> {
    int readahead() { return 0; }
    complex<T> interp(const complex<T> *pin, float mu, float phase) {
      return pin[0]*trig.expi(-phase);
    }
  private:
    trig16 trig;
  };  // nearest_sampler


  // LINEAR SAMPLER FOR CSTLN_RECEIVER

  template<typename T>
  struct linear_sampler : sampler_interface<T> {
    int readahead() { return 1; }

    complex<T> interp(const complex<T> *pin, float mu, float phase) {
      // Derotate pin[0] and pin[1]
      complex<T> s0 = pin[0]*trig.expi(-phase);
      complex<T> s1 = pin[1]*trig.expi(-(phase+freqw));
      // Interpolate linearly
      return s0*(1-mu) + s1*mu;
    }

    void update_freq(float _freqw) { freqw = _freqw; }

  private:
    trig16 trig;
    float freqw;
  };  // linear_sampler


  // FIR SAMPLER FOR CSTLN_RECEIVER

  template<typename T, typename Tc>
  struct fir_sampler : sampler_interface<T> {
    fir_sampler(int _ncoeffs, Tc *_coeffs, int _subsampling=1)
      : ncoeffs(_ncoeffs), coeffs(_coeffs), subsampling(_subsampling),
	shifted_coeffs(new complex<T>[ncoeffs]),
	update_freq_phase(0)
    {
    }

    int readahead() { return ncoeffs-1; }

    complex<T> interp(const complex<T> *pin, float mu, float phase) {
      // Apply FIR filter with subsampling
      complex<T> acc(0, 0);
      complex<T> *pc = shifted_coeffs + (int)((1-mu)*subsampling);
      complex<T> *pcend = shifted_coeffs + ncoeffs;
      if ( subsampling == 1 ) {
	// Special case for heavily oversampled signals,
	// where filtering is expensive.
	// gcc-4.9.2 can vectorize this form with NEON on ARM.
	while ( pc < pcend )
	  acc += (*pc++)*(*pin++);
      } else {
	// Not vectorized because the coefficients are not
	// guaranteed to be contiguous in memory.
	for ( ; pc<pcend; pc+=subsampling,++pin ) 
	  acc += (*pc)*(*pin);
      }
      // Derotate
      return trig.expi(-phase) * acc;
    }

    void update_freq(float freqw) {
      // Throttling: Update one coeff per 16 processed samples,
      // to keep the overhead of freq tracking below about 10%.
      update_freq_phase -= 128;  // chunk_size of cstln_receiver
      if ( update_freq_phase <= 0  ) {
	update_freq_phase = ncoeffs*16;
	do_update_freq(freqw);
      }
    }

  private:
    void do_update_freq(float freqw) {
      float f = freqw / subsampling;
      for ( int i=0; i<ncoeffs; ++i )
	shifted_coeffs[i] = trig.expi(-f*(i-ncoeffs/2)) * coeffs[i];
    }
    trig16 trig;
    int ncoeffs;
    Tc *coeffs;
    int subsampling;
    cf32 *shifted_coeffs;
    int update_freq_phase;
  };  // fir_sampler


  // CONSTELLATION RECEIVER

  // Linear interpolation: good enough for 1.2 samples/symbol,
  // but higher oversampling is recommended.

  template<typename T>
  struct cstln_receiver : runnable {
    sampler_interface<T> *sampler;
    cstln_lut<256> *cstln;
    unsigned long meas_decimation;      // Measurement rate
    float omega, min_omega, max_omega;  // Samples per symbol
    float freqw, min_freqw, max_freqw;  // Freq offs (65536 = 1 Hz)
    float pll_adjustment;
    bool allow_drift;                   // Follow carrier beyond safe limits
    static const unsigned int chunk_size = 128;
    float kest;
    
    cstln_receiver(scheduler *sch,
		   sampler_interface<T> *_sampler,
		   pipebuf< complex<T> > &_in,
		   pipebuf<softsymbol> &_out,
		   pipebuf<float> *_freq_out=NULL,
		   pipebuf<float> *_ss_out=NULL,
		   pipebuf<float> *_mer_out=NULL,
		   pipebuf<cf32> *_cstln_out=NULL)
      : runnable(sch, "Constellation receiver"),
	sampler(_sampler),
	cstln(NULL),
	meas_decimation(1048576),
	pll_adjustment(1.0),
	allow_drift(false),
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
      refresh_freq_tap();
    }

    void set_allow_drift(bool d) {
      allow_drift = d;
    }

    void update_freq_limits() {
      // Prevent PLL from crossing +-SR/n/2 and locking at +-SR/n.
      int n = 4;
      if ( cstln ) {
	switch ( cstln->nsymbols ) {
	case  2: n =  2; break;  // BPSK
	case  4: n =  4; break;  // QPSK
	case  8: n =  8; break;  // 8PSK
	case 16: n = 12; break;  // 16APSK
	case 32: n = 16; break;  // 32APSK
	default: n =  4; break;
	}
      }
      min_freqw = freqw - 65536/max_omega/n/2;
      max_freqw = freqw + 65536/max_omega/n/2;
    }
    
    void run() {
      if ( ! cstln ) fail("constellation not set");
      
      // Magic constants that work with the qa recordings.
      float freq_alpha = 0.04;
      float freq_beta = 0.0012 / omega * pll_adjustment;
      float gain_mu = 0.02 / (cstln_amp*cstln_amp) * 2;
      
      int max_meas = chunk_size/meas_decimation + 1;
      // Large margin on output_size because mu adjustments
      // can lead to more than chunk_size/min_omega symbols.
      while ( in.readable() >= chunk_size+sampler->readahead() &&
	      out.writable() >= chunk_size &&
	      ( !freq_out  || freq_out ->writable()>=max_meas ) &&
	      ( !ss_out    || ss_out   ->writable()>=max_meas ) &&
	      ( !mer_out   || mer_out  ->writable()>=max_meas ) &&
	      ( !cstln_out || cstln_out->writable()>=max_meas ) ) {
	
	sampler->update_freq(freqw);

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
	    sg = sampler->interp(pin, mu, phase);
	    s = sg * agc_gain;
	    
	    // Constellation look-up
	    cstln_lut<256>::result *cr = cstln->lookup(s.re, s.im);
	    *pout = cr->ss;
	    ++pout;
	    
	    // PLL
	    phase += cr->phase_error * freq_alpha;
	    freqw += cr->phase_error * freq_beta;
	    
	    // Modified Mueller and Müller
	    // mu[k]=real((c[k]-c[k-2])*conj(p[k-1])-(p[k]-p[k-2])*conj(c[k-1]))
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

	// Normalize phase so that it never exceeds 32 bits.
	// Max freqw is 2^31/65536/chunk_size = 256 Hz
	// (this may happen with leandvb --drift --decim).
	phase = fmodf(phase, 65536);

	if ( cstln_point ) {
	  
	  // Output the last interpolated PSK symbol, max once per chunk_size
	  if ( cstln_out )
	    cstln_out->write(s);
	
	  // AGC
	  // For APSK we must do AGC on the symbols, not the whole signal.
	  // TODO Use a better estimator at low SNR.
	  float insp = sg.re*sg.re + sg.im*sg.im;
	  est_insp = insp*kest + est_insp*(1-kest);
	  if ( est_insp )
	    agc_gain = cstln_amp / gen_sqrt(est_insp);
	  
	  // SS and MER
	  complex<float> ev(s.re-cstln_point->re, s.im-cstln_point->im);
	  float sig_power, ev_power;
	  if ( cstln->nsymbols == 2 ) {
	    // Special case for BPSK: Ignore quadrature component of noise.
	    // TBD Projection on I axis assumes BPSK at 45°
	    float sig_real = (cstln_point->re+cstln_point->im) * 0.707;
	    float ev_real = (ev.re+ev.im) * 0.707;
	    sig_power = sig_real * sig_real;
	    ev_power = ev_real * ev_real;
	  } else {
	    sig_power =
	      (int)cstln_point->re*cstln_point->re +
	      (int)cstln_point->im*cstln_point->im;
	    ev_power = ev.re*ev.re + ev.im*ev.im;
	  }
	  est_sp = sig_power*kest + est_sp*(1-kest);
	  est_ep = ev_power*kest + est_ep*(1-kest);

	}

	// This is best done periodically ouside the inner loop,
	// but will cause non-deterministic output.
	
	if ( ! allow_drift ) {
	  if ( freqw < min_freqw || freqw > max_freqw )
	    freqw = (max_freqw+min_freqw) / 2;
	}
      
	// Output measurements
	
	refresh_freq_tap();

	meas_count += pin-pin0;
	while ( meas_count >= meas_decimation ) {
	  meas_count -= meas_decimation;
	  if ( freq_out )
	    freq_out->write(freq_tap);
	  if ( ss_out )
	    ss_out->write(sqrtf(est_insp));
	  if ( mer_out )
	    mer_out->write(est_ep ? 10*logf(est_sp/est_ep)/logf(10) : 0);
	}
	
      }  // Work to do
    }
    
    float freq_tap;
    void refresh_freq_tap() {
      freq_tap = freqw / 65536;
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
    float phase;  // 65536=2pi
    // Signal estimation
    float est_sp;  // Estimated RMS signal power
    float est_ep;  // Estimated RMS error vector power
    unsigned long meas_count;
    pipewriter<float> *freq_out, *ss_out, *mer_out;
    pipewriter<cf32> *cstln_out;
  };
 
  
  // FAST QPSK RECEIVER

  // Optimized for u8 input, no AGC, uses phase information only.
  // Outputs hard symbols.

  template<typename T>
  struct fast_qpsk_receiver : runnable {
    typedef u8 hardsymbol;
    unsigned long meas_decimation;      // Measurement rate
    float omega, min_omega, max_omega;  // Samples per symbol
    signed long freqw, min_freqw, max_freqw;  // Freq offs (angle per sample)
    float pll_adjustment;
    bool allow_drift;                   // Follow carrier beyond safe limits
    static const unsigned int chunk_size = 128;
    
    fast_qpsk_receiver(scheduler *sch,
		       pipebuf< complex<T> > &_in,
		       pipebuf<hardsymbol> &_out,
		       pipebuf<float> *_freq_out=NULL,
		       pipebuf< complex<T> > *_cstln_out=NULL)
      : runnable(sch, "Fast QPSK receiver"),
	meas_decimation(1048576),
	pll_adjustment(1.0),
	allow_drift(false),
	in(_in), out(_out, chunk_size),
	mu(0), phase(0),
	meas_count(0)
    {
      set_omega(1);
      set_freq(0);
      freq_out = _freq_out ? new pipewriter<float>(*_freq_out) : NULL;
      cstln_out = _cstln_out ? new pipewriter< complex<T> >(*_cstln_out) : NULL;
      memset(hist, 0, sizeof(hist));
      init_lookup_tables();
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
      // Prevent PLL from locking at +-symbolrate/4.
      // TODO The +-SR/8 limit is suitable for QPSK only.
      min_freqw = freqw - 65536/max_omega/8;
      max_freqw = freqw + 65536/max_omega/8;
    }

    static const int RLUT_BITS = 8;
    static const int RLUT_ANGLES = 1 << RLUT_BITS;

    void run() {
      // Magic constants that work with the qa recordings.
      signed long freq_alpha = 0.04 * 65536;
      signed long freq_beta = 0.0012 * 256 * 65536 / omega * pll_adjustment;
      if ( ! freq_beta ) fail("Excessive oversampling");

      float gain_mu = 0.02 / (cstln_amp*cstln_amp) * 2;
      
      int max_meas = chunk_size/meas_decimation + 1;
      // Largin margin on output_size because mu adjustments
      // can lead to more than chunk_size/min_omega symbols.
      while ( in.readable() >= chunk_size+1 &&  // +1 for interpolation
	      out.writable() >= chunk_size &&
	      ( !freq_out  || freq_out ->writable()>=max_meas ) &&
	      ( !cstln_out || cstln_out->writable()>=max_meas ) ) {
	
	complex<T> *pin=in.rd(), *pin0=pin, *pend=pin+chunk_size;
	hardsymbol *pout=out.wr(), *pout0=pout;

	cu8 s;
	u_angle symbol_arg = 0;  // Exported for constellation viewer
	
	while ( pin < pend ) {
	  // Here mu is the time of the next symbol counted from 0 at pin.
	  if ( mu < 1 ) {
	    // Here 0<=mu<1 is the fractional time of the next symbol
	    // between pin and pin+1.

	    // Derotate and interpolate
#if 0  // Phase only (does not work)
	    // Careful with the float/signed/unsigned casts
	    u_angle a0 = fast_arg(pin[0]) - phase;
	    u_angle a1 = fast_arg(pin[1]) - (phase+freqw);
	    s_angle da = a1 - a0;
	    symbol_arg = a0 + (s_angle)(da*mu);
	    s = arg_to_symbol(symbol_arg);
#elif 1  // Linear by lookup-table. 1.2M on bench3bishs
	    polar *p0 = &lut_polar[pin[0].re][pin[0].im];
	    u_angle a0 = (u_angle)(p0->a-phase) >> (16-RLUT_BITS);
	    cu8 *p0r = &lut_rect[a0][p0->r>>1];
	    polar *p1 = &lut_polar[pin[1].re][pin[1].im];
	    u_angle a1 = (u_angle)(p1->a-(phase+freqw)) >> (16-RLUT_BITS);
	    cu8 *p1r = &lut_rect[a1][p1->r>>1];
	    s.re = (int)(p0r->re + (p1r->re-p0r->re)*mu);
	    s.im = (int)(p0r->im + (p1r->im-p0r->im)*mu);
	    symbol_arg = fast_arg(s);
#else  // Linear floating-point, for reference
	    float a0 = -(int)phase*M_PI/32768;
	    float cosa0=cosf(a0), sina0=sinf(a0);
	    complex<float>
	      p0r(((float)pin[0].re-128)*cosa0 - ((float)pin[0].im-128)*sina0,
		  ((float)pin[0].re-128)*sina0 + ((float)pin[0].im-128)*cosa0);
	    float a1 = -(int)(phase+freqw)*M_PI/32768;
	    float cosa1=cosf(a1), sina1=sinf(a1);
	    complex<float>
	      p1r(((float)pin[1].re-128)*cosa1 - ((float)pin[1].im-128)*sina1,
		  ((float)pin[1].re-128)*sina1 + ((float)pin[1].im-128)*cosa1);
	    s.re = (int)(128 + p0r.re + (p1r.re-p0r.re)*mu);
	    s.im = (int)(128 + p0r.im + (p1r.im-p0r.im)*mu);
	    symbol_arg = fast_arg(s);
#endif

	    int quadrant = symbol_arg >> 14;
	    static unsigned char quadrant_to_symbol[4] = { 0, 2, 3, 1 };
	    *pout = quadrant_to_symbol[quadrant];
	    ++pout;

	    // PLL
	    s_angle phase_error = (s_angle)(symbol_arg&16383) - 8192;
	    phase += (phase_error * freq_alpha + 32768) >> 16;
	    freqw += (phase_error * freq_beta + 32768*256) >> 24;
	    
	    // Modified Mueller and Müller
	    // mu[k]=real((c[k]-c[k-2])*conj(p[k-1])-(p[k]-p[k-2])*conj(c[k-1]))
	    //      =dot(c[k]-c[k-2],p[k-1]) - dot(p[k]-p[k-2],c[k-1])
	    // p = received signals
	    // c = decisions (constellation points)
	    hist[2] = hist[1];
	    hist[1] = hist[0];
#define HIST_FLOAT 0
#if HIST_FLOAT
	    hist[0].p.re = (float)s.re - 128;
	    hist[0].p.im = (float)s.im - 128;

	    cu8 cp = arg_to_symbol((symbol_arg&49152)+8192);
	    hist[0].c.re = (float)cp.re - 128;
	    hist[0].c.im = (float)cp.im - 128;

	    float muerr =
	      ( (hist[0].p.re-hist[2].p.re)*hist[1].c.re +
		(hist[0].p.im-hist[2].p.im)*hist[1].c.im ) -
	      ( (hist[0].c.re-hist[2].c.re)*hist[1].p.re +
		(hist[0].c.im-hist[2].c.im)*hist[1].p.im );
#else
	    hist[0].p = s;
	    hist[0].c = arg_to_symbol((symbol_arg&49152)+8192);

	    int muerr =
	      ( (signed char)(hist[0].p.re-hist[2].p.re)*((int)hist[1].c.re-128) +
		(signed char)(hist[0].p.im-hist[2].p.im)*((int)hist[1].c.im-128) ) -
	      ( (signed char)(hist[0].c.re-hist[2].c.re)*((int)hist[1].p.re-128) +
		(signed char)(hist[0].c.im-hist[2].c.im)*((int)hist[1].p.im-128) );
#endif
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

	if ( symbol_arg && cstln_out )
	  // Output the last interpolated PSK symbol, max once per chunk_size
	  cstln_out->write(s);
	
	// This is best done periodically ouside the inner loop,
	// but will cause non-deterministic output.
	
	if ( ! allow_drift ) {
	  if ( freqw < min_freqw || freqw > max_freqw )
	    freqw = (max_freqw+min_freqw) / 2;
	}
      
	// Output measurements
	
	meas_count += pin-pin0;
	while ( meas_count >= meas_decimation ) {
	  meas_count -= meas_decimation;
	  if ( freq_out )
	    freq_out->write((float)freqw / 65536);
	}
	
      }  // Work to do
  }
    
  private:

    struct polar { u_angle a; unsigned char r; } lut_polar[256][256];
    u_angle fast_arg(const cu8 &c) {
      // TBD read cu8 as u16 index, same endianness as in init()
      return lut_polar[c.re][c.im].a;
    }
    cu8 lut_rect[RLUT_ANGLES][256];
    cu8 lut_sincos[65536];
    cu8 arg_to_symbol(u_angle a) { return lut_sincos[a]; }
    void init_lookup_tables() {
      for ( int i=0; i<256; ++i )
	for ( int q=0; q<256; ++q ) {
	  // Don't cast float to unsigned directly
	  lut_polar[i][q].a = (s_angle)(atan2f(q-128,i-128)*65536/(2*M_PI));
	  lut_polar[i][q].r = (int)hypotf(i-128,q-128);
	}
      for ( unsigned long a=0; a<65536; ++a ) {
	float f = 2*M_PI * a / 65536;
	lut_sincos[a].re = 128 + cstln_amp*cosf(f);
	lut_sincos[a].im = 128 + cstln_amp*sinf(f);
      }
      for ( int a=0; a<RLUT_ANGLES; ++a )
	for ( int r=0; r<256; ++r ) {
	  lut_rect[a][r].re = (int)(128 + r*cos(2*M_PI*a/RLUT_ANGLES));
	  lut_rect[a][r].im = (int)(128 + r*sin(2*M_PI*a/RLUT_ANGLES));
	}
    }

    struct {
#if HIST_FLOAT
      complex<float> p;  // Received symbol
      complex<float> c;  // Matched constellation point
#else
      cu8 p;  // Received symbol
      cu8 c;  // Matched constellation point
#endif
    } hist[3];
    pipereader<cu8> in;
    pipewriter<hardsymbol> out;
    float mu;  // PSK time expressed in clock ticks. TBD fixed point.
    u_angle phase;
    unsigned long meas_count;
    pipewriter<float> *freq_out, *mer_out;
    pipewriter<cu8> *cstln_out;
  };  // fast_qpsk_receiver
  

  // CONSTELLATION TRANSMITTER

  // Maps symbols to I/Q points.

  template<typename Tout, int Zout>
  struct cstln_transmitter : runnable {
    cstln_lut<256> *cstln;
    cstln_transmitter(scheduler *sch,
		      pipebuf<u8> &_in, pipebuf< complex<Tout> > &_out)
      : runnable(sch, "cstln_transmitter"),
	in(_in), out(_out)
    {
    }
    void run() {
      if ( ! cstln ) fail("constellation not set");
      int count = min(in.readable(), out.writable());
      u8 *pin=in.rd(), *pend=pin+count;
      complex<Tout> *pout = out.wr();
      for ( ; pin<pend; ++pin,++pout ) {
	complex<signed char> *cp = &cstln->symbols[*pin];
	pout->re = Zout + cp->re;
	pout->im = Zout + cp->im;
      }
      in.read(count);
      out.written(count);
    }
  private:
    pipereader<u8> in;
    pipewriter< complex<Tout> > out;
  };  // cstln_transmitter


  // FREQUENCY SHIFTER

  // Resolution is sample_freq/65536.
  
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
  };  // rotator


  // SPECTRUM-BASED CNR ESTIMATOR

  // Assumes that the spectrum is as follows:
  //
  //  ---|--noise---|-roll-off-|---carrier+noise----|-roll-off-|---noise--|---
  //     |  (bw/2)  |   (bw)   |       (bw/2)       |   (bw)   |  (bw/2)  |
  //
  // Maximum roll-off 0.5

  template<typename T>
  struct cnr_fft : runnable {
    cnr_fft(scheduler *sch, pipebuf< complex<T> > &_in, pipebuf<float> &_out,
	    float _bandwidth, int nfft=4096)
      : runnable(sch, "cnr_fft"),
	bandwidth(_bandwidth), freq_tap(NULL), tap_multiplier(1),
	decimation(1048576), kavg(0.1),
	in(_in), out(_out),
	fft(nfft), avgpower(NULL), phase(0) {
      if ( bandwidth > 0.25 )
	fail("CNR estimator requires Fsampling > 4x Fsignal");
    }

    float bandwidth;
    float *freq_tap, tap_multiplier;    
    int decimation;
    float kavg;

    void run() {
      while ( in.readable()>=fft.n && out.writable()>=1 ) {
	phase += fft.n;
	if ( phase >= decimation ) {
	  phase -= decimation;
	  do_cnr();
	}
	in.read(fft.n);
      }
    }
    
  private:

    void do_cnr() {
      float center_freq = freq_tap ? *freq_tap * tap_multiplier : 0;
      int icf = floor(center_freq*fft.n+0.5);
      complex<T> data[fft.n];
      memcpy(data, in.rd(), fft.n*sizeof(data[0]));
      fft.inplace(data, true);
      T power[fft.n];
      for ( int i=0; i<fft.n; ++i )
	power[i] = data[i].re*data[i].re + data[i].im*data[i].im;
      if ( ! avgpower ) {
	// Initialize with first spectrum
	avgpower = new T[fft.n];
	memcpy(avgpower, power, fft.n*sizeof(avgpower[0]));
      }
      // Accumulate and low-pass filter
      for ( int i=0; i<fft.n; ++i )
	avgpower[i] = avgpower[i]*(1-kavg) + power[i]*kavg;
      
      int bwslots = (bandwidth/4) * fft.n;
      if ( ! bwslots ) return;
      // Measure carrier+noise in center band
      float c2plusn2 = avgslots(icf-bwslots, icf+bwslots);
      // Measure noise left and right of roll-off zones
      float n2 = ( avgslots(icf-bwslots*4, icf-bwslots*3) +
		   avgslots(icf+bwslots*3, icf+bwslots*4) ) / 2;
      float c2 = c2plusn2 - n2;
      float cnr = (c2>0 && n2>0) ? 10 * logf(c2/n2)/logf(10) : -50;
      out.write(cnr);
    }

    float avgslots(int i0, int i1) {  // i0 <= i1
      T s = 0;
      for ( int i=i0; i<=i1; ++i ) s += avgpower[i&(fft.n-1)];
      return s / (i1-i0+1);
    }
    
    pipereader< complex<T> > in;
    pipewriter< float > out;
    cfft_engine<T> fft;
    T *avgpower;
    int phase;
  };  // cnr_fft
  
}  // namespace

#endif  // LEANSDR_SDR_H
