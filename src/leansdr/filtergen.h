#ifndef LEANSDR_FILTERGEN_H
#define LEANSDR_FILTERGEN_H

#include <math.h>

namespace leansdr {
  namespace filtergen {
  
    template<typename T>
    void normalize_power(int n, T *coeffs, float gain=1) {
      float s2 = 0;
      for ( int i=0; i<n; ++i ) s2 = s2 + coeffs[i]*coeffs[i];  // TBD complex
      if ( s2 ) gain /= gen_sqrt(s2);
      for ( int i=0; i<n; ++i ) coeffs[i] = coeffs[i] * gain;
    }

    template<typename T>
    void normalize_dcgain(int n, T *coeffs, float gain=1) {
      float s = 0;
      for ( int i=0; i<n; ++i ) s = s + coeffs[i];
      if ( s ) gain /= s;
      for ( int i=0; i<n; ++i ) coeffs[i] = coeffs[i] * gain;
    }

    // Generate coefficients for a sinc filter.

    template<typename T>
    int lowpass(int order, float Fcut, T **coeffs, float gain=1) {
      int ncoeffs = order + 1;
      *coeffs = new T[ncoeffs];
      for ( int i=0; i<ncoeffs; ++i ) {
	float t = i - (ncoeffs-1)*0.5;
	float sinc = 2*Fcut * (t ? sin(2*M_PI*Fcut*t)/(2*M_PI*Fcut*t) : 1);
#if 0  // Hamming 
	float alpha = 25.0/46, beta = 21.0/46;
	float window = alpha - beta*cos(2*M_PI*i/order);
#else
	float window = 1;
#endif
	(*coeffs)[i] = sinc * window;
      }
      normalize_dcgain(ncoeffs, *coeffs, gain);
      return ncoeffs;
    }
    
    template<typename T>
    int root_raised_cosine (float Fm, float Fs, float rolloff, T **coeffs) {
      int ncoeffs = 0;
      *coeffs = new T[ncoeffs];
      fail("not implemented");
      return ncoeffs;
    }

  }  // namespace
}  // namespace

#endif  // LEANSDR_FILTERGEN_H
