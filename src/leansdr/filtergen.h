#ifndef LEANSDR_FILTERGEN_H
#define LEANSDR_FILTERGEN_H

#include <math.h>

namespace leansdr {
  namespace filtergen {
  
    template<typename T>
    void normalize_coeffs(int n, T *coeffs) {
      T s = 0;
      for ( int i=0; i<n; ++i ) s = s + coeffs[i];
      T k = (T)1.0 / s;
      for ( int i=0; i<n; ++i ) coeffs[i] = coeffs[i] * k;
    }

    // Generate coefficients for a sinc filter.

    template<typename T>
    int lowpass(int order, float bw, T **coeffs) {
      int ncoeffs = order + 1;
      *coeffs = new T[ncoeffs];
      for ( int i=0; i<ncoeffs; ++i ) {
	float t = i - (ncoeffs-1)*0.5;
	float sinc = 2*bw * (t ? sin(2*M_PI*bw*t)/(2*M_PI*bw*t) : 1);
#if 0  // Hamming 
	float alpha = 25.0/46, beta = 21.0/46;
	float window = alpha - beta*cos(2*M_PI*i/order);
#else
	float window = 1;
#endif
	(*coeffs)[i] = sinc * window;
      }
      normalize_coeffs(ncoeffs, *coeffs);
      return ncoeffs;
    }
    
  }  // namespace
}  // namespace

#endif  // LEANSDR_FILTERGEN_H
