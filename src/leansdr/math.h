#ifndef LEANSDR_MATH_H
#define LEANSDR_MATH_H

#include <math.h>

namespace leansdr {

  template<typename T>
  struct complex {
    T re, im;
    complex() { }
    complex(T x) : re(x), im(0) { }
    complex(T x, T y) : re(x), im(y) { }
    inline void operator +=(const complex<T> &x) { re+=x.re; im+=x.im; }
  };

  template<typename T>
  complex<T> operator +(const complex<T> &a, const complex<T> &b) {
    return complex<T>(a.re+b.re, a.im+b.im);
  }

  template<typename T>
  complex<T> operator *(const complex<T> &a, const complex<T> &b) {
    return complex<T>(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re);
  }

  template<typename T>
  complex<T> operator *(const complex<T> &a, const T &k) {
    return complex<T>(a.re*k, a.im*k);
  }

  template<typename T>
  complex<T> operator *(const T &k, const complex<T> &a) {
    return complex<T>(k*a.re, k*a.im);
  }
  
  // TBD Optimize with dedicated instructions
  int hamming_weight(unsigned char x) {
    static const int lut[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };
    return lut[x&15] + lut[x>>4];
  }
  int hamming_weight(unsigned short x) {
    return hamming_weight((unsigned char)x)
      +    hamming_weight((unsigned char)(x>>8));
  }
  int hamming_weight(unsigned long x) {
    return hamming_weight((unsigned short)x)
      +    hamming_weight((unsigned short)(x>>16));
  }
  int hamming_weight(unsigned long long x) {
    return hamming_weight((unsigned long)x)
      +    hamming_weight((unsigned long)(x>>32));
  }

  unsigned char parity(unsigned char x) {
    x ^= x>>4;
    return (0x6996 >> (x&15)) & 1;  // 16-entry look-up table
  }
  unsigned char parity(unsigned short x) {
    return parity((unsigned char)(x^(x>>8)));
  }
  unsigned char parity(unsigned long x) {
    return parity((unsigned short)(x^(x>>16)));
  }
  unsigned char parity(unsigned long long x) {
    return parity((unsigned long)(x^(x>>32)));
  }

  
  // Pre-computed sin/cos for 16-bit angles

  struct trig16 {
    complex<float> lut[65536];  // TBD static and shared
    trig16() {
      for ( int a=0; a<65536; ++a ) {
	float af = a * 2*M_PI / 65536;
	lut[a].re = cosf(af);
	lut[a].im = sinf(af);
      }
    }
    inline const complex<float> &expi(unsigned short a) const {
      return lut[a];
    }
    // a must fit in a signed long, otherwise behaviour is undefined
    inline const complex<float> &expi(float a) const {
      return expi((unsigned short)(signed short)(signed long)a);
    }
  };

}  // namespace

#endif  // LEANSDR_MATH_H
