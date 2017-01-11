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
  complex<T> operator *(const T &k, complex<T> &a) {
    return complex<T>(k*a.re, k*a.im);
  }
  
  // TBD Optimize with dedicated instructions
  int hamming_weight(unsigned char x) {
    static int lut[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };
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

}  // namespace

#endif  // LEANSDR_MATH_H
