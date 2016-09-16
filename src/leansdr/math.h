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
  complex<T> operator *(const complex<T> &a, const T &k) {
    return complex<T>(a.re*k, a.im*k);
  }

}  // namespace

#endif  // LEANSDR_MATH_H
