#ifndef	__OS_MATH_COMPAT_H__
#define	__OS_MATH_COMPAT_H__

#include <math.h>

/*
 * FreeBSD doesn't have these functions defined.
 */

#if defined(__FreeBSD__)
static inline void
sincos(double x, double *s, double *c)
{
	*s = sin(x);
	*c = cos(x);
}

static inline void
sincosf(float x, float *s, float *c)
{

	*s = sinf(x);
	*c = cosf(x);
}

static inline void
sincosl(long double x, long double *s, long double *c)
{

	*s = sinl(x);
	*c = cosl(x);
}
#endif	/* __FreeBSD__ */

#endif	/* __OS_MATH_COMPAT_H__ */
