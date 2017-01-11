#ifndef LEANSDR_RS_H
#define LEANSDR_RS_H

#include "leansdr/math.h"

namespace leansdr {

  // Finite group GF(2^N).

  // GF(2) is the ring ({0,1},+,*).
  // GF(2)[X] is the ring of polynomials with coefficients in GF(2).
  // P(X) is an irreducible polynomial of GF(2)[X].
  // N is the degree of P(x).
  // P is the bitfield representation of P(X), with degree 0 at LSB.
  // GF(2)[X]/(P) is GF(2)[X] modulo P(X).
  // (GF(2)[X]/(P), +) is a group with 2^N elements.
  // (GF(2)[X]/(P)*, *) is a group with 2^N-1 elements.
  // (GF(2)[X]/(P), +, *) is a field with 2^N elements, noted GF(2^N).
  // Te is a C++ integer type for representing elements of GF(2^N).
  //   "0" is 0
  //   "1" is 1
  //   "2" is X
  //   "3" is X+1
  //   "4" is X^2
  // Tp is a C++ integer type for representing P(X) (1 bit larger than Te).
  // ALPHA is a primitive element of GF(2^N). Usually "2"=[X] is chosen.

  template<typename Te, typename Tp, Tp P, int N, Te ALPHA>
  struct gf2x_p {
    gf2x_p() {
      if ( ALPHA != 2 ) fail("alpha!=2 not implemented");
      // Precompute log and exp tables.
      Tp alpha_i = 1;
      for ( Tp i=0; i<(1<<N); ++i ) {
	lut_exp[i] = alpha_i;
	lut_exp[((1<<N)-1)+i] = alpha_i;
	lut_log[alpha_i] = i;      
	alpha_i <<= 1;  // Multiply by alpha=[X] i.e. increase degrees
	if ( alpha_i & (1<<N) ) alpha_i ^= P;  // Modulo P iteratively
      }
    }
    static const Te alpha = ALPHA;
    inline Te add(Te x, Te y) { return x ^ y; }  // Addition modulo 2
    inline Te sub(Te x, Te y) { return x ^ y; }  // Subtraction modulo 2
    inline Te mul(Te x, Te y) {
      if ( !x || !y ) return 0;
      return lut_exp[lut_log[x] + lut_log[y]];
    }
    inline Te div(Te x, Te y) {
      //if ( ! y ) fail("div"); // TODO
      if ( ! x ) return 0;
      return lut_exp[lut_log[x] + ((1<<N)-1) - lut_log[y]];
    }
    inline Te inv(Te x) {
      //    if ( ! x ) fail("inv");
      return lut_exp[((1<<N)-1) - lut_log[x]];
    }
    inline Te exp(Te x) { return lut_exp[x]; }
    inline Te log(Te x) { return lut_log[x]; }
  private:
    Te lut_exp[(1<<N)*2];  // Wrap to avoid indexing modulo 2^N-1 
    Te lut_log[1<<N];
  };

  // Reed-Solomon for RS(204,188) shortened from RS(255,239).

  struct rs_engine {
    // EN 300 421, section 4.4.2, Field Generator Polynomial
    // p(X) = X^8 + X^4 + X^3 + X^2 + 1
    gf2x_p<unsigned char, unsigned short, 0x11d, 8, 2> gf;

    // RS-encoded messages are interpreted as coefficients in
    // GF(256) of a polynomial P.
    // The syndromes are synd[i] = P(alpha^i).
    // By convention coefficients are listed by decreasing degree here,
    // so we can evaluate syndromes of the shortened code without
    // prepending with 51 zeroes.
    bool syndromes(const u8 *poly, u8 *synd) {
      bool corrupted = false;
      for ( int i=0; i<16; ++i ) {
	synd[i] = eval_poly_rev(poly, 204, gf.exp(i));
	if ( synd[i] ) corrupted = true;
      }
      return corrupted;
    }
    u8 eval_poly_rev(const u8 *poly, int n, u8 x) {
      // poly[0]*x^(n-1) + .. + poly[n-1]*x^0 with Hörner method.
      u8 acc = 0;
      for ( int i=0; i<n; ++i ) acc = gf.add(gf.mul(acc,x), poly[i]);
      return acc;
    }

    // Evaluation with coefficients listed by increasing degree.
    u8 eval_poly(const u8 *poly, int deg, u8 x) {
      // poly[0]*x^0 + .. + poly[deg]*x^deg with Hörner method.
      u8 acc = 0;
      for ( ; deg>=0; --deg ) acc = gf.add(gf.mul(acc,x), poly[deg]);
      return acc;
    }

    // Try to fix errors in pout[].
    // If pin[] is provided, errors will be fixed in the original
    // message too and syndromes will be updated.

#define DEBUG_RS 0
  
    bool correct(u8 synd[16], u8 pout[188],
		 u8 pin[204]=NULL, int *bits_corrected=NULL) {
      // Berlekamp - Massey
      // http://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm#Code_sample
      u8 C[16] = { 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };  // Max degree is L
      u8 B[16] = { 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
      int L = 0;
      int m = 1;
      u8 b = 1;
      for ( int n=0; n<16; ++n ) {
	u8 d = synd[n];
	for ( int i=1; i<=L; ++i ) d ^= gf.mul(C[i], synd[n-i]);
	if ( ! d ) {
	  ++m;
	} else if ( 2*L <= n ) {
	  u8 T[16];
	  memcpy(T, C, sizeof(T));
	  for ( int i=0; i<16-m; ++i )
	    C[m+i] ^= gf.mul(d, gf.mul(gf.inv(b),B[i]));
	  L = n + 1 - L;
	  memcpy(B, T, sizeof(B));
	  b = d;
	  m = 1;
	} else {
	  for ( int i=0; i<16-m; ++i )
	    C[m+i] ^= gf.mul(d, gf.mul(gf.inv(b),B[i]));
	  ++m;
	}
      }
      // L is the number of errors
      // C of degree L is now the error locator polynomial (Lambda)
#if DEBUG_RS
      fprintf(stderr, "[L=%d  C=",L);
      for ( int i=0; i<16; ++i ) fprintf(stderr, " %d", C[i]);
      fprintf(stderr, "]\n");
      fprintf(stderr, "[S=");
      for ( int i=0; i<16; ++i ) fprintf(stderr, " %d", synd[i]);
      fprintf(stderr, "]\n");
#endif
    
      // Forney
      // http://en.wikipedia.org/wiki/Forney_algorithm  (2t=16)
    
      // Compute Omega
      u8 omega[16];
      memset(omega, 0, sizeof(omega));
      // TODO loops
      for ( int i=0; i<16; ++i )
	for ( int j=0; j<16; ++j )
	  if ( i+j < 16 ) omega[i+j] ^= gf.mul(synd[i], C[j]);
#if DEBUG_RS
      fprintf(stderr, "omega=");
      for ( int i=0; i<16; ++i ) fprintf(stderr, " %d", omega[i]);
      fprintf(stderr, "\n");
#endif
    
      // Compute Lambda'
      u8 Cprime[15];
      for ( int i=0; i<15; ++i )
	Cprime[i] = (i&1) ? 0 : C[i+1];
#if DEBUG_RS
      fprintf(stderr, "Cprime=");
      for ( int i=0; i<15; ++i ) fprintf(stderr, " %d", Cprime[i]);
      fprintf(stderr, "\n");
#endif
    
      // Find zeroes of C by exhaustive search?
      // TODO Chien method
      int roots_found = 0;
      for ( int i=0; i<255; ++i ) {
	u8 r = gf.exp(i);  // Candidate root alpha^0..alpha^254
	u8 v = eval_poly(C, L, r);
	if ( ! v ) {
	  // r is a root X_k^-1 of the error locator polynomial.
	  u8 xk = gf.inv(r);
	  int loc = (255-i) % 255;  // == log(xk)
#if DEBUG_RS
	  fprintf(stderr, "found root=%d, inv=%d, loc=%d\n", r, xk, loc);
#endif
	  if ( loc < 204 ) {
	    // Evaluate e_k
	    u8 num = gf.mul(xk, eval_poly(omega, L, r));
	    u8 den = eval_poly(Cprime, 14, r);
	    u8 e = gf.div(num, den);
	    // Subtract e from coefficient of degree loc.
	    // Note: Coeffients listed by decreasing degree in pin[] and pout[].
	    if ( bits_corrected ) *bits_corrected += hamming_weight(e);
	    if ( loc >= 16 ) pout[203-loc] ^= e;
	    if ( pin ) pin[203-loc] ^= e;
	  }
	  if ( ++roots_found == L ) break;
	}
      }
      if ( pin )
	return syndromes(pin, synd);
      else
	return false;
    }

  };

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
	if ( bitcount ) { *bitcount->wr()=nbits; bitcount->written(1); }
	if ( errcount ) { *errcount->wr()=nerrs; errcount->written(1); }
      }
    }
  private:
    pipereader< rspacket<Tbyte> > in;
    pipewriter<tspacket> out;
    pipewriter<int> *bitcount, *errcount;
  };

}  // namespace

#endif  // LEANSDR_RS_H
