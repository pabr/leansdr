#ifndef LEANSDR_VITERBI_H
#define LEANSDR_VITERBI_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace leansdr {

  // TS is an integer type for a least NSTATES+1 states.
  // NSTATES is the number of states (e.g. 2^(K-1)).
  // TUS is an integer type for uncoded symbols (branch identifiers).
  // NUS is the number of uncoded symbols.
  // TCS is an integer type for coded symbols (branch labels).
  // NCS is the number of coded symbols.
  // TP is a type for representing paths.
  // TPM, TBM are unsigned integer types for path/branch metrics.
  // TPM is at least as wide as TBM.

  template<typename TS, int NSTATES, typename TUS, int NUS, int NCS>
  struct trellis {
    static const int NOSTATE = NSTATES+1;

    struct state {
      struct branch {
	TS pred;         // Predecessor state or NOSTATE
	TUS us;          // Uncoded symbol
      } branches[NCS];   // Incoming branches indexed by coded symbol
    } states[NSTATES];
    
    trellis() {
      for ( TS s=0; s<NSTATES; ++s )
	for ( int cs=0; cs<NCS; ++cs )
	  states[s].branches[cs].pred = NOSTATE;
    }

    void init_convolutional(uint64_t G[]) {
      if ( NCS & (NCS-1) )  {
	fprintf(stderr, "NCS must be a power of 2\n");
	exit(1);
      }
      // Derive number of polynomials from NCS.
      int nG;
      for ( nG=1; (1<<nG)<NCS; ++nG ) ;

      for ( TS s=0; s<NSTATES; ++s ) {
	for ( TUS us=0; us<NUS; ++us ) {
	  // Run the convolutional encoder from state s with input us
	  uint64_t shiftreg = s | (us*NSTATES);
	  uint32_t cs = 0;
	  for ( int g=0; g<nG; ++g )
	    cs = (cs<<1) | parity(shiftreg&G[g]);
	  shiftreg /= NUS;  // Shift bits for 1 uncoded symbol
	  // [us] at state [s] emits [cs] and leads to state [shiftreg].
	  typename state::branch *b = &states[shiftreg].branches[cs];
	  if ( b->pred != NOSTATE ) {
	    fprintf(stderr, "Invalid convolutional code\n");
	    exit(1);
	  }
	  b->pred = s;
	  b->us = us;
	}
      }

    }

    void dump() {
      for ( TS s=0; s<NSTATES; ++s ) {
	fprintf(stderr, "State %02x:", s);
	for ( int cs=0; cs<NCS; ++cs ) {
	  typename state::branch *b = &states[s].branches[cs];
	  if ( b->pred == NOSTATE )
	    fprintf(stderr, "     - ");
	  else
	    fprintf(stderr, "   %02x+%x", b->pred, b->us);
	}
	fprintf(stderr, "\n");
      }
    }

  };

  template<typename TS, int NSTATES,
	   typename TUS, int NUS,
	   typename TCS, int NCS,
	   typename TBM, typename TPM,
	   typename TP>
  struct viterbi_dec {
    trellis<TS, NSTATES, TUS, NUS, NCS> *trell;

    struct state {
      TPM cost;    // Metric of best path leading to this state
      TP path;     // Best path leading to this state
    };
    typedef state statebank[NSTATES];
    state statebanks[2][NSTATES];
    statebank *states, *newstates;  // Alternate between banks

    viterbi_dec(trellis<TS, NSTATES, TUS, NUS, NCS> *_trellis) :
      trell(_trellis)
    {
      states = &statebanks[0];
      newstates = &statebanks[1];
      for ( TS s=0; s<NSTATES; ++s ) (*states)[s].cost = 0;
    }
    
    TUS update(TBM costs[NCS], TPM *quality=NULL) {
      TPM best_tpm = -1LL, best2_tpm = -1LL;
      TS best_state = 0;
      // Update all states
      for ( TS s=0; s<NSTATES; ++s ) {
	TPM best_m = -1LL;
	typename trellis<TS,NSTATES,TUS,NUS,NCS>::state::branch *best_b = NULL;
	// Select best branch
	for ( TCS cs=0; cs<NCS; ++cs ) {
	  typename trellis<TS,NSTATES,TUS,NUS,NCS>::state::branch *b =
	    &trell->states[s].branches[cs];
	  if ( b->pred == trell->NOSTATE ) continue;
	  TPM m = (*states)[b->pred].cost + costs[cs];
	  if ( m <= best_m ) {  // <= guarantees one match
	    best_m = m;
	    best_b = b;
	  }
	}
	(*newstates)[s].path = (*states)[best_b->pred].path;
	(*newstates)[s].path.append(best_b->us);
	(*newstates)[s].cost = best_m;
	// Select best and second-best states
	if ( best_m < best_tpm ) {
	  best_state = s;
	  best2_tpm = best_tpm;
	  best_tpm = best_m;
	} else if ( best_m < best2_tpm )
	  best2_tpm = best_m;
      }
      // Swap banks
      { statebank *tmp=states; states=newstates; newstates=tmp; }
      // Prevent overflow of path metrics
      for ( TS s=0; s<NSTATES; ++s ) (*states)[s].cost -= best_tpm;
#if 0
      // Observe that the min-max range remains bounded
      fprintf(stderr,"-%2d = [", best_tpm);
      for ( TS s=0; s<NSTATES; ++s ) fprintf(stderr," %d", (*states)[s].cost);
      fprintf(stderr," ]\n");
#endif
      // Return difference between best and second-best as quality metric.
      if ( quality ) *quality = best2_tpm - best_tpm;
      // Return uncoded symbol of best path
      return (*states)[best_state].path.read();
    }

    void dump() {
      fprintf(stderr, "[");
      for ( TS s=0; s<NSTATES; ++s )
	if ( states[s].cost )
	  fprintf(stderr, " %02x:%d", s, states[s].cost);
      fprintf(stderr, "\n");
    }

  };

  // Paths (sequences of uncoded symbols) represented as bitstreams.
  // NBITS is the number of bits per symbol.
  // DEPTH is the number of symbols stored in the path.
  // T is an unsigned integer type wider than NBITS*DEPTH.

  template<typename T, typename TUS, int NBITS, int DEPTH>
  struct bitpath {
    T val;
    bitpath() : val(0) { }
    void append(TUS us) { val = (val<<NBITS) | us; }
    TUS read() { return (val>>(DEPTH-1)*NBITS) & ((1<<NBITS)-1); }
  };


}  // namespace

#endif  // LEANSDR_VITERBI_H
