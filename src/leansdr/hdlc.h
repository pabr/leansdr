#ifndef LEANSDR_HDLC_H
#define LEANSDR_HDLC_H

#include "leansdr/framework.h"

namespace leansdr {

  // HDLC deframer
					  
  struct hdlc_dec {

    hdlc_dec(int _minframesize,  // Including CRC, excluding HDLC flags.
	     int _maxframesize,
	     bool _invert)
      : minframesize(_minframesize), maxframesize(_maxframesize),
	invertmask(_invert?0xff:0),
	framebuf(new u8[maxframesize]),
	debug(false)
    {
      reset();
    }
    
    void reset() { shiftreg=0; inframe=false; }

    void begin_frame() { framesize=0; crc16=crc16_init; }
      
    // Decode (*ppin)[count] as MSB-packed HDLC bitstream.
    // Return pointer to buffer[*pdatasize], or NULL if no valid frame.
    // Return number of discarded bytes in *discarded.
    // Return number of checksum errors in *fcs_errors.

    u8 *decode(u8 **ppin, int count,
	       int *pdatasize, int *discarded, int *fcs_errors) {
      *discarded = 0;
      *fcs_errors = 0;
      *pdatasize = -1;
      u8 *pin=*ppin, *pend=pin+count;
      for ( ; pin<pend; ++pin ) {
	if ( *pdatasize != -1 ) {
	  // The previous loop found a complete frame
	  *ppin = pin;
	  return framebuf;
	}
	u8 byte_in = *pin;
	byte_in ^= invertmask;
	for ( int bits=8; bits--; byte_in<<=1 ) {
	  u8 bit_in = byte_in & 128;
	  shiftreg = (shiftreg>>1) | bit_in;
	  if ( ! inframe ) {
	    if ( shiftreg == 0x7e ) {  // HDLC flag 01111110
	      inframe = true;
	      nbits_out = 0;
	      begin_frame();
	    }
	  } else {
	    if ( (shiftreg&0xfe) == 0x7c ) {  // 0111110x HDLC stuffing
	      // Unstuff this 0
	    } else if ( shiftreg == 0x7e ) {  // 01111110 HDLC flag
	      if ( nbits_out != 7 ) {
		if ( debug ) fprintf(stderr, "^");
		++*discarded;
	      }
	      nbits_out = 0;
	      // Checksum
	      crc16 ^= 0xffff;
	      if ( framesize<2 || framesize<minframesize ||
		   crc16!=crc16_check ) {
		if ( debug ) fprintf(stderr, "!");
		*discarded += framesize;
		// Do not report random noise as FCS errors
		if ( framesize >= minframesize ) ++*fcs_errors;
	      } else {
		if ( debug ) fprintf(stderr, "_");
		*pdatasize = framesize-2;
	      }
	      begin_frame();
	      // Keep processing up to 7 remaining bits from byte_in.
	      // Special cases 0111111 and 1111111 cannot affect *pdatasize.
	    } else if ( shiftreg == 0xfe ) {  // 11111110 HDLC invalid
	      if ( framesize && debug ) fprintf(stderr, "^");
	      *discarded += framesize;
	      inframe = 0;
	    } else {  // Data bit
	      byte_out = (byte_out>>1) | bit_in;  // HDLC is LSB first
	      ++nbits_out;
	      if ( nbits_out == 8 ) {
		if ( framesize == maxframesize )
		  ++*discarded;
		else {
		  framebuf[framesize++] = byte_out;
		  crc16_byte(byte_out);
		}
		nbits_out = 0;
	      }
	    }
	  }  // inframe
	}  // bits
      }
      *ppin = pin;
      return NULL;
    }

  private:
    // Config
    int minframesize, maxframesize;
    u8 invertmask;
    u8 *framebuf;   // [maxframesize]
    // State
    u8 shiftreg;    // Input bitstream
    bool inframe;   // Currently receiving a frame ?
    u8 byte_out;    // Accumulator for data bits
    int nbits_out;  // Number of data bits in byte_out
    int framesize;  // Number of bytes in framebuf, if inframe
    u16 crc16;      // CRC of framebuf[framesize]
    // CRC
    static const u16 crc16_init = 0xffff;
    static const u16 crc16_poly = 0x8408;  // 0x1021 MSB-first
    static const u16 crc16_check = 0x0f47;
    void crc16_byte(u8 data) {
      crc16 ^= data;
      for ( int bit=8; bit--; )
	crc16 = (crc16&1) ? (crc16>>1)^crc16_poly : (crc16>>1);
    }

  public:
    bool debug;
  };  // hdlc_dec


  // HDLC synchronizer with polarity detection
    
  struct hdlc_sync : runnable {
    hdlc_sync(scheduler *sch,
	      pipebuf<u8> &_in,   // Packed bits
	      pipebuf<u8> &_out,  // Bytes
	      int _minframesize,  // Including CRC, excluding HDLC flags.
	      int _maxframesize)
      : runnable(sch, "hdlc_sync"),
	chunk_size(maxframesize+2),
	in(_in), out(_out, _maxframesize+chunk_size),
	maxframesize(_maxframesize),
	current_sync(0), resync_phase(0), resync_period(32),
	header16(false)
    {
      for ( int s=0; s<NSYNCS; ++s )
	syncs[s] = new hdlc_dec(minframesize, maxframesize, s!=0);
    }

    void run() {
      // Note: hdlc_dec may already hold one frame ready for output.
      while ( in.readable() >= chunk_size &&
	      out.writable() >= maxframesize+chunk_size ) {
	if ( ! resync_phase ) {
	  // Once every resync_phase, try all decoders
	  int total_discarded[NSYNCS];
	  for ( int s=0; s<NSYNCS; ++s ) {
	    if ( s != current_sync ) syncs[s]->reset();
	    total_discarded[s] = 0;
	    for ( u8 *pin=in.rd(), *pend=pin+chunk_size; pin<pend; ) {
	      int datasize, discarded, fcs_errors;
	      u8 *f = syncs[s]->decode(&pin, pend-pin, &datasize,
				       &discarded, &fcs_errors);
	      total_discarded[s] += discarded;
	      if ( s==current_sync && f ) output_frame(f, datasize);
	    }
	  }
	  int best = current_sync;
	  for ( int s=0; s<NSYNCS; ++s )
	    if ( total_discarded[s] < total_discarded[best] ) best = s;
	  if ( best != current_sync ) {
	    if ( sch->debug ) fprintf(stderr, "[%d->%d]", current_sync, best);
	    syncs[current_sync]->debug = false;
	    current_sync = best;
	    syncs[current_sync]->debug = sch->debug;
	  }
	} else {  // resync_phase
	  for ( u8 *pin=in.rd(), *pend=pin+chunk_size; pin<pend; ) {
	    int datasize, discarded, fcs_errors;
	    u8 *f = syncs[current_sync]->decode(&pin, pend-pin, &datasize,
						&discarded, &fcs_errors);
	    if ( f ) output_frame(f, datasize);
	  }
	}  // resync_phase
	in.read(chunk_size);
	if ( ++resync_phase >= resync_period ) resync_phase = 0;
      }  // Work to do
    }

  private:
    void output_frame(u8 *f, int size) {
      if ( header16 )  {
	// Removed 16-bit CRC, add 16-bit prefix -> Still <= maxframesize.
	out.write(size >> 8);
	out.write(size & 255);
      }
      memcpy(out.wr(), f, size);
      out.written(size);
    }

    int chunk_size;
    pipereader<u8> in;
    pipewriter<u8> out;
    int minframesize, maxframesize;
    static const int NSYNCS = 2;  // Two possible polarities
    hdlc_dec *syncs[NSYNCS];
    int current_sync;
    int resync_phase;
  public:
    int resync_period;
    bool header16;  // Output length prefix
  };  // hdlc_sync

}  // namespace

#endif  // LEANSDR_HDLC_H
