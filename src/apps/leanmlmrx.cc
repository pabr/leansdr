// This file is part of LeanSDR Copyright (C) 2018 <pabr@pabr.org>.
// See the toplevel README for more information.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Multi-channel legacy decoder.

// Currently FM mono only.

// Multithreading is a quick hack; will be replaced with proper
// synchronization primitives.

#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include "fftw3.h"
#include <math.h>
#include <pthread.h>

// For control socket
#include <sys/socket.h>
#include <sys/un.h>

// For PMP
#include <sys/fcntl.h>
#include <sys/mman.h>

#define MLMRX_DEEMPH 1

void fatal(const char *s) { perror(s); exit(1); }
void fail(const char *s) { fprintf(stderr, "** leanmlmrx: %s\n", s); exit(1); }

// Max FM channels
static const int MAXCHANS = 201;  // 20 MHz, 100 kHz spacing

// Number of FFTs between thread synchronizations.
static const int wqsize = 1024;  // about 100 Hz for Fq=200kHz

struct ci16 { int16_t re; int16_t im; };

struct thread_data {
  struct runtime *run;
  pthread_t pth;
  struct job {
    ci16 *buf1;         // Leftover samples from previous chunk
    int nbuf1;
    ci16 *buf2;         // New samples
    volatile bool go;   // buf ready. Set by reader, cleared by fft.
    uint16_t *ph;       // Estimated phases [nchans]
    volatile bool done; // ph ready. Set by fft, cleared by joiner.
  } jobs[wqsize];
  int qread, qfft, qjoin;
  fftwf_complex *in, *out;
  fftwf_plan p;
  uint64_t nexecs;
};

static const int NTHREADS = 2;

inline void yield() {
  usleep(1000);
  //pthread_yield();
}

uint64_t nwaitr=0, nwaitf1=0, nwaitf2=0, nwaitj=0;

struct chan {
  double F;
  int ibin;           // Nearest left bin
  float bw[2][2][2];  // weight[bin0|1][line][col]
  bool enabled;
  uint16_t derot;
  // Runtime
  uint16_t prevph;    // Discriminator state
  float rms;          // For squelch
};

struct config {
  bool pmp;        // PMP input mode
  float Fs;        // Sample rate (Hz)
  float Fc;        // Center frequency (Hz)
  float Fq;        // Quadrature rate (Hz), 0=autoselect
  float maxdev;    // FM maximum deviation (Hz)
  float deemph;    // De-emphasis time constant (s)
  int N;           // FFT size for channelizer
  int nchans;
  chan chans[MAXCHANS];
  float squelch;   // RMS threshold (0..1, 0 = monitor)
  float Fau;       // Audio sample rate (Hz), 0=autoselect
  bool wav;        // Output wav header
  int fd_info;     // FD for aux output, or -1
  float info_rate; // Spectrum estimation rate (Hz)
  int spec_size;   // FFT size for spectrum display
  int spec_zoom;   // Spectrum zoom factor
  int fd_control;  // FD for control input, or -1
  
  config()
    : pmp(false), Fs(25.6e6), Fc(98e6), Fq(0),
      maxdev(75e3), deemph(50e-6),
      N(64), nchans(0), squelch(0),
      Fau(44100),
      wav(false), fd_info(-1), info_rate(1),
      spec_size(1024), spec_zoom(1),
      fd_control(-1)
  { };
};

struct spectrum_estimator {
  int size;
  int zoom;
  fftwf_complex *in, *out;
  fftwf_plan plan;
  float *avgpower;
  float smoothbw;
  spectrum_estimator(int _size)
    : size(_size), zoom(1), avgpower(NULL), smoothbw(0.5)
  {
    in  = (fftwf_complex*)fftwf_malloc(sizeof(*in )*size);
    out = (fftwf_complex*)fftwf_malloc(sizeof(*out)*size);
    plan = fftwf_plan_dft_1d(size, in, out, -1, FFTW_ESTIMATE);
  }
  void process(ci16 *buf) {
    for ( int i=0; i<size; ++i ) {
      in[i][0] = buf[i].re;
      in[i][1] = buf[i].im;
    }
    fftwf_execute(plan);
    float power[size];
    for ( int i=0; i<size; ++i )
      power[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
    if ( ! avgpower ) {
      avgpower = new float[size];
      memcpy(avgpower, power, sizeof(power));
    }
    for ( int i=0; i<size; ++i )
      avgpower[i] = power[i]*smoothbw + avgpower[i]*(1-smoothbw);
  }
  void output(FILE *f) {
    if ( zoom == 1 ) {
      // Use old format for backward compatibility
      fprintf(f, "SPECTRUM [");
      for ( int i=0; i<size; ++i ) {
	float db = 10 * log10f(avgpower[(size/2+i)&(size-1)]);
	fprintf(f, "%s%.1f", (i?",":""), db);
      }
      fprintf(f, "]\n");
    } else {
      int s = size / zoom;
      int pos = size/2 - s/2;
      fprintf(f, "SPECTRUMRANGE {\"size\":%d,\"pos\":%d,\"db\":[", size, pos);
      for ( int i=0; i<s; ++i ) {
	float db = 10 * log10f(avgpower[(size-s/2+i)&(size-1)]);
	fprintf(f, "%s%.1f", (i?",":""), db);
      }
      fprintf(f, "]}\n");
    }
  }
};

struct runtime {
  config *cfg;
  uint16_t lut_atan2[256][256];
  float lut_audioscale[MAXCHANS+1];
  
  pthread_t reader;
  thread_data threads[NTHREADS];

  // Input buffer
  int stride;
  int next_th;     // Worker thread for next FFT
  int nskip;       // Samples to skip before next FFT
  ci16 *leftover;  // Leftover samples for next FFT
  int nleftover;
  bool eof;        // Set by reader thread
  
  // Statistics
  ci16 iqmin, iqmax;
  uint64_t nsamples;

  // Aux output
  FILE *f_info;
  // Control input
  FILE *f_control;

  spectrum_estimator *spestim;
  uint64_t spestim_phase;
  
  runtime(config *_cfg) {
    cfg = _cfg;

    // Precompute atan2 table
    for ( int iy=0; iy<256; ++iy )
      for ( int ix=0; ix<256; ++ix ) {
	float a = atan2f((int8_t)(uint8_t)iy,(int8_t)(uint8_t)ix);
	lut_atan2[iy][ix] = (int16_t)(a*65536/(2*M_PI));
      }

    // Precompute audio scaling vs number of channels.
    // Assuming uncorrelated signals, we can scale by sqrt(N) instead of N.
    for ( int i=0; i<=MAXCHANS; ++i )
      lut_audioscale[i] = i ? 1/sqrtf(i) : 0;
    
    stride = floor(cfg->Fs/cfg->Fq + 0.5);
    if ( stride < cfg->N ) fail("FFT windows overlap");

    next_th = 0;
    nskip = 0;
    leftover = new ci16[stride];
    nleftover = 0;
    eof = false;
    iqmin.re = iqmin.im =  32767;
    iqmax.re = iqmax.im = -32768;
    nsamples = 0;
    
    // Aux output
    if ( cfg->fd_info < 0 ) {
      f_info = NULL;
      spestim = NULL;
    } else {
      f_info = fdopen(cfg->fd_info, "w");
      if ( ! f_info ) fatal("fdopen(fd_info)");
      spestim = new spectrum_estimator(cfg->spec_size);
      spestim->zoom = cfg->spec_zoom;
    }
    spestim_phase = 0;
    // Setup control input
    if ( cfg->fd_control < 0 )
      f_control = NULL;
    else {
      int flags = fcntl(cfg->fd_control, F_GETFL, 0);
      fcntl(cfg->fd_control, F_SETFL, flags|O_NONBLOCK);
      f_control = fdopen(cfg->fd_control, "r");
      if ( ! f_control ) fatal("fdopen(fd_control)");
    }
  }
};

void process_samples(runtime *run, ci16 *buf, int count) {
  config *cfg = run->cfg;

  // Spectrum estimation
  if ( run->f_info ) {
    run->spestim_phase += count;
    spectrum_estimator *spe = run->spestim;
    if ( run->spestim_phase>=cfg->Fs/cfg->info_rate && count>=spe->size ) {
      run->spestim_phase = 0;
      spe->process(buf);
      spe->output(run->f_info);
      // Also output generic information
      fprintf(run->f_info, "CHANNELS [");
      for ( int i=0; i<cfg->nchans; ++i )
	fprintf(run->f_info, "%s{\"freq\":\"%.5f\",\"enabled\":%s}",
		(i?",":""), cfg->chans[i].F/1e6,
		(cfg->chans[i].enabled?"true":"false"));
      fprintf(run->f_info, "]\n");
      fflush(run->f_info);
    }
  }
  
  run->nsamples += count;  

  // Inside a skipped segment ?
  if ( count < run->nskip ) {
    run->nskip -= count;
    return;
  }
  buf += run->nskip;
  count -= run->nskip;
  run->nskip = 0;

  while ( run->nleftover+count > cfg->N ) {
    // Enough samples for one FFT
    // in run->leftover[run->nleftover] followed by buf[count].
    {
      // Statistics on one sample per chunk
      ci16 *s = buf;
      if ( s->re < run->iqmin.re ) run->iqmin.re = s->re;
      if ( s->re > run->iqmax.re ) run->iqmax.re = s->re;
      if ( s->im < run->iqmin.re ) run->iqmin.im = s->im;
      if ( s->im > run->iqmax.re ) run->iqmax.im = s->im;
    }
    // Assign FFT job to next thread
    {
      thread_data *td = &run->threads[run->next_th];
      run->next_th = (run->next_th+1) % NTHREADS;
      // Add to queue
      thread_data::job *job = &td->jobs[td->qread];
      while ( job->go ) { ++nwaitr; yield(); }   // Busy
      job->buf1 = run->leftover;
      job->nbuf1 = run->nleftover;
      job->buf2 = buf;
      job->go = true;
      td->qread = (td->qread+1) % wqsize;
    }
    // Skip to end of chunk
    // TBD Allow overlapping here (careful if nleftover!=0)
    int nused = cfg->N - run->nleftover;
    buf += nused;
    count -= nused;
    run->nleftover = 0;
    // Skip to beginning of next chunk
    int skip = run->stride - cfg->N;
    if ( count > skip ) {
      buf += skip;
      count -= skip;
    } else {
      count = 0;
      run->nskip = skip - count;
      break;
    }
  }

  // Copy leftover samples, if any.
  // Cost of copying is negligible if batch size is much larger than N.
  memcpy(run->leftover+run->nleftover, buf, sizeof(ci16)*count);
  run->nleftover += count;
}

void poll_control(runtime *run) {
  if ( ! run->f_control ) return;  // Control channel not used
  char cmd[256];
  if ( ! fgets(cmd, sizeof(cmd), run->f_control) ) return;  // EWOULDBLOCK
  config *cfg = run->cfg;
  // Parse command on control channel
  int arg;
  if      ( sscanf(cmd,"MUTE %d",  &arg)==1 && arg>=0 && arg<cfg->nchans )
    cfg->chans[arg].enabled = false;
  else if ( sscanf(cmd,"UNMUTE %d",&arg)==1 && arg>=0 && arg<cfg->nchans )
    cfg->chans[arg].enabled = true;
  else if ( sscanf(cmd,"GET /MUTE=%d",  &arg)==1 && arg>=0 && arg<cfg->nchans )
    cfg->chans[arg].enabled = false;
  else if ( sscanf(cmd,"GET /UNMUTE=%d",&arg)==1 && arg>=0 && arg<cfg->nchans )
    cfg->chans[arg].enabled = true;
  else
    fprintf(stderr, "Ignoring unrecognized command '%s'\n", cmd);
}

// Reader thread (PMP streaming variant)

void thread_reader_pmp(runtime *run) {
  int fd = open("/dev/mem", O_RDONLY);
  static off_t memsize = (off_t)512 << 20;
  void *map = mmap(NULL, memsize, PROT_READ, MAP_SHARED, fd, 0);
  if ( map == MAP_FAILED ) fatal("mmap");
  char *physmem = (char*)map;
  
  struct {
    uint64_t magic;
    uint64_t physaddr;
    uint64_t size;
    uint64_t canary;
  } pointer;
  
  while ( fread(&pointer, sizeof(pointer), 1, stdin) == 1 ) {
    if ( pointer.magic != 0x504d5031 ) fatal("PMP: bad magic");
    char *buf = physmem + pointer.physaddr;
    if ( *(uint64_t*)buf == pointer.canary )
      process_samples(run, (ci16*)buf,  pointer.size/sizeof(ci16));
    else
      fprintf(stderr, "PMP: Buffer overrun\n");
    poll_control(run);
  }
}    

// Reader thread, pipe streaming variant

void thread_reader_stdin(runtime *run) {
  // Read alternately into two buffers each large enough to fill
  // the work queue, so that we won't overwrite data before it
  // has been processed.
  static const int bufsize = 1 << 20;
  ci16 *buf1 = new ci16[bufsize];
  ci16 *buf2 = new ci16[bufsize];
  
  while ( true ) {
    size_t nr1 = fread(buf1, sizeof(ci16), bufsize, stdin);
    if ( ! nr1 ) break;
    process_samples(run, buf1, nr1);
    size_t nr2 = fread(buf2, sizeof(ci16), bufsize, stdin);
    if ( ! nr2 ) break;
    process_samples(run, buf2, nr2);
    #if 0
    // Wait before overwriting the input buffer
    for ( int t=0; t<NTHREADS; ++t )
      for ( int i=0; i<wqsize; ++i )
	while ( run->threads[t].jobs[i].go ) { ++nwaitr; yield(); }
    #endif
    poll_control(run);
  }

  fprintf(stderr, "EOF\n");
  fprintf(stderr, "IQ range %d..%d + %d..%dw\n",
	  run->iqmin.re, run->iqmax.re, run->iqmin.im, run->iqmax.im);
}

void *thread_reader(void *arg) {
  runtime *run = (runtime*)arg;

  if ( run->cfg->pmp )
    thread_reader_pmp(run);
  else
    thread_reader_stdin(run);

#if 1
  float dur = run->nsamples / run->cfg->Fs;
  fprintf(stderr, "nsamples=%lld (%f s)\n", (long long)run->nsamples, dur);
  fprintf(stderr, "nwaitr=%lld (%.0f Hz)\n", (long long)nwaitr, nwaitr/dur);
  fprintf(stderr, "nwaitf1=%lld (%.0f Hz)\n", (long long)nwaitf1, nwaitf1/dur);
  fprintf(stderr, "nwaitf2=%lld (%.0f Hz)\n", (long long)nwaitf2, nwaitf2/dur);
  fprintf(stderr, "nwaitj=%lld (%.0f Hz)\n", (long long)nwaitj, nwaitj/dur);
  uint64_t nexecs = run->threads[0].nexecs + run->threads[1].nexecs;
  fprintf(stderr, "nfft=%lld+%lld=%lld (%.0f Hz)\n",
	  (long long)run->threads[0].nexecs,
	  (long long)run->threads[1].nexecs,
	  (long long)nexecs,
	  nexecs/dur);
#endif

  run->eof = true;
  return NULL;
}

// FFT worker threads

void *thread_fft(void *arg) {
  thread_data *td = (thread_data*)arg;
  runtime *run = td->run;
  config *cfg = run->cfg;
  
  while ( true ) {
    thread_data::job *job = &td->jobs[td->qfft];
    while ( ! job->go ) { ++nwaitf1; yield(); }  // buf is not ready
    fftwf_complex *pin = td->in;
#define IQSHIFT 0  // Use this to simulate 1-bit sampling
    // Use leftover samples
    ci16 *buf = job->buf1;
    for ( int n=job->nbuf1; n--; ++buf,++pin ) {
      (*pin)[0] = buf->re >> IQSHIFT;
      (*pin)[1] = buf->im >> IQSHIFT;
    }
    buf = job->buf2;
    for ( int n=cfg->N-job->nbuf1; n--; ++buf,++pin ) {
      (*pin)[0] = buf->re >> IQSHIFT;
      (*pin)[1] = buf->im >> IQSHIFT;
    }
    // Release buf
    job->go = false;
#if 0
      // Hamming
      for ( int i=0; i<N; ++i ) {
	float k = 0.54 - 0.46*cosf(2*M_PI*i/(N-1));
	in[i][0] *= k;
	in[i][1] *= k;
      }
#endif
    fftwf_execute(td->p);
    ++td->nexecs;
    while ( job->done ) { ++nwaitf2; yield(); }  // ph[] is busy
    // Compute phases
    uint16_t *pph = job->ph;
    for ( chan *ch=cfg->chans;
	  ch<cfg->chans+cfg->nchans;
	  ++ch,++pph ) {
      fftwf_complex *p = &td->out[ch->ibin];
      // Apply linear combination of nearest bins
      float d[2] = { 0, 0 };
      for ( int b=0; b<2; ++b )
	for ( int i=0; i<2; ++i )
	  d[i] += ch->bw[b][i][0]*p[b][0] + ch->bw[b][i][1]*p[b][1];
      // atan
#if 0
      int16_t sph = atan2f(d[1], d[0]) * 65536 / (2*M_PI);  // float->int
      uint16_t ph = sph;  //  uint -> int
#else
      while ( d[0]<-126 || d[0]>126 || d[1]<-126 || d[1]>126 ) {
	d[0] *= 0.5;
	d[1] *= 0.5;
      }
      uint8_t ix=(int8_t)d[0], iy=(int8_t)d[1];
      uint16_t ph = run->lut_atan2[iy][ix];
#endif
      *pph = ph;
    }
    job->done = true;
    td->qfft = (td->qfft+1) % wqsize;
  }
}

// WAV utilities

void fwbe32(FILE *f, uint32_t v) {
  fprintf(f, "%c%c%c%c", v>>24, (v>>16)&255, (v>>8)&255, v&255);
}
void fwle32(FILE *f, uint32_t v) {
  fprintf(f, "%c%c%c%c", v&255, (v>>8)&255, (v>>16)&255, (v>>24)&255);
}
void fwle16(FILE *f, uint16_t v) {
  fprintf(f, "%c%c", v&255, (v>>8)&255);
}
void write_wav_header(FILE *f, float Fau) {
  uint32_t nsamples = 0;
  fwbe32(f, 0x52494646);  // tag='RIFF'
  fwle32(f, nsamples ? 36+nsamples : 0); // chunk size
  fwbe32(f, 0x57415645);  // format='WAVE'
  fwbe32(f, 0x666d7420);  // id='fmt '
  fwle32(f, 16);          // fmt chunk size
  fwle16(f, 1);           // codec=PCM
  fwle16(f, 1);           // channels=1
  fwle32(f, Fau);         // sample rate
  fwle32(f, Fau);     // byte rate
  fwle16(f, 1);           // alignment
  fwle16(f, 8);           // bits per sample
  fwbe32(f, 0x64617461);  // id='data'
  fwle32(f, nsamples);    // data chunk size
}

int run(config &cfg) {
  int audiodecim;  // +decimation or -interpolation factor
  if ( ! cfg.Fq ) {
    if ( cfg.Fau ) {
      // Select smallest (sub)multiple of audio rate with enough bandwidth
      if ( cfg.Fau > 2*cfg.maxdev ) {
	audiodecim = -floor(cfg.Fau/(2*cfg.maxdev));
	cfg.Fq = cfg.Fau / (-audiodecim);
      } else {
	audiodecim = ceilf((2*cfg.maxdev)/cfg.Fau);
	cfg.Fq = cfg.Fau * audiodecim;
      }
    } else {
      // Demodulate consecutive blocs
      cfg.Fq = cfg.Fs / cfg.N;
      cfg.Fau = cfg.Fq;
      audiodecim = 1;
    }
  } else {
    if ( ! cfg.Fau ) {
      // Output audio without decimation
      cfg.Fau = cfg.Fq;
      audiodecim = 1;
    } else {
      // Check that audio decimation is integer
      audiodecim = floor(cfg.Fq/cfg.Fau + 0.5);
      if ( fabs(cfg.Fau*audiodecim-cfg.Fq) > 0.5 )
	fatal("Audio decimation ratio (Fq/Fa) must be an integer\n");
    }
  }
  
  if ( cfg.Fq < 2*cfg.maxdev )
    fprintf(stderr, "Warning: Fq too slow for FM deviation\n");

  fprintf(stderr, "IQ sample rate %.3f kHz\n", cfg.Fs/1000);
  fprintf(stderr, "Channel quadrature rate %.3f kHz\n", cfg.Fq/1000);
  fprintf(stderr, "Audio rate %.0f Hz\n", cfg.Fau);
  fprintf(stderr, "FFT cut-off %.0f kHz\n", cfg.Fs/cfg.N/1000);

  fprintf(stderr, "Realtime requires %.0f %d-point FFTs per second\n",
	  cfg.Fq, cfg.N);

  runtime run(&cfg);

  // Setup channels
  for ( struct chan *ch=cfg.chans; ch<cfg.chans+cfg.nchans; ++ch ) {
    float fbin = cfg.N * (ch->F-cfg.Fc) / cfg.Fs;
    int bin = floor(fbin);
    float Frel = fbin - bin;
    if ( Frel < 0.125 ) {  // Round to 0
      // Use low bin
      ch->bw[0][0][0] = 1; ch->bw[0][0][1] = 0; 
      ch->bw[0][1][0] = 0; ch->bw[0][1][1] = 1;
      ch->bw[1][0][0] = 0; ch->bw[1][0][1] = 0; 
      ch->bw[1][1][0] = 0; ch->bw[1][1][1] = 0;
    } else if ( Frel < 0.375 ) {  // Rount to 0.25
      // Rotate -45x3, +135
      ch->bw[0][0][0] =  0.707; ch->bw[0][0][1] =  0.707; 
      ch->bw[0][1][0] = -0.070; ch->bw[0][1][1] =  0.707;
      ch->bw[1][0][0] = -0.2;   ch->bw[1][0][1] = -0.2; 
      ch->bw[1][1][0] =  0.2;   ch->bw[1][1][1] = -0.2;
    } else if ( Frel < 0.625 ) {  // Round to 0.5
      // Rotate -90, +90
      ch->bw[0][0][0] =  0; ch->bw[0][0][1] =  1; 
      ch->bw[0][1][0] = -1; ch->bw[0][1][1] =  0; 
      ch->bw[1][0][0] =  0; ch->bw[1][0][1] = -1; 
      ch->bw[1][1][0] =  1; ch->bw[1][1][1] =  0;
    } else if ( Frel < 0.875 ) {  // Round to 0.75
      // Rotate -135, +45x3
      ch->bw[0][0][0] = -0.2;   ch->bw[0][0][1] =  0.2; 
      ch->bw[0][1][0] = -0.2;   ch->bw[0][1][1] = -0.2;
      ch->bw[1][0][0] =  0.707; ch->bw[1][0][1] = -0.707; 
      ch->bw[1][1][0] =  0.707; ch->bw[1][1][1] =  0.707;
    } else { // Round to 1
      // Use high bin
      ch->bw[0][0][0] = 0; ch->bw[0][0][1] = 0; 
      ch->bw[0][1][0] = 0; ch->bw[0][1][1] = 0;
      ch->bw[1][0][0] = 1; ch->bw[1][0][1] = 0; 
      ch->bw[1][1][0] = 0; ch->bw[1][1][1] = 1;
    }

    // Scale for fast atan2 lookup
    for ( int b=0; b<2; ++b )
      for ( int i=0; i<2; ++i )
	for ( int j=0; j<2; ++j )
	  ch->bw[b][i][j] *= 8.0 * 128 / 2048 / cfg.N;
	    
    ch->ibin = (cfg.N+bin) % cfg.N;
    float derot = 2*M_PI * (ch->F-cfg.Fc) * run.stride/cfg.Fs;
    while ( derot >  M_PI ) derot -= 2*M_PI;
    while ( derot < -M_PI ) derot += 2*M_PI;
    ch->derot = (int16_t) (derot * 65536 / (2*M_PI));
    fprintf(stderr, "  channel fbin=%.2f ibins=%d,%d derot=%f\n",
	    fbin, ch->ibin, ch->ibin+1, derot);
    ch->prevph = 0;
    ch->rms = 1;
  }

  // Spawn FFT threads
  for ( thread_data *td=run.threads; td<run.threads+NTHREADS; ++td ) {
    td->run = &run;
    for ( int i=0; i<wqsize; ++i ) {
      td->jobs[i].ph = new uint16_t[cfg.nchans];
      td->jobs[i].go = false;
      td->jobs[i].done = false;
    }
    td->qread = td->qfft = td->qjoin = 0;
    td->in  = (fftwf_complex*) fftwf_malloc(sizeof(*td->in) *cfg.N);
    td->out = (fftwf_complex*) fftwf_malloc(sizeof(*td->out)*cfg.N);
    td->p = fftwf_plan_dft_1d(cfg.N, td->in, td->out, -1, FFTW_ESTIMATE);
    if ( pthread_create(&td->pth, NULL, thread_fft, td) )
      fatal("pthread_create");
    td->nexecs = 0;
  }

  // Spawn the reader/dispatcher thread
  if ( pthread_create(&run.reader, NULL, thread_reader, &run) )
    fatal("pthread_create");
  
#if 0
  {
    float Ftest = 4.25;
    fprintf(stderr, "Carrier %f between bins\n", Ftest);
    thread_data *td = &run.threads[0];
    for ( int i=0; i<cfg.N; ++i ) {
      td->in[i][0]=cosf(2*M_PI*Ftest*i/cfg.N);
      td->in[i][1]=sinf(2*M_PI*Ftest*i/cfg.N);
    }
    fftwf_execute(td->p);
    for ( int i=0; i<10; ++i )a
      fprintf(stderr, "%3d %f %f\n", i, td->out[i][0], td->out[i][1]);
  }
  exit(0);
#endif

  // Demodulate and output

  fprintf(stderr, "De-emphasis: %.0f us\n", cfg.deemph*1e6);
  float alpha_deemph = 1 / (cfg.Fq*cfg.deemph);

  float t_squelch = 0.1;  // Squelch response time (s)
  float alpha_squelch = 1 / (cfg.Fau*t_squelch);
  
  uint16_t audioclock = 0;

  float discr_gain = cfg.Fq/65536/(2*cfg.maxdev); // Scale from 16-bit angles
  discr_gain *= 0.75;  // Allow some roll-off
  discr_gain *= 256;   // Scale to 8-bit audio out

  float deemph = 0;  // Filter state

  if ( cfg.wav ) { write_wav_header(stdout, cfg.Fau); fflush(stdout); }

  static const int audiobatch = 8192;
  int8_t audiobuf[audiobatch], *audioptr=audiobuf;
		  
  while ( ! run.eof ) {
    for ( thread_data *td=run.threads; td<run.threads+NTHREADS; ++td ) {
      thread_data::job *job = &td->jobs[td->qjoin];
      while ( !job->done && !run.eof ) { ++nwaitj; yield(); }  // Not ready
      if ( run.eof ) break;
      // Compute one audio sample
      {
	float audio = 0;  // Sum of channels
	int nactive = 0;  // Unmuted and not squelched
	uint16_t *pph = job->ph;
	for ( chan *ch=cfg.chans; ch<cfg.chans+cfg.nchans; ++ch,++pph ) {
	  if ( ! ch->enabled ) continue;
	  int16_t dph = *pph - ch->prevph - ch->derot;  // uint modulo -> int
	  ch->prevph = *pph;
	  float dev = dph;
	  if ( cfg.squelch ) {
	    ch->rms = ch->rms*(1-alpha_squelch)
	      + (dev*dev/(32768*32768))*alpha_squelch;
	    if ( ch->rms > 1-cfg.squelch ) continue;
	  }
	  audio += dev;
	  ++nactive;
	}

	deemph = deemph*(1-alpha_deemph) + audio*alpha_deemph;
	audio = deemph;

	if ( audiodecim < 0 ) {
	  // Upsample (repeat)
	  audio *= run.lut_audioscale[nactive];
	  int8_t au = audio * discr_gain;
	  if ( cfg.wav ) au ^= 128;  // Make unsigned
	  for ( int repeat=-audiodecim; --repeat>=0; ) {
	    *audioptr = au;
	    if ( ++audioptr == audiobuf+audiobatch ) {
	      int nw = write(1, audiobuf, sizeof(audiobuf));
	      if ( nw !=  sizeof(audiobuf) ) fatal("write");
	      audioptr = audiobuf;
	    }
	  }
	} else {
	  // Downsample (decimate)
	  if ( ++audioclock == audiodecim ) {
	    audioclock = 0;
	    audio *= run.lut_audioscale[nactive];
	    int8_t au = audio * discr_gain;
	    if ( cfg.wav ) au ^= 128;  // Make unsigned
	    *audioptr = au;
	    if ( ++audioptr == audiobuf+audiobatch ) {
	      int nw = write(1, audiobuf, sizeof(audiobuf));
	      if ( nw !=  sizeof(audiobuf) ) fatal("write");
	      audioptr = audiobuf;
	    }
	  }
	}
      }
      job->done = false;
      td->qjoin = (td->qjoin+1) % wqsize;
    } // threads
  } // main loop

  for ( thread_data *td=run.threads; td<run.threads+NTHREADS; ++td ) {
    fftwf_destroy_plan(td->p);
    //fftwf_free(in);  // double-free bug ?
    //fftwf_free(out);
  }

  fprintf(stderr, "Exiting.\n");

  return 0;
}

// CLI

void usage(const char *name, FILE *f, int c, const char *info=NULL) {
  fprintf(f, "Usage: %s [options] CHANNEL ...  < IQ  > RAWAUDIO\n", name);
  fprintf(f, "Read int16 I/Q from stdin, demodulate multiple FM channels,\n"
	  "write int8 mono audio to stdout.\n");
  fprintf
    (f,
     "\nOptions:\n"
     "  --pmp                Input by reference to /dev/mem\n"
     "  --fs FLOAT           IQ sampling rate (Hz)\n"
     "  --fc FLOAT           Center RF frequency (Hz)\n"
     "  -N INT               FFT size\n"
     "  --fq FLOAT           Quadrature rate (Hz)\n"
     "  --maxdev FLOAT       FM deviation (Hz)\n"
     "  --deemph FLOAT       De-emphasis time constant (s)\n"
     "  --fa FLOAT           Audio sampling rate\n"
     "  --wav                Output WAV header\n"
     "  --fd-info INT        Output aux info to this file descriptor\n"
     "  --info-rate FLOAT    Aux info refresh rate (Hz)\n"
     "  --fd-control INT     Read MUTE/UNMUTE requests from this FD\n"
     "\n"
     "Channel syntax:\n"
     "  FreqMHz          Single channel\n"
     "  Min:Step:Max     Multiple channels with fixed separation\n"
     "  (...)            Same as above, muted initially\n"
     );
  if ( info ) fprintf(f, "Error while processing '%s'\n", info);
  exit(c);
}
  
void add_chan(config *cfg, double fMHz, bool enabled) {
  float f = fMHz * 1e6;
  if ( cfg->nchans == MAXCHANS ) fail("Too many channels");
  struct chan *ch = &cfg->chans[cfg->nchans];
  ch->F = f;
  ch->enabled = enabled;
  ++cfg->nchans;
}

#ifndef VERSION
#define VERSION "undefined"
#endif

int main(int argc, char *argv[]) {
  config cfg;
  
  for ( int i=1; i<argc; ++i ) {
    double fmin, fmax, fstep;
    int nchars;
    if      ( ! strcmp(argv[i],"--fs") && i+1<argc )
      cfg.Fs = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--fc") && i+1<argc )
      cfg.Fc = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--pmp") )
      cfg.pmp = true;
    else if ( ! strcmp(argv[i],"-N") && i+1<argc )
      cfg.N = atoi(argv[++i]);
    else if ( ! strcmp(argv[i],"--maxdev") && i+1<argc )
      cfg.maxdev = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--fq") && i+1<argc )
      cfg.Fq = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--deemph") && i+1<argc )
      cfg.deemph = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--squelch") && i+1<argc )
      cfg.squelch = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--fa") && i+1<argc )
      cfg.Fau = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--wav") )
      cfg.wav = true;
    else if ( ! strcmp(argv[i],"--fd-info") && i+1<argc )
      cfg.fd_info = atoi(argv[++i]);
    else if ( ! strcmp(argv[i],"--info-rate") && i+1<argc )
      cfg.info_rate = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--spectrum-size") && i+1<argc )
      cfg.spec_size = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--spectrum-zoom") && i+1<argc )
      cfg.spec_zoom = atof(argv[++i]);
    else if ( ! strcmp(argv[i],"--fd-control") && i+1<argc )
      cfg.fd_control = atoi(argv[++i]);
    else if ( !strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "--version") ) {
      printf("%s\n", VERSION);
      exit(0);
    }
    else if ( argv[i][0] == '-' )
      fail(argv[i]);
    else if ( sscanf(argv[i], "%lf:%lf:%lf", &fmin,&fstep,&fmax) == 3 ) {
      for ( double f=fmin; f<fmax+fstep/2; f+=fstep )
	add_chan(&cfg, f, true);
    }
    else if ( sscanf(argv[i], "(%lf:%lf:%lf)", &fmin,&fstep,&fmax) == 3 ) {
      for ( double f=fmin; f<fmax+fstep/2; f+=fstep )
	add_chan(&cfg, f, false);
    }
    else if ( sscanf(argv[i], "%lf", &fmin) == 1 ) {
      add_chan(&cfg, fmin, true);
    }
    else if ( sscanf(argv[i], "(%lf)", &fmin) == 1 ) {
      add_chan(&cfg, fmin, false);
    }
    else
      usage(argv[0], stderr, 1, argv[i]);
  }

  if ( ! cfg.nchans ) fprintf(stderr, "Warning: no channel specified\n");

  return run(cfg);
}
