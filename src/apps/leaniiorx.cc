// This file is part of LeanSDR (c) <pabr@pabr.org>.
// See the toplevel README for more information.

// LeanSDR frontend for AD936x via libiio.
// Largely derived from libiio documentation and examples.

#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <iio.h>

#define LEANIIORX_PHYS_MEM_POINTERS 1

void fatal(const char *s) { perror(s); exit(1); }

void iiofatal(int err, const char *info) {
  char msg[256];
  iio_strerror(err, msg, sizeof(msg));
  fprintf(stderr, "** IIO: %s failed with %s\n", info, msg);
  exit(1);
}

void fail(const char *s) {
  fprintf(stderr, "** leaniiorx: %s\n", s);
  exit(1);
}

// Config

struct config {
  double Flo;   // LO frequency (Hz)
  float Fm;     // Modulation rate / Analog bandwidth (Hz)
  float Fs;     // Sampling rate (Hz)
  int nbufs;    // 0 for default
  int bufsize;
  bool pmp;
  bool verbose, debug;
  config() : Flo(2449e6), Fm(0), Fs(2.5e6), bufsize(64*1024),
	     pmp(false), verbose(false), debug(false) { }
};

#if LEANIIORX_PHYS_MEM_POINTERS

// PMP support

static const int PMP_CACHE_SIZE = 64;

static struct pmp_cache_entry {
  uint64_t phys;
  void *virt;
} pmp_cache[64];
static int pmp_cache_used = 0;
char *cma_map = NULL;

uint64_t pmp_find(void *virt, size_t size) {
  pmp_cache_entry *pmp;
  for ( pmp=pmp_cache; pmp<pmp_cache+pmp_cache_used; ++pmp )
    if ( pmp->virt == virt ) return pmp->phys;
  if ( pmp_cache_used == PMP_CACHE_SIZE ) fatal("PMP: Too many buffers");
  // Not found. Scan mem for matching content.
  // Hardcoded for PlutoSDR v0.22. TBD auto-detect instead.
  off_t cmabase = 0x0e400000;
  int cmasize = 256<<20;
  if ( ! cma_map ) {
    fprintf(stderr, "Mapping CMA...\n");
    int fd = open("/dev/mem", O_RDONLY);
    if ( fd < 0 ) fatal("/dev/mem");
    void *map = mmap(0, cmasize, PROT_READ, MAP_SHARED, fd, cmabase);
    if ( map == MAP_FAILED ) fatal("mmap");
    cma_map = (char*)map;
    close(fd);
  }
  // Look for match on first 1 KiB of content
  size_t cmpsize = (size<1024) ? size : 1024;
  // The first buffer is usually at cmabase+0x10000
  int stride = (size<65536) ? size : 65536;
  for ( char *p=cma_map; p+size<cma_map+cmasize; p+=stride )
    if ( ! memcmp(p, virt, cmpsize) ) {
      uint64_t phys = cmabase + (p-cma_map);
      fprintf(stderr, "Found buffer %p[%zd] at 0x%08lx (cma+0x%08lx)\n",
	      virt, size,
	      (unsigned long)phys,
	      (unsigned long)(p-cma_map));
      pmp->virt = virt;
      pmp->phys = phys;
      ++pmp_cache_used;
      return phys;
    }
  fatal("PMP: Buffer not found");
  return 0;  // Avoid compiler warning
}

#endif  // LEANIIORX_PHYS_MEM_POINTERS

// Main loop

int stream(config &cfg,
	   struct iio_buffer *rxbuf,
	   struct iio_channel *rxichan) {
  while ( 1 ) {
    ssize_t nr = iio_buffer_refill(rxbuf);
    if ( nr < 0 ) iiofatal(nr, "iio_buffer_refill");
    if ( cfg.debug ) fprintf(stderr, ".");
    void *buf = iio_buffer_first(rxbuf, rxichan);
    if ( ! cfg.pmp ) {
      // Simply write samples to stdout
      ssize_t nw = write(1, buf, nr);
      if ( nw < 0 ) fatal("write");
      if ( ! nw ) break;
      if ( nw != nr ) fail("partial write");
    } else {
#if LEANIIORX_PHYS_MEM_POINTERS
      // Output pointer to phys mem
      struct {
	uint64_t magic;
	uint64_t physaddr;
	uint64_t size;
	uint64_t canary;
      } pointer;
      pointer.magic = 0x504d5031;  // PMP1
      pointer.physaddr = pmp_find(buf, nr);
      pointer.size = (size_t)nr;
      pointer.canary = *(volatile uint64_t*)buf;
      int nw = write(1, &pointer, sizeof(pointer));
      if ( nw < 0 ) fatal("write(pmp)");
      if ( nw != sizeof(pointer) ) fail("partial write(pmp)");
      // As a substitute for flow control, delay for half the duration.
      uint64_t buffer_time = 1000000 * (nr/4) / cfg.Fs;
      usleep(buffer_time/2);
#else
      fatal("pmp support not compiled");
#endif  // LEANIIORX_PHYS_MEM_POINTERS
    }
  }
  return 0;
}

// Device setup

void leaniio_devattr_write(const struct iio_device *dev, const char *attr,
			   const void *src, size_t len) {
  int res = iio_device_attr_write_raw(dev, attr, src, len);
  if ( res < 0 ) iiofatal(res, attr);
}

void leaniio_devattr_write(const struct iio_device *dev, const char *attr,
			   bool val) {
  int res = iio_device_attr_write_bool(dev, attr, val);
  if ( res < 0 ) iiofatal(res, attr);
}

void leaniio_chattr_write(const struct iio_channel *chn, const char *attr,
			  const char *src) {
  int res = iio_channel_attr_write(chn, attr, src);
  if ( res < 0 ) iiofatal(res, attr);
}

void leaniio_chattr_write(const struct iio_channel *chn, const char *attr,
			  bool val) {
  int res = iio_channel_attr_write_bool(chn, attr, val);
  if ( res < 0 ) iiofatal(res, attr);
}

void leaniio_chattr_write(const struct iio_channel *chn, const char *attr,
			  long long val) {
  int res = iio_channel_attr_write_longlong(chn, attr, val);
  if ( res < 0 ) iiofatal(res, attr);
}

// Simple quarter-band decimator 0.1125 .. 0.125.
// Length must be 16*[1..8].

static const char *fir_rx_dec4 =
  "RX 3 GAIN 0 DEC 4\n"
  "-334  \n -335  \n -52   \n 565   \n 1293  \n 1610  \n 938   \n -929  \n "
  "-3436 \n -5221 \n -4569 \n -231  \n 7774  \n 17826 \n 27121 \n 32709 \n "
  "32709 \n 27121 \n 17826 \n 7774  \n -231  \n -4569 \n -5221 \n -3436 \n "
  "-929  \n 938   \n 1610  \n 1293  \n 565   \n -52   \n -335  \n -334  \n ";

// Dummy interpolator in case we need to configure a compatible
// sampling rate on the TX chain.

static const char *fir_tx_int4 =
  "TX 3 GAIN 0 INT 4\n"
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n "
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n ";

int run(config &cfg) {

  // AD936x RX internals:
  //   Downconverter: 70..6000 MHz
  //   Analog baseband filter: 200..56000 kHz 
  //   ADC: 25..640 MHz
  //   Half-band filters, decimation 1..12
  //   FIR: 128 taps, decimation 1..4
  //
  // 61444000 to 25000000 Hz: Supported natively by the ADC.
  // 25000000 to  2083333 Hz: The IIO driver enables HB filters automatically.
  //  2083333 to   520833 Hz: We must enable FIR with decimation.
  //   520833 to    65104 Hz: We must enable decimation in the FPGA. 

  long long bbFs = cfg.Fs;
  const char *fir = NULL;
  if ( bbFs > 61440000 )
    fail("Requested sampling rate is too high");
  else if ( bbFs < 2083333 ) {
    if ( cfg.verbose ) fprintf(stderr, "Using AD936x FIR decimator /4\n");
    fir = fir_rx_dec4;
    if ( bbFs < 520833 )
      fail("Requested sampling rate needs FPGA decimation (not implemented)");
  }
    
  struct iio_context *ctx = iio_create_default_context();
  if ( ! ctx ) fatal("iio_create_default_context");
  if ( iio_context_get_devices_count(ctx) < 1 ) fail("No device found");

  static const char *phydevname = "ad9361-phy";
  struct iio_device *phydev = iio_context_find_device(ctx, phydevname);
  if ( ! phydev ) fatal(phydevname);
      
  static const char *phychanname = "voltage0";  
  struct iio_channel *phychan =
    iio_device_find_channel(phydev, phychanname, false);
  if ( ! phychan ) fatal(phychanname);

  leaniio_chattr_write(phychan, "rf_port_select", "A_BALANCED");
  
  if ( ! cfg.Fm ) {
    if ( cfg.verbose ) fprintf(stderr, "No analog bandpass filtering.\n");
    cfg.Fm = cfg.Fs;
  }
  if ( cfg.verbose )
    fprintf(stderr, "Setting RF bandwidth %.0f kHz\n", cfg.Fm/1e3);
  if ( cfg.Fm < 200e3 )
    fprintf(stderr, "Warning: Minimum RF bandwidth is 200 kHz\n");
  leaniio_chattr_write(phychan, "rf_bandwidth", (long long)cfg.Fm);

  // Select a safe rate regardless of current fir state.
  leaniio_chattr_write(phychan, "in_voltage_sampling_frequency",
		       2500000LL);
  leaniio_devattr_write(phydev, "in_voltage_filter_fir_en", false);
  // leaniio_chattr_write(phychan, "filter_fir_en", false);
  if ( fir ) {
    leaniio_devattr_write(phydev, "filter_fir_config", fir, strlen(fir));
    leaniio_devattr_write(phydev, "in_voltage_filter_fir_en", true);
    // leaniio_chattr_write(phychan, "filter_fir_en", true);
  }

  if ( cfg.verbose )
    fprintf(stderr, "Setting sampling rate %.0f kHz\n", bbFs/1e3);
  if ( iio_channel_attr_write_longlong
       (phychan, "in_voltage_sampling_frequency", bbFs) < 0 ) {
    fprintf(stderr, "Failed to set in_voltage_sampling_frequency.\n");
    fprintf(stderr, "Possible cause is no BBPLL for RX+TX rates.\n");
    fprintf(stderr, "Trying again with matching TX interpolation.\n");
    leaniio_chattr_write(phychan, "out_voltage_sampling_frequency",
			 2500000LL);
    leaniio_devattr_write(phydev, "out_voltage_filter_fir_en", false);
    if ( fir ) {
      leaniio_devattr_write(phydev, "filter_fir_config",
			    fir_tx_int4, strlen(fir_tx_int4));
      leaniio_devattr_write(phydev, "out_voltage_filter_fir_en", true);
    }
    leaniio_chattr_write(phychan, "in_voltage_sampling_frequency", bbFs);
  }
  
  static const char *lochanname = "altvoltage0";  
  struct iio_channel *lochan =
    iio_device_find_channel(phydev, lochanname, true);
  if ( ! lochan ) fail(lochanname);
  if ( cfg.verbose )
    fprintf(stderr, "Tuning to %.6f MHz\n", cfg.Flo/1e6);
  leaniio_chattr_write(lochan, "frequency", (long long)cfg.Flo);
  
  static const char *rxdevname = "cf-ad9361-lpc";
  struct iio_device *rxdev = iio_context_find_device(ctx, rxdevname);
  if ( ! rxdev ) fatal(rxdevname);
  
  struct iio_channel
    *rxichan = iio_device_find_channel(rxdev, "voltage0", false),
    *rxqchan = iio_device_find_channel(rxdev, "voltage1", false);
  if ( !rxichan || !rxqchan ) fatal("RX I/Q");

  iio_channel_enable(rxichan);
  iio_channel_enable(rxqchan);

  if ( cfg.nbufs ) {
    if ( cfg.verbose )
      fprintf(stderr, "Allocating %d buffers\n", cfg.nbufs);
    int res = iio_device_set_kernel_buffers_count(rxdev, cfg.nbufs);
    if ( res < 0 ) iiofatal(res, "iio_device_set_kernel_buffers_count");
  }
  
  if ( cfg.verbose )
    fprintf(stderr, "Allocating %d samples per buffer\n", cfg.bufsize);
  struct iio_buffer *rxbuf = iio_device_create_buffer(rxdev,cfg.bufsize,false);
  if ( ! rxbuf ) fail("iio_device_create_buffer");
  if ( iio_buffer_step(rxbuf) != 4 ) fail("Unsupported buffer layout");

  return stream(cfg, rxbuf, rxichan);
}

// CLI

void usage(const char *name, FILE *f, int c, const char *info=NULL) {
  fprintf(f, "Usage: %s [options]  > IQ\n", name);
  fprintf(f, "Capture from libiio device, write int16 I/Q to stdout.\n");
  fprintf(f,
	  "\nOptions:\n"
	  "  -f FLOAT             RF center frequency (Hz)\n"
	  "  --bw FLOAT           Filter bandwidth (Hz)\n"
	  "  -s FLOAT             Sampling rate (Hz)\n"
	  "  --nbufs INT          Number of buffers\n"
	  "  --bufsize INT        Buffer size (samples)\n"
	  "  --pmp                Output by reference to /dev/mem\n"
	  "  -v                   Verbose at startup and exit\n"
	  "  -d                   Verbose while running\n"
	  "  --version            Display version and exit\n"
	  );
  if ( info ) fprintf(f, "Error while processing '%s'\n", info);
  exit(c);
}

#ifndef VERSION
#define VERSION "undefined"
#endif

int main(int argc, char *argv[]) {
  config cfg;

  for ( int i=1; i<argc; ++i ) {
    if      ( !strcmp(argv[i], "-f") && i+1<argc )
      cfg.Flo = atof(argv[++i]);
    else if ( !strcmp(argv[i], "-s") && i+1<argc )
      cfg.Fs = atof(argv[++i]);
    else if ( !strcmp(argv[i], "--bw") && i+1<argc )
      cfg.Fm = atof(argv[++i]);
    else if ( !strcmp(argv[i], "--nbufs") && i+1<argc )
      cfg.nbufs = strtoul(argv[++i], NULL, 0);
    else if ( !strcmp(argv[i], "--bufsize") && i+1<argc )
      cfg.bufsize = strtoul(argv[++i], NULL, 0);
    else if ( !strcmp(argv[i], "--pmp") )
      cfg.pmp = true;
    else if ( !strcmp(argv[i], "-v") )
      cfg.verbose = true;
    else if ( !strcmp(argv[i], "-d") )
      cfg.debug = true;
    else if ( !strcmp(argv[i], "-h") )
      usage(argv[0], stdout, 0);
    else if ( ! strcmp(argv[i], "--version") ) {
      printf("%s\n", VERSION);
      exit(0);
    }
    else
      usage(argv[0], stderr, 1, argv[i]);
  }

  return run(cfg);
}
