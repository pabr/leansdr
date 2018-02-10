// This file is part of LeanSDR (c) <pabr@pabr.org>.
// See the toplevel README for more information.

// LeanSDR backend for AD936x via libiio.
// Largely derived from libiio documentation and examples.

#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <iio.h>

void fatal(const char *s) { perror(s); exit(1); }

void iiofatal(int err, const char *info) {
  char msg[256];
  iio_strerror(err, msg, sizeof(msg));
  fprintf(stderr, "** IIO: %s failed with %s\n", info, msg);
  exit(1);
}

void fail(const char *s) {
  fprintf(stderr, "** leaniiotx: %s\n", s);
  exit(1);
}

struct config {
  double Flo;   // LO frequency (Hz)
  float Fm;     // Modulation rate / Analog bandwidth (Hz)
  float Fs;     // Sampling rate (Hz)
  int nbufs;    // 0 for default
  int bufsize;
  bool verbose, debug;
  config() : Flo(2449e6), Fm(0), Fs(2.5e6), bufsize(64*1024),
	     verbose(false), debug(false) { }
};

int stream(config &cfg,
	   struct iio_buffer *txbuf,
	   struct iio_channel *txichan) {
  while ( 1 ) {
    ssize_t nw = iio_buffer_push(txbuf);
    if ( nw <= 0 ) iiofatal(nw, "iio_buffer_push");
    if ( cfg.debug ) fprintf(stderr, ".");
    void *buf = iio_buffer_first(txbuf, txichan);
    // Fill from stdin
    for ( int i=0; i<nw; ) {
      ssize_t nr = read(0, (char*)buf+i, nw-i);
      if ( nr < 0 ) fatal("read");
      if ( ! nr ) return 0;
      i += nr;
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

// Simple x4 interpolator 0.1125 .. 0.125.
// Length must be 16*[1..8].

// Note: The TX FIR filter wants at least 64 coefs, otherwise
// it silently fails to transmit anything (seen on v0.24).

static const char *fir_tx_int4 =
  "TX 3 GAIN 0 INT 4\n"
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n "
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n "
  "-334  \n -335  \n -52   \n 565   \n 1293  \n 1610  \n 938   \n -929  \n "
  "-3436 \n -5221 \n -4569 \n -231  \n 7774  \n 17826 \n 27121 \n 32709 \n "
  "32709 \n 27121 \n 17826 \n 7774  \n 231   \n -4569 \n -5221 \n -3436 \n "
  "-929  \n 938   \n 1610  \n 1293  \n 565   \n -52   \n -335  \n -334  \n "
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n "
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n ";

// Dummy decimator in case we need to configure a compatible
// sampling rate on the RX chain.

static const char *fir_rx_dec4 =
  "RX 3 GAIN 0 DEC 4\n"
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n "
  "0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n 0 \n ";

int run(config &cfg) {

  // AD936x TX internals:
  //   FIR: 128 taps, interpolation 1..4
  //   Half-band filters, interpolation 1..12
  //   DAC: 320 MHz
  //   Analog baseband filter: 200..40000 kHz 
  //   Upconverter: 70..6000 MHz
  //
  // 61444000 to 25000000 Hz: Supported natively.
  // 25000000 to  2083333 Hz: The IIO driver enables HB interpolators.
  //  2083333 to   520833 Hz: We must enable FIR with interpolation.
  //   520833 to    65104 Hz: We must enable interpolation in the FPGA. 

  long long bbFs = cfg.Fs;
  const char *fir = NULL;
  if ( bbFs > 61440000 )
    fail("Requested sample rate is too high");
  else if ( bbFs < 2083333 ) {
    if ( cfg.verbose ) fprintf(stderr, "Using AD936x FIR interpolator x4\n");
    fir = fir_tx_int4;
    if ( bbFs < 520833 )
      fail("Requested sample rate needs FPGA interpolation (not implemented)");
  }
    
  struct iio_context *ctx = iio_create_default_context();
  if ( ! ctx ) fail("iio_create_default_context");
  if ( iio_context_get_devices_count(ctx) < 1 ) fail("No device found");

  static const char *phydevname = "ad9361-phy";
  struct iio_device *phydev = iio_context_find_device(ctx, phydevname);
  if ( ! phydev ) fatal(phydevname);
  
  static const char *phychanname = "voltage0";  
  struct iio_channel *phychan =
    iio_device_find_channel(phydev, phychanname, true);
  if ( ! phychan ) fatal(phychanname);

  leaniio_chattr_write(phychan, "rf_port_select", "A");
  
  if ( ! cfg.Fm ) {
    if ( cfg.verbose ) fprintf(stderr, "No analog bandpass filtering.\n");
    cfg.Fm = cfg.Fs;
  }
  if ( cfg.verbose )
    fprintf(stderr, "Setting RF bandwidth %.0f kHz\n", cfg.Fm/1e3);
  if ( cfg.Fm < 200e3 )
    fprintf(stderr, "Warning: Minimum RF bandwidth is 200 kHz\n");
  leaniio_chattr_write(phychan, "rf_bandwidth", (long long)cfg.Fm);

  // Try to alter only the TX sampling rate.
  // The driver apparently changes RX too anyway.
  
  // Select safe rate regardless of current fir state.
  leaniio_chattr_write(phychan, "out_voltage_sampling_frequency",
		       2500000LL);
  leaniio_devattr_write(phydev, "out_voltage_filter_fir_en", false);
  // leaniio_chattr_write(phychan, "filter_fir_en", false);
  if ( fir ) {
    leaniio_devattr_write(phydev, "filter_fir_config", fir, strlen(fir));
    leaniio_devattr_write(phydev, "out_voltage_filter_fir_en", true);
    // leaniio_chattr_write(phychan, "filter_fir_en", true);
  }

  if ( cfg.verbose )
    fprintf(stderr, "Setting sample rate %.0f kHz\n", cfg.Fs/1e3);
  if ( iio_channel_attr_write_longlong
       (phychan, "out_voltage_sampling_frequency", bbFs) < 0 ) {
    fprintf(stderr, "Failed to set out_voltage_sampling_frequency.\n");
    fprintf(stderr, "Possible cause is no BBPLL for RX+TX rates.\n");
    fprintf(stderr, "Trying again with matching RX decimation.\n");
    leaniio_chattr_write(phychan, "in_voltage_sampling_frequency",
			 2500000LL);
    leaniio_devattr_write(phydev, "in_voltage_filter_fir_en", false);
    if ( fir ) {
      leaniio_devattr_write(phydev, "filter_fir_config",
			    fir_rx_dec4, strlen(fir_rx_dec4));
      leaniio_devattr_write(phydev, "in_voltage_filter_fir_en", true);
    }
    leaniio_chattr_write(phychan, "out_voltage_sampling_frequency", bbFs);
  }
  
  static const char *lochanname = "altvoltage1";  
  struct iio_channel *lochan =
    iio_device_find_channel(phydev, lochanname, true);
  if ( ! lochan ) fail(lochanname);
  if ( cfg.verbose )
    fprintf(stderr, "Upconverting to %.6f MHz\n", cfg.Flo/1e6);
  leaniio_chattr_write(lochan, "frequency", (long long)cfg.Flo);

  static const char *txdevname = "cf-ad9361-dds-core-lpc";
  struct iio_device *txdev = iio_context_find_device(ctx, txdevname);
  if ( ! txdev ) fatal(txdevname);

  struct iio_channel
    *txichan = iio_device_find_channel(txdev, "voltage0", true),
    *txqchan = iio_device_find_channel(txdev, "voltage1", true);
  if ( !txichan || !txqchan ) fatal("TX I/Q");

  iio_channel_enable(txichan);
  iio_channel_enable(txqchan);

  if ( cfg.nbufs ) {
    if ( cfg.verbose )
      fprintf(stderr, "Allocating %d buffers\n", cfg.nbufs);
    if ( iio_device_set_kernel_buffers_count(txdev, cfg.nbufs) )
      fail("iio_device_set_kernel_buffers_count");
  }
  
  if ( cfg.verbose )
    fprintf(stderr, "Allocating %d samples per buffer\n", cfg.bufsize);
  struct iio_buffer *txbuf = iio_device_create_buffer(txdev,cfg.bufsize,false);
  if ( ! txbuf ) fail("iio_device_create_buffer");
  if ( iio_buffer_step(txbuf) != 4 ) fail("Unsupported buffer layout");

  int res = stream(cfg, txbuf, txichan);
  if ( cfg.verbose ) fprintf(stderr, "Shutting down TX channels\n");
  iio_channel_disable(txichan);
  iio_channel_disable(txqchan);
  iio_context_destroy(ctx);
  return res;
}

// CLI

void usage(const char *name, FILE *f, int c, const char *info=NULL) {
  fprintf(f, "Usage: %s [options]  < IQ\n", name);
  fprintf(f, "Read int16 I/Q from stdin, transmit via libiio device.\n");
  fprintf
    (f,
     "\nOptions:\n"
     "  --nbufs INT          Number of buffers\n"
     "  --bufsize INT        Buffer size (samples)\n"
     "  -s FLOAT             Sampling rate (Hz)\n"
     "  --bw FLOAT           Filter bandwidth (Hz)\n"
     "  -f FLOAT             RF center frequency (Hz)\n"
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
    if      ( !strcmp(argv[i], "--nbufs") && i+1<argc )
      cfg.nbufs = strtoul(argv[++i], NULL, 0);
    else if ( !strcmp(argv[i], "--bufsize") && i+1<argc )
      cfg.bufsize = strtoul(argv[++i], NULL, 0);
    else if ( !strcmp(argv[i], "-s") && i+1<argc )
      cfg.Fs = atof(argv[++i]);
    else if ( !strcmp(argv[i], "--bw") && i+1<argc )
      cfg.Fm = atof(argv[++i]);
    else if ( !strcmp(argv[i], "-f") && i+1<argc )
      cfg.Flo = atof(argv[++i]);
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
