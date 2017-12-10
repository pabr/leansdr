HEAD
  * leandvb: Spectrum output, JSON syntax.
  * leandvb, leandvbtx: DVB-S with non-standard constellations (with --viterbi).
  * leandvb: Support all DVB-S code rates with --viterbi.
  * leandvbtx: Support all DVB-S code rates.
  * New viterbi_sync with simplified metrics (may degrade QPSK sensitivity).

2017-07-23 v1.2.0
  * Fixed --hs mode on x86_64.
  * Added --s8, --s16, --u16.
  * Improved timing recovery at low SNR (small performance penalty).
  * Added --hdlc mode (compatible with some satellite modems).
  * Added benchmarking/test suite (test/leandvb_bench.sh)
  * Added tools for benchmarking: leantsgen, leandvbtx, leanchansim.
  * Added LOCKTIME output with --fd-info.
  * Added --hs mode (high throughput for raw u8 QPSK 1/2).
  * Added --drift.
  * Filters follow --drift.
  * Added --float-scale.
  * Documented --fd-const.
  * Renamed --probesize (was --maxsend).
  * Added --viterbi.
  * Added VBER (bit error rate after deconvolution, detected by RS).
  * Added FFT-based CNR estimator.

2016-09-16 v1.1.0
  * Support for DVB-S2 constellations (not FEC).
  * Fixed --derotate.
  * Fixed excess deconvolution errors with FEC7/8.
  * Added simple baseband filter (not RRC).
  * Added --fastlock mode for low SR and off-line processing.
  * Added leandvb_vt100ui.sh with constellation plot.

2016-09-04 v1.0.0
  * Development moved to git and github.
  * leandvb is now distributed as part of leansdr.
  * Support for all DVB-S code rates.
  * Status output for third-party UIs: lock, MER, frequency offset.
  * Software AGC is always enabled. rtl-sdr HW AGC is not recommended.
  * stderr is quiet by default. Use -v -d for troubleshooting.
  * Deconvolution is now algebraic (instead of look-up table).
  * Added leansdrscan for cycling through DVB settings.
  * Added leansdrcat for debugging real-time behaviour.

2016-02-29 Preview release of leandvb
  * Code rate 1/2 only.
  * Hard-decision deconvolution (without Viterbi).
