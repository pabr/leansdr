HEAD
  * Added --float-scale..
  * Documented --fd-const.
  * Renamed --probesize (was --maxsend).

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
