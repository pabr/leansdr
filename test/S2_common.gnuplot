set xrange [-3:20]
set xtics 1
set xlabel "Es/N0 (dB)"

set logscale y
set yrange [0.8e-6:1]
set ylabel "BER before BCH decoding"
set format y "%.1tE%T"

set grid

# Error performance targets = BCH t parameter vs LDPC frame size

set label " QPSK 1/4 target"    at -2.35,750e-6 rotate by 90 point pt 5
set label " QPSK 1/2 target"    at  1.00,375e-6 rotate by 90 point pt 5
set label " QPSK 9/10 target"   at  6.42,137e-6 rotate by 90 point pt 5
set label " 8PSK 2/3 target"    at  6.62,232e-6 rotate by 90 point pt 5
set label " 8PSK 9/10 target"   at 10.98,137e-6 rotate by 90 point pt 5
#set label " 16APSK 2/3 target"  at  8.97,232e-6 rotate by 90 point pt 5
#set label " 16APSK 9/10 target" at 13.13,137e-6 rotate by 90 point pt 5
#set label " 32APSK 3/4 target"  at 12.73,248e-6 rotate by 90 point pt 5
#set label " 32APSK 9/10 target" at 16.05,137e-6 rotate by 90 point pt 5
