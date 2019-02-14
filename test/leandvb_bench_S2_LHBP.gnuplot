load "S2_common.gnuplot"

set title "leandvb: DVB-S2 error performance with external LDPC decoder"

#set term wxt size 800,480
#set term qt size 800,480

set term gif medium size 1280,768
set output "leandvb_bench_S2_results.gif"

eps=1e-6
plot "leandvb_bench_results.txt" index "S2_LHBP_QPSK_1/4."      using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "QPSK 1/4",    \
     "leandvb_bench_results.txt" index "S2_LHBP_QPSK_1/2."      using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "QPSK 1/2",    \
     "leandvb_bench_results.txt" index "S2_LHBP_QPSK_2/3."      using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "QPSK 2/3",    \
     "leandvb_bench_results.txt" index "S2_LHBP_QPSK_9/10."     using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "QPSK 9/10",   \
     "leandvb_bench_results.txt" index "S2_LHBP_8PSK_2/3."      using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "8PSK 2/3",    \
     "leandvb_bench_results.txt" index "S2_LHBP_8PSK_9/10."     using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "8PSK 9/10",   \
     "leandvb_bench_results.txt" index "S2_LHBP_16APSK_2/3."    using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "16APSK 2/3",  \
     "leandvb_bench_results.txt" index "S2_LHBP_16APSK_9/10."   using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "16APSK 9/10", \
     "leandvb_bench_results.txt" index "S2_LHBP_32APSK_3/4."    using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "32APSK 3/4",  \
     "leandvb_bench_results.txt" index "S2_LHBP_32APSK_9/10."   using 4:($8+$9+2*eps)/2:($8+eps):($9+eps) with yerrorlines title "32APSK 9/10", \
     0.5*erfc(sqrt(10**(x/10))) with lines title "Uncoded QPSK"

set term qt size 2048, 1024
replot
