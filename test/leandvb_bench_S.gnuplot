set title "leandvb: DVB-S error performance vs samples-per-symbol (sps) and demodulation options"

set xrange [2:25]
set xtics 1
set xlabel "Eb/N0 (dB)"

set logscale y
set yrange [1e-6:1e-1]
set ylabel "VBER (BER before RS decoding)"
set format y "%.1tE%T"

set grid

#set term wxt size 800,480
#set term qt size 800,480

set term gif medium size 1280,768
set output "leandvb_bench_results.gif"

set object rectangle from 0,1e-7 to 4.5,2e-4  behind  fillcolor rgbcolor "black"  fillstyle solid 0.1 noborder
set label "DVB-S\nerror performance\nrequirements\n(vber 2e-4 at 4.5 dB)" at 2.8,1.1e-5 center rotate by 90

plot "leandvb_bench_results.txt" index "1.2sps-hs."                using 4:($8+$9)/2:8:9 with yerrorlines title "1.2sps-hs (HamTV on Raspberry Pi 2)"  , \
     "leandvb_bench_results.txt" index "2.4sps-hs."                using 4:($8+$9)/2:8:9 with yerrorlines title "2.4sps-hs"                            , \
     "leandvb_bench_results.txt" index "4.2sps-hs."                using 4:($8+$9)/2:8:9 with yerrorlines title "4.2sps-hs"                            , \
     "leandvb_bench_results.txt" index "1.2sps."                   using 4:($8+$9)/2:8:9 with yerrorlines title "1.2sps (HamTV on tablet)"             , \
     "leandvb_bench_results.txt" index "4.2sps."                   using 4:($8+$9)/2:8:9 with yerrorlines title "4.2sps"                               , \
     "leandvb_bench_results.txt" index "4.2sps-rrc."               using 4:($8+$9)/2:8:9 with yerrorlines title "4.2sps-rrc"                      , \
     "leandvb_bench_results.txt" index "1.2sps-viterbi."           using 4:($8+$9)/2:8:9 with yerrorlines title "1.2sps-viterbi (HamTV on desktop)"    , \
     0.5*erfc(sqrt(10**(x/10)))          with lines title "Uncoded QPSK"                                                   , \
     "leandvb_bench_results.txt" index "1.2sps-viterbi-rrc."       using 4:($8+$9)/2:8:9 with yerrorlines title "1.2sps-viterbi-rrc"              , \
     "leandvb_bench_results.txt" index "2.4sps-viterbi-rrc."       using 4:($8+$9)/2:8:9 with yerrorlines title "2.4sps-viterbi-rrc"              , \
     "leandvb_bench_results.txt" index "4sps-viterbi-rrc."         using 4:($8+$9)/2:8:9 with yerrorlines title "4sps-viterbi-rrc"                , \
     "leandvb_bench_results.txt" index "4.2sps-viterbi-rrc."       using 4:($8+$9)/2:8:9 with yerrorlines title "4.2sps-viterbi-rrc"              , \
     "leandvb_bench_results.txt" index "8.2sps."                   using 4:($8+$9)/2:8:9 with yerrorlines title "8.2sps"                               , \
     "leandvb_bench_results.txt" index "8sps-viterbi-rrc."         using 4:($8+$9)/2:8:9 with yerrorlines title "8sps-viterbi-rrc"                , \
     "leandvb_bench_results.txt" index "32sps-viterbi-rrc."        using 4:($8+$9)/2:8:9 with yerrorlines title "32sps-viterbi-rrc"                , \
     "leandvb_bench_results.txt" index "test."                     using 4:($8+$9)/2:8:9 with yerrorlines title "test"                                 , \
     "leandvb_bench_results.txt" index "satmodem4200-60sps."       using 4:($8+$9)/2:8:9 with yerrorlines title "satmodem4200-60sps"                   , \
     2e-4 with lines title "Quasi-error-free threshold"

set term qt size 2048, 1024
replot