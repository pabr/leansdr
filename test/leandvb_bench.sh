#!/bin/sh

BINDIR=../src/apps
RESULTS=leandvb_bench_results.txt

echo "# leansdr_commit platform sampling_ratio RXSNR CNR SS MER VBERMIN VBERMAX"  > $RESULTS

COMMIT="git:$(git rev-parse --short HEAD)"

case $(cat /proc/cpuinfo) in
    *i3-2357M*) PLATFORM="laptop_i3_1300M" ;;
    *i7-4770R*) PLATFORM="desktop_i7_3200M" ;;
    *ARMv7*)    PLATFORM="raspberry_pi_2" ;;
    *)          PLATFORM="other" ;;
esac

run() {
    echo $*  1>&1
    local ratio="$1"
    local symbrate=1000000
    local samprate=$(echo "$symbrate*$ratio" | bc)
    local noise="$2"  # RMS total noise
    local ibnoise=$(echo "$noise/sqrt($ratio)" | bc -l)  # RMS in-band noise
    local rxsnr=$(echo "20 * l(75/$ibnoise)/l(10)" | bc -l)
    rxsnr=$(printf "%.2f" $rxsnr)
    local flags="$3"

    local cnr=
    [ $samprate -gt $((3*symbrate)) ]  &&  cnr="--cnr"
    local ou8=
    [ "$flags" = "--u8 --hs" ]  &&  ou8="--ou8"

    cmd="$BINDIR/leantsgen  -c 2000                                       \
	|  $BINDIR/leandvbtx  -f $ratio                                   \
	|  $BINDIR/leanchansim  --awgn $noise  --deterministic  $ou8      \
	|  $BINDIR/leandvb  --f32  -f $samprate  --sr $symbrate  --anf 0  \
                            $cnr --gui --fd-info 2 $flags  2>&1  > /dev/null"
    echo $cmd  1>&2
    eval $cmd | (
	local cnr ss mer vbermin vbermax 
	local success=0
	while read info arg; do
	    case $info in
		LOCK)
		    if [ $arg = 0 ]; then
			vbermin=1000000
			vbermax=0
			cnr=0
			ss=0
			mer=0
		    fi
		    ;;
		VBER)
		    V=$(echo "$arg*1000000" | bc -l);
		    V=${V%.*}
		    [ $V -lt $vbermin ] && vbermin=$V;
		    [ $V -gt $vbermax ] && vbermax=$V;
		    ;;
		CNR) cnr=$arg ;;
		SS) ss=$arg ;;
		MER) mer=$arg ;;
		LOCKTIME) if [ $arg -ge 1000 ]; then success=1; break; fi ;;
	    esac
	done;
	if [ "$success" = 1 ]; then
	    vbmin=$(echo "$vbermin * 10^-6" | bc -l)
	    vbmax=$(echo "$vbermax * 10^-6" | bc -l)
	    echo "$COMMIT $PLATFORM $(echo $ratio | bc -l) $rxsnr $cnr $ss $mer $vbmin $vbmax"  |  tee -a $RESULTS
	fi
    )
}

index() {
    echo
    echo "== $* =="
    echo
    echo  >> $RESULTS
    echo  >> $RESULTS
    echo "# $*."  >> $RESULTS
}

######################################################################

index "1.2sps-hs"
run 6/5 5   "--u8 --hs"
run 6/5 6   "--u8 --hs"
run 6/5 8   "--u8 --hs"
run 6/5 10  "--u8 --hs"
run 6/5 12  "--u8 --hs"
run 6/5 14  "--u8 --hs"
# run 6/5 16  "--u8 --hs"
# run 6/5 18  "--u8 --hs"
# run 6/5 20  "--u8 --hs"

index "2.4sps-hs"
run 12/5  9  "--u8 --hs"
run 12/5 10  "--u8 --hs"
run 12/5 12  "--u8 --hs"
run 12/5 15  "--u8 --hs"
run 12/5 20  "--u8 --hs"
run 12/5 25  "--u8 --hs"
# run 12/5 30  "--u8 --hs"

index "4.2sps-hs"
run 21/5 12  "--u8 --hs"
run 21/5 15  "--u8 --hs"
run 21/5 18  "--u8 --hs"
run 21/5 20  "--u8 --hs"
run 21/5 22  "--u8 --hs"
run 21/5 25  "--u8 --hs"
# run 12/5 30  "--u8 --hs"

index "1.2sps"
run 6/5 8
run 6/5 9
run 6/5 10
run 6/5 11
run 6/5 12
run 6/5 13
run 6/5 14
 run 6/5 15

if false; then
index "4sps-viterbi-resample"
run 4 60 "--viterbi --resample"
run 4 65 "--viterbi --resample"
run 4 70 "--viterbi --resample"
run 4 75 "--viterbi --resample"
run 4 80 "--viterbi --resample"
fi

index "4.2sps"
run 21/5 15
run 21/5 18
run 21/5 20
run 21/5 22
run 21/5 25
run 21/5 30
# run 21/5 35
# run 21/5 40

if false; then
index "8.2sps"
run 41/5 10
run 41/5 15
run 41/5 20
run 41/5 22
run 41/5 25
run 41/5 30
fi

index "4.2sps-resample"
run 21/5 18 "--resample"
run 21/5 20 "--resample"
run 21/5 22 "--resample"
run 21/5 24 "--resample"
run 21/5 26 "--resample"
run 21/5 28 "--resample"
run 21/5 30 "--resample"
run 21/5 32 "--resample"
run 21/5 34 "--resample"

index "1.2sps-viterbi"
run 6/5 26.0 "--viterbi"
run 6/5 27.0 "--viterbi"
run 6/5 28.0 "--viterbi"


if false; then
index "1.2sps-viterbi-resample"
run 6/5 24.0 "--viterbi --resample"
run 6/5 24.1 "--viterbi --resample"
run 6/5 24.2 "--viterbi --resample"
fi

if false; then
index "2.4sps-viterbi-resample"
run 12/5 46 "--viterbi --resample"
run 12/5 47 "--viterbi --resample"
run 12/5 48 "--viterbi --resample"
run 12/5 49.0 "--viterbi --resample"
run 12/5 49.2 "--viterbi --resample"
run 12/5 49.4 "--viterbi --resample"
run 12/5 49.6 "--viterbi --resample"
run 12/5 49.8 "--viterbi --resample"
run 12/5 50.0 "--viterbi --resample"
fi

index "4.2sps-viterbi-resample"
#run 21/5 50 "--viterbi --resample"
#run 21/5 55 "--viterbi --resample"
run 21/5 60 "--viterbi --resample"
run 21/5 65 "--viterbi --resample"
run 21/5 70 "--viterbi --resample"
run 21/5 75 "--viterbi --resample"
run 21/5 80 "--viterbi --resample"
run 21/5 85 "--viterbi --resample"
run 21/5 90 "--viterbi --resample"

if false; then
index "8sps-viterbi-resample"
run 8  75 "--viterbi --resample"
run 8  80 "--viterbi --resample"
run 8  85 "--viterbi --resample"
run 8  90 "--viterbi --resample"
run 8  95 "--viterbi --resample"
run 8 100 "--viterbi --resample"
run 8 105 "--viterbi --resample"
run 8 110 "--viterbi --resample"
run 8 115 "--viterbi --resample"
fi

######################################################################

gnuplot -e 'load "leandvb_bench.gnuplot"; pause -1'
