#!/bin/sh

BINDIR=../src/apps
RESULTS=leandvb_bench_results.txt

echo "# leansdr_commit platform sampling_ratio RXSNR CNR SS MER VBERMIN VBERMAX"  > $RESULTS

COMMIT="git:$(git rev-parse --short HEAD)"

case $(cat /proc/cpuinfo) in
    *i3-2357M*) PLATFORM="laptop_i3_1300M" ;;
    *i5-3337U*) PLATFORM="laptop_i5_1800M" ;;
    *i7-4770R*) PLATFORM="desktop_i7_3200M" ;;
    *ARMv7*)    PLATFORM="raspberry_pi_2" ;;
    *)          PLATFORM="other" ;;
esac

run() {
    echo $*  1>&1
    local ratio="$1"
    local snrtarget="$2"
    local flags="$3"

    local symbrate=1000000
    local samprate=$(echo "$symbrate*$ratio" | bc)
    local noisepow sigpow float_scale=0
    if [ "$flags" = "--u8 --hs" ]; then
	# In HS mode we expect the receiver gain to be set
	# so that the u8 modulation amplitude matches cstln_amp.
	sigpow="37.5"
	noisepow=$(echo "$sigpow - $snrtarget" | bc -l)
    else
	# Otherwise, set arbitrary fixed noise floor
	# and adjust display scale accordingly.
	sigpow="$snrtarget"
	noisepow=0
	float_scale=$(echo "10*sqrt($ratio)" | bc -l)
    fi
    local rxsnr=$(echo "$sigpow - $noisepow" | bc -l)
    rxsnr=$(printf "%.2f" $rxsnr)

    local cnr=
    [ $samprate -gt $((3*symbrate)) ]  &&  cnr="--cnr"
    local ou8=
    [ "$flags" = "--u8 --hs" ]  &&  ou8="--ou8"

    # How long to try to lock onto the signal.
    local NPACKETS=3000
    # Measurement window. Should contain at least 10/BER bits.
    local MINPACKETS=1000

    cmd="$BINDIR/leantsgen  -c $NPACKETS                                  \
	|  $BINDIR/leandvbtx  -f $ratio  --power $sigpow  --agc      \
	|  $BINDIR/leanchansim  --awgn $noisepow  --deterministic  $ou8      \
	|  $BINDIR/leandvb  --f32  --float-scale $float_scale  -f $samprate  --sr $symbrate  --anf 0  \
                            $cnr  --gui  --fd-info 2  $flags  2>&1  > /dev/null"
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
		LOCKTIME) if [ $arg -ge $MINPACKETS ]; then success=1; break; fi ;;
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

SELECTED_SERIES="$1"

test_series() {
    local name="$1"
    local ratio="$2"
    local cnrvals="$3"
    local options="$4"

    [ -n "$SELECTED_SERIES" -a "$name" != "$SELECTED_SERIES" ] && return;

    index "$name"
    for cnr in $cnrvals; do
	run $ratio $cnr "$options"
    done
}

test_series "1.2sps-hs" 6/5 "20 19 18 17 16 15 14 13 12 11 10" "--u8 --hs"
test_series "2.4sps-hs" 12/5 "20 18 16 14 12 10" "--u8 --hs"
test_series "4.2sps-hs" 21/5 "20 18 16 14 12 10" "--u8 --hs"
test_series "1.2sps" 6/5 "22 21 20 19 18 17 16 15"
test_series "4sps-viterbi-rrc" 4 "6.5 6.0 5.5 5.0 4.5" "--viterbi --sampler rrc"
test_series "4.2sps" 21/5 "20 19 18 17 16 15 14"
test_series "8.2sps" 41/5 "21 20 19 18"
test_series "4.2sps-rrc" 21/5 "16 15 14 13 12 11 10" "--sampler rrc"
test_series "1.2sps-viterbi" 6/5 "12 11 10.5 10 9.5 9 8.5" "--viterbi"
test_series "1.2sps-viterbi-rrc" 6/5 "10 9 8.5 8 7 6 5 4" "--viterbi --sampler rrc"
test_series "2.4sps-viterbi-rrc" 12/5 "8 7 6 5.8 5.6 5.4 5.2 5.0 4.8" "--viterbi --sampler rrc"
test_series "4.2sps-viterbi-rrc" 21/5 "6 5 4.8 4.6 4.5 4.4 4.3 4.2 4.0 3.8" "--viterbi --sampler rrc"
test_series "8sps-viterbi-rrc" 8 "6 5 4.8 4.6 4.5 4.4 4.3 4.2 4.0 3.8" "--viterbi --sampler rrc"
test_series "32sps-viterbi-rrc" 32 "6 5 4.8 4.6 4.5 4.4 4.3 4.2 4.0 3.8" "--viterbi --sampler rrc"

test_series "satmodem4200-60sps" 60 "6 5.2 5 4.8 4.6 4.4 4.2 4.0 3.8" "--viterbi --sampler rrc"

######################################################################

gnuplot -e 'load "leandvb_bench.gnuplot"; pause -1'
