#!/bin/sh

BINDIR=../src/apps
RESULTS=leandvb_bench_results.txt

echo "# leansdr_commit platform sampling_ratio TXSNR CNR SS MER VBERMIN VBERMAX"  > $RESULTS

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
    local txflags="$3"
    local rxflags="$4"

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
	float_scale=$(echo "sqrt($ratio)/e(l(10)*$snrtarget/20)" | bc -l)
    fi
    local rxsnr=$(echo "$sigpow - $noisepow" | bc -l)
    rxsnr=$(printf "%.2f" $rxsnr)

    local cnr=
    [ $samprate -gt $((3*symbrate)) ]  &&  cnr="--cnr"
    local ou8=
    [ "$flags" = "--u8 --hs" ]  &&  ou8="--ou8"

    # How long to try to lock onto the signal.
    local NPACKETS=3000   # DVB-S
    local NPACKETS=20000  # DVB-S2 32APSK
    # Measurement window. Should contain at least 10/BER bits.
    local MINPACKETS=1000

    cmd="$BINDIR/leantsgen  -c $NPACKETS                                  \
	|  $BINDIR/leandvbtx  -f $ratio  --power $sigpow  --agc  $txflags     \
	|  $BINDIR/leanchansim  --awgn $noisepow  --deterministic  $ou8      \
	|  $BINDIR/leandvb  --f32  --float-scale $float_scale  -f $samprate  --sr $symbrate  --anf 0  \
                            $cnr  --gui  --fd-info 2  $rxflags  2>&1  > /dev/null"
    echo $cmd  1>&2
    eval $cmd | (
	local cnr ss mer vbermin vbermax 
	local success=0
	while read info arg; do
	    case $info in
		LOCK|FRAMELOCK)
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

begin_datablock() {
    echo
    echo "== $* =="
    echo
    echo  >> $RESULTS
    echo  >> $RESULTS
    echo "# $*."  >> $RESULTS
}

end_datablock() {
    # gnuplot falls through empty data blocks, ignoring names.
    # Append an invalid entry to ensure the block is never empty.
    echo "0"  >> $RESULTS
}

######################################################################

SELECTED_SERIES="$1"

test_series() {
    local name="$1"
    local ratio="$2"
    local cnrvals="$3"
    local txoptions="$4"
    local rxoptions="$5"

    if [ -n "$SELECTED_SERIES" ]; then
	case "$name" in
	    $SELECTED_SERIES*) ;;
	    *) return;
	esac
    fi

    begin_datablock "$name"
    for cnr in $cnrvals; do
	run $ratio $cnr "$txoptions" "$rxoptions"
    done
    end_datablock
}

# DVB-S functional test (huge SNR margin, 4 samples/symbol)

STX="--standard DVB-S"
SRX="--standard DVB-S --sampler rrc"

s_ft() {
    local args="$1"
    local name="$2"
    local snr=$(echo "$3 + 10" | bc -l)
    test_series "S_FT_$name" 4 $snr "$STX $args" "$SRX $args"
}

s_ft  "--cr 1/2" "QPSK_1/2"     4.5
s_ft  "--cr 2/3" "QPSK_2/3"     5.0
s_ft  "--cr 3/4" "QPSK_3/4"     5.5
s_ft  "--cr 5/6" "QPSK_5/6"     6.0
s_ft  "--cr 7/8" "QPSK_7/8"     6.4

# DVB-S2 functional test (huge SNR margin, 4 samples/symbol)

S2TX="--standard DVB-S2 --pilots"
S2RX="--standard DVB-S2 --sampler rrc"

s2_ft() {
    local modcod="$1"
    local name="$2"
    local snr=$(echo "$3 + 20" | bc -l)
    test_series "S2_FT_NF_$name" 4 $snr "$S2TX --modcod $modcod " "$S2RX"
    case "$name" in
	*9/10) return ;;
    esac
    test_series "S2_FT_SF_$name" 4 $snr "$S2TX --modcod $modcod --shortframes" "$S2RX"
}

s2_ft  1 "QPSK_1/4"     -2.35
s2_ft  2 "QPSK_1/3"     -1.24
s2_ft  3 "QPSK_2/5"     -0.30
s2_ft  4 "QPSK_1/2"      1.00
s2_ft  5 "QPSK_3/5"      2.23
s2_ft  6 "QPSK_2/3"      3.10
s2_ft  7 "QPSK_3/4"      4.03
s2_ft  8 "QPSK_4/5"      4.68
s2_ft  9 "QPSK_5/6"      5.18
s2_ft 10 "QPSK_8/9"      6.20
s2_ft 11 "QPSK_9/10"     6.42
s2_ft 12 "8PSK_3/5"      5.50
s2_ft 13 "8PSK_2/3"      6.62
s2_ft 14 "8PSK_3/4"      7.91
s2_ft 15 "8PSK_5/6"      9.35
s2_ft 16 "8PSK_8/9"     10.69
s2_ft 17 "8PSK_9/10"    10.98
s2_ft 18 "16APSK_2/3"    8.97
s2_ft 19 "16APSK_3/4"   10.21
s2_ft 20 "16APSK_4/5"   11.03
s2_ft 21 "16APSK_5/6"   11.61
s2_ft 22 "16APSK_8/9"   12.89
s2_ft 23 "16APSK_9/10"  13.13
s2_ft 24 "32APSK_3/4"   12.73
s2_ft 25 "32APSK_4/5"   13.64
s2_ft 26 "32APSK_5/6"   14.28
s2_ft 27 "32APSK_8/9"   15.69
s2_ft 28 "32APSK_9/10"  16.05

# DVB-S error performance

test_series "1.2sps-hs" 6/5 "20 19 18 17 16 15 14 13 12 11 10" "" "--u8 --hs"
test_series "2.4sps-hs" 12/5 "20 18 16 14 12 10" "" "--u8 --hs"
test_series "4.2sps-hs" 21/5 "20 18 16 14 12 10" "" "--u8 --hs"
test_series "1.2sps" 6/5 "22 21 20 19 18 17 16 15" "" ""
test_series "4sps-viterbi-rrc" 4 "6.5 6.0 5.5 5.0 4.5" "" "--viterbi --sampler rrc"
test_series "4.2sps" 21/5 "20 19 18 17 16 15 14" "" ""
test_series "8.2sps" 41/5 "21 20 19 18" "" ""
test_series "4.2sps-rrc" 21/5 "16 15 14 13 12 11 10" "" "--sampler rrc"
test_series "1.2sps-viterbi" 6/5 "12 11 10.5 10 9.5 9 8.5" "" "--viterbi"
test_series "1.2sps-viterbi-rrc" 6/5 "10 9 8.5 8 7 6 5 4" "" "--viterbi --sampler rrc"
test_series "2.4sps-viterbi-rrc" 12/5 "8 7 6 5.8 5.6 5.4 5.2 5.0 4.8" "" "--viterbi --sampler rrc"
test_series "4.2sps-viterbi-rrc" 21/5 "6 5 4.8 4.6 4.5 4.4 4.3 4.2 4.0 3.8" "" "--viterbi --sampler rrc"
test_series "8sps-viterbi-rrc" 8 "6 5 4.8 4.6 4.5 4.4 4.3 4.2 4.0 3.8" "" "--viterbi --sampler rrc"
test_series "32sps-viterbi-rrc" 32 "6 5 4.8 4.6 4.5 4.4 4.3 4.2 4.0 3.8" "" "--viterbi --sampler rrc"

test_series "satmodem4200-60sps" 60 "6 5.2 5 4.8 4.6 4.4 4.2 4.0 3.8" "" "--viterbi --sampler rrc"

# DVB-S2 error performance with built-in LDPC bit-flipping

BF=100
S2LDPC="--ldpc-bf $BF"

test_series "S2_BF_QPSK_1/4"    2 "10.5 10 9.99 9.98 9.97 9.9 9 8.8 8.6 8.57 8.54 8.5"        "$S2TX --modcod 1 " "$S2RX $S2LDPC"
test_series "S2_BF_QPSK_1/2"    2 "10 9.95 9.94 9.93 9.92 9.91 9.9 " "$S2TX --modcod 4 " "$S2RX $S2LDPC"
test_series "S2_BF_QPSK_9/10"   2 "12.5 11 10.5 10.3 10.24 10.23 10.22 10.21 10.2"       "$S2TX --modcod 11" "$S2RX $S2LDPC"
test_series "S2_BF_8PSK_2/3"    2 "16 15.7 15.5 15.4 15.37 15.34 15.3"      "$S2TX --modcod 13" "$S2RX $S2LDPC"
test_series "S2_BF_8PSK_9/10"   2 "17 16 15.7  15.64 15.635 15.63 15.6 "      "$S2TX --modcod 17" "$S2RX $S2LDPC"
test_series "S2_BF_16APSK_2/3"  2 "20 19 18.5 18.495 18.49"    "$S2TX --modcod 18" "$S2RX $S2LDPC"
test_series "S2_BF_16APSK_9/10" 2 "20 19 18.5 18.4 18.35 18.3 18.2" "$S2TX --modcod 23" "$S2RX $S2LDPC"
test_series "S2_BF_32APSK_3/4"  2 "30 28 27 26.7 26.69 26.68 26.67 26.5"      "$S2TX --modcod 24" "$S2RX $S2LDPC"
test_series "S2_BF_32APSK_9/10" 2 "30 29 28 27 26 25.5 25"      "$S2TX --modcod 28" "$S2RX $S2LDPC"


# DVB-S2 error performance with external LDPC decoder

S2LDPC="--ldpc-helper ./ldpc_tool"

test_series "S2_LHBP_QPSK_1/4"    2 "7 6 5 4.75"        "$S2TX --modcod 1 " "$S2RX $S2LDPC"
test_series "S2_LHBP_QPSK_1/2"    2 "7 6 5.5 5 4.75 4.5 4.25" "$S2TX --modcod 4 " "$S2RX $S2LDPC"
test_series "S2_LHBP_QPSK_2/3"    2 "7 6 5.5 5.25 5" "$S2TX --modcod 6 " "$S2RX $S2LDPC"
test_series "S2_LHBP_QPSK_9/10"   2 "9 8 7.75 7.5"       "$S2TX --modcod 11" "$S2RX $S2LDPC"
test_series "S2_LHBP_8PSK_2/3"    2 "12 11 10.9 10.8 10 9 8"      "$S2TX --modcod 13" "$S2RX $S2LDPC"
test_series "S2_LHBP_8PSK_9/10"   2 "13 12.9 12.8 12 11.5 11"      "$S2TX --modcod 17" "$S2RX $S2LDPC"
test_series "S2_LHBP_16APSK_2/3"  2 "14 13.75 13.5"    "$S2TX --modcod 18" "$S2RX $S2LDPC"
test_series "S2_LHBP_16APSK_9/10" 2 "15 14.5 14"    "$S2TX --modcod 23" "$S2RX $S2LDPC"
test_series "S2_LHBP_32APSK_3/4"  2 "16 15 14"      "$S2TX --modcod 24" "$S2RX $S2LDPC"
test_series "S2_LHBP_32APSK_9/10" 2 "18 17"      "$S2TX --modcod 28" "$S2RX $S2LDPC"

######################################################################

echo "Results written to $RESULTS."  1>&2
