#!/bin/sh

if [ -z "$1" ]; then
  echo "Usage: $0 /path/to/leandvb --fd-info 2 --fd-const 2 [leandvb options]  < IQ  > TS"  1>&2
  exit 1
fi

printat() {
    row=$1; col=$2; fmt=$3; shift 3
    printf "\033[${row};${col}f$fmt" $*  1>&2
}

print_symbols() {
    char="$1"
    set $2
    shift
    for iq in $*; do
	echo $iq  |  (
	    IFS="," read i q
	    [ $i -lt -128 ] && i=-128
	    [ $i -gt  127 ] && i=127
	    [ $q -lt -128 ] && q=-128
	    [ $q -gt  127 ] && q=127
	    printat $((11-$q/12)) $((50+$i/6)) "$char"
	    r=$((11-$q/12))
	    c=$((50+$i/6))
#	    printat $r $c "$char"
	)
    done
}

constel_decim=50

constel_count=0
print_constel() {
    if [ $constel_count -lt $constel_decim ]; then
	constel_count=$((constel_count+1))
	return
    fi
    constel_count=0

    y=1
    while [ $y -le 21 ]; do
	#printat $y 50 "|"
	printat $y 29 "                     |                     "
	y=$((y+1))
    done
    printat 11 29 "---------------------+---------------------"
    print_symbols "x" "$2"
    print_symbols "#" "$1"
}

printf "\033[2J\033[?25l"  1>&2

printat 2 3 "leandvb VT100 UI example"

(exec $*  2>&1  1>&4  |
    (   while read tag val; do
	    case "$tag" in
		SS)   printat 4 3 "SS: %3.0f" "$val" ;;
		MER)  printat 6 3 "MER: %4.1f dB" "$val" ;;
		FREQ) printat 8 3 "Offset: %+7.0f Hz" "$val" ;;
		LOCK) 
		    case "$val" in
			0) printat 10 3 "SYNCHRONIZING..." ;;
			1) printat 10 3 "LOCKED          " ;;
		    esac ;;
		CR) printat 12 3 "FEC: %s" "$val" ;;
		SR) printat 14 3 "SR: %7.0f Hz" "$val" ;;
		CONST) constel="$val" ;;
		SYMBOLS) print_constel "$constel" "$val" ;;
		*) printat 20 1 "$tag $val"
	    esac
	done)
)  4>&1

printat 30 1 "\033[?25h"
