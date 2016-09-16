#!/bin/sh

if [ -z "$1" ]; then
  echo "Usage: $0 /path/to/leandvb --fd-info 2 [leandvb options]  < IQ  > TS"  1>&2
  exit 1
fi

echo -e "\nleandvb text UI example" 1>&2

(exec $*  2>&1  1>&4  |
    (   while read tag val; do
	    case "$tag" in
		LOCK) 
		    case "$val" in
			0) lock="[SEARCH]" ;;
			1) lock="[LOCKED]" ;;
		    esac ;;
		MER) mer=$(printf "[MER %4.1f dB]" "$val") ;;
		SS) ss=$(printf "[SS %3.0f]" "$val") ;;
		FREQ) freq=$(printf "[Offset %+7.0f Hz]" "$val") ;;
		CR) cr=$(printf "[FEC %s]" "$val") ;;
		SR) sr=$(printf "[SR %7.0f Hz]" "$val") ;;
		*)  echo -e "\n$tag $val"  1>&2 ;;
	    esac
	    echo -ne "\r$ss $freq $mer $lock $sr $cr"  1>&2
	done)
)  4>&1
echo  1>&2