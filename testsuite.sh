#!/bin/zsh

dataset=/scratch/data/moni/coronavirus.fa
pattern=/tmp/pattern

function genDNA {
	chars=ACGT
	num="$1"
	for i in $(seq 1 "$num") ; do
		echo -n ${chars[1 + $RANDOM % ${#chars[@]} ]}
	done
	echo
}


while [ 1 ]; do
	while ! make; do sleep 10; done


	truncate -s0 $pattern
	for i in $(seq 1 100); do 
		echo "> Seq $i " >> $pattern
		genDNA 100 >> $pattern
	done

	/home/niki/code/pfp/moni/debug/test/src/ms -f "$dataset" -p $pattern
	/home/niki/code/pfp/moni/build/_deps/pfp_thresholds-build/test/src/sdsl_matching_statistics -f "$dataset" -p $pattern
	# /home/niki/code/pfp/moni/build/_deps/pfp_thresholds-build/test/src/matching_statistics -f "$dataset" -p $pattern

	if ! diff -q $pattern.sdsl.lengths $pattern.lengths; then
		echo "SDSL"
		cat $pattern.sdsl.lengths
		echo "PHONI"
		cat $pattern.lengths
		break
	fi
done

