#!/bin/zsh

rdataset=/scratch/data/moni/rrepaired/coronavirus.fa
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

	set -e
	# if [[ ! -e $rdataset.bwt ]]; then
		/home/niki/code/pfp/moni/debug/moni -f $rdataset
		rm $rdataset.R $rdataset.C $rdataset.slp
		[[ ! -f $rdataset.plain ]] && ~/extractfasta.py $rdataset > $rdataset.plain
		cd /home/niki/code/pfp/rrepair/build/ && /home/niki/code/pfp/rrepair/build/rrepair $rdataset.plain 10 100 10 100000000
		#/home/niki/code/pfp/moni/build/_deps/bigrepair-src/bigrepair $rdataset.plain
		/home/niki/code/pfp/moni/build/_deps/shaped_slp-build/SlpEncBuild -e SelfShapedSlp_SdSd_Sd -f rrepair -i $rdataset.plain -o $rdataset.slp
	# fi
	

	while ! make; do sleep 10; done


	truncate -s0 $pattern
	for i in $(seq 1 100); do 
		echo "> Seq $i " >> $pattern
		genDNA 100 >> $pattern
	done

	/home/niki/code/pfp/moni/debug/test/src/ms -f "$dataset" -p $pattern
	/home/niki/code/pfp/moni/build/_deps/pfp_thresholds-build/test/src/sdsl_matching_statistics -f "$dataset" -p $pattern
	# /home/niki/code/pfp/moni/build/_deps/pfp_thresholds-build/test/src/matching_statistics -f "$dataset" -p $pattern
	
	cp $pattern ${pattern}2
	/home/niki/code/pfp/moni/debug/test/src/ms -f "$rdataset" -p ${pattern}2

	if ! diff -q $pattern.sdsl.lengths $pattern.lengths; then
		echo "SDSL"
		cat $pattern.sdsl.lengths
		echo "PHONI"
		cat $pattern.lengths
		break
	fi

	if ! diff -q ${pattern}2.lengths $pattern.lengths; then
		echo "PHONI R"
		cat ${pattern}2.lengths
		echo "PHONI"
		cat $pattern.lengths
		break
	fi

done

