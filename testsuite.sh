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
		/home/niki/code/pfp/moni/debug/test/src/bwt2rlbwt $rdataset
		rm $rdataset.R $rdataset.C $rdataset.slp
		[[ ! -f $rdataset.plain ]] && ~/extractfasta.py $rdataset > $rdataset.plain
		cd /home/niki/code/pfp/rrepair/build/ && /home/niki/code/pfp/rrepair/build/rrepair $rdataset.plain 10 100 10 100000000
		#/home/niki/code/pfp/moni/build/_deps/bigrepair-src/bigrepair $rdataset.plain
		/home/niki/code/pfp/moni/build/_deps/shaped_slp-build/SlpEncBuild -e SelfShapedSlp_SdSd_Sd -f rrepair -i $rdataset.plain -o $rdataset.slp
	# fi
	

	cd /home/niki/code/pfp/moni/debug
	while ! make; do sleep 10; done


	truncate -s0 $pattern
	for i in $(seq 1 100); do 
		echo "> Seq $i " >> $pattern
		genDNA 100 >> $pattern
	done

	[[ -d ${pattern}.dir ]] && rm -r ${pattern}.dir
	mkdir -p ${pattern}.dir
	/home/niki/code/pfp/moni/splitpattern.py $pattern ${pattern}.dir
	[[ -d ${pattern}2.dir ]] && rm -r ${pattern}2.dir
	cp -a ${pattern}.dir ${pattern}2.dir

	set -x
	/home/niki/code/pfp/moni/debug/test/src/build_phoni -f "$dataset"
	/home/niki/code/pfp/moni/debug/test/src/phoni -f "$dataset" -p "$pattern"
	/home/niki/code/pfp/moni/build/_deps/pfp_thresholds-build/test/src/sdsl_matching_statistics -f "$dataset" -p $pattern
	set +x
	# /home/niki/code/pfp/moni/build/_deps/pfp_thresholds-build/test/src/matching_statistics -f "$dataset" -p $pattern
	
	set -x
	/home/niki/code/pfp/moni/debug/test/src/build_phoni -f "$rdataset"
	/home/niki/code/pfp/moni/debug/test/src/phoni -f "$rdataset" -p "$pattern"2
	set +x

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

