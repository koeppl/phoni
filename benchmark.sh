#!/usr/bin/env zsh

dataset_folder='/s/nut/a/homes/dominik.koeppl/fast/data'
datasets=(chr19.1.fa chr19.16.fa chr19.32.fa chr19.64.fa chr19.100.fa chr19.256.fa  chr19.512.fa) # chr19.1000.fa) #  chr19.100.fa  chr19.128.fa  chr19.16.fa  chr19.1.fa  chr19.256.fa  chr19.512.fa)
#datasets=(chr19.64.fa) #chr19.16.fa chr19.100.fa chr19.128.fa chr19.256.fa  chr19.512.fa chr19.1000.fa) # chr19.256.fa  chr19.512.fa) #  chr19.100.fa  chr19.128.fa  chr19.16.fa  chr19.1.fa  chr19.256.fa  chr19.512.fa)
#datasets=(chr19.100.fa  chr19.128.fa chr19.256.fa  chr19.512.fa chr19.1000.fa)
pattern=$dataset_folder/10_samples.fa

log_folder=$HOME/fast/log/

moni_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/no_thresholds
rrepair_prg=/s/nut/a/homes/dominik.koeppl/fast/code/rrepair/build/rrepair
bigrepair_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/_deps/bigrepair-src/bigrepair
slpenc_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/_deps/shaped_slp-build/SlpEncBuild
#ms_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/test/src/ms
ms_prg=/s/nut/a/homes/dominik.koeppl/msrelease
readlog_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/readlog.sh

indexedms_index_prg=/s/nut/a/homes/dominik.koeppl/fast/code/indexed_ms/fast_ms/bin/dump_cst.x
indexedms_query_prg=/s/nut/a/homes/dominik.koeppl/fast/code/indexed_ms/fast_ms/bin/matching_stats_parallel.x


indexedms_pprg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/splitpattern.py
indexedms_tprg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/catfasta.py
pattern_dir=$dataset_folder/samplesdir

alias Time='/usr/bin/time --format="Wall Time: %e\nMax Memory: %M"'
set -e
test -e $pattern

if [[ ! -d $pattern_dir ]]; then
	echo "writing each pattern to $pattern_dir"
	mkdir -p $pattern_dir
	$indexedms_pprg $pattern $pattern_dir
fi

for filename in $datasets; do
	dataset=$dataset_folder/$filename
	[[ -r $dataset.raw ]] && continue
	test -e $dataset
	echo "creating raw dataset $dataset.raw"
	$indexedms_tprg $dataset > $dataset.raw
done

 ##################
 ## START indexedms
 ##################
for filename in $datasets; do
	dataset=$dataset_folder/$filename
	test -e $dataset

	basestats="RESULT file=${filename} "

	if [[ ! -r $dataset.raw.fwd.stree ]]; then

		logFile=$log_folder/$filename.indexms.const.log
		stats="$basestats type=indexedconst "

		set -x
		Time $indexedms_index_prg -s_path $dataset.raw > "$logFile" 2>&1
		set +x

		echo -n "$stats"
		echo -n "fwdsize=$(stat --format="%s" $dataset.raw.fwd.stree) revsize=$(stat --format="%s" $dataset.raw.rev.stree) "
		$readlog_prg $logFile
	fi
 
 for patternseq in $pattern_dir/[0-9]; do
	patternnumber=$(basename $patternseq)
 	logFile=$log_folder/$filename.indexms.query.${patternnumber}log
 	stats="$basestats type=indexedquery pattern=${patternnumber} "
 
 	set -x
 	Time $indexedms_query_prg -s_path $dataset.raw -t_path $patternseq -load_cst 1 > "$logFile" 2>&1
 	set +x
 
 	echo -n "$stats"
 	$readlog_prg $logFile
 done

done
 ##################
 ## END indexedms
 ##################




for filename in $datasets; do
	dataset=$dataset_folder/$filename
	test -e $dataset

	basestats="RESULT file=${filename} "


#	if [[ ! -e $dataset.bwt ]]; then 
		stats="$basestats type=baseconstruction "
		logFile=$log_folder/$filename.constr.log
		set -x
		Time $moni_prg -f $dataset > "$logFile" 2>&1
		set +x
		echo -n "$stats"
		echo -n "bwtsize=$(stat --format="%s" $dataset.bwt) "
		echo -n "ssasize=$(stat --format="%s" $dataset.ssa) "
		echo -n "esasize=$(stat --format="%s" $dataset.esa) "
		$readlog_prg $logFile
#	fi
	continue

	_basestats=$basestats
	rawdataset=$dataset.raw
	for rrepair_round in $(seq 0 2); do
		basestats="$_basestats rrepair=$rrepair_round "
		if [[ $rrepair_round -eq 0 ]]; then
			logFile=$log_folder/$filename.repair.${rrepair_round}.log
			stats="$basestats type=repair "

			set -x
			Time $bigrepair_prg $rawdataset > "$logFile" 2>&1
			set +x

			echo -n "$stats"
			echo -n "Csize=$(stat --format="%s" $rawdataset.C) Rsize=$(stat --format="%s" $rawdataset.R) "
			$readlog_prg $logFile

			logFile=$log_folder/$filename.grammar.${rrepair_round}.log
			stats="$basestats type=grammar "
			set -x
			Time $slpenc_prg -e SelfShapedSlp_SdSd_Sd -f Bigrepair -i $rawdataset -o $dataset.slp > "$logFile" 2>&1
			set +x
			echo -n "$stats"
			echo -n "slpsize=$(stat --format="%s" $dataset.slp) "
			$readlog_prg $logFile
		else 
			logFile=$log_folder/$filename.repair.${rrepair_round}.log
			stats="$basestats type=repair "

			cd $(dirname $rrepair_prg)
			set -x
			Time $rrepair_prg $rawdataset 10 100 $rrepair_round 100000000 > "$logFile" 2>&1
			set +x

			echo -n "$stats"
			echo -n "Csize=$(stat --format="%s" $rawdataset.C) Rsize=$(stat --format="%s" $rawdataset.R) "
			$readlog_prg $logFile

			logFile=$log_folder/$filename.grammar.${rrepair_round}.log
			stats="$basestats type=grammar "
			set -x
			Time $slpenc_prg -e SelfShapedSlp_SdSd_Sd -f rrepair -i $rawdataset -o $dataset.slp > "$logFile" 2>&1
			set +x
			echo -n "$stats"
			echo -n "slpsize=$(stat --format="%s" $dataset.slp) "
			$readlog_prg $logFile
		fi
		logFile=$log_folder/$filename.phoni.${rrepair_round}.log
		stats="$basestats type=phoni "

		set -x
		Time $ms_prg -f "$dataset" -p ${pattern}  > "$logFile" 2>&1
		set +x
		echo -n "$stats"
		$readlog_prg $logFile
		cp ${pattern}.lengths $log_folder/$filename.phoni.${rrepair_round}.lengths
		cp ${pattern}.pointers  $log_folder/$filename.phoni.${rrepair_round}.pointers
		if [[ $rrepair_round -gt 0 ]]; then
			echo -n "$basestats type=mscheck "
			echo "check=\"$(diff -q $log_folder/$filename.phoni.${rrepair_round}.lengths $log_folder/$filename.phoni.0.lengths)\""
		fi
	done


done




echo "finished"
