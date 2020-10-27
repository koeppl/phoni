#!/usr/bin/env zsh

dataset_folder='/s/nut/a/homes/dominik.koeppl/data'
datasets=(chr19.1.fa chr19.16.fa chr19.100.fa chr19.1000.fa chr19.256.fa  chr19.512.fa) #  chr19.100.fa  chr19.128.fa  chr19.16.fa  chr19.1.fa  chr19.256.fa  chr19.512.fa)
pattern=$dataset_folder/chr19.10.fa
tmp_folder=$HOME/tmp/

moni_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/no_thresholds
rrepair_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/rrepair/build/rrepair
bigrepair_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/_deps/bigrepair-src/bigrepair
slpenc_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/_deps/shaped_slp-build/SlpEncBuild
ms_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/build/test/src/ms
readlog_prg=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/readlog.sh
alias Time='/usr/bin/time --format="Wall Time: %e\nMax Memory: %M"'
set -e
test -e $pattern

for filename in $datasets; do
	dataset=$dataset_folder/$filename
	test -e $dataset

	basestats="RESULT file=${filename} "
	stats="$basestats type=baseconstruction "

	logFile=$tmp_folder/$filename.constr.log
	set -x
	Time $moni_prg -f $dataset > "$logFile" 2>&1
	set +x
	echo -n "$stats"
	echo -n "bwtsize=$(stat --format="%s" $dataset.bwt) "
	echo -n "ssasize=$(stat --format="%s" $dataset.ssa) "
	echo -n "esasize=$(stat --format="%s" $dataset.esa) "
	$readlog_prg $logFile

	_basestats=$basestats
	for rrepair_round in $(seq 0 2); do
		basestats="$_basestats rrepair=$rrepair_round "
		if [[ $rrepair_round -eq 0 ]]; then
			logFile=$tmp_folder/$filename.repair.${rrepair_round}.log
			stats="$basestats type=repair "

			set -x
			Time $bigrepair_prg $dataset > "$logFile" 2>&1
			set +x
			echo -n "$stats"
			echo -n "Csize=$(stat --format="%s" $dataset.C) Rsize=$(stat --format="%s" $dataset.R) "
			$readlog_prg $logFile

			logFile=$tmp_folder/$filename.grammar.${rrepair_round}.log
			stats="$basestats type=grammar "
			set -x
			Time $slpenc_prg -e SelfShapedSlp_SdSd_Sd -f Bigrepair -i $dataset -o $dataset.slp > "$logFile" 2>&1
			set +x
			echo -n "$stats"
			echo -n "slpsize=$(stat --format="%s" $dataset.slp) "
			$readlog_prg $logFile
		else 
			logFile=$tmp_folder/$filename.repair.${rrepair_round}.log
			stats="$basestats type=repair "

			cd $(dirname $rrepair_prg)
			set -x
			Time $rrepair_prg $dataset 10 100 $rrepair_round 100000000 > "$logFile" 2>&1
			set +x
			echo -n "$stats"
			$readlog_prg $logFile

			logFile=$tmp_folder/$filename.grammar.${rrepair_round}.log
			stats="$basestats type=grammar "
			set -x
			Time $slpenc_prg -e SelfShapedSlp_SdSd_Sd -f rrepair -i $dataset -o $dataset.slp > "$logFile" 2>&1
			set +x
			echo -n "$stats"
			$readlog_prg $logFile
		fi
		logFile=$tmp_folder/$filename.phoni.${rrepair_round}.log
		stats="$basestats type=phoni "

		set -x
		Time $ms_prg -f "$dataset" -p ${pattern}  > "$logFile" 2>&1
		set +x
		echo -n "$stats"
		$readlog_prg $logFile
	done
done
