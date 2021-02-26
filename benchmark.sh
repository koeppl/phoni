#!/usr/bin/env zsh
# Evaluation script
# Computes an output that is parsable with sqlplot (https://github.com/koeppl/sqlplot)
# Either hardcode PATTERN_FILE for your `hostname` below, or run the script with the pattern filename as the first argument

function die {
	echo "$1" >&2
	exit 1
}

# TODO: the text datasets we want to index: 
datasets=(chr19.1.fa chr19.16.fa chr19.32.fa chr19.64.fa chr19.100.fa chr19.128.fa chr19.256.fa  chr19.512.fa chr19.1000.fa) 

# TODO: add here your hostname and change the variables such that the programs for the evaluation can be found
if [[ $(hostname) = "shannon" ]]; then

PHONI_ROOTDIR=/home/dkoeppl/code/moni
PHONI_BULIDDIR=$PHONI_ROOTDIR/build
DATASET_DIR=/home/dkoeppl/data
INDEXEDMS_BUILDDIR=/s/nut/a/homes/dominik.koeppl/fast/code/indexed_ms/fast_ms/bin/
RREPAIR_BULIDDIR=/home/dkoeppl/code/rrepair/build/
PATTERN_FILE=$DATASET_DIR/chr19.10.fa
LOG_DIR=$HOME/log/

elif [[ $(hostname) = "elm" ]]; then

PHONI_ROOTDIR=/s/nut/a/homes/dominik.koeppl/code/pfp/moni/ 
PHONI_BULIDDIR=$PHONI_ROOTDIR/build
DATASET_DIR='/s/nut/a/homes/dominik.koeppl/fast/data'
INDEXEDMS_BUILDDIR=/s/nut/a/homes/dominik.koeppl/fast/code/indexed_ms/fast_ms/bin/
RREPAIR_BULIDDIR=/s/nut/a/homes/dominik.koeppl/fast/code/rrepair/build/
PATTERN_FILE=$DATASET_DIR/10_samples.fa
LOG_DIR=$HOME/fast/log/

SOLCA_DIR=TODO

elif [[ $(hostname) = "eris" ]]; then

PHONI_ROOTDIR=$HOME/code/phoni
PHONI_BULIDDIR=$PHONI_ROOTDIR/build
DATASET_DIR='/scratch/data/moni'
# INDEXEDMS_BUILDDIR=/s/nut/a/homes/dominik.koeppl/fast/code/indexed_ms/fast_ms/bin/
# RREPAIR_BULIDDIR=/s/nut/a/homes/dominik.koeppl/fast/code/rrepair/build/
PATTERN_FILE=$DATASET_DIR/yeast.fasta
datasets=(coronavirus.fa)
LOG_DIR=/scratch/log

SOLCA_DIR=$HOME/src/solca

else 
	die "need setting for host $hostname"
fi

[[ -n $1 ]] && PATTERN_FILE="$1"

mkdir -p $LOG_DIR || die "cannot create $LOG_DIR for logging"

base_prg=$PHONI_BULIDDIR/no_thresholds
thresholds_prg=$PHONI_BULIDDIR/thresholds
rrepair_prg=$RREPAIR_BULIDDIR/rrepair
bigrepair_prg=$PHONI_BULIDDIR/_deps/bigrepair-src/bigrepair
slpenc_prg=$PHONI_BULIDDIR/_deps/shaped_slp-build/SlpEncBuild
readlog_prg=$PHONI_ROOTDIR/readlog.sh

solca_prg=$SOLCA_DIR/src/compress
[[ ! -x $solca_prg ]] && die "solca not found at $solca_prg"

indexedms_index_prg=$INDEXEDMS_BUILDDIR/dump_cst.x
indexedms_query_prg=$INDEXEDMS_BUILDDIR/matching_stats_parallel.x

indexedms_pprg=$PHONI_ROOTDIR/splitpattern.py
indexedms_tprg=$PHONI_ROOTDIR/catfasta.py

pattern_dir=${PATTERN_FILE}.dir

monibuild_prg=$PHONI_BULIDDIR/test/src/rlebwt_ms_build
moni_prg=$PHONI_BULIDDIR/test/src/rlems
phoni_prg=$PHONI_BULIDDIR/test/src/phoni
phonibuild_prg=$PHONI_BULIDDIR/test/src/build_phoni
rlbwt_prg=$PHONI_BULIDDIR/test/src/bwt2rlbwt



alias Time='/usr/bin/time --format="Wall Time: %e\nMax Memory: %M"'

set -e
for filename in $datasets; do
	dataset=$DATASET_DIR/$filename
	test -e $dataset
done
test -e $PATTERN_FILE

if [[ ! -d $pattern_dir ]]; then
	echo "writing each pattern to $pattern_dir"
	mkdir -p $pattern_dir
	$indexedms_pprg $PATTERN_FILE $pattern_dir
fi

for filename in $datasets; do
	dataset=$DATASET_DIR/$filename
	[[ -r $dataset.raw ]] && continue
	test -e $dataset
	echo "creating raw dataset $dataset.raw"
	$indexedms_tprg $dataset > $dataset.raw
done



for filename in $datasets; do
	dataset=$DATASET_DIR/$filename
	test -e $dataset
	rawdataset=$dataset.raw

	basestats="RESULT patternfile=${PATTERN_FILE} file=${filename} "

	if [[ ! -e $dataset.bwt ]]; then 
		stats="$basestats type=baseconstruction "
		logFile=$LOG_DIR/$filename.constr.log
		set -x
		Time $base_prg -f $dataset > "$logFile" 2>&1
		set +x
		echo -n "$stats"
		echo -n "bwtsize=$(stat --format="%s" $dataset.bwt) "
		echo -n "ssasize=$(stat --format="%s" $dataset.ssa) "
		echo -n "esasize=$(stat --format="%s" $dataset.esa) "
		$readlog_prg $logFile
	fi

	if [[ ! -e $dataset.bwt.heads ]]; then 
		stats="$basestats type=rlebwt "
		logFile=$LOG_DIR/$filename.rlbwt.log
		set -x
		Time $rlbwt_prg $dataset > "$logFile" 2>&1
		set +x
		echo -n "$stats"
		echo -n "headssize=$(stat --format="%s" $dataset.bwt.heads) "
		echo -n "lensize=$(stat --format="%s" $dataset.bwt.len) "
		$readlog_prg $logFile
	fi

	_basestats=$basestats

	
#  #############
#  ## BEGIN MONI
#  #############
#  
#  if [[ ! -e $dataset.thr ]]; then
#  	stats="$basestats type=thresholds "
#  	logFile=$LOG_DIR/$filename.thresholds.log
#  	set -x
#  	Time $thresholds_prg -f $dataset > "$logFile" 2>&1
#  	set +x
#  	echo -n "$stats"
#  	echo -n "thressize=$(stat --format="%s" $dataset.thr) "
#  	echo -n "thrpossize=$(stat --format="%s" $dataset.thr_pos) "
#  	$readlog_prg $logFile
#  fi
#
#  if [[ ! -e $dataset.ms ]]; then
#  stats="$basestats type=monibuild "
#  logFile=$LOG_DIR/$filename.monibuild.log
#  set -x
#  Time $monibuild_prg -f $dataset -p $PATTERN_FILE > "$logFile" 2>&1
#  set +x
#  echo -n "$stats"
#  $readlog_prg $logFile
#  fi
#  
#  stats="$basestats type=moni "
#  logFile=$LOG_DIR/$filename.moni.log
#  set -x
#  Time $moni_prg -f $dataset -p $PATTERN_FILE > "$logFile" 2>&1
#  set +x
#  echo -n "$stats"
#  $readlog_prg $logFile
#  
#  ###########
#  ## END MONI
#  ###########
 
 ##############
 ## BEGIN SOLCA
 ##############
 
 if [[ ! -e $rawdataset.solca ]]; then
	 stats="$basestats type=solca "
	 logFile=$LOG_DIR/$filename.solca.log
	 set -x
	 Time $solca_prg -i ${rawdataset} -o ${rawdataset}.solca -n 1 > "$logFile" 2>&1
	 set +x
	 echo -n "$stats"
	 echo -n "size=$(stat --format="%s" ${rawdataset}.solca) "
	 $readlog_prg $logFile
 fi

 ############
 ## END SOLCA
 ############

 basestats="$_basestats"

if [[ ! -e ${rawdataset}.C ]]; then
 logFile=$LOG_DIR/$filename.repair.log
 stats="$basestats type=repair "

 set -x
 Time $bigrepair_prg $rawdataset > "$logFile" 2>&1
 set +x

 echo -n "$stats"
 echo -n "Csize=$(stat --format="%s" $rawdataset.C) Rsize=$(stat --format="%s" $rawdataset.R) "
 $readlog_prg $logFile
fi


if [[ ! -e $dataset.phoni ]]; then
	logFile=$LOG_DIR/$filename.phonibuild.log
	stats="$basestats type=phonibuild "
	set -x
	Time $phonibuild_prg -f "$dataset" > "$logFile" 2>&1
	set +x
	echo -n "$stats"
	echo -n "phonisize=$(stat --format="%s" $dataset.phoni) "
	$readlog_prg $logFile
fi

for phoni_param in 'none' 'solca' 'solcanaive' 'naive'; do

	if [[ $phoni_param = 'solca' ]] || [[ $phoni_param = 'naivesolca' ]]; then
		grammarinfile=${rawdataset}.solca
		grammartype=solca
	else
		grammarinfile=${rawdataset}
		grammartype=Bigrepair
	fi

	if [[ $phoni_param = 'none' ]] || [[ $phoni_param = 'solca' ]]; then
		grammarencoding='SelfShapedSlp_SdSd_Sd'
		phoniexecparam=''
	else
		grammarencoding='PlainSlp_FblcFblc'
		phoniexecparam='-g naive'
	fi

	logFile=$LOG_DIR/$filename.grammar.${phoni_param}.log
	stats="$basestats type=grammar param=${phoni_param} "
	set -x
	Time $slpenc_prg -e ${grammarencoding}  -f ${grammartype} -i ${grammarinfile} -o $dataset.slp > "$logFile" 2>&1
	set +x
	echo -n "$stats"
	echo -n "slpsize=$(stat --format="%s" $dataset.slp) "
	$readlog_prg $logFile

     logFile=$LOG_DIR/$filename.phoni.log
     stats="$basestats type=phoni param=${phoni_param} "
     set -x
     Time $phoni_prg -f "$dataset" $(echo ${phoniexecparam}) -p ${PATTERN_FILE}  > "$logFile" 2>&1
     set +x
     echo -n "$stats"
     $readlog_prg $logFile
     cp ${PATTERN_FILE}.lengths $LOG_DIR/$filename.phoni.${phoni_param}.lengths
     cp ${PATTERN_FILE}.pointers  $LOG_DIR/$filename.phoni.${phoni_param}.pointers

     if [[ $phoni_param != 'none' ]]; then
	     echo -n "$basestats type=mscheck param=${phoni_param} "
	     echo "check=\"$(diff -q $LOG_DIR/$filename.phoni.none.lengths $LOG_DIR/$filename.phoni.${phoni_param}.lengths)\""
     fi
 done # for phoni_param


done # for dataset



##################
## START indexedms
##################
basestats="$_basestats"
for filename in $datasets; do
	dataset=$DATASET_DIR/$filename
	test -e $dataset


	if [[ ! -r $dataset.raw.fwd.stree ]]; then

		logFile=$LOG_DIR/$filename.indexms.const.log
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
	logFile=$LOG_DIR/$filename.indexms.query.${patternnumber}log
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


echo "finished"
