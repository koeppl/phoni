#!/usr/bin/env zsh

function die {
    echo "$1" >&2
    exit 1
}
[[ ! -r $1 ]] && die "cannot read log file $1"

logFile=$1

replacements=(
'^malloc_count ### exiting.* peak: \([0-9.,]\+\),.*' 'mem'
'^Wall Time:\s\+\([0-9]\+\)' 'time'
'^Max Memory:\s\+\([0-9]\+\)' 'mrs'
)
replacements_it=1
while [[ $replacements_it -lt $#replacements ]]; do
    ((subst_it=replacements_it+1))
    pattern=$replacements[$replacements_it]
    value=$(grep "$pattern" "$logFile" | sed "s@$pattern@\1@")
    [[ -n $value ]] && echo -n "$replacements[$subst_it]=$value "
    ((replacements_it+=2))
done

replacements=(
'^Time backwardsearch:\s\+\([0-9]\+\)' 'timebackwardsearch'
'^Time lce:\s\+\([0-9]\+\)' 'timelce'
'^Count lce:\s\+\([0-9]\+\)' 'countlce'
'^Count lce skips:\s\+\([0-9]\+\)' 'countlceskips'
)

replacements_it=1
while [[ $replacements_it -lt $#replacements ]]; do
    ((subst_it=replacements_it+1))
    pattern=$replacements[$replacements_it]
    value=$(grep "$pattern" "$logFile" | sed "s@$pattern@\1@" | tr '\n' '+')
	value=$(echo "print($value 0)" | python3)
    [[ -n $value ]] && echo -n "$replacements[$subst_it]=$value "
    ((replacements_it+=2))
done
echo ""


