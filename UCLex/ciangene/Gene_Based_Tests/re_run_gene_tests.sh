#!/bin/bash
shopt -s expand_aliases
source ~/.bashrc

Dir=$(pwd)
for dir in $(find . -name "*_kin*" )
do
	for f in {1..9}
	do 
		if [ ! -f $dir/regress${f} ]; then
		echo $dir/regress${f}
		cd $dir
		runSh partition_${f}.sh
		cd $Dir
		fi
	done
done
