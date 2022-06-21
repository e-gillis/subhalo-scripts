#!/usr/bin/bash

jobs_printout () {
	for dir in peri100/m??/; do
		echo fully processed in $dir:
		ls $dir/??????_??????/done | wc -l
	done
}

cd ~/scratch/prod/bary
echo Baryonic Jobs:
jobs_printout

cd ~/scratch/prod/dm
echo Dark Matter Jobs
jobs_printout

