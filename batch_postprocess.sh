#!/usr/bin/bash

MAX_JOBS=28

echo PID $$

for dir in ??????_*/; do

        RUNNING=$(ps -u gillis | grep postprocessing | wc -l)
        while [[ $MAX_JOBS -le $RUNNING ]]; do
                sleep 10s
                RUNNING=$(ps -u gillis | grep postprocessing | wc -l)
        done 

        echo processing $dir
        # postprocessing.py $dir >> $dir"postprocessing_log.txt" 2>&1 &
        postprocessing_class.py $dir $1 --silent &
        
done

