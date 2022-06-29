#!/usr/bin/bash

MAX_JOBS=10

echo PID $$

for dir in ??????_??????/; do

        RUNNING=$(ps -u gillis | grep postprocessing | wc -l)
        while [[ $MAX_JOBS -le $RUNNING ]]; do
                sleep 150s
                RUNNING=$(ps -u gillis | grep postprocessing | wc -l)
        done 

        echo processing $dir
        postprocessing_class.py $dir $1 --silent &
        
done

