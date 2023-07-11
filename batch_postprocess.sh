#!/usr/bin/bash

MAX_JOBS=20

INFILE="files_to_process.dat"
OUTFILE="processed.dat"
LOGFILE="postprocessing_log"

touch $LOGFILE
echo "PID $$ Started at $(date +"%R %d %B")" >> $LOGFILE
echo "Maximum Jobs $MAX_JOBS" >> $LOGFILE


# rm files_to_process.dat
touch processed.dat

# subhalo_update_list.py

JOBS_LEFT=$(cat $INFILE | wc -l)

while (( $JOBS_LEFT > 0 )); do

    JOBS_RUNNING=$(ps -u gillis | grep postprocessing | wc -l)
    if (( $JOBS_RUNNING >= $MAX_JOBS )); then
            sleep 240s
            JOBS_RUNNING=$(ps -u gillis | grep postprocessing | wc -l)
    else 
        DIR=$(head -n 1 $INFILE)
        
        for f in $DIR*z??.nemo; do s2a $f $f".dat"; done
        echo $(date +"%R %d %B"): Submitting job to process: $DIR >> $LOGFILE
        postprocessing_class.py $DIR $1 --silent --remove_datfiles >> $LOGFILE 2>&1 &
        
        echo $DIR >> $OUTFILE
        sed -i -e "1d" $INFILE
    fi
    
    JOBS_LEFT=$(cat $INFILE | wc -l)

done

