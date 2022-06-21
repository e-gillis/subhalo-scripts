#!/usr/bin/bash

MAX_JOBS=$1
INFILE="jobs_to_run.dat"
OUTFILE="finished.dat"
LOGFILE="evolution_logfile"
JOBS_LEFT=$(cat $INFILE | wc -l)


touch $LOGFILE
echo "PID $$ Started at $(date +"%R %d %B")" >> $LOGFILE
echo "Maximum Jobs $MAX_JOBS" >> $LOGFILE

while (( $JOBS_LEFT > 0 )); do

    JOBS_RUNNING=$(ps -u gillis | grep job_commands.sh | wc -l)
    
    if (( $JOBS_RUNNING >= $MAX_JOBS )); then
        echo $(date +"%R %d %B"): Not submitting job >> $LOGFILE
        sleep 5m
    else
        # Take directory to process
        DIR=$(head -n 1 $INFILE)
        echo $(date +"%R %d %B"): Submitting job to process: $DIR >> $LOGFILE
        # Move and run the commands
        cd $DIR
        ./job_commands.sh & 
        cd ..
        # Record that the job has been run
        echo $DIR >> $OUTFILE
        # Remove the job from the infile
        sed -i -e "1d" $INFILE
    fi
    
    JOBS_LEFT=$(cat $INFILE | wc -l)
    
done

