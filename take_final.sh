#!/usr/bin/bash

for dir in ??????_??????/; do
    echo processing $dir
    # Convert restart if it exists
    for f in $dir*z07.nemo; do
        s2a $dir*restart.nemo $dir"z07_final_frame.dat" 
    done
    
    # Get final snap from z05
    for f in $dir*z05.nemo; do
        snaptrim in=$f out=$dir"z05_final_frame.nemo" times=#11
        s2a $dir"z05_final_frame.nemo" $dir"z05_final_frame.dat" 
        rm $dir"z05_final_frame.nemo"
    done

done
