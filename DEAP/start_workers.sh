#!/bin/bash

# Array to hold worker PIDs
declare -a workerPIDs
MCRROOT="/usr/pubsw/packages/matlab/R2024b"

LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
export LD_LIBRARY_PATH;
#fstem_imaging=JA
fstem_imaging=RNI

#DATA_DIR="./data"
DATA_DIR="./data_v5.0"

# Function to start workers
start_workers() {
    local num_workers=$1
    for (( i=1; i<=num_workers; i++ ))
    do
        #./test_worker/test_worker "$i" &
        ./FEMA_DEAP_worker/FEMA_DEAP_worker $fstem_imaging voxel "$i" $num_workers $DATA_DIR/FEMA_cache $DATA_DIR/jobs &
        workerPIDs+=($!)  # Store the PID of the background process
    done
}

# Function to kill all workers
kill_workers() {
    echo "Killing all worker processes..."
    for pid in "${workerPIDs[@]}"
    do
        kill "$pid"
    done
}

# Function to handle SIGINT (Ctrl+C)
trap_ctrlc() {
    echo "Ctrl-C caught...performing clean up"
    
    kill_workers
    
    echo "Exiting script"
    exit 2
}

# Trap SIGINT (Ctrl+C) and execute trap_ctrlc function
trap "trap_ctrlc" 2

# Main script execution
read -p "Enter the number of workers to start: " num_workers
start_workers "$num_workers"

echo "Workers started. Press Ctrl+C to stop them."

# Wait indefinitely until Ctrl+C is pressed
while true
do
    sleep 1
done
