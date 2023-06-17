#!/usr/bin/env bash

# usage: pinCpuCores.sh prog <args>

# THIS IS ONLY MEANT TO WORK ON DERECHO
# USE ON OTHER SYSTEMS AT YOUR OWN RISK

let newNumThreads=$OMP_NUM_THREADS/$PMI_LOCAL_SIZE
let minThread=$newNumThreads*$PMI_LOCAL_RANK
let maxThread=$minThread+$newNumThreads-1
export OMP_NUM_THREADS=$newNumThreads

#echo "I am rank $PMI_LOCAL_RANK of $PMI_LOCAL_SIZE on this node using cpus $minThread to $maxThread"

# call the arguments to this script
taskset --cpu-list $minThread-$maxThread $1 "${@:2}"

