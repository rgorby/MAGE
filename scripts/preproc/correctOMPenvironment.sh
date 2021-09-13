#!/bin/sh
# usage: correctOMPenvironment.sh NODEFILE prog <args>

# get my host and figure out how many other jobs are on this node
myHost=`hostname -f`
numJobsOnMyNode=`grep -c $myHost $1`

# split the CPUs on this node between me and the other jobs
#echo "myHost = $myHost"
#echo "numJobsOnMyNode = $numJobsOnMyNode"
let newNumThreads=$OMP_NUM_THREADS/$numJobsOnMyNode
export OMP_NUM_THREADS=$newNumThreads
#echo "Setting OMP_NUM_THREADS = $OMP_NUM_THREADS"

# call the arguments to this script
$2 "${@:3}"

