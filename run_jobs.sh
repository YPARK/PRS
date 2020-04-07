#!/bin/bash -l

jobs=$1

jobname=$(echo $jobs | awk '{ gsub("/","_"); print }')
command=$(zless $jobs | head -n $SGE_TASK_ID | tail -n 1)
mkdir -p log/rscript/${jobname}/
log=log/rscript/${jobname}/$SGE_TASK_ID.log
[ -f $log ] && rm -f $log

source /broad/software/scripts/useuse >> $log

# reuse Anaconda3  >> $log
reuse GCC-5.2  >> $log
reuse OpenblasR  >> $log
reuse R-3.5  >> $log
reuse Python-3.4  >> $log
reuse .openblas-0.2.8  >> $log
reuse .boost-1.60.0  >> $log
reuse .zlib-1.2.8  >> $log
reuse GSL  >> $log
reuse .htslib-1.3.2 >> $log

export LD_LIBRARY_PATH=${LIBRARY_PATH}:${LD_LIBRARY_PATH}

printf "[%s] Running ... \n$command" "$(date)" >> $log
printf "\n\n" >> $log

exec $command >> $log 2>&1 || exit 1

## save some disk space if it is done successfully
[ -f $log ] && rm $log
printf "[%s] Running ... \n$command" "$(date)" >> $log
printf "[%s] \nDone\n\n" "$(date)" >> $log

