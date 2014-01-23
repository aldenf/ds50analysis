#!/bin/bash

#runlist=( 5370 5372 5373 5375 5376 5378 5379 5381 5382 5384 5385 5387 5388 5390 5391 5393 5396 5398 5399 5401 5402 5404 5405 5412 5413 5420 5421 5423 5424 5440 5441 5444 5445 5449 5450 5452 5453 5455 5456 5460 5461 5464 5468 5469 5471 5472 5473 5475 5478 5481 5484 5485 5488 5489 )

runlist=( 5370 5372 5373 5375 5376 5379 5382 5384 5385 5387 5388 5390 5391 5393 5396 5398 5399 5401 5402 5404 5405 5407 5408 5413 5421 5423 5424 5426 5427 5432 5433 5435 5436 5438 5440 5441 5444 5449 5450 5463 5464 5468 5469 5471 5472 5473 5479 5481 5484 5488 5489 5495 5496 5503 5504 5506 5507 5509 5510 5511 5513 5514 5518 5519 )

#runlist=( 5370 )

njobs=32

exe=s1frac_studies
#exe="/ds50/app/user/aldenf/ds50analysis/darkart/analysis/s1frac_studies"
#exe="/home/aldenf/ds50analysis/darkart/analysis/analysis_richard_exec"

outdir="s1frac_results"    # DIRECTORY in which to store all outputs
#outdir="daq_rate_results"
outbase="result" # base name of output files

datadir="/ds50/data/test_processing/darkart_release3"



################################################################

lbase=$(dirname $exe)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lbase

mkdir -p $outdir

nfiles=${#runlist[@]} 

# determine the number of runs per job
filesperjob=0
if [ $(( nfiles%njobs )) -eq 0 ]
then 
    filesperjob=$(( nfiles/njobs ))
else 
    filesperjob=$(( nfiles/njobs+1 ))
fi


# build an array of all file names
filesarray=""
for run in "${runlist[@]}" ; do
    filename="${datadir}/Run$(printf %06d $run)/Run$(printf %06d $run).root"
    filesarray="$filesarray $filename"
done

# create a bunch of text files with lists of processed ROOT file paths
filenum=0
counter=0
txtfile=""
filelist=""
for i in $filesarray; do
    if [ $(( counter%filesperjob )) -eq 0 ]
    then
        # create the next text file with raw data file names    
        filenum3=$(printf "%03d\n" $filenum)
        txtfile="${outdir}/${outbase}_${filenum3}.txt"
        > $txtfile #clear the file first
        filelist="${filelist} ${txtfile}"
        (( filenum++ ))
    fi
    echo "$i" >> $txtfile
    (( counter++ ))
done


# submit the condor jobs
echo "Submitting $nfiles files with $filenum jobs"
outfiles=""
jobnum=0
jobs=""

for i in $filelist; do

jobnum3=$(printf "%03d\n" $jobnum)
outfile=${outdir}/${outbase}_${jobnum3}.root
cmd_file=${outdir}/${outbase}_${jobnum3}.cmd
err_file=${outdir}/${outbase}_${jobnum3}.err
out_file=${outdir}/${outbase}_${jobnum3}.out
log_file=${outdir}/${outbase}_${jobnum3}.log

cat > $cmd_file <<EOF
 
universe = vanilla
executable = $exe
getenv = True
output = $out_file
error = $err_file
log = $log_file
arguments = $i $outfile
notify_user = aldenf@physics.ucla.edu
notification = Always

queue
EOF

echo "Submitting $cmd_file ..."
condor_submit $cmd_file

(( jobnum++ ))

done
