#!/bin/bash

#runlist=( 5370 5372 5373 5375 5376 5378 5379 5381 5382 5384 5385 5387 5388 5390 5391 5393 5396 5398 5399 )
#runlist=( 5370 5372 5373 5375 5376 5379 5382 5384 5385 5387 5388 5390 5391 5393 5396 5398 5399 5401 5402 5404 5405 5407 5408 5413 5421 5423 5424 5426 5427 5432 5433 5435 5436 5438 5440 5441 5444 5449 5450 5463 5464 5468 5469 5471 5472 5473 5479 5481 5484 5488 5489 5495 5496 5503 5504 5506 5507 5509 5510 5511 5513 5514 5516 5518 5519 5521 5522 5524 5525 5529 5530 5531 5533 5534 5536 5539 5540 5542 5545 5547 5548 5549 5551 5552 )
#runlist=( 5370 )

# 2-day campaign (Oct - Dec 2013)
#runlist=( 5369 5370 5372 5373 5375 5376 5378 5379 5381 5382 5384 5385 5387 5388 5390 5391 5393 5396 5398 5399 5401 5402 5404 5405 5407 5408 5412 5413 5417 5418 5421 5423 5424 5426 5427 5432 5433 5435 5436 5438 5440 5441 5444 5445 5449 5450 5452 5453 5455 5456 5460 5461 5463 5464 5468 5469 5471 5472 5473 5475 5476 5478 5479 5481 5484 5485 5488 5489 5495 5496 5500 5501 5503 5504 5506 5507 5509 5510 5511 5513 5514 5516 5519 5521 5522 5524 5525 5529 5530 5533 5534 5536 5539 5540 5542 5544 5547 5548 5549 5551 5552 )

# n-day campaign (Jan - Feb 2013)
runlist=( 6553 6554 6562 6563 6574 6575 6578 7013 7014 7018 7019 7043 7045 7047 7050 7051 7052 7054 7057 7072 7074 7080 7081 7082 7084 7086 7094 7095 7096 7097 7099 7101 7102 7104 7106 7107 7108 7110 7112 7113 7115 7116 7118 7120 7128 7129 7130 7131 7133 7134 7135 7137 7138 7139 7141 7143 7144 7146 7148 7149 7150 7152 )

njobs=16

exe=analysis_jan2014note
#exe=s1frac_studies

outdir="analysis_jan2014note_results"    # DIRECTORY in which to store all outputs
#outdir="s1frac_tdrift_slices_beforecut_results_v4"

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
logfile=${outdir}/${outbase}_${jobnum3}.log

command="./$exe $i $outfile"
echo $command
$command > $logfile &

(( jobnum++ ))

done

echo "Done!"