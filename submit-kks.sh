#!/bin/bash
#PBS -N kks_v003_k4_without_noise
#PBS -l ncpus=1
#PBS -M yusuke.seguchi@mat.eng.osaka-u.ac.jp
#PBS -m abe
#PBS -o /home/seguchi/log/
#PBS -e /home/seguchi/log/

ROOT_DIR="$PBS_O_WORKDIR"
FILE_POT="kks.f90"
FILE_DATA="mkdir.sh"
FILE_INPUT="kks"

job_number=`echo $PBS_JOBID|sed 's/.mpr01//'`
JOB_DIR="$PBS_O_WORKDIR/job_${job_number}-$PBS_JOBNAME"
echo "Job $job_number started at `date`"
echo "$job_number $PBS_JOBNAME $HOSTNAME $PBS_O_WORKDIR" >> /home/seguchi/jobs.dat

mkdir -p ${JOB_DIR}

cp ${ROOT_DIR}/${FILE_DATA} ${JOB_DIR}
cp ${ROOT_DIR}/${FILE_POT} ${JOB_DIR}
cp ${ROOT_DIR}/${FILE_INPUT} ${JOB_DIR}


cd ${JOB_DIR}


sh ${JOB_DIR}/${FILE_DATA}
./${FILE_INPUT} 

echo "Job $job_number ended at `date`"
