#!/bin/bash
#SBATCH --job-name a1671
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p defq
#SBATCH --mem-per-cpu 2500

module add c3ddb/gaussian/09.d01
which g09
echo 'SLURM_JOB_NAME': $SLURM_JOB_NAME
echo 'SLURM_JOB_NAME': $SLURM_JOB_ID

WorkDir=/scratch/users/duminda/$SLURM_JOB_NAME
SubmitDir=`pwd`

GAUSS_SCRDIR=/scratch/users/duminda/g09/$SLURM_JOB_NAME-$SLURM_JOB_ID
export  GAUSS_SCRDIR

mkdir -p $GAUSS_SCRDIR
mkdir -p $WorkDir

cd  $WorkDir
. $g09root/g09/bsd/g09.profile

cp $SubmitDir/input.gjf .

g09 < input.gjf > input.log
formchk  check.chk check.fchk
cp * $SubmitDir/

rm -rf $GAUSS_SCRDIR
rm -rf $WorkDir

