#!/bin/bash -l
#SBATCH -p defq
#SBATCH -J a5
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0-5:00:00
#SBATCH --mem-per-cpu 4500

module add c3ddb/gaussian/09.d01
which g09

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

WorkDir=/scratch/users/duminda/$SLURM_JOB_NAME-$SLURM_JOB_ID
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

