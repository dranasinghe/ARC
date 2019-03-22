#!/bin/bash
#SBATCH -p normal
#SBATCH -J a1658
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2048

export PATH=/opt/molpro/molprop_2015_1_linux_x86_64_i8/bin:$PATH

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

WorkDir=`pwd`
cd
source .bashrc
sdir=/scratch/duminda/$SLURM_JOB_ID
mkdir -p $sdir
export TMPDIR=$sdir
cd $WorkDir

molpro input.in

rm -rf $sdir

