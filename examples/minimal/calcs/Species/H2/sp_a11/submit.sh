#! /bin/bash -l

#$ -N a11
#$ -l long
#$ -l harpertown
#$ -l h_rt=5:00:00
#$ -pe singlenode 6
#$ -l h=!node60.cluster
#$ -cwd
#$ -o out.txt
#$ -e err.txt

export PATH=/opt/molpro2012/molprop_2012_1_Linux_x86_64_i8/bin:$PATH

sdir=/scratch/duminda
mkdir -p /scratch/duminda/qlscratch

molpro -d $sdir -n 6 input.in
