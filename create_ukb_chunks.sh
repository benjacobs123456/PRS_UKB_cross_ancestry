#! /bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:1000


cd /data/home/hmy117/PRS_UKB_cross_ancestry/
module load R/3.6.1
Rscript chunk.R ${SGE_TASK_ID}

