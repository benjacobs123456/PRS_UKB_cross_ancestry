#! /bin/bash
#$ -pe smp 20
#$ -l h_vmem=20G
#$ -l highmem
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:7


cd /data/home/hmy117/PRS_UKB_cross_ancestry/

rs=(0.001 0.01 0.05 0.1 0.2 0.4 0.6)
index=$((${SGE_TASK_ID}-1))
r=${rs[$index]}

echo clumping r2 is $r

~/PRSice_linux \
--base ms_prs_snps.tsv \
--chr CHR \
--A1 effect_allele \
--A2 other_allele \
--stat beta \
--beta \
--snp SNP \
--bp BP \
--pvalue pval \
--thread $NSLOTS \
--binary-target T \
--score sum \
--target /data/scratch/hmy117/ms_chr_# \
--pheno ms_pheno.tsv \
--ld /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr# \
--clump-r2 $r \
--missing CENTER \
--bar-levels 1e-8,5e-8,1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,0.01,5e-2,0.1,0.5,1 \
--fastscore \
--prevalence 0.002 \
--cov ms_covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--print-snp \
--out /data/scratch/hmy117/whole_cohort_mhc_output_prs_r$r

# After first run there are some duplicates > re-run with extract
~/PRSice_linux \
--base ms_prs_snps.tsv \
--chr CHR \
--A1 effect_allele \
--A2 other_allele \
--stat beta \
--beta \
--snp SNP \
--bp BP \
--pvalue pval \
--thread $NSLOTS \
--binary-target T \
--score sum \
--target /data/scratch/hmy117/ms_chr_# \
--pheno ms_pheno.tsv \
--ld /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr# \
--clump-r2 $r \
--missing CENTER \
--bar-levels 1e-8,5e-8,1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,0.01,5e-2,0.1,0.5,1 \
--fastscore \
--prevalence 0.002 \
--cov ms_covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--extract /data/scratch/hmy117/whole_cohort_mhc_output_prs_r$r\.valid \
--print-snp \
--out /data/scratch/hmy117/whole_cohort_mhc_output_prs_r$r

