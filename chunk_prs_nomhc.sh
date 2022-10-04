#! /bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:1000


cd /data/home/hmy117/PRS_UKB_cross_ancestry/

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
--keep prs_subset_iter_${SGE_TASK_ID}\.tsv \
--target /data/scratch/hmy117/ms_chr_# \
--pheno ms_pheno.tsv \
--ld /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr# \
--clump-r2 0.4 \
--perm 1000 \
--missing CENTER \
--bar-levels 5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,0.5,1 \
--fastscore \
--prevalence 0.002 \
--cov ms_covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--print-snp \
--x-range 6:25000000-35000000 \
--out /data/scratch/hmy117/output_prs_nomhc_${SGE_TASK_ID}

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
--keep prs_subset_iter_${SGE_TASK_ID}\.tsv \
--target /data/scratch/hmy117/ms_chr_# \
--extract /data/scratch/hmy117/output_prs_nomhc_${SGE_TASK_ID}\.valid \
--pheno ms_pheno.tsv \
--ld /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr# \
--clump-r2 0.4 \
--perm 1000 \
--missing CENTER \
--bar-levels 5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,0.5,1 \
--fastscore \
--prevalence 0.002 \
--cov ms_covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--print-snp \
--x-range 6:25000000-35000000 \
--out /data/scratch/hmy117/output_prs_nomhc_${SGE_TASK_ID}

