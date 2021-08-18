# PRS distributions in UKB MS

We're going to generate polygenic risk scores for MS using GWAS from European populations. We're going to apply these scores in UKB and demonstrate that they perform poorly for people of non-European ancestries.

## 0. Filter UKB imputed genotype files in PLINK2

````unix
#! /bin/bash
#$ -pe smp 20
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -t 8
echo ${SGE_TASK_ID}

# Move into the temporary scratch space where your data now is
cd $TMPDIR

module load plink/2.0-20170920
plink2 \
--bgen /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/ukb_imp_chr${SGE_TASK_ID}\_v3.bgen \
--sample /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/ukb43101_imp_chr1_v3_s487314.sample \
--threads $NSLOTS \
--out temp_chr${SGE_TASK_ID}
echo "converted bgen to pgen"

plink2 \
--pfile temp_chr${SGE_TASK_ID} \
--mach-r2-filter 0.3 \
--threads $NSLOTS \
--out chr_${SGE_TASK_ID} \
--make-pgen \
--maf 0.01 \
--hwe 0.00001 \
--mind 0.1 \
--geno 0.1

echo "applied filters to pgen"

cp $TMPDIR/chr_${SGE_TASK_ID}\.pgen /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files
cp $TMPDIR/chr_${SGE_TASK_ID}\.pvar /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files
cp $TMPDIR/chr_${SGE_TASK_ID}\.psam /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files
echo "moved pgen back to normal directory"
````

## 1. Generate polygenic scores  
````R
library(dplyr)
library(readr)

# read in MS chip discovery GWAS stats
df = read_table2("/data/Wolfson-UKBB-Dobson/mschip_discovery_gwas/discovery_metav3.0.meta")

# read in UKB SNPs from QCd plink2 files  
ukb_snps = data.frame()
for(i in c(1:22)){
  snps = read_table2(paste0("/data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files/chr_",i,".pvar"))
  ukb_snps <<- bind_rows(ukb_snps,snps)
}

# filter MS GWAS to SNPs in post-QC UKB dataset
df = df %>% mutate(SNPID = paste0(CHR,":",BP))
ukb_snps = ukb_snps %>% mutate(SNPID = paste0(`#CHROM`,":",POS))
df = df %>% filter(SNPID %in% ukb_snps$SNPID)

# filter to P < 5e-8
df = df %>% filter(P<5e-8)

# filter to biallelic SNVs
df = df %>% filter(nchar(A1)==1 & nchar(A2)==1)

# remove overlapping positions
df = df %>% distinct(SNPID,.keep_all=TRUE)

# filter to rsIDs only
df = df %>% filter(grepl("rs",SNP))

# prepare file for score and clumping
df = df %>% mutate(beta = log(OR)) %>% select(SNP,A1,beta,P)

# write to file
write_tsv(df,"ms_prs_snps")
````

# Clump SNPs per chromosome in plink using EUR 1kg samples
````unix
cd /data/Wolfson-UKBB-Dobson/PRS_cross_ancestry/
module load plink

rm ms_prs_snps_clumped
for i in {1..22};
  do
    plink --bfile /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr$i \
  --clump ms_prs_snps \
  --clump-r2 0.1 \
  --out clumped_snps_chr$i
    awk 'NR>1{print $3}' clumped_snps_chr$i\.clumped >> ms_prs_snps_clumped
  done
````

## 2. Apply polygenic scores to UKB
````unix
cd /data/Wolfson-UKBB-Dobson/PRS_cross_ancestry/
module load plink/2.0-20200328

for i in {1..22};
  do
    plink2 --pfile /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files/chr_$i \
    --score ms_prs_snps list-variants variance-standardize cols=+scoresums \
    --extract ms_prs_snps_clumped \
    --out ms_prs_scores_chr$i
  done
````

### PCA with 1kg samples and UKB
First download 1kg reference samples
````unix
#! /bin/bash
#$ -pe smp 40
#$ -l h_vmem=10G
#$ -l highmem
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117

# wget 1kg
wget ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr${SGE_TASK_ID}\.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
````
Then merge with UKB PLINK files
````unix
#! /bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -o /data/scratch/hmy117
#$ -t 1:22

cd /data/scratch/hmy117
module load plink
# convert 1kg to plink
plink --vcf ALL.chr${SGE_TASK_ID}\.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--make-bed \
--biallelic-only \
--out plink_1kg_chr${SGE_TASK_ID}

# join 1kg and plink
plink --bfile plink_1kg_chr${SGE_TASK_ID} \
--bmerge /data/Wolfson-UKBB-Dobson/raw_ukb_genotype_files/chr${SGE_TASK_ID} \
--out combined_plink_chr${SGE_TASK_ID}

plink --bfile plink_1kg_chr${SGE_TASK_ID} \
--exclude combined_plink_chr${SGE_TASK_ID}\.missnp \
--make-bed \
--out fixed_1kg_chr${SGE_TASK_ID}

# join 1kg and plink
plink --bfile fixed_1kg_chr${SGE_TASK_ID} \
--bmerge /data/Wolfson-UKBB-Dobson/raw_ukb_genotype_files/chr${SGE_TASK_ID} \
--out combined_plink_chr${SGE_TASK_ID}

# prune
plink --bfile combined_plink_chr${SGE_TASK_ID} \
--indep-pairwise 1000 50 0.5 \
--geno 0.01 \
--maf 0.05 \
--hwe 1e-5 \
--out pruned_snps_chr${SGE_TASK_ID}

plink --bfile combined_plink_chr${SGE_TASK_ID} \
--extract pruned_snps_chr${SGE_TASK_ID}\.prune.in \
--make-bed \
--out combined_pruned_genotypes_for_pca_chr${SGE_TASK_ID}
````
### Combine genotypes and compute PCs
````unix
#! /bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -o /data/scratch/hmy117


cd /data/scratch/hmy117

module load singularity
module load plink
module load R

rm merge_filelist
for i in {2..22};
  do
    echo "combined_pruned_genotypes_for_pca_chr$i.bed combined_pruned_genotypes_for_pca_chr$i.bim combined_pruned_genotypes_for_pca_chr$i.fam" >> merge_filelist
  done

plink --bfile combined_pruned_genotypes_for_pca_chr1 \
--merge-list merge_filelist \
--make-bed \
--out combined_genotypes_for_pca

singularity exec /data/containers/flashpca/debian-flashpca.simg flashpca --bfile combined_genotypes_for_pca
````

Now we can explore the data in R.

````R
library(dplyr)
library(readr)
library(ggplot2)
setwd("/data/Wolfson-UKBB-Dobson/PRS_cross_ancestry/")

# read in scores per chrom
prs_scores = data.frame()
for(i in c(1:22)){
  if( file.exists(paste0("ms_prs_scores_chr",i,".sscore"))){
    df = read_tsv(paste0("ms_prs_scores_chr",i,".sscore"))
    prs_scores <<- bind_rows(prs_scores,df)
  }
}

prs_scores = as.tbl(prs_scores)

# collate scores across chroms
prs_scores = prs_scores %>% group_by(IID) %>%
summarise(PRS = sum(SCORE1_SUM))

# remove missing indivs
prs_scores = prs_scores %>% filter(IID>=0)

# read in phenotype data
source_of_report_data = read_tsv("../ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G35 (multiple sclerosis).0.0` = col_character(),
  EID = col_double()
  ))

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))

# read in remainder of phenotype data
ukb_pheno = read_tsv("../ukb_pheno_911/ukb_pheno_final_MS_1301")
ukb_pheno = ukb_pheno %>% select(EID,`Age at recruitment.0.0`,Sex.0.0,`Ethnic background.0.0`,contains("rincipal components"),contains("ethnic grouping"))

# combine
ukb_pheno = ukb_pheno %>% left_join(source_of_report_data,by="EID")

# tabulate MS cases vs self-reported ethnicity
ukb_pheno = ukb_pheno %>% mutate(genetically_british = ifelse(!is.na(`Genetic ethnic grouping.0.0`),1,0))
table(ukb_pheno$MS_status,ukb_pheno$genetically_british)
````
Tabulating the genetic ancestry (Genetic ethnic grouping) vs MS status, we can see that ~300 people from non-European genetic backgrounds have MS

MS status | Non-british | British
--------- | ----------- | -------
Control   | 78340       | 407493
Case      | 303         | 2102

````R
# combine with PRS
ukb_pheno = ukb_pheno %>% left_join(prs_scores %>% rename(EID=IID),by="EID")

# remove NAs
ukb_pheno = ukb_pheno %>% filter(!is.na(MS_status))
````
Plotting the PRS distribution in the whole cohort shows a nice separation of cases and controls.
````R
png("prs_distro_overall.png",height=6,width=6,res=300,units="in")
ggplot(ukb_pheno,aes(PRS,fill=factor(MS_status)))+geom_density(alpha=0.5)+labs(fill="MS case/control status")+scale_fill_brewer(palette="Set2",labels=c("Control","MS"))
dev.off()
````
!(https://github.com/benjacobs123456/PRS_UKB_cross_ancestry/blob/main/prs_distro_overall.png)

The separation of cases and controls is weakened among people of non-British (i.e. non-European) genetic ancestry.
````R
ukb_pheno = ukb_pheno %>% mutate(genetically_british = ifelse(genetically_british==0,"Genetically non-British","Genetically British"))

png("prs_distro_ancestry.png",height=6,width=6,res=300,units="in")
ggplot(ukb_pheno,aes(PRS,fill=factor(MS_status)))+geom_density(alpha=0.5)+labs(fill="MS case/control status")+scale_fill_brewer(palette="Set2",labels=c("Control","MS"))+facet_wrap(~genetically_british)
dev.off()
````
!(https://github.com/benjacobs123456/PRS_UKB_cross_ancestry/blob/main/prs_distro_ancestry.png)
Calculating variance explained (Nagelkerke's Pseudo-R2) demonstrates this quantitatively.
````R
library(rcompanion)

null_model = glm(data=ukb_pheno,MS_status~`Age at recruitment.0.0`+Sex.0.0,family=binomial(link="logit"))
full_model = glm(data=ukb_pheno,MS_status~`Age at recruitment.0.0`+Sex.0.0+PRS,family=binomial(link="logit"))

nagelkerke(full_model,null_model)$Pseudo.R.squared.for.model.vs.null
