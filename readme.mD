# MS-PRS in a South Asian-ancestry cohort

This pipeline was used to generate PRS in UKB participants of EUR ancestry, and compare the performance with a similarly-constructed PRS in Genes & Health. The pipeline used for Genes & Health is gh_prs_pipeline.mD[https://github.com/benjacobs123456/PRS_UKB_cross_ancestry/blob/main/gh_prs_pipeline.mD].

## Filter UKB imputed genotype files in PLINK2

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

## Variant & Individual QC  
````R
library(dplyr)
library(readr)
setwd("/data/home/hmy117/PRS_UKB_cross_ancestry/")
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

# filter to biallelic SNVs
df = df %>% filter(nchar(A1)==1 & nchar(A2)==1)

# remove overlapping positions
df = df %>% distinct(SNPID,.keep_all=TRUE)

# remove duplicated rsids
df = df %>% distinct(SNP,.keep_all=TRUE)

# filter to rsIDs only
df = df %>% filter(grepl("rs",SNP))

# filter NAs
df = df %>% filter(!is.na(P))

# filter out palindromic
df = df %>%
filter(!(A1 == "T" & A2 =="A") & !(A1=="A" & A2=="T")) %>%
filter(!(A1 == "C" & A2 =="G") & !(A1=="G" & A2=="C"))

# prepare file for score and clumping
df = df %>% mutate(beta = log(OR)) %>% select(-N)

df = df %>% mutate("z_score" = qnorm(1-(P/2))) # se as beta/z
df = df %>% mutate("se" = sqrt((beta/z_score)^2))
df = df %>% rename('effect_allele' = A1, 'other_allele' = A2, 'pval'=P)

# filter P <= 0.5
df = df %>% filter(pval <= 0.5)

# write to file
write_tsv(df,"ms_prs_snps.tsv")

# prepare pheno file

# read in phenotype data
source_of_report_data = read_tsv("/data/Wolfson-UKBB-Dobson/ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G35 (multiple sclerosis).0.0` = col_character(),
  EID = col_double()
  ))

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))

source_of_report_data = source_of_report_data %>% mutate("FID" = EID, "IID" = EID) %>%
select(FID,IID,MS_status)

# exclude participants who have withdrawn
withdrawn = read_tsv("/data/Wolfson-UKBB-Dobson/helper_progs_and_key/excluded_indivs",col_names=FALSE)
source_of_report_data = source_of_report_data %>% filter(!IID %in% withdrawn$X1)

# exclude related individuals
ms = source_of_report_data %>% filter(MS_status==1)
kin = read_table2("/data/Wolfson-UKBB-Dobson/helper_progs_and_key/ukb43101_rel_s488282.dat")
kin_to_excl = kin %>% filter(Kinship>0.0884) %>% filter(!ID1 %in% withdrawn$X1) %>% filter(!ID2 %in% withdrawn$X1)

kin_to_excl = kin_to_excl %>%
mutate(ID1_MS = ifelse(ID1 %in% ms$IID,"yes","no")) %>%
mutate(exclusion_id = ifelse(ID1_MS =="yes",ID2,ID1))
source_of_report_data = source_of_report_data %>% filter(!IID %in% kin_to_excl$exclusion_id)

# exclude non-EUR
covar = read_tsv("ukb_pheno_21_01_22.tsv")
ethnicity = covar %>% select(EID,`Genetic ethnic grouping.0.0`)
eur = ethnicity %>% filter(`Genetic ethnic grouping.0.0` == "Caucasian")
source_of_report_data = source_of_report_data %>% filter(IID %in% eur$EID)

# create covar file
covar = covar %>% select(
  EID,
  EID,
  `Age at recruitment.0.0`,
  Sex.0.0,
  `Genetic principal components.0.1`,
  `Genetic principal components.0.2`,
  `Genetic principal components.0.3`,
  `Genetic principal components.0.4`
  )

colnames(covar) = c("IID","age","sex","PC1","PC2","PC3","PC4")
covar = covar %>% filter(!IID %in% withdrawn$X1)
covar = covar %>% mutate("FID" = IID) %>% dplyr::select(FID,colnames(covar))
covar = covar %>% filter(IID %in% source_of_report_data$IID)

write_tsv(source_of_report_data,"ms_pheno.tsv")
write_tsv(covar,"ms_covars.tsv")

````

## Convert genotypes to plink1 format

````unix
#! /bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:22
echo ${SGE_TASK_ID}

cd /data/home/hmy117/PRS_UKB_cross_ancestry/
module load plink/2.0-20200328
plink2 \
--pfile /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files/chr_${SGE_TASK_ID} \
--make-bed \
--out /data/scratch/hmy117/ms_chr_${SGE_TASK_ID} \
--keep ms_pheno.tsv \
--extract ms_prs_snps.tsv
````


## Run PRSice in UKB in chunks (1000 iterations of case-control sizes identical to Genes & Health)
````unix
qsub create_ukb_chunks.sh
qsub chunk_prs_nomhc.sh
qsub chunk_prs_mhc.sh
````


## Run PRSIce for whole cohort
````unix
qsub prs_all.sh
qsub prs_all_mhc.sh
````


# liftover to hg38
````unix
cd /data/home/hmy117/PRS_UKB_cross_ancestry/

awk 'NR>1{print "chr"$1,$3-1,$3,$2}' ./exports_0410/whole_cohort_nomhc_output_prs_r0.4.snp > gh_nomhc_hg38
/data/Wolfson-UKBB-Dobson/liftover/liftOver gh_nomhc_hg38 /data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz gh_nomhc_hg19 unlifted

awk 'NR>1{print "chr"$1,$3-1,$3,$2}' ./exports_0410/whole_cohort_mhc_output_prs_r0.05.snp > gh_mhc_hg38
/data/Wolfson-UKBB-Dobson/liftover/liftOver gh_mhc_hg38 /data/Wolfson-UKBB-Dobson/liftover/hg38ToHg19.over.chain.gz gh_mhc_hg19 unlifted
````

## Process results

````R
library(tidyverse)
setwd("/data/home/hmy117/PRS_UKB_cross_ancestry/")

# first process for whole cohort
outputs = list()
for(r in c("0.001","0.01","0.05","0.1","0.2","0.4","0.6")){
  output_nomhc = read_table2(paste0("/data/scratch/hmy117/whole_cohort_nomhc_output_prs_r",r,".summary"))
  output_mhc = read_table2(paste0("/data/scratch/hmy117/whole_cohort_mhc_output_prs_r",r,".summary"))
  output_nomhc$mhc = "No_MHC"
  output_mhc$mhc = "MHC"
  output = bind_rows(output_mhc,output_nomhc)
  output$clump_r2 = r   
  outputs[[length(outputs)+1]] = output
}
outputs = do.call("bind_rows",outputs)
write_csv(outputs,"whole_cohort_prs_results.csv")
whole_cohort_res = outputs


# process iterations
outputs = list()
for(i in c(1:1000)){
  outputs[[i]] = read_tsv(paste0("/data/scratch/hmy117/output_prs_",i,".summary"))
}
outputs = do.call("bind_rows",outputs)
outputs$mhc = "MHC"

n_extreme = outputs %>% filter(PRS.R2.adj<0.0108414) %>% nrow
n_total = outputs %>% nrow
pval = (n_extreme + 1) / (n_total + 1)

# repeat without mhc
outputs_nomhc = list()
for(i in c(1:1000)){
  outputs_nomhc[[i]] = read_tsv(paste0("/data/scratch/hmy117/output_prs_nomhc_",i,".summary"),col_types = "ccddddddddddd")
}
outputs_nomhc = do.call("bind_rows",outputs_nomhc)
outputs_nomhc$mhc = "No_MHC"
n_extreme = outputs_nomhc %>% filter(PRS.R2.adj<0.0146176) %>% nrow
n_total = outputs_nomhc %>% nrow
pval = (n_extreme + 1) / (n_total + 1)

# add in case n
outputs_nomhc$n_case = 42
outputs$n_case = 42

# add in whole cohort stats
whole_cohort_res$n_case = 2091

# get best score for MHC and non-MHC in whole cohort
whole_cohort_res = whole_cohort_res %>% group_by(mhc) %>% slice_max(PRS.R2.adj)

# combine
combo = bind_rows(outputs,outputs_nomhc,whole_cohort_res)
combo %>% group_by(mhc,n_case) %>% summarise(mean(PRS.R2.adj))

plot_dat = combo %>% group_by(mhc,n_case) %>% summarise(mean_r2 = mean(PRS.R2.adj),lower_ci = quantile(PRS.R2.adj,0.025),upper_ci = quantile(PRS.R2.adj,0.975))
plot_dat$cohort = "UKB"

# add in GH estimates
gh_prs = read_csv("./exports_0410/prs_summaries.csv") %>%
  filter(MHC %in% c("nomhc","mhc")) %>%
  group_by(MHC) %>%
  slice_max(PRS.R2.adj) %>%
  mutate(n_case = 42, cohort = "G&H", lower_ci=PRS.R2.adj, upper_ci=PRS.R2.adj) %>%
  dplyr::rename("mean_r2" = PRS.R2.adj,"mhc" = MHC) %>%
  dplyr::select(mhc,n_case,mean_r2,lower_ci,upper_ci,cohort) %>%
  mutate(mhc = ifelse(mhc == "mhc","MHC","No_MHC"))

plot_dat = plot_dat %>% bind_rows(gh_prs)

plot_dat = plot_dat %>%
  mutate(cohort = ifelse(cohort == "UKB" & n_case == 42,"UKB subset (EUR)",
  ifelse(cohort == "UKB" & n_case == 2091,"UKB all (EUR)","G&H (SAS)")))

p=ggplot(plot_dat,aes(mean_r2*100,cohort,col=mhc))+
  geom_point(size=3,position=ggstance::position_dodgev(height=0.3))+
  geom_errorbarh(mapping=aes(xmin=lower_ci*100,xmax=upper_ci*100,y=cohort),height=0,
    position = ggstance::position_dodgev(height=0.3))+
  theme_minimal()+
  scale_color_manual(values = c("purple","orange"),labels=c("MHC","No MHC"))+
  labs(x="% liability to MS \n explained by PRS",y="PRS",color = "MHC included/excluded")+
  scale_x_continuous(limits=c(0,10))+
  geom_vline(xintercept=0,alpha=0.1)+
  geom_vline(xintercept=1.08,alpha=0.4,color="purple")+
  geom_vline(xintercept=1.46,alpha=0.4,color="orange")

png("prs_comparison.png",res=300,units="in",width=6,height=3)
p
dev.off()

# get total n snps
snps = list()
for(i in c(1:22)){
  snps[[i]] = read_tsv(paste0("/data/scratch/hmy117/ms_chr_",i,".bim"),col_names=F)
}
snps = do.call("bind_rows",snps)


# compare snps used in PRS
gh_prs_snps_nomhc = read_table2("./exports_0410/whole_cohort_nomhc_output_prs_r0.4.snp",col_types = "cccdc")
gh_prs_snps_mhc = read_table2("./exports_0410/whole_cohort_mhc_output_prs_r0.05.snp",col_types = "cccdc")

# read in lifted coordinates
gh_prs_snps_nomhc_hg19 = read_table2("gh_nomhc_hg19",col_types = "cccc",col_names = F)
gh_prs_snps_mhc_hg19 = read_table2("gh_mhc_hg19",col_types = "cccc",col_names = F)

# combine
gh_prs_snps_nomhc = gh_prs_snps_nomhc %>%
  left_join(gh_prs_snps_nomhc_hg19 %>%
      dplyr::rename("SNP" = X4),
      by="SNP") %>%
      dplyr::select(-X1,-X2) %>%
      mutate(hg19_coord = paste0(CHR,":",X3))

gh_prs_snps_mhc = gh_prs_snps_mhc %>%
  left_join(gh_prs_snps_mhc_hg19 %>%
      dplyr::rename("SNP" = X4),
      by="SNP") %>%
      dplyr::select(-X1,-X2) %>%
      mutate(hg19_coord = paste0(CHR,":",X3))


snps = snps %>% mutate(chrpos =paste0(X1,":",X4))
gh_prs_snps_mhc %>% filter(hg19_coord %in% snps$chrpos)
gh_prs_snps_nomhc %>% filter(hg19_coord %in% snps$chrpos)

````
