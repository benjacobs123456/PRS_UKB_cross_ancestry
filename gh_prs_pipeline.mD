# MS-PRS in a South Asian-ancestry cohort

This pipeline was used to generate PRS in UKB participants of EUR ancestry, and compare the performance with a similarly-constructed PRS in Genes & Health. Here is the pipeline used for the Genes & Health cohort.

###############################
# Step 1 - GWAS & LD QC       #
###############################

- Convert GWAS & LD reference to hg38
- Basic QC of GWAS sum stats

# liftover to hg38
````unix
cd /genesandhealth/red/BenJacobs/ms_prs/data/
awk 'NR>1{print "chr"$1,$2-1,$2,$3}' discovery_metav3.0.meta > ms_gwas_hg19
liftOver ms_gwas_hg19 hg19ToHg38.over.chain.gz ms_gwas_hg38 unlifted

awk '{print "chr"$1,$4-1,$4,$2}' /genesandhealth/library-green/sanger/LD_reference_panel/1kg/EUR503/plink_1000G_EUR503_SNPs_geno0.05_MAF0.001_biallelic.bim > 1kg_eur_ref_hg19.bed

liftOver 1kg_eur_ref_hg19.bed hg19ToHg38.over.chain.gz 1kg_eur_ref_hg38 unlifted
````

# filter in R
````R
library(tidyverse)
setwd("/genesandhealth/red/BenJacobs/ms_prs/data/")
# ms gwas
ms_hg38 = read_table("ms_gwas_hg38",col_names=FALSE)
ms_hg19 = read_table("discovery_metav3.0.meta")
ms_hg19 = ms_hg19 %>% filter(SNP %in% ms_hg38$X4) %>% filter(!is.na(P))
ms_lifted_hg38 = ms_hg19 %>% dplyr::select(-BP) %>% left_join(ms_hg38 %>% dplyr::select(X3,X4) %>% rename(BP = X3, SNP = X4),by="SNP")
ms_lifted_hg38 = ms_lifted_hg38 %>% dplyr::select(-N) %>% dplyr::select(CHR,BP,SNP,A1,A2,OR,P) %>% mutate(beta = log(OR))
ms_lifted_hg38 = ms_lifted_hg38 %>% mutate(SNP = paste0(CHR,":",BP))

# filter to biallelic SNVs
ms_lifted_hg38 = ms_lifted_hg38 %>% filter(nchar(A1)==1 & nchar(A2)==1)

# remove overlapping positions
ms_lifted_hg38 = ms_lifted_hg38 %>% distinct(SNP,.keep_all=TRUE)

# filter NAs
ms_lifted_hg38 = ms_lifted_hg38 %>% filter(!is.na(P))

# filter out palindromic SNPs
ms_lifted_hg38 = ms_lifted_hg38 %>%
filter(!(A1 == "T" & A2 =="A") & !(A1=="A" & A2=="T")) %>%
filter(!(A1 == "C" & A2 =="G") & !(A1=="G" & A2=="C"))

# calculate beta and se
ms_lifted_hg38 = ms_lifted_hg38 %>% mutate("z_score" = qnorm(1-(P/2))) # se as beta/z
ms_lifted_hg38 = ms_lifted_hg38 %>% mutate("se" = sqrt((beta/z_score)^2))
ms_lifted_hg38 = ms_lifted_hg38 %>% rename('effect_allele' = A1, 'other_allele' = A2, 'pval'=P)
ms_lifted_hg38 = ms_lifted_hg38 %>% filter(pval <= 0.5)

# write to file

write_tsv(ms_lifted_hg38,"ms_gwas_sumstats_hg38.tsv")

# 1kg
kg_hg38 = read_table("1kg_eur_ref_hg38",col_names=FALSE,col_types = cols_only(
  X1 = col_character(),
  X2 = col_double(),
  X3 = col_double(),
  X4 = col_character()
)
)
kg_hg19 = read_table("/genesandhealth/library-green/sanger/LD_reference_panel/1kg/EUR503/plink_1000G_EUR503_SNPs_geno0.05_MAF0.001_biallelic.bim",col_names=FALSE,col_types = cols_only(
  X1 = col_double(),
  X2 = col_character(),
  X3 = col_double(),
  X4 = col_double(),
  X5 = col_character(),
  X6 = col_character()
))
kg_hg19 = kg_hg19 %>% filter(X2 %in% kg_hg38$X4)

kg_lifted_hg38 = kg_hg19 %>%
rename("SNP" = X2) %>%
left_join(kg_hg38 %>% rename("BP" = X3, "SNP" = X4) %>% select(BP,SNP),by="SNP" ) %>%
select(-X4) %>% rename("X4" = BP) %>% mutate(SNP = paste0(X1,":",X4)) %>% select(X1,SNP,X3,X4,X5,X6)
write_tsv(kg_hg19,"filtered_snps.tsv")
write_tsv(kg_lifted_hg38,"EUR_1kg_hg38.bim",col_names=FALSE)
````

# filter 1kg plink files to lifted snps
````unix
plink --bfile /genesandhealth/library-green/sanger/LD_reference_panel/1kg/EUR503/plink_1000G_EUR503_SNPs_geno0.05_MAF0.001_biallelic \
--extract filtered_snps.tsv \
--make-bed \
--out EUR_1kg

# convert to hg38
plink --bed EUR_1kg.bed \
--bim EUR_1kg_hg38.bim \
--fam EUR_1kg.fam \
--make-bed \
--out EUR_1kg_hg38_bfile

awk '{print $2}' EUR_1kg_hg38_bfile.bim | uniq -d > duplicate_snps

# deduplicate 1kg reference SNPs
plink --bfile EUR_1kg_hg38_bfile \
--exclude duplicate_snps \
--make-bed \
--out EUR_1kg_hg38_bfile_unique
````

###############################
# Step 2 - Target QC	      #
###############################
# Prior to this QC extensive QC has been run by Sanger

```unix
cd /genesandhealth/red/BenJacobs/ms_prs/outputs
rm good_vars_maf_0.01_info_0.7
# filter by INFO > 0.7 and MAF > 0.01
for i in {1..22}
do
	awk '{if($7>0.01 && $10 > 0.7) print $1,$7,$9,$10}' /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/topmed-r2_merged_version03/chr$i\_variant_info.tab >> good_vars_maf_0.01_info_0.7
echo "finishied copying list of good vars for chrom $i\ "
done

cd /genesandhealth/red/BenJacobs/ms_prs/outputs
# convert from BGEN to PLINK1 binary and perform basic QC

# copy big pgen
cp /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/topmed-r2_merged_version03/chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2.8bit.sorted* /home/ivm/

cd /genesandhealth/red/BenJacobs/ms_prs/outputs
plink2 --pfile /home/ivm/chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2.8bit.sorted \
--extract good_vars_maf_0.01_info_0.7 \
--out /home/ivm/chr_all_filtered \
--maf 0.01 \
--hwe 1e-10 \
--geno 0.1 \
--mind 0.1 \
--chr 1-22 \
--make-pgen
````

# convert IIDs
````R
# modify psam file
geno = read_tsv("/home/ivm/chr_all_filtered.psam")
geno = geno %>%
mutate(old_iid = IID) %>%
mutate(old_fid = 1) %>%
separate(IID,sep="_",into=c("IID","FID","other")) %>%
select(old_fid,old_iid,FID,IID)
write_tsv(geno,"./outputs/new_iids")
````

# rename variants as chr:pos and filter to MS PRS SNPs
# convert to plink1
# filter to MS GWAS SNPs
# rename IIDs
# rename SNPs as chr:pos
````unix
plink2 --pfile /home/ivm/chr_all_filtered \
--set-all-var-ids @:# \
--update-ids ./outputs/new_iids \
--extract ms_gwas_sumstats_hg38.tsv \
--out /home/ivm/post_qc_all_chrs_for_prs \
--make-bed
````



###############################
# Step 3 - Pheno & Covars     #
###############################

````R
library(readr)
library(dplyr)
library(tidyr)
setwd("/genesandhealth/red/BenJacobs/ms_prs/")

# pheno file
ms = read_tsv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version005_2022_06/1stoccurrence_3digitICD10_1SNOMEDto1ICD10/2022_06_version005_3digitICD10_1to1.txt")
ms = ms %>% select(PseudoNHSNumber,ICD10__G35) %>% filter(ICD10__G35==1)

link_file = read_tsv("/genesandhealth/library-red/genesandhealth/2022_05_12_pseudoNHS_oragene_withmissing_DEIDENTIFIED.txt")
link_file = link_file %>%
  mutate(MS_status = ifelse(PseudoNHS %in% ms$PseudoNHSNumber,"1","0")) %>%
  mutate(FID=OrageneID) %>%
  select(FID,OrageneID,MS_status) %>%
  rename("IID" = OrageneID) %>%
distinct(IID,.keep_all=T) %>%
select(-FID)
write_tsv(link_file,"./outputs/ms_pheno.tsv")

# modify covar file
covars = read_tsv("/genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/GNH.44190.noEthnicOutliers.covariates.20PCs.withS1QST_diabetes.txt")
covars = covars %>% select(-2) %>% rename("FID"=OrageneID) %>% mutate(IID = FID) %>% select(FID,IID,3:ncol(covars))
covars = covars %>% select(1:3,6,10:13) %>% rename("sex" = S1QST_Gender,"age" = AgeAtRecruitment) %>% distinct(IID,.keep_all=T) %>%
select(-FID)
write_tsv(covars,"./outputs/covars.tsv")
````

###############################
# Step 4 - PRSice 	      #
###############################

# run PRSice

````unix
cd /genesandhealth/red/BenJacobs/ms_prs/
for r in 0.001 0.01 0.05 0.1 0.2 0.4 0.6
do

	# After first run there are some duplicates > re-run with extract
	PRSice \
--base ./data/ms_gwas_sumstats_hg38.tsv \
--chr CHR \
--A1 effect_allele \
--A2 other_allele \
--stat beta \
--beta \
--snp SNP \
--bp BP \
--pvalue pval \
--binary-target T \
--score sum \
--target /home/ivm/post_qc_all_chrs_for_prs \
--ignore-fid \
--pheno ./outputs/ms_pheno.tsv \
--ld ./data/EUR_1kg_hg38_bfile_unique \
--clump-r2 $r \
--missing CENTER \
--bar-levels 5e-8,1e-8,5e-7,1e-7,5e-6,1e-6,5e-5,1e-5,5e-4,1e-4,5e-3,1e-3,5e-2,1e-2,0.5,0.1,1 \
--fastscore \
--prevalence 0.002 \
--cov ./outputs/covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--x-range chr6:25000000-35000000 \
--print-snp \
--extract /home/ivm/whole_cohort_nomhc_output_prs_r$r\.valid \
--thread 64 \
--out /home/ivm/whole_cohort_nomhc_output_prs_r$r &
done
````


# repeat with MHC
````unix
cd /genesandhealth/red/BenJacobs/ms_prs/
for r in 0.001 0.01 0.05 0.1 0.2 0.4 0.6
do

	# After first run there are some duplicates > re-run with extract
	PRSice \
--base ./data/ms_gwas_sumstats_hg38.tsv \
--chr CHR \
--A1 effect_allele \
--A2 other_allele \
--stat beta \
--beta \
--snp SNP \
--bp BP \
--pvalue pval \
--binary-target T \
--score sum \
--target /home/ivm/post_qc_all_chrs_for_prs \
--ignore-fid \
--pheno ./outputs/ms_pheno.tsv \
--ld ./data/EUR_1kg_hg38_bfile_unique \
--clump-r2 $r \
--missing CENTER \
--bar-levels 5e-8,1e-8,5e-7,1e-7,5e-6,1e-6,5e-5,1e-5,5e-4,1e-4,5e-3,1e-3,5e-2,1e-2,0.5,0.1,1 \
--fastscore \
--prevalence 0.002 \
--cov ./outputs/covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--print-snp \
--thread 64 \
--extract /home/ivm/whole_cohort_mhc_output_prs_r$r\.valid \
--out /home/ivm/whole_cohort_mhc_output_prs_r$r &
done


# repeat with JUST MHC
cd /genesandhealth/red/BenJacobs/ms_prs/
for r in 0.001 0.01 0.05 0.1 0.2 0.4 0.6
do
	awk '{if($2==6 && $3>25000000 && $3<35000000) print}' /home/ivm/whole_cohort_mhc_output_prs_r$r\.valid > just_mhc_snps_r$r

	# After first run there are some duplicates > re-run with extract
	PRSice \
--base ./data/ms_gwas_sumstats_hg38.tsv \
--chr CHR \
--A1 effect_allele \
--A2 other_allele \
--stat beta \
--beta \
--snp SNP \
--bp BP \
--pvalue pval \
--binary-target T \
--score sum \
--target /home/ivm/post_qc_all_chrs_for_prs \
--ignore-fid \
--pheno ./outputs/ms_pheno.tsv \
--ld ./data/EUR_1kg_hg38_bfile_unique \
--clump-r2 $r \
--missing CENTER \
--bar-levels 5e-8,1e-8,5e-7,1e-7,5e-6,1e-6,5e-5,1e-5,5e-4,1e-4,5e-3,1e-3,5e-2,1e-2,0.5,0.1,1 \
--fastscore \
--prevalence 0.002 \
--cov ./outputs/covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--print-snp \
--thread 64 \
--extract just_mhc_snps_r$r \
--out /home/ivm/whole_cohort_JUST_mhc_output_prs_r$r &
done
````

# sense check in PLINK
````unix
cd /genesandhealth/red/BenJacobs/ms_prs/
awk '{print $3,$4,$8}' ./data/ms_gwas_sumstats_hg38.tsv > prs_snps
plink2 \
--bfile /home/ivm/post_qc_all_chrs_for_prs \
--score prs_snps \
--extract just_mhc_snps_r0.001 \
--exclude plink2.rmdup.mismatch

--clump-r2 $r \
--missing CENTER \
--bar-levels 5e-8,1e-8,5e-7,1e-7,5e-6,1e-6,5e-5,1e-5,5e-4,1e-4,5e-3,1e-3,5e-2,1e-2,0.5,0.1,1 \
--fastscore \
--prevalence 0.002 \
--cov ./outputs/covars.tsv \
--cov-factor sex \
--ultra \
--seed 123456 \
--print-snp \
--thread 64 \
--extract just_mhc_snps_r$r \
--out /home/ivm/whole_cohort_JUST_mhc_output_prs_r$r &
done
````



############################################
Analyse outputs

############################################
````R
library(tidyverse)
library(Hmisc)
library(ggplot2)
setwd("/genesandhealth/red/BenJacobs/ms_prs/")

# find best prs
prs_outputs = list()
for(r in c(0.001,0.01,0.05,0.1,0.2,0.4,0.6)){
  for(mhc in c("mhc","nomhc","JUST_mhc")){
    x = read_table(paste0("/home/ivm/whole_cohort_",mhc,"_output_prs_r",r,".prsice"))  
    x$MHC = mhc
    x$clumping_r = r
    prs_outputs[[length(prs_outputs)+1]] = x
  }
}
prs_outputs = do.call("bind_rows",prs_outputs)


p0=ggplot(prs_outputs %>% mutate(MHC = ifelse(MHC == "mhc","MHC PRS",ifelse(MHC=="nomhc","Non-MHC PRS","Just MHC"))),
       aes(factor(Threshold),
           -log10(P),
           fill=factor(clumping_r)))+
  geom_col(position = position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+
  labs(x="P value threshold",y="-log10(P)",fill=Clumping~R^2)+
  facet_wrap(~MHC)+
  theme(axis.text.x = element_text(angle=45))
p0

save_plot = function(x,height=6,width=6){
  filename_to_write = paste0("./outputs/exports/",x,".png")
  png(filename_to_write,res=300,units="in",height=height,width=width)
  print(eval(parse(text=x)))
  dev.off()
}

save_plot("p0",5,8)
prs_outputs %>%
  group_by(MHC) %>%
  slice_min(P)

prs_outputs2 = list()
for(r in c(0.001,0.01,0.05,0.1,0.2,0.4,0.6)){
  for(mhc in c("mhc","nomhc","JUST_mhc")){
    x = read_table(paste0("/home/ivm/whole_cohort_",mhc,"_output_prs_r",r,".summary"))  
    x$MHC = mhc
    x$clumping_r = r
    prs_outputs2[[length(prs_outputs2)+1]] = x
  }
}
prs_outputs2 = do.call("bind_rows",prs_outputs2)
write_csv(prs_outputs2,"./outputs/exports/prs_summaries.csv")
prs_outputs2 %>%
  group_by(MHC) %>%
  slice_min(P) %>%
  select(PRS.R2.adj)

prs_outputs2 %>%
  group_by(MHC) %>%
  slice_min(order_by = P)%>%
  select(Num_SNP,P,PRS.R2.adj,clumping_r,Threshold)


# read in best prs
prs_mhc = read_table("/home/ivm/whole_cohort_mhc_output_prs_r0.05.best")
prs_mhc_snps = read_table("/home/ivm/whole_cohort_mhc_output_prs_r0.05.snp",col_types = "dcddd")
prs_mhc_snps %>% filter(CHR==6 & BP>25000000 & BP < 35000000) %>% arrange(P)
prs_mhc = prs_mhc %>% rename("PRS_MHC" = PRS)
prs_nomhc = read_table("/home/ivm/whole_cohort_nomhc_output_prs_r0.4.best")
prs_nomhc_snps = read_table("/home/ivm/whole_cohort_nomhc_output_prs_r0.4.snp",col_types = "dcddd")
prs_nomhc_snps %>% filter(CHR==6 & BP>25000000 & BP < 35000000) %>% arrange(P)
prs_nomhc = prs_nomhc %>% rename("PRS_No_MHC" = PRS)


# read in pheno data
pheno = read_table("./outputs/ms_pheno.tsv")

# read in covars
covars = read_table("./outputs/covars.tsv")

# combine
pheno = pheno %>%
  filter(IID %in% covars$IID) %>%
  left_join(covars,by=c("IID")) %>%
  left_join(prs_mhc,by=c("IID")) %>%
  left_join(prs_nomhc,by=c("IID")) %>%
  na.omit()

# read in plink sense check
just_mhc_plink = read_tsv("plink2.sscore")
just_mhc_plink = just_mhc_plink %>% select(IID,SCORE1_AVG)
just_mhc_plink = just_mhc_plink %>% left_join(pheno,by="IID")
ggplot(just_mhc_plink,aes(SCORE1_AVG,fill=factor(MS_status)))+geom_density(alpha=0.5)

model1 = glm(data = just_mhc_plink,
             MS_status ~ age + sex + PC1 + PC2 + PC3 + PC4,
             family=binomial(link="logit"))
model2 = glm(data = just_mhc_plink,
             MS_status ~ age + sex + PC1 + PC2 + PC3 + PC4 + SCORE1_AVG,
             family=binomial(link="logit"))

rcompanion::nagelkerke(model2,model1)
# basic demographics
covar_full = read_tsv("/genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/GNH.44190.noEthnicOutliers.covariates.20PCs.withS1QST_diabetes.txt")
covar_full = covar_full %>% select(OrageneID,inferred_elghPCA) %>% rename(IID = OrageneID)
pheno = pheno %>% left_join(covar_full,by="IID")

summarise_cat = function(x){
  pheno %>%
    group_by(MS_status) %>%
    count(.data[[x]]) %>%
    mutate(prop = n/sum(n)) %>%
    mutate(value = paste0(n," (",round(prop*100,1),"%)"))%>%
    select(1,5) %>%
    mutate(var = x)
}
x1=pheno %>% count(MS_status) %>%
  mutate(var = "n") %>%
  select(1,3,2) %>%
  rename("value" = n) %>%
  mutate(value = as.character(value))
x2=summarise_cat("sex")
x3=summarise_cat("inferred_elghPCA")
x4=pheno %>% group_by(MS_status) %>% summarise(mean(age),sd(age)) %>%
  mutate(value = paste0(round(`mean(age)`,1)," (",round(`sd(age)`,1),")")) %>%
  mutate(var = "age") %>%
  select(MS_status,value,var)

demogs = bind_rows(x1,x2,x3,x4)

write_csv(demogs,"./outputs/exports/demographics.csv")

# normalise PRS
normalise = function(x){
  (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
}
pheno = pheno %>%
  mutate(PRS_MHC_Z = normalise(PRS_MHC)) %>%
  mutate(PRS_No_MHC_Z = normalise(PRS_No_MHC))

# manual nagelkerke
pheno %>% glimpse
null_model = glm(data=pheno,
                 MS_status~age+
                   factor(sex)+PC1+PC2+PC3+PC4,
                 family=binomial(link="logit"))

hlaprs_model = glm(data=pheno,
                 MS_status~age+
                   factor(sex)+PC1+PC2+PC3+PC4+
                   PRS_MHC_Z,
                 family=binomial(link="logit"))

nohla_prs_model = glm(data=pheno,
                   MS_status~age+
                     factor(sex)+PC1+PC2+PC3+PC4+
                     PRS_No_MHC_Z,
                   family=binomial(link="logit"))

library(rcompanion)
print(nagelkerke(hlaprs_model,null_model))
print(nagelkerke(nohla_prs_model,null_model))
print(nagelkerke(hlaprs_model,nohla_prs_model))


# basic histograms
p1=ggplot(
  pheno,
  aes(PRS_MHC_Z,fill=factor(MS_status))
)+
  geom_density(alpha=0.5)+
  theme_minimal()+
  scale_fill_brewer(palette="Paired",labels=c("Controls","MS"))+
  labs(x="PRS (with MHC)",fill="MS status")+
  ggtitle("MHC PRS")

p2=ggplot(
  pheno,
  aes(PRS_No_MHC,fill=factor(MS_status))
)+
  geom_density(alpha=0.5)+
  theme_minimal()+
  scale_fill_brewer(palette="Paired",labels=c("Controls","MS"))+
  labs(x="PRS (without MHC)",fill="MS status")+
  ggtitle("Non-MHC PRS")

# quartile plot
pheno$prs_quartile = cut2(pheno$PRS_MHC,g=4)
levels(pheno$prs_quartile) = c(1:4)
prs_quartile_model = glm(data=pheno,
                         MS_status~age+
                           factor(sex)+PC1+PC2+PC3+PC4+
                           prs_quartile,
                         family=binomial(link="logit"))

tbl = data.frame(summary(prs_quartile_model)$coefficients[-c(1:7),])
tbl$quartile=c(2:4)
tbl$or=exp(tbl$Estimate)
tbl$lower_ci=exp(tbl$Estimate-1.96*tbl$Std..Error)
tbl$upper_ci=exp(tbl$Estimate+1.96*tbl$Std..Error)
plot_dat = tbl %>% select(quartile,Estimate,Std..Error)
rownames(plot_dat) = NULL
plot_dat = rbind(plot_dat,c(1,0,0))
plot_dat$mhc = "MHC"

# repeat without MHC
pheno$prs_quartile_nomhc = cut2(pheno$PRS_No_MHC,g=4)
levels(pheno$prs_quartile_nomhc) = c(1:4)
prs_quartile_model = glm(data=pheno,
                         MS_status~age+
                           factor(sex)+PC1+PC2+PC3+PC4+
                           prs_quartile_nomhc,
                         family=binomial(link="logit"))

tbl = data.frame(summary(prs_quartile_model)$coefficients[-c(1:7),])
tbl$quartile=c(2:4)
tbl$or=exp(tbl$Estimate)
tbl$lower_ci=exp(tbl$Estimate-1.96*tbl$Std..Error)
tbl$upper_ci=exp(tbl$Estimate+1.96*tbl$Std..Error)
plot_dat2 = tbl %>% select(quartile,Estimate,Std..Error)
rownames(plot_dat2) = NULL
plot_dat2 = rbind(plot_dat2,c(1,0,0))
plot_dat2$mhc = "No MHC"

plot_dat = plot_dat %>% bind_rows(plot_dat2)

p3=ggplot(plot_dat,aes(factor(quartile),Estimate,fill=mhc))+
  geom_errorbar(aes(x=factor(quartile),ymin=Estimate - 1.96*`Std..Error`,ymax=Estimate + 1.96*`Std..Error`,width=0.2),position = position_dodge(width=0.5))+
  geom_point(size=3,shape=22,color="black",position = position_dodge(width=0.5))+
  theme_minimal()+
  labs(x="PRS quartile",y="Log(OR) for MS (vs lowest quartile)",fill="PRS with/out MHC")+
  geom_hline(yintercept = 0,alpha=0.3)+
  scale_fill_brewer(palette="Paired")


# roc analysis
library(ROCR)
make_prediction = function(x){
  predict(x,newdata=pheno,type="response")
}

preds = data.frame("hlaprs_model"=make_prediction(hlaprs_model),
                   "null_model"=make_prediction(null_model),
                   "prs_quartile_mhc"=pheno$prs_quartile,
                   "MS_status"=pheno$MS_status)

preds_summary = preds %>%
  group_by(prs_quartile_mhc) %>%
  summarise(
    "prs" = mean(hlaprs_model,na.rm=TRUE),
    "null" = mean(null_model,na.rm=TRUE)) %>%
  na.omit()

obs_risk = pheno %>%
  group_by(prs_quartile) %>%
  count(MS_status) %>%
  mutate(prop = n/sum(n)) %>%
  filter(MS_status==1) %>%
  select(prs_quartile,prop) %>%
  rename("prs_quartile_mhc" = prs_quartile)

pheno %>%
  group_by(prs_quartile) %>%
  count(MS_status) %>%
  mutate(prop = n/sum(n), total = sum(n), percent = prop * 100)

preds_summary = preds_summary %>% left_join(obs_risk,by="prs_quartile_mhc")

library(reshape2)
pred_df = melt(preds_summary,id="prs_quartile_mhc")
pred_df$prs_quartile = factor(pred_df$prs_quartile)

p4 = ggplot(pred_df,
                    aes(prs_quartile,value,col=variable,group=variable))+
  geom_point(col="black",shape=22,alpha=0.3)+
  geom_line()+
  labs(col="Source of risk estimate",x="PRS quartile",y="MS risk (probability scale)")+
  scale_color_brewer(palette="Set2",labels=c("MHC PRS model","Null model (Age, Sex, PCs)","Observed risk"))+
  theme_minimal()+
  scale_y_continuous(limits=c(0,0.002))+
  ggtitle("MHC PRS")


# repeat for non-MHC
preds = data.frame("nohla_prs_model"=make_prediction(nohla_prs_model),
                   "null_model"=make_prediction(null_model),
                   "prs_quartile"=pheno$prs_quartile_nomhc,
                   "MS_status"=pheno$MS_status)

preds_summary = preds %>%
  group_by(prs_quartile) %>%
  summarise(
    "prs" = mean(nohla_prs_model,na.rm=TRUE),
    "null" = mean(null_model,na.rm=TRUE)) %>%
  na.omit()

pheno %>%
  group_by(prs_quartile_nomhc) %>%
  count(MS_status) %>%
  mutate(prop = n/sum(n), total = sum(n), percent = prop * 100)

obs_risk = pheno %>%
  group_by(prs_quartile) %>%
  count(MS_status) %>%
  mutate(prop = n/sum(n)) %>%
  filter(MS_status==1) %>%
  select(prs_quartile,prop)

preds_summary = preds_summary %>% left_join(obs_risk,by="prs_quartile")

library(reshape2)
pred_df = melt(preds_summary,id="prs_quartile")
pred_df$prs_quartile = factor(pred_df$prs_quartile)

p5 = ggplot(pred_df,
            aes(prs_quartile,value,col=variable,group=variable))+
  geom_point(col="black",shape=22,alpha=0.3)+
  geom_line()+
  labs(col="Source of risk estimate",x="PRS quartile",y="MS risk (probability scale)")+
  scale_color_brewer(palette="Set2",labels=c("Non-MHC PRS model","Null model (Age, Sex, PCs)","Observed risk"))+
  theme_minimal()+
  scale_y_continuous(limits=c(0,0.002))+
  ggtitle("Non-MHC PRS")


# auc
preds = data.frame("nohla_prs_model"=make_prediction(nohla_prs_model),
                   "hla_prs_model"=make_prediction(hlaprs_model),
                   "null_model"=make_prediction(null_model),
                   "MS_status"=pheno$MS_status)
preds = preds %>% na.omit()
preds$MS=factor(preds$MS_status)
predictions = prediction(list(preds$hla_prs_model,preds$nohla_prs_model,preds$null_model),
                         list(preds$MS,preds$MS,preds$MS))
roc.perf = ROCR::performance(predictions, measure = "tpr", x.measure = "fpr")

hla_prs_model_auc = data.frame(x=roc.perf@x.values[[1]],y=roc.perf@y.values[[1]],model="HLA PRS model")
nohla_prs_model_auc = data.frame(x=roc.perf@x.values[[2]],y=roc.perf@y.values[[2]],model="Non-HLA PRS model")
null_model_auc = data.frame(x=roc.perf@x.values[[3]],y=roc.perf@y.values[[3]],model="Null model")

df = bind_rows(hla_prs_model_auc,nohla_prs_model_auc,null_model_auc)
auc.perf = performance(predictions, measure = "auc")
auc.perf@y.values

p6=ggplot(df,aes(x,y,col=model))+
  geom_line()+
  scale_color_brewer(palette="Set2")+
  labs(x="False Positive Rate",y="True Positive Rate",col="Model")+
  theme_classic()+
  geom_abline()+
  annotate("text",x=0.8,y=0.4,label="AUC",hjust=0)+
  annotate("text",x=1,y=0.3,label=paste0("MHC PRS model: ",round(auc.perf@y.values[[1]],3)),hjust=1)+
  annotate("text",x=1,y=0.2,label=paste0("Non-MHC PRS model: ",round(auc.perf@y.values[[2]],3)),hjust=1)+
  annotate("text",x=1,y=0.1,label=paste0("Null model: ",round(auc.perf@y.values[[3]],3)),hjust=1)

png("./outputs/exports/plots.png",res=300,height=8,width=12,units="in")
gridExtra::grid.arrange(p1+ggtitle("A"),p2+ggtitle("B"),p3+ggtitle("C"),p6+ggtitle("D"),p4+ggtitle("E"),p5+ggtitle("F"))
dev.off()



# compare PCs

compare_pcs = function(x){
ms = pheno[pheno$MS_status==1,][[x]]
cont = pheno[pheno$MS_status==0,][[x]]
t.test(ms,cont)$p.value
}
sapply(
  c("PC1","PC2","PC3","PC4"),
  compare_pcs
)


# correlation with HLA
hla_calls = read_tsv("../hla/outputs/hla_calls_sas_per_allele_counts.tsv")
pheno = pheno %>% left_join(hla_calls %>% rename("IID" = sample.id) %>% mutate(IID = as.numeric(IID)),by="IID")

# look at association with MS
overall_results = list()
hla_colnames = colnames(pheno)[20:ncol(pheno)]
for(i in c(1:length(hla_colnames))){
  message("Doing allele ",i," of ",length(hla_colnames))
  hla_allele_to_test = hla_colnames[i]
  af = sum(pheno[[hla_allele_to_test]],na.rm=TRUE)/(nrow(pheno)*2)
  if(af <= 0.05){
    message("allele is rare. skipping")
    results = c(hla_allele_to_test,NA,NA,NA,NA) %>%
      t() %>%
      data.frame()
    colnames(results) = c("HLA_allele","Beta","SE","Z","P")
  } else {
    message("allele is common. Processing")
    model = glm(
    data=pheno,
    MS_status ~ age + sex + PC1 + PC2 + PC3 + PC4 +
      pheno[[hla_allele_to_test]],
    family=binomial(link="logit")
  )

  results = c(hla_allele_to_test,summary(model)$coefficients[nrow(summary(model)$coefficients),]) %>%
    t() %>%
    data.frame()
  colnames(results) = c("HLA_allele","Beta","SE","Z","P")
  }

  overall_results[[i]] = results
}

overall_results = do.call("bind_rows",overall_results)
overall_results$P = as.numeric(overall_results$P)
overall_results$Beta = as.numeric(overall_results$Beta)
overall_results$SE = as.numeric(overall_results$SE)
overall_results$fdr = p.adjust(overall_results$P,method="fdr")
overall_results = overall_results %>% na.omit()
write_tsv(overall_results,paste0("./outputs/",pheno_to_test,"_hla_association_results.tsv"))

overall_results %>% arrange(P)
plot_dat = overall_results %>% arrange(desc(Beta))
plot_dat$HLA_allele = factor(plot_dat$HLA_allele,ordered = TRUE,levels = plot_dat$HLA_allele)

p=ggplot(plot_dat,
         aes(Beta,HLA_allele))+
  geom_point()+
  geom_errorbarh(mapping = aes(xmin = Beta - 1.96*SE, xmax = Beta + 1.96*SE,height=0.3))+
  ggtitle("MS")

png(filename = paste0("./outputs/",pheno_to_test,"_hla_association_results.png"),res=300,units="in",height=8,width=8)
print(p)
dev.off()
}


# power calcs
library(tidyverse)

control_n = 40490
case_n = 42

get_power = function(z_diff){
  n_perm = 1000
  pvals = list()
  for(i in c(1:n_perm)){
    message(i)
    message("Z diff:",z_diff)
    control_prs = rnorm(n = control_n,mean = 0,sd=1)
    case_prs = rnorm(n = case_n,mean = z_diff,sd=1)
    df = data.frame(PRS = case_prs,cc_status = 1) %>%
      bind_rows(
        data.frame(PRS = control_prs,cc_status = 0)    
      )
    model = glm(data=df,cc_status ~ PRS,family=binomial(link="logit"))
    pval = summary(model)$coefficients[2,4]
    pvals[[i]] = pval
  }
  pvals = unlist(pvals)
  power = sum(pvals<0.05) / length(pvals)
  return(power)
}


z_diffs = seq(0,3,by=0.5)

powers = sapply(z_diffs,get_power)

df = data.frame(Z = z_diffs,power = powers)
ggplot(df,aes(Z,power))+geom_point()+theme_minimal()+geom_smooth(se=F)
````
