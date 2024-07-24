sink("Druggable_gene_MR.txt")
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(stringr)
library(ieugwasr)
library(readxl) 
library(biomaRt)

# Get gene list -----------------------------------------------------------
#From Gaziano et al. - wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01310-z/MediaObjects/41591_2021_1310_MOESM3_ESM.xlsx
#Get ENSEMBL IDs from UniProt ID - do this on hpcapp01 as need internet access
#dat<-read_excel("CRC_TWAS/41591_2021_1310_MOESM3_ESM.xlsx",sheet=3,skip=11)
#ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
#genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"), values=dat$ENSEMBL_GENE_ID, mart=ensembl)
#genes$cisstart <- genes$start_position-1000000
#genes$cisend <- genes$start_position+1000000
#fwrite(genes,"CRC_TWAS/druggable_genes_with_cis_pos.csv")
genes<-fread("CRC_TWAS/druggable_genes_with_cis_pos.csv")


# GTEx data ---------------------------------------------------------------

as<-fread("v8_eQTL_european_associations/Adipose_Subcutaneous.allpairs.txt.gz")
av<-fread("v8_eQTL_european_associations/Adipose_Visceral_Omentum.allpairs.txt.gz")
cs<-fread("CRC_TWAS/bris_GTEx/Colon_Sigmoid.allpairs.csv")
ct<-fread("CRC_TWAS/bris_GTEx/Colon_Transverse.allpairs.csv")
ly<-fread("CRC_TWAS/bris_GTEx/Cells_EBV-transformed_lymphocytes.allpairs.csv")
wb<-fread("CRC_TWAS/bris_GTEx/Whole_Blood.allpairs.csv")
as$id<-paste(as$gene_id,"Adipose_Subcutaneous",sep=";")
av$id<-paste(av$gene_id,"Adipose_Visceral_Omentum",sep=";")
ly$id<-paste(ly$phenotype_id,"Cells_EBV-transformed_lymphocytes",sep=";")
ct$id<-paste(ct$phenotype_id,"Colon_Transverse",sep=";")
cs$id<-paste(cs$phenotype_id,"Colon_Sigmoid",sep=";")
wb$id<-paste(wb$phenotype_id,"Whole_Blood",sep=";")

as$X__index_level_0__<-0:(nrow(as)-1)
av$X__index_level_0__<-0:(nrow(av)-1)
colnames(ly)[1]<-"gene_id"
colnames(cs)[1]<-"gene_id"
colnames(ct)[1]<-"gene_id"
colnames(wb)[1]<-"gene_id"
eqtl<-rbind(as,av,ly,cs,ct,wb)
eqtl$gene<-substr(eqtl$gene_id,1,15)
eqtl<-eqtl[eqtl$gene %in% genes$ensembl_gene_id,]
#Find SNP rsIDs
snps<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)
eqtl<-eqtl %>%
  mutate(chr = str_remove_all(chr, "chr"))
# Limit to cis SNPs -------------------------------------------------------
exp<-merge(eqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp<-exp[exp$chr==exp$chromosome_name,]
exp<-exp[exp$variant_pos>=exp$cisstart,]
exp<-exp[exp$variant_pos<=exp$cisend,]
# MR ----------------------------------------------------------------------
exp1<-exp[exp$gene %in% unique(exp$gene)[1:200],]
exp2<-exp[exp$gene %in% unique(exp$gene)[201:400],]
exp3<-exp[exp$gene %in% unique(exp$gene)[401:600],]
exp4<-exp[exp$gene %in% unique(exp$gene)[601:800],]
exp5<-exp[exp$gene %in% unique(exp$gene)[801:1000],]
exp6<-exp[exp$gene %in% unique(exp$gene)[1001:1100],]
exp7<-exp[exp$gene %in% unique(exp$gene)[1101:length(unique(exp$gene))],]

exposure_dat1<-format_data(
  exp1,
  type = "exposure",
  phenotype_col="id",
  header = TRUE,
  id_col = "id",
  snp_col = "rs_id_dbSNP151_GRCh38p7",
  beta_col = "slope",
 se_col = "slope_se",
  effect_allele_col = "alt",
  eaf_col = "maf",
  other_allele_col = "ref",
  chr_col = "chr",
  pval_col = "pval",
  pos_col = "variant_pos",
  log_pval = FALSE
)  
exposure_dat2<-format_data(
  exp2,
  type = "exposure",
  header = TRUE,
  phenotype_col="id",
  id_col = "id",
  snp_col = "rs_id_dbSNP151_GRCh38p7",
  beta_col = "slope",
  eaf_col = "maf",
  effect_allele_col = "alt",
  se_col = "slope_se",
  other_allele_col = "ref",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "variant_pos",
  log_pval = FALSE
)  
exposure_dat3<-format_data(
  exp3,
  type = "exposure",
  header = TRUE,
  phenotype_col="id",
  snp_col = "rs_id_dbSNP151_GRCh38p7",
  beta_col = "slope",
  id_col = "id",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  chr_col = "chr",
  log_pval = FALSE,
  pos_col = "variant_pos"
)  
exposure_dat4<-format_data(
  type = "exposure",
  exp4,
  header = TRUE,
  phenotype_col="id",
  id_col = "id",
  snp_col = "rs_id_dbSNP151_GRCh38p7",
beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "variant_pos",
  log_pval = FALSE
)  
exposure_dat5<-format_data(
  exp5,
  header = TRUE,
  type = "exposure",
  phenotype_col="id",
  id_col = "id",
 snp_col = "rs_id_dbSNP151_GRCh38p7",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "variant_pos",
  log_pval = FALSE
)  
exposure_dat6<-format_data(
  exp6,
  type = "exposure",
  header = TRUE,
  phenotype_col="id",
  id_col = "id",
  snp_col = "rs_id_dbSNP151_GRCh38p7",
 se_col = "slope_se",
  beta_col = "slope",
  eaf_col = "maf",
  effect_allele_col = "alt",
other_allele_col = "ref",
pval_col = "pval",
 chr_col = "chr",
 pos_col = "variant_pos",
  log_pval = FALSE
)  
exposure_dat7<-format_data(
  exp7,
  type = "exposure",
  header = TRUE,
  phenotype_col="id",
  id_col = "id",
  snp_col = "rs_id_dbSNP151_GRCh38p7",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  pval_col = "pval",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  chr_col = "chr",
  pos_col = "variant_pos",
  log_pval = FALSE
)  
exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4,exposure_dat5,exposure_dat6,exposure_dat7)
exposure_dat$pval.exposure<-as.numeric(exposure_dat$pval.exposure)
exposure_dat<-exposure_dat[exposure_dat$pval.exposure<5*10^-8,]
exposure_dat$rsid<-exposure_dat$SNP
exposure_dat$id<-exposure_dat$id.exposure
exposure_dat$pval<-exposure_dat$pval.exposure

#clump
#Save file to clump server that has internet access
#remotes::install_github("MRCIEU/genetics.binaRies")
fwrite(exposure_dat,"CRC_TWAS/exposure_dat_temp_druggable.csv")
exposure_dat<-fread("CRC_TWAS/exposure_dat_temp_druggable.csv")
exposure_dat<-exposure_dat[!is.na(exposure_dat$mr_keep.exposure),]
kgsnps<-fread("1000GenomesReferenceFiles/EUR.bim")
exposure_dat<-exposure_dat[exposure_dat$rsid %in% kgsnps$V2,]
exposure_dat<- ld_clump(
  exposure_dat,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "1000GenomesReferenceFiles/EUR",
  pop="EUR")
fwrite(exposure_dat,"CRC_TWAS/clumped_exposure_data_druggable.csv")
exposure_dat<-fread("CRC_TWAS/clumped_exposure_data_druggable.csv")

# r2 and F stat -----------------------------------------------------------

exposure_dat$samplesize.exposure<-NA
exposure_dat$samplesize.exposure[grepl("Adipose_Subcutaneous",exposure_dat$exposure)]<-581
exposure_dat$samplesize.exposure[grepl("Adipose_Visceral_Omentum",exposure_dat$exposure)]<-469
exposure_dat$samplesize.exposure[grepl("Cells_EBV-transformed_lymphocytes",exposure_dat$exposure)]<-187
exposure_dat$samplesize.exposure[grepl("Colon_Sigmoid",exposure_dat$exposure)]<-318
exposure_dat$samplesize.exposure[grepl("Colon_Transverse",exposure_dat$exposure)]<-368
exposure_dat$samplesize.exposure[grepl("Whole_Blood",exposure_dat$exposure)]<-670

exposure_dat$num <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
exposure_dat$den <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure) + ((exposure_dat$se.exposure^2)*2*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure))
exposure_dat$pve <- exposure_dat$num/exposure_dat$den              
exposure_dat$F=((exposure_dat$pve)*(exposure_dat$samplesize.exposure-2))/(1-exposure_dat$pve)
f_r2<-dplyr::select(exposure_dat,id,pve,F)
f_r2<-aggregate(.~id, f_r2, sum)
colnames(f_r2)<-c("id","r2","F_stat")
exposure_dat<-merge(exposure_dat,f_r2)
f_r2<-dplyr::select(exposure_dat,exposure,r2,F_stat)
f_r2<-f_r2[!duplicated(f_r2),]
f_r2<-f_r2 %>% tidyr::separate(exposure, c("Gene", "Tissue"),sep=";")
fwrite(f_r2,"CRC_TWAS/f_stat_and_r2_druggable.csv")


# Running MR --------------------------------------------------------------
exposure_dat <- read_exposure_data(
  filename = "CRC_TWAS/clumped_exposure_data_druggable.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval",
  id_col = "id.exposure",
  phenotype_col="exposure")



results<-data.frame()

#overall
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Overall"
results<-rbind(results,res)

#Colon
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Colon"
results<-rbind(results,res)

#Distal
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Distal"
results<-rbind(results,res)

#Female
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/female_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Female"
results<-rbind(results,res)

#Male
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/male_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Male"
results<-rbind(results,res)

#Proximal
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Proximal"
results<-rbind(results,res)

#Rectal
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat<-harmonise_data(exposure_dat,outcome_dat)
res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
res$subtype<-"Rectal"
results<-rbind(results,res)

fwrite(results,"CRC_TWAS/MR_results_druggable_genome.csv")
sink()
