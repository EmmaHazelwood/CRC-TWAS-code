sink("coloc_abf_risk_factors.txt")
library(coloc)
library(dplyr)
library(gwasglue)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(plinkbinr)


# Get analyses needed -----------------------------------------------------
res<-fread("CRC_TWAS/Coloc_results_abf_druggable.csv")
res<-res %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
length(unique(res$Gene))
list<-res[res$PP.H4.abf>0.8,]
list$gene<-substr(list$Gene,1,15)

# Read in exposure data ---------------------------------------------------
as<-fread("v8_eQTL_european_associations/Adipose_Subcutaneous.allpairs.txt.gz")
av<-fread("v8_eQTL_european_associations/Adipose_Visceral_Omentum.allpairs.txt.gz")
ly<-fread("CRC_TWAS/bris_GTEx/Cells_EBV-transformed_lymphocytes.allpairs.csv")
cs<-fread("CRC_TWAS/bris_GTEx/Colon_Sigmoid.allpairs.csv")
ct<-fread("CRC_TWAS/bris_GTEx/Colon_Transverse.allpairs.csv")
wb<-fread("CRC_TWAS/bris_GTEx/Whole_Blood.allpairs.csv")

as$X__index_level_0__<-0:(nrow(as)-1)
av$X__index_level_0__<-0:(nrow(av)-1)
colnames(ly)[1]<-"gene_id"
colnames(cs)[1]<-"gene_id"
colnames(ct)[1]<-"gene_id"
colnames(wb)[1]<-"gene_id"

as$gene<-substr(as$gene_id,1,15)
as$id<-paste(as$gene,"Adipose_Subcutaneous",sep=";")
av$gene<-substr(av$gene_id,1,15)
av$id<-paste(av$gene,"Adipose_Visceral_Omentum",sep=";")
ly$gene<-substr(ly$gene_id,1,15)
ly$id<-paste(ly$gene,"Cells_EBV-transformed_lymphocytes",sep=";")
cs$gene<-substr(cs$gene_id,1,15)
cs$id<-paste(cs$gene,"Colon_Sigmoid",sep=";")
ct$gene<-substr(ct$gene_id,1,15)
ct$id<-paste(ct$gene,"Colon_Transverse",sep=";")
wb$gene<-substr(wb$gene_id,1,15)
wb$id<-paste(wb$gene,"Whole_Blood",sep=";")

as$samplesize<-581
av$samplesize<-469
ly$samplesize<-187
cs$samplesize<-318
ct$samplesize<-368
wb$samplesize<-670

eqtl<-rbind(as,av,ly,cs,ct,wb)

snps<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)
eqtl<-eqtl %>%
  mutate(chr = str_remove_all(chr, "chr"))


# Limit to cis SNPs -------------------------------------------------------

eqtl<-eqtl[eqtl$gene %in% list$gene,]

genes<-fread("CRC_TWAS/druggable_genes_with_cis_pos.csv")
exp<-merge(eqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp<-exp[exp$chr==exp$chromosome_name,]
exp<-exp[exp$variant_pos>=exp$cisstart,]
exp<-exp[exp$variant_pos<=exp$cisend,]


# Perform coloc -----------------------------------------------------------

list$id<-paste(list$gene,list$Tissue,sep=";")
list<-dplyr::select(list,id)
list<-unique(list)

results<-data.frame()
for (a in list$id){
  print(a)
  dat1<-exp[exp$id==a,]
  
  #Format exposure data
  exposure_dat<-format_data(
    dat1,
    type = "exposure",
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
    log_pval = FALSE,
    samplesize_col = "samplesize"
  )  
  
  #BMI
  outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "data/CRC_risk_factors/ieu-b-40.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  samplesize_col="samplesize.exposure"
  )

  
  snps<-harmonise_data(exposure_dat,outcome_dat)
  snps<-snps[!is.na(snps$beta.exposure),]
  snps<-snps[!is.na(snps$beta.outcome),]
  snps<-snps[!is.na(snps$se.outcome),]
  snps<-snps[!is.na(snps$eaf.outcome),]
  snps<-snps[!is.na(snps$eaf.exposure),]
  
  df1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = snps$samplesize.exposure)
  df2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = snps$samplesize.outcome)
  
  res=coloc.abf(df1,df2)
  
  df3 <- data.frame(t(res$summary))
  df3$id<-a
  df3$risk_factor<-"BMI"
  results<-rbind(results,df3)
 
  #WHR
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = "data/CRC_risk_factors/fat-distn.giant.ukbb.meta-analysis.whr.combined_formatted.txt",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "Effect",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value"
  )
  outcome_dat$samplesize.outcome<-697734
  outcome_dat$p_one_sided<-outcome_dat$pval.outcome/2
  outcome_dat$z<-qnorm(outcome_dat$p_one_sided)
  outcome_dat$se.outcome<-as.numeric(abs(outcome_dat$beta.outcome/outcome_dat$z))
  
  snps<-harmonise_data(exposure_dat,outcome_dat)
  snps<-snps[!is.na(snps$beta.exposure),]
  snps<-snps[!is.na(snps$beta.outcome),]
  snps<-snps[!is.na(snps$se.outcome),]
  snps<-snps[!snps$se.outcome==0,]
  snps<-snps[!is.na(snps$eaf.outcome),]
  snps<-snps[!is.na(snps$eaf.exposure),]
  
  df1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = snps$samplesize.exposure)
  df2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = snps$samplesize.outcome)
  
  res=coloc.abf(df1,df2)
  
  df3 <- data.frame(t(res$summary))
  df3$id<-a
  df3$risk_factor<-"WHR"
  results<-rbind(results,df3)
  
  #Alcohol consumption
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = "data/CRC_risk_factors/DRINKS_PER_WEEK_GWAS.txt",
    sep = "\t",
    snp_col = "MarkerName",
    beta_col = "Beta",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "EAF_A1",
    se_col = "SE"
  )
  outcome_dat$samplesize.outcome<-414343
  outcome_dat$p_one_sided<-outcome_dat$pval.outcome/2
  outcome_dat$z<-qnorm(outcome_dat$p_one_sided)
  outcome_dat$se.outcome<-as.numeric(abs(outcome_dat$beta.outcome/outcome_dat$z))
  
  snps<-harmonise_data(exposure_dat,outcome_dat)
  snps<-snps[!is.na(snps$beta.exposure),]
  snps<-snps[!is.na(snps$beta.outcome),]
  snps<-snps[!is.na(snps$se.outcome),]
  snps<-snps[!snps$se.outcome==0,]
  snps<-snps[!is.na(snps$eaf.outcome),]
  snps<-snps[!is.na(snps$eaf.exposure),]
  
  df1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = snps$samplesize.exposure)
  df2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = snps$samplesize.outcome)
  
  res=coloc.abf(df1,df2)
  
  df3 <- data.frame(t(res$summary))
  df3$id<-a
  df3$risk_factor<-"Alcohol"
  results<-rbind(results,df3)
  
  #Tobacco use
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = "data/CRC_risk_factors/ieu-b-4877.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    samplesize_col="samplesize.exposure"
  )
  
  
  snps<-harmonise_data(exposure_dat,outcome_dat)
  snps<-snps[!is.na(snps$beta.exposure),]
  snps<-snps[!is.na(snps$beta.outcome),]
  snps<-snps[!is.na(snps$se.outcome),]
  snps<-snps[!is.na(snps$eaf.outcome),]
  snps<-snps[!is.na(snps$eaf.exposure),]
  
  df1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = snps$samplesize.exposure)
  df2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = snps$samplesize.outcome)
  
  res=coloc.abf(df1,df2)
  
  df3 <- data.frame(t(res$summary))
  df3$id<-a
  df3$risk_factor<-"Tobacco"
  results<-rbind(results,df3)

  
}


#Sex-specific not needed as none in druggable

# Save file ---------------------------------------------------------------

fwrite(results,"CRC_TWAS/Risk_factor_druggable_colocalization_results.csv")

print("fin")

sink()
