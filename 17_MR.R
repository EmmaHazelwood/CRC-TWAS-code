sink("MR.txt")
#module load R/4.2.0

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(stringr)
library(ieugwasr)

# Get gene list -----------------------------------------------------------

splice<-fread("CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
multi<-fread("CRC_TWAS/SMultiXcan_Expression_All_Strong.csv")
jti<-fread("CRC_TWAS/JTI_All_Strong.csv")
jti<-jti[!jti$subtype=="",]

splice<-dplyr::select(splice,gene,t_i_best,subtype)
multi<-dplyr::select(multi,gene,t_i_best,subtype)
jti<-dplyr::select(jti,gene,tissue,subtype)
list<-rbind(splice,multi,jti,use.names=F)
list<-list[!duplicated(list),]

list$gene<-substr(list$gene,1,15)

list$subtype<-tolower(list$subtype)

list$gene_id<-paste(list$gene,list$t_i_best,sep=";")
list$id<-paste(list$gene_id,list$subtype,sep=";")

list<-list[!list$gene=="",]
list<-list[!list$t_i_best=="",]
list<-list[!list$subtype=="",]
list<-list[!duplicated(list),]


# Get eQTL data -----------------------------------------------------------


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

as$tissue<-"Adipose_Subcutaneous"
av$tissue<-"Adipose_Visceral_Omentum"
ly$tissue<-"Cells_EBV-transformed_lymphocytes"
cs$tissue<-"Colon_Sigmoid"
ct$tissue<-"Colon_Transverse"
wb$tissue<-"Whole_Blood"

eqtl<-rbind(as,av,ly,cs,ct,wb)
eqtl$gene<-substr(eqtl$gene_id,1,15)
eqtl$id<-paste(eqtl$gene,eqtl$tissue,sep=";")

eqtl<-eqtl[eqtl$id %in% list$gene_id,]

#Find SNP rsIDs
snps<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)
eqtl<-eqtl %>%
  mutate(chr = str_remove_all(chr, "chr"))


# Limit to cis SNPs -------------------------------------------------------

#Find gene coding region for each gene (note: do this on hpcapp01 as need internet access)
library("biomaRt")                                                                                                                   
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"), values=unique(list$gene), mart=ensembl)
genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000
fwrite(genes,"CRC_TWAS/genes_with_cis_pos.csv")

#Read file back in on interactive job
genes<-fread("CRC_TWAS/genes_with_cis_pos.csv")
exp<-merge(eqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp<-exp[exp$chr==exp$chromosome_name,]
exp<-exp[exp$variant_pos>=exp$cisstart,]
exp<-exp[exp$variant_pos<=exp$cisend,]


# MR ----------------------------------------------------------------------

exposure_dat<-format_data(
  exp,
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
  log_pval = FALSE
)  
exposure_dat$pval.exposure<-as.numeric(exposure_dat$pval.exposure)
exposure_dat<-exposure_dat[exposure_dat$pval.exposure<5*10^-8,]
exposure_dat$rsid<-exposure_dat$SNP
exposure_dat$id<-exposure_dat$id.exposure
exposure_dat$pval<-exposure_dat$pval.exposure

#clump
#Save file to clump server that has internet access
remotes::install_github("MRCIEU/genetics.binaRies")
exposure_dat<-exposure_dat[!is.na(exposure_dat$gene.exposure),]
fwrite(exposure_dat,"CRC_TWAS/exposure_dat_temp.csv")

exposure_dat<-fread("CRC_TWAS/exposure_dat_temp.csv")
exposure_dat<- ld_clump(
  exposure_dat,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "1000GenomesReferenceFiles/EUR",
  pop="EUR")
fwrite(exposure_dat,"CRC_TWAS/clumped_exposure_dat_temp.csv")

# r2 and F stat -----------------------------------------------------------
exposure_dat<-fread("CRC_TWAS/clumped_exposure_dat_temp.csv")
                    
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
fwrite(f_r2,"CRC_TWAS/f_stat_and_r2.csv")


# Running MR --------------------------------------------------------------
exp <- read_exposure_data(
  filename = "CRC_TWAS/clumped_exposure_dat_temp.csv",
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

list<-list[!duplicated(list$id),]
results<-data.frame()

#Filter for SNPs present in GECCO as breaks otherwise
ov<-fread("gecco/annotated/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
co<-fread("gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
di<-fread("gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
fe<-fread("gecco/annotated/female_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
ma<-fread("gecco/annotated/male_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
pr<-fread("gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
re<-fread("gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")

snps<-unique(c(ov$SNP,co$SNP,di$SNP,fe$SNP,ma$SNP,pr$SNP,re$SNP))

exp<-exp[exp$SNP %in% snps,]

results<-data.frame()

for (a in list$id){
  print(a)
  mr<-list[list$id==a,]
  dat1<-exp[exp$exposure %in% mr$gene_id,]
  if(nrow(dat1)>0){
  exposure_dat<-format_data(
    dat1,
    type = "exposure",
    header = TRUE,
    phenotype_col="exposure",
    id_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    log_pval = FALSE
  )  
  exposure_dat$exposure<-a
  b<-mr$subtype[1]
  print(b)
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = paste("gecco/annotated/",b,"_CRC_GWAS_noUKBio_summary_stats_annotated.txt",sep=""),
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  dat<-harmonise_data(exposure_dat,outcome_dat)
  
  if(length(which(dat$mr_keep=="TRUE"))>0){
  res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
  res$subtype<-b
  results<-rbind(results,res)
}
}}
print("fin")
fwrite(results,"CRC_TWAS/MR_results.csv")
sink()
