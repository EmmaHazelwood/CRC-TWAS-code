sink("coloc_abf_druggable.txt")
library(coloc)
library(dplyr)
library(gwasglue)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(plinkbinr)
library(readxl)
library(biomaRt)


# Get ones needed ---------------------------------------------------------

dat<-read_excel("CRC_TWAS/41591_2021_1310_MOESM3_ESM.xlsx",sheet=3,skip=11)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"), values=dat$ENSEMBL_GENE_ID, mart=ensembl)
genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000

# Get SNPs ----------------------------------------------------------------

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

as$id<-paste(as$gene_id,"Adipose_Subcutaneous",sep=";")
av$id<-paste(av$gene_id,"Adipose_Visceral_Omentum",sep=";")
ly$id<-paste(ly$gene_id,"Cells_EBV-transformed_lymphocytes",sep=";")
cs$id<-paste(cs$gene_id,"Colon_Sigmoid",sep=";")
ct$id<-paste(ct$gene_id,"Colon_Transverse",sep=";")
wb$id<-paste(wb$gene_id,"Whole_Blood",sep=";")

as$samplesize<-581
av$samplesize<-469
ly$samplesize<-187
cs$samplesize<-318
ct$samplesize<-368
wb$samplesize<-670

eqtl<-rbind(as,av,ly,cs,ct,wb)
eqtl$gene<-substr(eqtl$gene_id,1,15)

eqtl<-eqtl[eqtl$gene %in% genes$ensembl_gene_id,]

snps<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)
eqtl<-eqtl %>%
  mutate(chr = str_remove_all(chr, "chr"))

exp<-merge(eqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp<-exp[exp$chr==exp$chromosome_name,]
exp<-exp[exp$variant_pos>=exp$cisstart,]
exp<-exp[exp$variant_pos<=exp$cisend,]

subtype<-c("proximal","distal","colon","rectal","overall","male","female")
samplesize<-c(57515,55978,57249,71835,98715,50622,48530)
samplesize_df<-data.frame(subtype,samplesize)

exp<-exp[!exp$chromosome_name=="X",]

# Run coloc ---------------------------------------------------------------
results<-data.frame()
for (a in unique(exp$id)){
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
  exposure_dat<-exposure_dat[!exposure_dat$beta.exposure=="0",]
  
  for (b in c("overall","colon","distal","female","male","proximal","rectal")){
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
    outcome_dat$samplesize.outcome<-samplesize_df$samplesize[samplesize_df$subtype==b]
    snps<-harmonise_data(exposure_dat,outcome_dat)
    snps<-snps[!is.na(snps$beta.exposure),]
    snps<-snps[!is.na(snps$beta.outcome),]
    
    df1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = snps$samplesize.exposure)
    df2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "cc", N = snps$samplesize.outcome)
    
    res=coloc.abf(df1,df2)
    
    df3 <- data.frame(t(res$summary))
    df3$id<-a
    df3$subtype<-b
    results<-rbind(results,df3)
    
  }
}


fwrite(results,"CRC_TWAS/Coloc_results_abf_druggable.csv")

print("fin")

sink()