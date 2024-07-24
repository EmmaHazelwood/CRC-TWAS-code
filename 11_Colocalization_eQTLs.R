sink("coloc_abf.txt")
library(coloc)
library(dplyr)
library(gwasglue)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(plinkbinr)


# Get ones needed ---------------------------------------------------------

splice<-fread("SMultiXcan_Splicing_Expression_All_Strong.csv")
multi<-fread("SMultiXcan_Expression_All_Strong.csv")
jti<-fread("JTI_All_Strong.csv")
genes<-c(splice$gene,substr(multi$gene,1,15),jti$gene)
genes<-unique(genes)

splice<-select(splice,gene,t_i_best,subtype)
multi<-select(multi,gene,t_i_best,subtype)
jti<-select(jti,gene,tissue,subtype)
list<-rbind(splice,multi,jti,use.names=F)
list<-list[!duplicated(list),]

list$gene<-substr(list$gene,1,15)

list$subtype<-tolower(list$subtype)


# Get SNPs ----------------------------------------------------------------

as<-fread("v8_eQTL_european_associations/Adipose_Subcutaneous.allpairs.txt.gz")
av<-fread("v8_eQTL_european_associations/Adipose_Visceral_Omentum.allpairs.txt.gz")
ly<-fread("bris_GTEx/Cells_EBV-transformed_lymphocytes.allpairs.csv")
cs<-fread("bris_GTEx/Colon_Sigmoid.allpairs.csv")
ct<-fread("bris_GTEx/Colon_Transverse.allpairs.csv")
wb<-fread("bris_GTEx/Whole_Blood.allpairs.csv")

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

eqtl<-eqtl[eqtl$gene %in% list$gene,]

genes<-fread("genes_with_cis_pos.csv")
exp<-merge(eqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp<-exp[exp$chr==exp$chromosome_name,]
exp<-exp[exp$variant_pos>=exp$cisstart,]
exp<-exp[exp$variant_pos<=exp$cisend,]

subtype<-c("proximal","distal","colon","rectal","overall","male","female")
samplesize<-c(57515,55978,57249,71835,98715,50622,48530)
samplesize_df<-data.frame(subtype,samplesize)

# Run coloc ---------------------------------------------------------------
list$subtype<-tolower(list$subtype)
list$gene_id<-paste(list$gene,list$t_i_best,sep=";")
list$id<-paste(list$gene_id,list$subtype,sep=";")

list<-list[!list$gene=="",]
list<-list[!list$t_i_best=="",]
list<-list[!list$subtype=="",]
list<-list[!duplicated(list),]

results<-data.frame()
for (a in unique(list$id)){
  print(a)
  mr<-list[list$id==a,]
  dat1<-exp[exp$id %in% mr$gene_id,]
  
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

print("fin")

fwrite(results,"Coloc_results_abf.csv")


sink()