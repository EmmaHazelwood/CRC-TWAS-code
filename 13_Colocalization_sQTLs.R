sink("coloc_abf_splice.txt")
library(coloc)
library(dplyr)
library(gwasglue)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(plinkbinr)


# Get ones needed ---------------------------------------------------------

splice<-fread("CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
splice$splice_id<-paste(splice$chr,splice$start,splice$end,sep=":")
splice<-select(splice,splice_id,subtype,gene,t_i_best)
splice2<-select(splice,splice_id)

as<-splice2[splice$t_i_best=="Adipose_Subcutaneous",]
av<-splice2[splice$t_i_best=="Adipose_Visceral_Omentum",]
cs<-splice2[splice$t_i_best=="Colon_Sigmoid",]
ct<-splice2[splice$t_i_best=="Colon_Transverse",]
ly<-splice2[splice$t_i_best=="Cells_EBV-transformed_lymphocytes",]
wb<-splice2[splice$t_i_best=="Whole_Blood",] #Note: none of the splice events pass Bonferroni in blood


fwrite(as,"splice_as",sep="\t",quote = F)
fwrite(av,"splice_av",sep="\t",quote = F)
fwrite(cs,"splice_cs",sep="\t",quote = F)
fwrite(ct,"splice_ct",sep="\t",quote = F)
fwrite(ly,"splice_ly",sep="\t",quote = F)
fwrite(wb,"splice_wb",sep="\t",quote = F)

# Get SNPs ----------------------------------------------------------------

as<-fread("/home/sw20203/Splicing_GTEx/Adipose_Subcutaneous.v8.sqtl_allpairs_filtered.txt")
av<-fread("/home/sw20203/Splicing_GTEx/Adipose_Visceral_Omentum.v8.sqtl_allpairs_filtered.txt")
ly<-fread("/home/sw20203/Splicing_GTEx/Cells_EBV-transformed_lymphocytes.v8.sqtl_allpairs_filtered.txt")
cs<-fread("/home/sw20203/Splicing_GTEx/Colon_Sigmoid.v8.sqtl_allpairs_filtered.txt")
ct<-fread("/home/sw20203/Splicing_GTEx/Colon_Transverse.v8.sqtl_allpairs_filtered.txt")
wb<-fread("/home/sw20203/Splicing_GTEx/Whole_Blood.v8.sqtl_allpairs_filtered.txt")


as$tissue<-"Adipose_Subcutaneous"
av$tissue<-"Adipose_Visceral_Omentum"
ly$tissue<-"Cells_EBV-transformed_lymphocytes"
cs$tissue<-"Colon_Sigmoid"
ct$tissue<-"Colon_Transverse"
wb$tissue<-"Whole_Blood"

eqtl<-rbind(as,av,ly,cs,ct,wb)
colnames(eqtl)<-c("splice_id","id","gene","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nomminal","slope","slope_se","tissue")

eqtl$gene<-substr(eqtl$gene,1,15)

snps<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)

eqtl$list<-paste(eqtl$splice_id,eqtl$tissue,sep=":")

eqtl<-eqtl %>%
  mutate(chr = str_remove_all(chr, "chr"))

genes<-fread("CRC_TWAS/genes_with_cis_pos.csv")
exp<-merge(eqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp<-exp[exp$chr==exp$chromosome_name,]
exp<-exp[exp$variant_pos>=exp$cisstart,]
exp<-exp[exp$variant_pos<=exp$cisend,]

exp$samplesize<-NA
exp$samplesize[exp$tissue=="Adipose_Subcutaneous"]<-581
exp$samplesize[exp$tissue=="Adipose_Visceral_Omentum"]<-469
exp$samplesize[exp$tissue=="Cells_EBV-transformed_lymphocytes"]<-187
exp$samplesize[exp$tissue=="Colon_Sigmoid"]<-318
exp$samplesize[exp$tissue=="Colon_Transverse"]<-368
exp$samplesize[exp$tissue=="Whole_Blood"]<-670

exp<-exp[!exp$chr=="X",]

subtype<-c("proximal","distal","colon","rectal","overall","male","female")
samplesize<-c(57515,55978,57249,71835,98715,50622,48530)
samplesize_df<-data.frame(subtype,samplesize)


# Run coloc ---------------------------------------------------------------
splice$id<-paste(splice$splice_id,splice$t_i_best,sep=":")
splice$list<-paste(splice$id,splice$subtype,sep=";")


fwrite(exp,"coloc_data_splice.csv")
exp<-fread("coloc_data_splice.csv")

results<-data.frame()
for (a in unique(splice$list)){
  print(a)
  mr<-splice[splice$list==a,]
  dat1<-exp[exp$list %in% mr$id,]
  if(nrow(dat1)>0){
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
    log_pval = FALSE
  )  
  
  b<-tolower(mr$subtype[1])
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
    snps<-snps[!is.na(snps$pval.exposure),]
    
    df1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = snps$samplesize.exposure)
    df2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "cc", N = snps$samplesize.outcome)
    
    res=coloc.abf(df1,df2)
    
    df3 <- data.frame(t(res$summary))
    df3$id<-a
    df3<-df3%>%tidyr::separate(id,c("id","subtype"),sep=";")
    results<-rbind(results,df3)
    
  }
  
}

fwrite(results,"CRC_TWAS/Coloc_results_abf_splice.csv")

print("fin")

sink()
