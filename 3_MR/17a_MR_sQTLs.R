sink("MR.txt")
#module load R/4.2.0

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(stringr)
library(ieugwasr)

# Get gene list -----------------------------------------------------------
paste("Get gene list")

splice<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")

splice$splice_id<-paste(splice$chr,splice$start,splice$end,sep=":")
splice$id<-paste(splice$splice_id,splice$t_i_best,sep=";")

list<-dplyr::select(splice,splice_id,id,t_i_best,subtype,gene)
list$subtype<-tolower(list$subtype)

list<-list[!list$splice_id=="",]
list<-list[!list$t_i_best=="",]
list<-list[!list$subtype=="",]
list<-list[!duplicated(list),]


# Get sQTL data -----------------------------------------------------------
paste("Get sQTL data")

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

snps<-fread("data/GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)

eqtl$list<-paste(eqtl$splice_id,eqtl$tissue,sep=";")


eqtl$gene<-substr(eqtl$gene,1,15)
eqtl$id<-paste(eqtl$gene,eqtl$tissue,sep=";")
eqtl <- as.data.frame(eqtl)

eqtl<-eqtl[eqtl$list %in% list$id,]

#Find SNP rsIDs
snps<-fread("data/GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
sqtl<-merge(eqtl,snps,by=c("variant_id","chr","variant_pos","ref","alt","num_alt_per_site","rs_id_dbSNP151_GRCh38p7","variant_id_b37"),all.x=T)

sqtl<-sqtl[sqtl$splice_id %in% list$splice_id,]

sqtl<-sqtl %>%
  mutate(chromosome = str_remove_all(chr, "chr"))

# Limit to cis SNPs -------------------------------------------------------
paste("Limit to cis SNPs")

#Find gene coding region for each gene (note: do this on hpcapp01 as need internet access)
library("biomaRt")

# Set up the BioMart connection
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define the genes to query
gene_list <- unique(list$gene)

# Retrieve data from BioMart
genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  values = gene_list,
  mart = ensembl,
  useCache=FALSE
)

genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000
fwrite(genes,"data/CRC_TWAS/genes_with_cis_pos_splice.csv")

#Read file back in on interactive job
genes<-fread("data/CRC_TWAS/genes_with_cis_pos_splice.csv")
exp1<-merge(sqtl,genes,by.x="gene",by.y="ensembl_gene_id")
exp1<-exp1[exp1$chromosome==exp1$chromosome_name,]
exp1<-exp1[exp1$variant_pos>=exp1$cisstart,]
exp1<-exp1[exp1$variant_pos<=exp1$cisend,]
exp1<-data.frame(exp1)
exp1$id<-paste(exp1$splice_id,exp1$tissue,sep=";")


# MR ----------------------------------------------------------------------
paste("MR")

exposure_dat<-format_data(
  exp1,
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
fwrite(exposure_dat,"data/CRC_TWAS/exposure_dat_temp.csv")

exposure_dat<-fread("data/CRC_TWAS/exposure_dat_temp.csv")
exposure_dat<- ld_clump(
  exposure_dat,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "data/1000GenomesReferenceFiles/EUR",
  pop="EUR")
fwrite(exposure_dat,"data/CRC_TWAS/clumped_exposure_dat_temp_sqtls.csv")

# r2 and F stat -----------------------------------------------------------
paste("r2 and F stat")

exposure_dat<-fread("data/CRC_TWAS/clumped_exposure_dat_temp_sqtls.csv")

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
f_r2<-f_r2 %>% tidyr::separate(exposure, c("Gene", "Tissue","Splicing event"),sep=";")
fwrite(f_r2,"data/CRC_TWAS/f_stat_and_r2_sqtls.csv")


# Running MR --------------------------------------------------------------
paste("Running MR")
no_outcome<-list()
no_harmonise<-list()
steiger<-list()

exp <- read_exposure_data(
  filename = "data/CRC_TWAS/clumped_exposure_dat_temp_sqtls.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval",
  id_col = "id.exposure",
  phenotype_col="exposure"
  )

list<-list[!duplicated(list$id),]
results<-data.frame()

#Filter for SNPs present in GECCO as breaks otherwise
ov<-fread("data/gecco/annotated/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
co<-fread("data/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
di<-fread("data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
fe<-fread("data/gecco/annotated/female_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
ma<-fread("data/gecco/annotated/male_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
pr<-fread("data/gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
re<-fread("data/gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")

snps<-unique(c(ov$SNP,co$SNP,di$SNP,fe$SNP,ma$SNP,pr$SNP,re$SNP))

exp<-exp[exp$SNP %in% snps,]

#Cases and controls for Steiger filtering
Ns<-data.frame(matrix(ncol=3,nrow=7))
colnames(Ns)<-c("CRC","ncase","ncontrol")
Ns$CRC<-c("overall","colon","rectal","proximal","distal","male","female")
Ns$ncase<-c(98715,28736,14150,14416,12879,28271,24594)
Ns$ncontrol<-c(52775,43099,43099,43099,43099,22351,23936)

#Sample sizes
exp$samplesize.exposure<-NA
exp$samplesize.exposure[grepl("Adipose_Subcutaneous",exp$exposure)]<-581
exp$samplesize.exposure[grepl("Adipose_Visceral_Omentum",exp$exposure)]<-469
exp$samplesize.exposure[grepl("Cells_EBV-transformed_lymphocytes",exp$exposure)]<-187
exp$samplesize.exposure[grepl("Colon_Sigmoid",exp$exposure)]<-318
exp$samplesize.exposure[grepl("Colon_Transverse",exp$exposure)]<-368
exp$samplesize.exposure[grepl("Whole_Blood",exp$exposure)]<-670


results<-data.frame()

for (a in list$id){
  print(a)
  mr<-list[list$id==a,]
  dat1<-exp[exp$exposure %in% mr$id,]
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
      samplesize_col = "samplesize.exposure",
      log_pval = FALSE
    )  
    exposure_dat$exposure<-a
    b<-mr$subtype[1]
    print(b)
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = paste("data/gecco/annotated/",b,"_CRC_GWAS_noUKBio_summary_stats_annotated.txt",sep=""),
      sep = " ",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      eaf_col = "Freq1"
    )
    
    no_outcome<-c(no_outcome,exposure_dat$exposure[!exposure_dat$SNP %in% outcome_dat$SNP])
    
    dat<-harmonise_data(exposure_dat,outcome_dat)
    
    no_harmonise<-c(no_harmonise,exposure_dat$exposure[!exposure_dat$exposure %in% dat$exposure])
    
    if(length(which(dat$mr_keep=="TRUE"))>0){
      
      #Steiger filtering  
      dat$ncase.outcome<-Ns$ncase[which(Ns$CRC==b)]
      dat$ncontrol.outcome<-Ns$ncontrol[which(Ns$CRC==b)]
      dat$samplesize.outcome<-dat$ncase.outcome+dat$ncontrol.outcome
      dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
      dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
      dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                      af=dat$eaf.outcome,
                                      ncase=dat$ncase.outcome,
                                      ncontrol=dat$ncontrol.outcome,
                                      prevalence=dat$prevalence.outcome)
      dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
      dat$eaf.exposure[is.na(dat$eaf.exposure)]<-dat$eaf.outcome[is.na(dat$eaf.exposure)]
      dat <- steiger_filtering(dat)
      steiger<-c(steiger,dat$exposure[dat$steiger_dir==FALSE])
      
      dat<-dat[dat$steiger_dir==TRUE,]
      
      if(length(which(dat$mr_keep=="TRUE"))>0){
        
      res <- mr(dat,method_list = c("mr_wald_ratio", "mr_ivw"))
      res$subtype<-b
      res$Steiger<-dat$direction
      results<-rbind(results,res)
    }}
  }}
print("fin")
fwrite(results,"data/CRC_TWAS/MR_results_sqtls.csv")

fwrite(list(no_outcome),"data/CRC_TWAS/no_outcome_spl.csv")
fwrite(list(no_harmonise),"data/CRC_TWAS/no_harmonise_spl.csv")
fwrite(list(steiger),"data/CRC_TWAS/steiger_spl.csv")

sink()
