sink("finemapping.txt")
library(susieR)
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(coloc)
library(gwasglue)
library(stringr)
library(plinkbinr)

events<-fread("CRC_TWAS/Colocalization_splicing_results_strong.csv")

events<-events%>%tidyr::separate(id,c("chr","start","end","tissue"),sep=":")
events$id<-paste(events$chr,events$start,events$end,sep=":")

events$start<-as.numeric(events$start)
events$end<-as.numeric(events$end)

events$window_start<-events$start-100000
events$window_end<-events$end+100000

events$subtype<-tolower(events$subtype)


events<- events %>%
  mutate(chr = str_remove_all(chr, "chr"))

as<-fread("/home/sw20203/Splicing_GTEx/Adipose_Subcutaneous.v8.sqtl_allpairs_filtered.txt")
av<-fread("/home/sw20203/Splicing_GTEx/Adipose_Visceral_Omentum.v8.sqtl_allpairs_filtered.txt")
ly<-fread("/home/sw20203/Splicing_GTEx/Cells_EBV-transformed_lymphocytes.v8.sqtl_allpairs_filtered.txt")
cs<-fread("/home/sw20203/Splicing_GTEx/Colon_Sigmoid.v8.sqtl_allpairs_filtered.txt")
ct<-fread("/home/sw20203/Splicing_GTEx/Colon_Transverse.v8.sqtl_allpairs_filtered.txt")


as$tissue<-"Adipose_Subcutaneous"
av$tissue<-"Adipose_Visceral_Omentum"
ly$tissue<-"Cells_EBV-transformed_lymphocytes"
cs$tissue<-"Colon_Sigmoid"
ct$tissue<-"Colon_Transverse"

eqtl<-rbind(as,av,ly,cs,ct)
colnames(eqtl)<-c("splice_id","id","gene","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nomminal","slope","slope_se","tissue")

eqtl$gene<-substr(eqtl$gene,1,15)

snps<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)

eqtl$samplesize[eqtl$tissue=="Adipose_Subcutaneous"]<-581
eqtl$samplesize[eqtl$tissue=="Adipose_Visceral_Omentum"]<-469
eqtl$samplesize[eqtl$tissue=="Cells_EBV-transformed_lymphocytes"]<-187
eqtl$samplesize[eqtl$tissue=="Colon_Sigmoid"]<-318
eqtl$samplesize[eqtl$tissue=="Colon_Transverse"]<-368

results<-data.frame()

for (a in unique(events$id)){
  print(a)
  mr<-events[events$id==a,]
  
  for (c in unique(mr$tissue)){
    print(c)
    sqtl<-eqtl[eqtl$tissue==c,]
    sqtl<-sqtl[sqtl$splice_id==a,]
    
    sqtl<-sqtl[!is.na(sqtl$slope),]
    sqtl<-sqtl[!is.na(sqtl$slope_se),]
    
    sqtl$SNP_alleles_1<-paste(sqtl$rs_id_dbSNP151_GRCh38p7,sqtl$ref,sqtl$alt,sep="_")
    sqtl$SNP_alleles_2<-paste(sqtl$rs_id_dbSNP151_GRCh38p7,sqtl$alt,sqtl$ref,sep="_")
    
    ld<-ld_matrix_local(sqtl$rs_id_dbSNP151_GRCh38p7, with_alleles = TRUE, bfile="UK_Biobank_Reference/mergey",plink_bin=get_plink_exe())
    
    ld<-ld[which(rownames(ld) %in% c(sqtl$SNP_alleles_1,sqtl$SNP_alleles_1)),which(rownames(ld) %in% c(sqtl$SNP_alleles_1,sqtl$SNP_alleles_1))]
    colnames(ld)<-substr(colnames(ld),1,nchar(colnames(ld))-4)
    rownames(ld)<-substr(rownames(ld),1,nchar(rownames(ld))-4)
    
    sqtl<-sqtl[which(sqtl$rs_id_dbSNP151_GRCh38p7 %in% rownames(ld)),]
    ld<-ld[match(sqtl$rs_id_dbSNP151_GRCh38p7,rownames(ld)),]
    ld<-ld[,match(sqtl$rs_id_dbSNP151_GRCh38p7,colnames(ld))]
    sqtl<-sqtl[match(rownames(ld),sqtl$rs_id_dbSNP151_GRCh38p7),]
    
    zscore=sqtl$slope/sqtl$slope_se
    estimate_s_rss(z=zscore, R=ld, n=sqtl$samplesize[1])
    
    res<-susie_rss(R=ld,n=sqtl$samplesize[1],bhat=sqtl$slope,shat=sqtl$slope_se)
    
    dataset <- list(beta=sqtl$slope, varbeta = sqtl$slope_se^2, MAF=sqtl$maf, type = "quant", snps=rownames(ld),LD=ld,N=sqtl$samplesize[1])
    
    df<-runsusie(dataset,N=sqtl$samplesize[1],estimate_prior_variance=FALSE)
    
    results2<-data.frame(matrix(ncol=3,nrow=length(res$sets$cs$L1)))
    
    if (nrow(results2)>0){
    results2[1]<-a
    results2[2]<-c
    results2[3]<-rownames(ld)[res$sets$cs$L1]
    
    results<-rbind(results,results2)
    }
    
  }
  
  
}


fwrite(results,"CRC_TWAS/Finemapping_results.csv")


# For missing events ------------------------------------------------------

results<-fread("CRC_TWAS/Finemapping_results.csv")

colnames(results)<-c("id","Tissue","SNP")

results$id_tissue<-paste(results$id,results$Tissue,sep="_")
events$id_tissue<-paste(events$id,events$tissue,sep="_")

missing<-events[!events$id_tissue %in% results$id_tissue,]

results<-data.frame()

for (a in unique(missing$id)){
  print(a)
  mr<-missing[missing$id==a,]
  
  for (c in unique(mr$tissue)){
    print(c)
    sqtl<-eqtl[eqtl$tissue==c,]
    sqtl<-sqtl[sqtl$splice_id==a,]
    
    sqtl<-sqtl[!is.na(sqtl$slope),]
    sqtl<-sqtl[!is.na(sqtl$slope_se),]
    
    sqtl$SNP_alleles_1<-paste(sqtl$rs_id_dbSNP151_GRCh38p7,sqtl$ref,sqtl$alt,sep="_")
    sqtl$SNP_alleles_2<-paste(sqtl$rs_id_dbSNP151_GRCh38p7,sqtl$alt,sqtl$ref,sep="_")
    
    ld<-ld_matrix_local(sqtl$rs_id_dbSNP151_GRCh38p7, with_alleles = TRUE, bfile="UK_Biobank_Reference/mergey",plink_bin=get_plink_exe())
    
    ld<-ld[which(rownames(ld) %in% c(sqtl$SNP_alleles_1,sqtl$SNP_alleles_1)),which(rownames(ld) %in% c(sqtl$SNP_alleles_1,sqtl$SNP_alleles_1))]
    colnames(ld)<-substr(colnames(ld),1,nchar(colnames(ld))-4)
    rownames(ld)<-substr(rownames(ld),1,nchar(rownames(ld))-4)
    
    sqtl<-sqtl[which(sqtl$rs_id_dbSNP151_GRCh38p7 %in% rownames(ld)),]
    ld<-ld[match(sqtl$rs_id_dbSNP151_GRCh38p7,rownames(ld)),]
    ld<-ld[,match(sqtl$rs_id_dbSNP151_GRCh38p7,colnames(ld))]
    sqtl<-sqtl[match(rownames(ld),sqtl$rs_id_dbSNP151_GRCh38p7),]
    
    all<-eqtl[eqtl$tissue==c,]
    all<-all[all$splice_id==a,]
    
    min_p<-min(all$pval_nomminal)
    log10<-round(log10(min_p))
    new_log10<-log10+2
    new_min_p<-min_p*10^abs(new_log10)
    
    sqtl<-sqtl[sqtl$pval_nomminal<new_min_p,]
    
    if (nrow(sqtl)>0){
      results<-rbind(results,sqtl)
    }
    
  }
  
  
}


results<-dplyr::select(results,splice_id,tissue,rs_id_dbSNP151_GRCh38p7)

results2<-fread("CRC_TWAS/Finemapping_results.csv")
colnames(results2)<-c("id","Tissue","SNP")

results<-rbind(results2,results,use.names=FALSE)

fwrite(results,"CRC_TWAS/Finemapping_results.csv")

print("fin")




sink()
