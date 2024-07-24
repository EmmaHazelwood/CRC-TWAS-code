sink("finemapping_CRC.txt")

library(data.table)
library(dplyr)
library(stringr)

res<-fread("/CRC_TWAS/Finemapping_results.csv")
colnames(res)<-c("splice_id","tissue","SNP")

events<-fread("/CRC_TWAS/Colocalization_splicing_results_strong.csv")

events<-events%>%tidyr::separate(id,c("chr","start","end","tissue"),sep=":")
events$id<-paste(events$chr,events$start,events$end,sep=":")
events$subtype<-tolower(events$subtype)

events2<-select(events,id,tissue,subtype)

res<-merge(res,events2,by.x=c("splice_id","tissue"),by.y=c("id","tissue"),allow.cartesian=TRUE)

overall<-fread("/gecco/annotated/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
colon<-fread("/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
distal<-fread("/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
female<-fread("/gecco/annotated/female_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
male<-fread("/gecco/annotated/male_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
proximal<-fread("/gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
rectal<-fread("/gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")

results<-data.frame()
for (a in unique(res$subtype)){
  print(a)
  crc<-get(a)
  crc<-select(crc,SNP,P.value)
  crc$subtype<-a
  results<-rbind(results,crc)
}

res<-merge(res,results,by=c("SNP","subtype"),all.x=T)

events$start<-as.numeric(events$start)
events$end<-as.numeric(events$end)

events$window_start<-events$start-100000
events$window_end<-events$end+100000

events<- events %>%
  mutate(chr = str_remove_all(chr, "chr"))

events<-select(events,id,chr,window_start,window_end)

list<-select(res,splice_id,subtype)
list<-distinct(list)

list<-merge(list,events,by.x="splice_id",by.y="id",allow.cartesian=TRUE)
list<-distinct(list)

results<-data.frame()
for (b in 1:nrow(list)){
  print(b)
  crc<-get(list$subtype[b])
  crc<-crc[crc$Chromosome==list$chr[b],]
  crc<-crc[crc$Position>=list$window_start[b],]
  crc<-crc[crc$Position<=list$window_end[b],]
  
  crc_p<-data.frame(matrix(ncol=0,nrow=1))
  
  crc_p$subtype<-list$subtype[b]
  crc_p$min_p<-min(crc$P.value)
  crc_p$splice_id<-list$splice_id[b]
  
  results<-rbind(results,crc_p)
   
}

res<-merge(res,results,by=c("subtype","splice_id"))

res$log10<-round(log10(res$min_p))
res$new_log10<-res$log10+2
res$new_min_p<-res$min_p*10^abs(res$new_log10)

res$crc<-"no"
res$crc[res$P.value<=res$new_min_p]<-"yes"


fwrite(res,"/CRC_TWAS/Finemapping_results_CRC_GWAS_P.csv")
res<-res[res$crc=="yes"]

res<-select(res,splice_id,SNP)
res<-distinct(res)

fwrite(res,"/CRC_TWAS/SNP_list_for_SpliceAI.csv")

res<-fread("/CRC_TWAS/SNP_list_for_SpliceAI.csv")

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

snps<-fread("/GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
eqtl<-merge(eqtl,snps,by="variant_id",all.x=T)

eqtl<-select(eqtl,rs_id_dbSNP151_GRCh38p7,splice_id,gene,variant_id,ref,alt,slope)

res<-merge(res,eqtl,by.x=c("splice_id","SNP"),by.y=c("splice_id","rs_id_dbSNP151_GRCh38p7"),all.x=T)

res<-distinct(res)

fwrite(res,"/CRC_TWAS/SNPs_for_SpliceAI.csv")


print("fin")




sink()

