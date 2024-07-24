library(data.table)
library(dplyr)

mr_res<-fread("CRC_TWAS/MR_results.csv")

mr_res2<-mr_res %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
mr_res3<-mr_res2[mr_res2$pval<(0.05/(length(unique(mr_res$exposure))*7)),]
length(unique(mr_res3$Gene))
length(unique(mr_res2$Gene))


#Filter for ones with effects in TWAS
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

list$subtype<-tolower(list$subtype)
list$gene_id<-paste(list$gene,list$t_i_best,sep=";")
list$id<-paste(list$gene_id,list$subtype,sep=";")

list<-list[!list$gene=="",]
list<-list[!list$t_i_best=="",]
list<-list[!list$subtype=="",]
list<-list[!duplicated(list),]

mr_res2$Gene<-substr(mr_res2$Gene,1,15)
mr_res2$subtype<-tolower(mr_res2$subtype)
mr_res2$id<-paste(mr_res2$Gene,mr_res2$Tissue,mr_res2$subtype,sep=";")

mr_res2<-mr_res2[mr_res2$id %in% list$id,]

length(unique(mr_res2$Gene))

mr_res3<-mr_res2[mr_res2$pval<(0.05/(length(unique(mr_res2$id))*7)),]

length(unique(mr_res3$Gene))

fwrite(mr_res3,"CRC_TWAS/MR_results_Strong.csv")
