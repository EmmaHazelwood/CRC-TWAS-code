library(data.table)
library(dplyr)

mr_res<-fread("data/CRC_TWAS/MR_results.csv")

mr_res2<-mr_res %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
mr_res3<-mr_res2[mr_res2$pval<(0.05/(length(unique(mr_res$exposure))*7)),]
length(unique(mr_res3$Gene))
length(unique(mr_res2$Gene))


#Filter for ones with effects in TWAS
splice<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
multi<-fread("data/CRC_TWAS/SMultiXcan_Expression_All_Strong.csv")
jti<-fread("data/CRC_TWAS/JTI_All_Strong.csv")
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

fwrite(mr_res3,"data/CRC_TWAS/MR_results_Strong.csv")


# Splicing genes ----------------------------------------------------------
mr_res_sqtl<-fread("data/CRC_TWAS/MR_results_sqtls.csv")
mr_res_sqtl2<-mr_res_sqtl %>% tidyr::separate(exposure,c("id","Tissue"),sep=";")


splice<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")

splice$id<-paste(splice$chr,splice$start,splice$end,sep=":")
length(unique(splice$gene_name))
splice$id_tissue<-paste(splice$id,splice$t_i_best,sep=";")
splice<-dplyr::select(splice,id,gene,id_tissue)
splice<-splice[!duplicated(splice$id_tissue),]

mr_res_sqtl2<-merge(mr_res_sqtl2,splice,by.x="id.exposure",by.y="id_tissue",allow.cartesian=TRUE,all.x=TRUE)
mr_res_sqtl2<-mr_res_sqtl2[!is.na(mr_res_sqtl2$b),]
mr_res_sqtl2$id_subtype<-paste(mr_res_sqtl2$id.exposure,mr_res_sqtl2$subtype,sep=";")
mr_res_sqtl_sig<-mr_res_sqtl2[mr_res_sqtl2$pval<(0.05/(length(unique(mr_res_sqtl2$id_subtype)))),]

length(unique(mr_res_sqtl2$id.x))
length(unique(mr_res_sqtl_sig$id.x))
length(unique(mr_res_sqtl_sig$gene))

