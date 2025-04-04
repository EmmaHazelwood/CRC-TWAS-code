library(dplyr)
library(data.table)

splice<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
multi<-fread("data/CRC_TWAS/SMultiXcan_Expression_All_Strong.csv")
jti<-fread("data/CRC_TWAS/JTI_All_Strong.csv")

splice<-dplyr::select(splice,gene,t_i_best,subtype,z_mean,pvalue)
splice$subtype<-tolower(splice$subtype)
splice$analysis<-"splice"
splice<-splice[!splice$gene=="",]
splice<-splice[!splice$subtype=="",]
multi<-dplyr::select(multi,gene,t_i_best,subtype,z_mean,pvalue)
multi$subtype<-tolower(multi$subtype)
multi$analysis<-"multi"
multi$gene<-substr(multi$gene,1,15)
jti<-dplyr::select(jti,gene,tissue,subtype,zscore,pvalue)

jti$analysis<-"jti"

twas<-rbind(splice,multi,jti,use.names=FALSE)
twas<-distinct(twas)
twas<-twas[!twas$subtype=="",]
twas<-twas[!twas$analysis=="splice",]

res<-fread("data/CRC_TWAS/Coloc_results_abf.csv")

res<-res %>% tidyr::separate(id,c("Gene","Tissue","subtype"),sep=";")

mr_res<-fread("data/CRC_TWAS/MR_results.csv")

mr_res<-mr_res %>% tidyr::separate(exposure,c("Gene","Tissue","subtype"),sep=";")

both<-merge(res,mr_res,by=c("Gene","Tissue","subtype"),all.x=T)


all<-merge(twas,both,by.x=c("gene","t_i_best","subtype"),by.y=c("Gene","Tissue","subtype"),all.x=T)

ov<-all[all$subtype=="overall",]

library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(ov$gene), mart=ensembl)

ov<-merge(ov,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T)

ov<-ov[ov$PP.H4.abf>=0.5,]
ov_up<-ov[ov$z_mean>0,]
ov_down<-ov[ov$z_mean<0,]


length(unique(ov_up$hgnc_symbol))
length(unique(ov_down$hgnc_symbol))





