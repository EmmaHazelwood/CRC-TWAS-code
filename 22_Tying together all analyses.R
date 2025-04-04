library(data.table)
library(dplyr)
library(biomaRt)

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

res<-fread("data/CRC_TWAS/Coloc_results_abf.csv")

res<-res %>% tidyr::separate(id,c("Gene","Tissue","subtype"),sep=";")

length(unique(sig$Gene))

mr_res<-fread("data/CRC_TWAS/MR_results.csv")

mr_res<-mr_res %>% tidyr::separate(exposure,c("Gene","Tissue","subtype"),sep=";")

both<-merge(res,mr_res,by=c("Gene","Tissue","subtype"),all.x=T)

all<-merge(twas,both,by.x=c("gene","t_i_best","subtype"),by.y=c("Gene","Tissue","subtype"),all.x=T)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",host = "http://www.ensembl.org")                                                                         

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(all$gene), mart=ensembl)

all<-merge(all,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T)

#Manually annotate ones which ensembl missed
unique(all$gene[all$hgnc_symbol==""])
all$hgnc_symbol[all$gene=="ENSG00000257298"]<-"Novel Transcript, Sense Intronic To LIMA1"
all$hgnc_symbol[all$gene=="ENSG00000261888"]<-"Novel Transcript"
all$hgnc_symbol[all$gene=="ENSG00000262003"]<-"Uncharacterized LOC101927727"
all$hgnc_symbol[all$gene=="ENSG00000271993"]<-"Novel Transcript, Antisense To LRRFIP2"
all$hgnc_symbol[all$gene=="ENSG00000272334"]<-"Novel Transcript"
all$hgnc_symbol[all$gene=="ENSG00000272368"]<-"Novel Transcript, Antisense To CERS5"
all$hgnc_symbol[all$gene=="ENSG00000273619"]<-"Novel Transcript, Antisense To RPS21"
all$hgnc_symbol[all$gene=="ENSG00000274370"]<-"Novel Transcript"
all$hgnc_symbol[all$gene=="ENSG00000275437"]<-"Novel Transcript, Sense Intronic To CABLES2"

fwrite(all,"data/CRC_TWAS/all_3_together.csv")

all<-all[all$PP.H4.abf>=0.8,]
all<-all[all$pval<(0.05/nrow(mr_res)) | is.na(all$pval),]


fwrite(all,"data/CRC_TWAS/Consistent_all_3.csv")



# List to send  ------------------------------------------------------
#all$z_mean[all$analysis=="splice"]<-"splice"

#all<-dplyr::select(all,hgnc_symbol,gene,z_mean,t_i_best,subtype)
#all<-all[order(all$hgnc_symbol),]


#fwrite(all,"data/CRC_TWAS/Gene_list.csv")


# Druggable ---------------------------------------------------------------
drug<-fread("data/CRC_TWAS/Druggable_both.csv")

list<-unique(c(all$gene,drug$`Ensembl ID`))
length(intersect(all$gene,drug$`Ensembl ID`))


# Splicing ----------------------------------------------------------------

mr_res_sqtl2 #from MR results script

res #from colocalisation results script

res$id_subtype.y<-gsub(";","",res$id_subtype.y)

mr_res_sqtl2$id.exposure<-gsub(";",":",mr_res_sqtl2$id.exposure)
mr_res_sqtl2$merge_col<-paste(mr_res_sqtl2$id.exposure,mr_res_sqtl2$subtype,sep=";")  

res$merge_col<-paste(res$id_subtype.y,tolower(res$subtype),sep=";")  

mr_col_spl<-merge(res,mr_res_sqtl2,by="merge_col",all=T,allow.cartesian=FALSE)
mr_col_spl<-distinct(mr_col_spl)

mr_sig_col_spl<-mr_col_spl[mr_col_spl$pval<0.001190476,]

length(unique(mr_sig_col_spl$id[mr_sig_col_spl$PP.H4.abf>0.8]))
length(unique(mr_sig_col_spl$gene.x[mr_sig_col_spl$PP.H4.abf>0.8]))
