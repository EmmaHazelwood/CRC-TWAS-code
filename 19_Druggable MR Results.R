library(data.table)
library(dplyr)
library(readxl)

firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

res<-fread("Results/MR_results_druggable_genome.csv")
res2<-res %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
res2$gene<-substr(res2$Gene,1,15)
res2<-res2[res2$gene %in% druggable_list,]
res3<-res2[res2$pval<(0.05/(length(unique(res2$gene)))),]

length(unique(res3$Gene))
length(unique(res2$Gene))

drugs<-read_excel("CRC_TWAS/41591_2021_1310_MOESM3_ESM.xlsx",sheet=4,skip=27)
drugs<-drugs %>% tidyr::separate(DRUG_IN_VITRO_ACTIVITY,c("Drug_ID","Drug_name"),sep=":")
drugs<-dplyr::select(drugs,ENSEMBL_GENE_ID,Drug_ID,Drug_name,UNIPROT)

drugs2<-read_excel("CRC_TWAS/41573_2017_BFnrd2016230_MOESM16_ESM.xlsx")

res3$Gene<-substr(res3$Gene,1,15)
res_drugs<-merge(res3,drugs,by.x="Gene",by.y="ENSEMBL_GENE_ID",all.x=T)
length(unique(res_drugs$Gene))

res_drugs<-merge(res_drugs,drugs2,by.x="UNIPROT",by.y="ACCESSION",all.x=T,incomparables = NA)
res_drugs<-res_drugs[!duplicated(res_drugs),]

res_drugs$OR<-format(round(exp(res_drugs$b),digits=2),nsmall=2)
res_drugs$LI<-format(round(exp(res_drugs$b-1.96*res_drugs$se),digits=2),nsmall=2)
res_drugs$UI<-format(round(exp(res_drugs$b+1.96*res_drugs$se),digits=2),nsmall=2)
res_drugs$CI<-paste(res_drugs$LI," to ",res_drugs$UI,sep="")
res_drugs<-dplyr::select(res_drugs,Gene,UNIPROT,Tissue,subtype,method,Drug_name,MECHANISM_OF_ACTION,OR,CI,pval)

res_drugs$Drug_name<-firstup(res_drugs$Drug_name)
res_drugs$MECHANISM_OF_ACTION<-firstup(res_drugs$MECHANISM_OF_ACTION)

colnames(res_drugs)<-c("Ensembl ID","Uniprot ID","Tissue","Subtype","Method","Drug name","Mechanism of action","OR","95% Confidence interval","P value")

res_drugs<-res_drugs[!duplicated(res_drugs),]

fwrite(res_drugs,"Results/Druggable_genes_strong.csv")

res<-fread("Results/MR_results_druggable_genome.csv")
res3<-res %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
drugs<-read_excel("CRC_TWAS/41591_2021_1310_MOESM3_ESM.xlsx",sheet=4,skip=27)
drugs<-drugs %>% tidyr::separate(DRUG_IN_VITRO_ACTIVITY,c("Drug_ID","Drug_name"),sep=":")
drugs<-dplyr::select(drugs,ENSEMBL_GENE_ID,Drug_ID,Drug_name,UNIPROT)

drugs2<-read_excel("CRC_TWAS/41573_2017_BFnrd2016230_MOESM16_ESM.xlsx")

res3$Gene<-substr(res3$Gene,1,15)
res_drugs<-merge(res3,drugs,by.x="Gene",by.y="ENSEMBL_GENE_ID",all.x=T)
length(unique(res_drugs$Gene))

res_drugs<-merge(res_drugs,drugs2,by.x="UNIPROT",by.y="ACCESSION",all.x=T,incomparables = NA)
res_drugs<-res_drugs[!duplicated(res_drugs),]

res_drugs$OR<-format(round(exp(res_drugs$b),digits=2),nsmall=2)
res_drugs$LI<-format(round(exp(res_drugs$b-1.96*res_drugs$se),digits=2),nsmall=2)
res_drugs$UI<-format(round(exp(res_drugs$b+1.96*res_drugs$se),digits=2),nsmall=2)
res_drugs$CI<-paste(res_drugs$LI," to ",res_drugs$UI,sep="")
res_drugs<-dplyr::select(res_drugs,Gene,UNIPROT,Tissue,subtype,method,Drug_name,MECHANISM_OF_ACTION,OR,CI,pval)

res_drugs$Drug_name<-firstup(res_drugs$Drug_name)
res_drugs$MECHANISM_OF_ACTION<-firstup(res_drugs$MECHANISM_OF_ACTION)

colnames(res_drugs)<-c("Ensembl ID","Uniprot ID","Tissue","Subtype","Method","Drug name","Mechanism of action","OR","95% Confidence interval","P value")

res_drugs<-res_drugs[!duplicated(res_drugs),]


res<-fread("Results/2_Colocalization/Coloc_results_abf_druggable.csv")
res<-res %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
res$`Ensembl ID`<-substr(res$Gene,1,15)
res$Subtype<-firstup(res$subtype)

sig<-res[res$PP.H4.abf>0.8,]

res2<-merge(res,res_drugs,by=c("Ensembl ID","Subtype","Tissue"))
res2$col<-paste(res2$`Ensembl ID`,res2$Subtype,res2$Tissue)

res$col<-paste(res$`Ensembl ID`,res$Subtype,res$Tissue)

res3<-res[res$col %in% res2$col,]
res3<-res3[,-(ncol(res3))]

res<-rbind(sig,res3)

res<-merge(res,res_drugs,by=c("Ensembl ID","Subtype","Tissue"),all.x=T)

sig<-res[res$PP.H4.abf>0.8,]
sig<-sig[sig$`P value`<=(0.05/(556)) | is.na(sig$`P value`),]

sig<-distinct(sig)
sig<-sig[!is.na(sig$`Ensembl ID`),]

fwrite(sig,"Results/3_MR/Druggable_both.csv")
