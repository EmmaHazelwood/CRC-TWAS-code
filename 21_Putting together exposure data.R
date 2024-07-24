library(data.table)
library(dplyr)

exp1<-fread("CRC_TWAS/clumped_exposure_dat_temp.csv")
exp2<-fread("CRC_TWAS/clumped_exposure_data_druggable.csv")

exp<-rbind(exp1,exp2)
exp<-exp%>%tidyr::separate(exposure,c("Gene.version","tissue"),sep=";")
exp$id<-paste(exp$gene.exposure,exp$tissue,sep=";")
exp<-select(exp,id,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exp)<-c("Exposure","SNP","Effect_allele","Other_allele","EAF","Beta","SE","P-value")

fwrite(exp,"CRC_TWAS/Exposure_data.csv")
