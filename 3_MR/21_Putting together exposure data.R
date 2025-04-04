library(data.table)
library(dplyr)

exp1<-fread("exposure_dat_temp_exp.csv")
exp2<-fread("clumped_exposure_data_druggable.csv")
exp3<-fread("clumped_exposure_dat_temp_sqtls.csv")


exp<-rbind(exp1,exp2)
exp<-exp%>%tidyr::separate(exposure,c("Gene.version","tissue"),sep=";")
exp3<-exp3%>%tidyr::separate(exposure,c("Gene.version","tissue"),sep=";")
exp3$gene.exposure<-paste(exp3$Gene.version,exp3$gene.exposure,sep=", ")

exp<-rbind(exp,exp3)

exp$id<-paste(exp$gene.exposure,exp$tissue,sep=";")
exp<-select(exp,id,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exp)<-c("Exposure","SNP","Effect_allele","Other_allele","EAF","Beta","SE","P-value")

fwrite(exp,"Exposure_data.csv")
