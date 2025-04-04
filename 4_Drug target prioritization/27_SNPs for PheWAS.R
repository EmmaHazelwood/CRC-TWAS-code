library(data.table)
library(dplyr)

df<-fread("working/data/CRC_TWAS/Exposure_data.csv")

genes_twas<-fread("working/data/CRC_TWAS/Consistent_all_3.csv")
genes_twas$Exposure<-paste(genes_twas$gene,genes_twas$t_i_best,sep=";")
genes_drug<-fread("Results/3_MR/Druggable_both.csv")
genes_drug$Exposure<-paste(genes_drug$`Ensembl ID`,genes_drug$Tissue,sep=";")

genes<-unique(c(genes_twas$Exposure,genes_drug$Exposure))
genes2<-unique(c(genes_twas$gene,genes_drug$`Ensembl ID`))

df$Exposure2<-df$Exposure
df<-df %>% tidyr::separate(Exposure2,c("Gene","Tissue"),sep = ";")

snps<-df[df$Exposure %in% genes,]

snps<-snps[-c(347,1255,1256,1258)]

length(unique(snps$Gene))
length(unique(snps$SNP))

snplist<-unique(snps$SNP)
snplist<-snplist[order(snplist)]

#Check any missing
phewas<-list.files("Results/4_Drug target prioritization/pheWAS/")

phewas<-gsub(".csv","",phewas)


setdiff(phewas,snplist)
setdiff(snplist,phewas)


