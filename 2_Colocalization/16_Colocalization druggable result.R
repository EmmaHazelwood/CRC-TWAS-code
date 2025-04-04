library(data.table)
library(dplyr)

res<-fread("Results/2_Colocalization/Coloc_results_abf_druggable.csv")

res<-res %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")

length(unique(res$Gene))

res$gene<-substr(res$Gene,1,15)

res2<-res[res$gene %in% druggable_list,]

length(unique(res2$Gene))


sig<-res[res$PP.H4.abf>0.8,]
length(unique(sig$Gene))



sig<-res2[res2$PP.H4.abf>0.8,]
length(unique(sig$Gene))
