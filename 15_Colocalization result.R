library(data.table)
library(dplyr)

firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

res<-fread("CRC_TWAS/Coloc_results_abf.csv")

res<-res %>% tidyr::separate(id,c("Gene","Tissue","subtype"),sep=";")
res$subtype<-firstup(res$subtype)

sig<-res[res$PP.H4.abf>0.8,]
length(unique(sig$Gene))

mr_res<-fread("CRC_TWAS/MR_results_Strong.csv")

both<-merge(sig,mr_res,by=c("Gene","Tissue","subtype"))

length(unique(both$Gene))

length(unique(mr_res$Gene))
length(unique(res$Gene))

length(intersect(unique(mr_res$Gene),unique(res$Gene)))

setdiff(unique(mr_res$Gene),unique(res$Gene))

all<-merge(mr_res,res,by=c("Gene","Tissue","subtype"),all.x=T)

#Find gene-tissue-subtypes missing from MultiXcan
splice<-fread("CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
splice$analysis<-"splice"
multi<-fread("CRC_TWAS/SMultiXcan_Expression_All_Strong.csv")
multi$analysis<-"multi"
jti<-fread("CRC_TWAS/JTI_All_Strong.csv")
jti$analysis<-"jti"
jti<-jti[!jti$subtype=="",]
jti$subtype[grepl("female",jti$file)]<-"female"
jti$subtype<-firstup(jti$subtype)
genes<-c(splice$gene,substr(multi$gene,1,15),jti$gene)
genes<-unique(genes)

splice<-dplyr::select(splice,gene,t_i_best,subtype,analysis)
multi<-dplyr::select(multi,gene,t_i_best,subtype,analysis)
jti<-dplyr::select(jti,gene,tissue,subtype,analysis)
list<-rbind(splice,multi,jti,use.names=F)
list<-list[!duplicated(list),]
list<-list[!list$gene=="",]

list$gene<-substr(list$gene,1,15)

list2<-merge(list,res,by.x=c("gene","t_i_best","subtype"),by.y=c("Gene","Tissue","subtype"),all.x=T)


#Looking at consistent directions across all 3 (note: excluding splicing TWAS from this as don't know direction of effect on expression)
splice<-fread("CRC_TWAS/SMultiXcan_Splicing_Expression_All.csv")
multi<-fread("CRC_TWAS/SMultiXcan_Expression_All.csv")
jti<-fread("CRC_TWAS/JTI_All.csv")
jti$subtype[grepl("female",jti$file)]<-"female"

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

druggable_list<-unique(twas$gene[twas$pvalue<0.05])

twas$id<-paste(twas$gene,twas$t_i_best,twas$subtype,sep=";")
twas_exp<-twas[!twas$analysis=="splice",]




mr_res_sig<-fread("CRC_TWAS/MR_results_Strong.csv")

mr_res_sig<-dplyr::select(mr_res_sig,Gene,Tissue,b,pval,se,subtype)

twas_mr<-merge(twas,mr_res_sig,by.x=c("gene","t_i_best","subtype"),by.y=c("Gene","Tissue","subtype"))

twas_mr$twas_mr_consistent<-"No"
twas_mr$twas_mr_consistent[twas_mr$z_mean<0 & twas_mr$b<0]<-"Yes"
twas_mr$twas_mr_consistent[twas_mr$z_mean>0 & twas_mr$b>0]<-"Yes"
twas_mr$twas_mr_consistent[twas_mr$analysis=="splice"]<-"NA"

sig<-dplyr::select(sig,Gene,Tissue,subtype,PP.H4.abf)

twas_mr_coloc<-merge(twas_mr,sig,by.x=c("gene","t_i_best","subtype"),by.y=c("Gene","Tissue","subtype"))

nrow(twas_mr_coloc)
table(twas_mr_coloc$twas_mr_consistent)

length(unique(twas_mr_coloc$gene[twas_mr_coloc$twas_mr_consistent=="Yes"]))
       
twas_exp<-twas_mr_coloc[!twas_mr_coloc$analysis=="splice",]
#Consistent in both TWAS where present
list<-unique(twas_exp$id[duplicated(twas_exp$id)])

dup<-twas_exp[twas_exp$id %in% list,]
dup<-dplyr::select(dup,id,analysis,z_mean)

same<-data.frame(matrix(ncol=2,nrow=length(unique(dup$id))))
colnames(same)<-c("id","consistent")

same$id<-unique(dup$id)
same$consistent<-"No"

for (b in unique(dup$id)){
  df<-dup[dup$id==b]
  if(df$z_mean[1]>0 & df$z_mean[2]>0){
    same$consistent[same$id==b]<-"Yes"
  }
  if(df$z_mean[1]<0 & df$z_mean[2]<0){
    same$consistent[same$id==b]<-"Yes"
  }
}


list<-data.frame(list)
list<-list%>%tidyr::separate(list,c("Gene","Tissue","subtype"),sep=";")
length(unique(list$Gene))
twas_mr_coloc$OR<-format(round(exp(twas_mr_coloc$b),digits=2),nsmall=2)
twas_mr_coloc$LI<-format(round(exp(twas_mr_coloc$b-1.96*twas_mr_coloc$se),digits=2),nsmall=2)
twas_mr_coloc$UI<-format(round(exp(twas_mr_coloc$b+1.96*twas_mr_coloc$se),digits=2),nsmall=2)
twas_mr_coloc$CI<-paste(twas_mr_coloc$LI," to ",twas_mr_coloc$UI,sep="")

library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","external_gene_name"), values=twas_mr_coloc$gene, mart=ensembl)
twas_mr_coloc<-merge(twas_mr_coloc,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T)

#Clean up file to save
twas_mr_coloc<-dplyr::select(twas_mr_coloc,gene,external_gene_name,t_i_best,subtype,z_mean,pvalue,analysis,PP.H4.abf,OR,CI,pval,twas_mr_consistent)
colnames(twas_mr_coloc)<-c("Ensembl ID","Gene","Tissue","Subtype","Z-score TWAS","P value TWAS","TWAS analysis","Posterior probability H4 colocalization","Odds ratio MR","95% confidence interval MR","P value MR","Consistent in TWAS and MR?")


firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

twas_mr_coloc$Subtype<-firstup(twas_mr_coloc$Subtype)

twas_mr_coloc$`TWAS analysis`[twas_mr_coloc$`TWAS analysis`=="jti"]<-"JTI"
twas_mr_coloc$`TWAS analysis`[twas_mr_coloc$`TWAS analysis`=="multi"]<-"S-MultiXcan Expression"
twas_mr_coloc$`TWAS analysis`[twas_mr_coloc$`TWAS analysis`=="splice"]<-"S-MultiXcan Splicing"

list<-unique(twas_mr_coloc$`Ensembl ID`[twas_mr_coloc$Gene==""])
twas_mr_coloc$Gene[twas_mr_coloc$`Ensembl ID`=="ENSG00000261888"]<-"Novel transcript ENSG00000261888"
twas_mr_coloc$Gene[twas_mr_coloc$`Ensembl ID`=="ENSG00000262003"]<-"LOC101927727"
twas_mr_coloc$Gene[twas_mr_coloc$`Ensembl ID`=="ENSG00000272334"]<-"Novel transcript ENSG00000272334"
twas_mr_coloc$Gene[twas_mr_coloc$`Ensembl ID`=="ENSG00000273619"]<-"Novel Transcript ENSG00000273619, Antisense To RPS21"
twas_mr_coloc$Gene[twas_mr_coloc$`Ensembl ID`=="ENSG00000275437"]<-"Novel Transcript ENSG00000275437, Sense Intronic To CABLES2"




fwrite(twas_mr_coloc,"Results/Consistent_all_3.csv")




# Splicing events ---------------------------------------------------------

res<-fread("CRC_TWAS/Coloc_results_abf_splice.csv")

res_sig<-res[res$PP.H4.abf>0.8,]

res$id_subtype<-paste(res$id,res$subtype,sep=";")

splice<-fread("CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
length(unique(splice$gene_name))
length(unique(splice$gene))
splice$id<-paste(splice$chr,splice$start,splice$end,splice$t_i_best,sep=":")
splice<-dplyr::select(splice,id,gene)
splice<-splice[!duplicated(splice$id),]
splice$id_subtype<-paste(splice$id,splice$subtype,sep=";")

res<-merge(res,splice,by="id",allow.cartesian=TRUE,all.y=TRUE)
res_sig<-res[res$PP.H4.abf>0.8,]


length(unique(res$id))
length(unique(res$gene))

length(unique(res_sig$gene_name))
length(unique(res_sig$gene))


twas_exp<-res_sig[!res_sig$analysis=="splice",]



fwrite(res_sig,"CRC_TWAS/Coloc_results_abf_splice_strong.csv")
