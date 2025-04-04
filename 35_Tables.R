library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)

firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

setwd("Results/")

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", " x 10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

missing<-fread("missing_with_gene_names_all.txt")

# Gene names --------------------------------------------------------------

tb_2<-fread("Consistent_all_3.csv")
hgnc_names<-distinct(dplyr::select(tb_2,gene,hgnc_symbol))


# Table 1 -----------------------------------------------------------------
mr_res<-fread("data/CRC_TWAS/MR_results.csv")
mr_res2<-mr_res %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")

res<-fread("data/CRC_TWAS/Coloc_results_abf.csv")
res<-res %>% tidyr::separate(id,c("Gene","Tissue","subtype"),sep=";")
res$id<-paste(res$Gene,res$Tissue,res$subtype,sep=";")
mr_res2$id<-paste(mr_res2$id.exposure,mr_res2$subtype,sep=";")
both<-merge(res,mr_res2,by="id",all.x=TRUE)

both<-both[both$PP.H4.abf>0.8,]
expression_TWAS<-both[is.na(both$pval) | both$pval<(0.05/(length(unique(mr_res$exposure))*7)),]

#Splicing specifically
res<-fread("data/CRC_TWAS/Coloc_results_abf_splice.csv")
res_sig<-res[res$PP.H4.abf>0.8,]
res$id_subtype<-paste(res$id,res$subtype,sep=";")
splice<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
splice$id<-paste(splice$chr,splice$start,splice$end,splice$t_i_best,sep=":")
splice<-dplyr::select(splice,id,gene,gene_name,z_mean,pvalue)
splice<-splice[!duplicated(splice$id),]
splice$id_subtype<-paste(splice$id,splice$subtype,sep=";")

res<-merge(res,splice,by="id",allow.cartesian=TRUE,all.y=TRUE)
splice_TWAS<-res[res$PP.H4.abf>0.8,]

#druggable
res<-fread("Results/2_Colocalization/Coloc_results_abf_druggable.csv")
res<-res %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
res$gene<-substr(res$Gene,1,15)

#Limit to those with nominal significance in TWAS
#Read in full TWAS results
splice<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All.csv")
multi<-fread("data/CRC_TWAS/SMultiXcan_Expression_All.csv")
jti<-fread("data/CRC_TWAS/JTI_All.csv")

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

#Get list of P<0.05
druggable_list<-unique(twas$gene[twas$pvalue<0.05])
fwrite(list(druggable_list),"Druggable_list.csv",quote = F,row.names = F)

#Filter for those with P<0.05 in TWAS
druggable<-res[res$gene %in% druggable_list,]

splice_strong<-fread("data/CRC_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
multi_strong<-fread("data/CRC_TWAS/SMultiXcan_Expression_All_Strong.csv")
jti_strong<-fread("data/CRC_TWAS/JTI_All_Strong.csv")

splice_strong<-dplyr::select(splice_strong,gene,t_i_best,subtype,z_mean,pvalue)
splice_strong$subtype<-tolower(splice_strong$subtype)
splice_strong$analysis<-"splice_strong"
splice_strong<-splice_strong[!splice_strong$gene=="",]
splice_strong<-splice_strong[!splice_strong$subtype=="",]
multi_strong<-dplyr::select(multi_strong,gene,t_i_best,subtype,z_mean,pvalue)
multi_strong$subtype<-tolower(multi_strong$subtype)
multi_strong$analysis<-"multi_strong"
multi_strong$gene<-substr(multi_strong$gene,1,15)
jti_strong<-dplyr::select(jti_strong,gene,tissue,subtype,zscore,pvalue)

jti_strong$analysis<-"jti_strong"

twas_strong<-rbind(splice_strong,multi_strong,jti_strong,use.names=FALSE)

both<-merge(both,twas_strong,by.x=c("Gene.x","Tissue.x","subtype.x"),by.y=c("gene","t_i_best","subtype"),all.x=TRUE,allow.cartesian=FALSE)
both<-merge(both,hgnc_names,by.x="Gene.x",by.y="gene")

both<-dplyr::select(both,hgnc_symbol,Gene.x,Tissue.x,subtype.x,analysis,z_mean,pvalue,PP.H4.abf,b,pval)
both$b<-exp(both$b)
colnames(both)<-c("Gene", "Ensembl ID","Tissue","Subtype","TWAS analysis","TWAS Z-score","TWAS P-value","H4 posterior probability","MR odds ratio","MR P-value")
both$`TWAS analysis`[both$`TWAS analysis`=="multi_strong"]<-"Expression MultiXcan"
both$`TWAS analysis`[both$`TWAS analysis`=="splice_strong"]<-"Splicing MultiXcan"
both$`TWAS analysis`[both$`TWAS analysis`=="jti_strong"]<-"Expression JTI"
both$Subtype<-firstup(both$Subtype)

both$`TWAS Z-score`<-format(round(both$`TWAS Z-score`,digits=2),nsmall=2)
both$`H4 posterior probability`<-format(round(both$`H4 posterior probability`,digits=2),nsmall=2)
both$`MR odds ratio`<-format(round(both$`MR odds ratio`,digits=2),nsmall=2)
both$`TWAS P-value`<-format(both$`TWAS P-value`,digits=3)
both$`TWAS P-value`<-changeSciNot(both$`TWAS P-value`)
both$`MR P-value`<-format(both$`MR P-value`,digits=3)
both$`MR P-value`<-changeSciNot(both$`MR P-value`)
both$`MR P-value`[which(both$`MR P-value`=="      NA")]<-NA
both$`TWAS P-value`[which(both$`TWAS P-value`=="      NA")]<-NA
both$`MR odds ratio`[which(both$`MR odds ratio`=="  NA")]<-NA

fwrite(both,"Table 1 r output.csv")

#Druggable genome
#druggable_list
druggable_list<-unique(twas$gene[twas$pvalue<0.05])

#Read in MR results
res<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Analysis/Results/MR_results_druggable_genome.csv")
res2<-res %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
res2$gene<-substr(res2$Gene,1,15)
res2<-res2[res2$gene %in% druggable_list,]
druggable_mr_res<-res2

#druggable
res<-fread("Results/2_Colocalization/Coloc_results_abf_druggable.csv")
res<-res %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
res$gene<-substr(res$Gene,1,15)
res2<-res[res$gene %in% druggable_list,]
sig<-res2[res2$PP.H4.abf>0.8,]

druggable_mr_res$subtype<-tolower(druggable_mr_res$subtype)
sig<-merge(sig,druggable_mr_res,by=c("gene","subtype","Tissue"),all.x=TRUE)
sig<-sig[is.na(sig$pval) | sig$pval<(0.05/380),]

library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",host = "http://www.ensembl.org")                                                                         
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=sig$gene, mart=ensembl)

sig<-merge(sig,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T)
sig$analysis<-"Druggable genome"
sig$z_mean<-NA
sig$pvalue<-NA

sig<-dplyr::select(sig,hgnc_symbol,gene,Tissue,subtype,analysis,z_mean,pvalue,PP.H4.abf,b,pval)
sig$b<-exp(sig$b)
colnames(sig)<-c("Gene", "Ensembl ID","Tissue","Subtype","TWAS analysis","TWAS Z-score","TWAS P-value","H4 posterior probability","MR odds ratio","MR P-value")
sig$Subtype<-firstup(sig$Subtype)
sig$`H4 posterior probability`<-format(round(sig$`H4 posterior probability`,digits=2),nsmall=2)
sig$`MR odds ratio`<-format(round(sig$`MR odds ratio`,digits=2),nsmall=2)
sig$`MR P-value`<-format(sig$`MR P-value`,digits=3)
sig$`MR P-value`<-changeSciNot(sig$`MR P-value`)
sig$`MR P-value`[which(sig$`MR P-value`=="      NA")]<-NA
sig$`MR odds ratio`[which(sig$`MR odds ratio`=="  NA")]<-NA

fwrite(sig,"Table 1 r output druggable.csv")

# Table 2 -----------------------------------------------------------------

splice_TWAS<-splice_TWAS %>% tidyr::separate(id,c("chr","bp start","bp end","Tissue"),sep=":")
splice_TWAS<-merge(splice_TWAS,hgnc_names,by="gene",all.x=TRUE)

splice_TWAS<-dplyr::select(splice_TWAS,hgnc_symbol,gene,Tissue,subtype,z_mean,pvalue,PP.H4.abf)
colnames(splice_TWAS)<-c("Gene","Ensembl ID","Tissue","Subtype","TWAS Z-score","TWAS P-value","H4 posterior probability")

splice_TWAS$`TWAS Z-score`<-format(round(splice_TWAS$`TWAS Z-score`,digits=2),nsmall=2)
splice_TWAS$`H4 posterior probability`<-format(round(splice_TWAS$`H4 posterior probability`,digits=2),nsmall=2)

splice_TWAS$`TWAS P-value`<-format(splice_TWAS$`TWAS P-value`,digits=3)
splice_TWAS$`TWAS P-value`<-changeSciNot(splice_TWAS$`TWAS P-value`)
splice_TWAS$`TWAS P-value`[which(splice_TWAS$`TWAS P-value`=="      NA")]<-NA

fwrite(splice_TWAS,"Table 2 r output.csv")

# Supplementary table 1 ---------------------------------------------------
st_1a<-readxl::read_excel("gwas.xlsx")
st_1a<-distinct(st_1a)

# Supplementary table 2 -----------------------------------------------------------------

tb_2<-fread("Consistent_all_3.csv")
tb_2$OR<-exp(tb_2$b)
tb_2$LCI<-format(round(exp(tb_2$b - 1.96*tb_2$se),digits=2),nsmall=2)
tb_2$UCI<-format(round(exp(tb_2$b + 1.96*tb_2$se),digits=2),nsmall=2)
tb_2$CI<-paste(tb_2$LCI, " to ",tb_2$UCI,sep="")
tb_2$CI[which(tb_2$CI=="NA to NA")]<-NA
tb_2$PP.H0.abf<-format(round(tb_2$PP.H0.abf,digits=3),nsmall=3)
tb_2$PP.H1.abf<-format(round(tb_2$PP.H1.abf,digits=3),nsmall=3)
tb_2$PP.H2.abf<-format(round(tb_2$PP.H2.abf,digits=3),nsmall=3)
tb_2$PP.H3.abf<-format(round(tb_2$PP.H3.abf,digits=3),nsmall=3)
tb_2$PP.H4.abf<-format(round(tb_2$PP.H4.abf,digits=3),nsmall=3)
tb_2$z_mean<-format(round(tb_2$z_mean,digits=2),nsmall=2)
tb_2$pvalue<-format(tb_2$pvalue,digits=3)
tb_2$pvalue<-changeSciNot(tb_2$pvalue)
tb_2$pval<-format(tb_2$pval,digits=3)
tb_2$pval<-changeSciNot(tb_2$pval)
tb_2$pval[which(tb_2$pval=="      NA")]<-NA
tb_2$OR<-format(round(tb_2$OR,digits=2),nsmall=2)
tb_2$t_i_best<-sub("_"," ",tb_2$t_i_best)
tb_2$subtype<-str_to_title(tb_2$subtype)
tb_2$analysis[which(tb_2$analysis=="multi")]<-"S-MutliXcan eQTLs"
tb_2$analysis[which(tb_2$analysis=="splice")]<-"S-MutliXcan sQTLs"
tb_2$analysis[which(tb_2$analysis=="jti")]<-"JTI eQTLs"
tb_2<-dplyr::select(tb_2,gene,hgnc_symbol,t_i_best,subtype,analysis,z_mean,pvalue,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,OR,CI,pval)

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(tb_2$gene), mart=ensembl)
missing<-fread("missing_with_gene_names_all.txt")
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)

tb_2<-merge(tb_2,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T)
tb_2<-dplyr::select(tb_2,gene,hgnc_symbol.y,t_i_best,subtype,analysis,z_mean,pvalue,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,OR,CI,pval)

colnames(tb_2)<-c("Ensembl gene ID","Gene","Tissue","Subtype","TWAS method","TWAS Z mean","TWAS P-value","H0","H1","H2","H3","H4","MR odds ratio","MR 95% confidence interval","MR P-value")
tb_2<-distinct(tb_2)
tb_2$`MR 95% confidence interval`[tb_2$`MR 95% confidence interval`=="  NA to   NA"]<-"NA"
tb_2$`MR P-value`[is.na(tb_2$`MR P-value`)]<-"NA"

# Table 1 -----------------------------------------------------------------

genes<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Analysis/Results/All_genes_with_loci.csv")
loci<-dplyr::select(genes,gene,Locus)
genes<-tb_2
genes<-merge(genes,loci,by.x="Ensembl gene ID",by.y="gene",all.x=T,allow.cartesian=TRUE)
fwrite(genes,"genes_with_loci.csv")


# Supplementary table 3 ---------------------------------------------------
st_1<-fread("1_TWAS/SMultiXcan_Splicing_Expression_All.csv")

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(st_1$gene), mart=ensembl)
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)

st_1<-merge(st_1,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T)
st_1$gene[st_1$gene_name=="intron_11_111356433_111357445"]<-"ENSG00000110777"
st_1$gene[st_1$gene_name=="intron_11_61802462_61802801"]<-"ENSG00000149485"
st_1$gene[st_1$gene_name=="intron_12_50130694_50130863"]<-"ENSG00000139624"
st_1$gene[st_1$gene_name=="intron_9_89418118_89424476"]<-"ENSG00000187764"
st_1$gene[st_1$gene_name=="intron_9_89418118_89424476"]<-"ENSG00000187764"

st_1$hgnc_symbol[st_1$gene=="ENSG00000213753"]<-"CENPBD2P"

st_1<-st_1[order(st_1$pvalue),]

st_1<-dplyr::select(st_1,gene,hgnc_symbol,indel,chr,start,end,t_i_best,subtype,z_mean,pvalue)

colnames(st_1)<-c("Ensembl gene ID","Gene","Splice event type","Splice event chromosome","Splice event start position","Splice event end position","Best tissue","Subtype","Z mean","P-value")

st_1$`Z mean`<-round(st_1$`Z mean`,digits=2)
st_1$`P-value`<-format(st_1$`P-value`,digits=3)

st_1<-st_1 %>%
  mutate(`Splice event chromosome` = str_remove_all(`Splice event chromosome`, "chr"))

st_1<-distinct(st_1)

# Supplementary table 4 ---------------------------------------------------
st_2<-fread("1_TWAS/SMultiXcan_Expression_All.csv")

st_2<-st_2[order(st_2$pvalue),]
st_2$gene<-substr(st_2$gene,1,15)

st_2<-dplyr::select(st_2,gene,gene_name,t_i_best,subtype,z_mean,pvalue)

st_2$gene_name[st_2$gene=="ENSG00000150750"]<-"POU2AF2"

colnames(st_2)<-c("Ensembl gene ID","Gene","Best tissue","Subtype","Z mean","P-value")

st_2$`Z mean`<-round(st_2$`Z mean`,digits=2)
st_2$`P-value`<-format(st_2$`P-value`,digits=3)
st_2<-distinct(st_2)

# Supplementary table 5 ---------------------------------------------------
st_3<-fread("1_TWAS/JTI_All.csv")

st_3<-st_3[order(st_3$pvalue),]

st_3<-dplyr::select(st_3,gene,gene_name,tissue,subtype,effect_size,pvalue)

st_3$gene_name[st_3$gene=="ENSG00000150750"]<-"POU2AF2"

colnames(st_3)<-c("Ensembl gene ID","Gene","Tissue","Subtype","Effect size","P-value")

st_3$`Effect size`<-round(st_3$`Effect size`,digits=2)
st_3$`P-value`<-format(st_3$`P-value`,digits=3)
st_3<-distinct(st_3)

# Supplementary table 6 ---------------------------------------------------
files<-list.files("1_TWAS/PrediXcan/SPrediXcan_results/")
files<-paste("1_TWAS/PrediXcan/SPrediXcan_results/",files,sep="")
files<-files[grepl("CRC_splicing.csv",files)]
st_4<-data.frame()
for (a in files){
  df<-fread(a)
  a<-sub("1_TWAS/PrediXcan/SPrediXcan_results/","",a)
  a<-sub(".csv","",a)
  df$analysis<-a
  st_4<-rbind(st_4,df)
}
st_4<-st_4[,-c(4,7,8,9)]
st_4<-distinct(st_4)

# Supplementary table 7 ---------------------------------------------------
files<-list.files("1_TWAS/PrediXcan/SPrediXcan_results/")
files<-paste("1_TWAS/PrediXcan/SPrediXcan_results/",files,sep="")
files<-files[grepl("CRC.csv",files)]
st_5<-data.frame()
for (a in files){
  df<-fread(a)
  a<-sub("1_TWAS/PrediXcan/SPrediXcan_results/","",a)
  a<-sub(".csv","",a)
  df$analysis<-a
  st_5<-rbind(st_5,df)
}
st_5<-st_5[,-c(4,7,8,9)]
st_5<-distinct(st_5)


st_6<-fread("3_MR/MR_results.csv")
st_6<-st_6 %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(st_6$Gene), mart=ensembl)
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)

st_6<-merge(st_6,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)

st_6$hgnc_symbol[st_6$Gene=="ENSG00000213753"]<-"CENPBD2P"

st_6$OR<-exp(st_6$b)

st_6<-dplyr::select(st_6,Gene,hgnc_symbol,Tissue,subtype,OR,pval)

st_6$OR<-round(st_6$OR,digits=2)
st_6$pval<-format(st_6$pval,digits=3)

st_6<-st_6[order(st_6$pval),]

colnames(st_6)<-c("Ensembl gene ID","Gene","Tissue","Subtype","Odds ratio","P-value")
st_6<-distinct(st_6)

st_7<-fread("3_MR/f_stat_and_r2.csv")
nsnps<-fread("3_MR/MR_results.csv")
nsnps<-nsnps %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
nsnps<-dplyr::select(nsnps,Gene,Tissue,nsnp)
st_7<-merge(st_7,nsnps,by=c("Gene","Tissue"))
st_7$r2<-round(st_7$r2,digits=2)
st_7$F_stat<-round(st_7$F_stat,digits=2)
st_7<-distinct(st_7)
min(st_7$F_stat)
max(st_7$F_stat)
median(st_7$F_stat)
st_7<-merge(st_7,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
st_7<-dplyr::select(st_7,Gene,hgnc_symbol,Tissue,r2,F_stat,nsnp)
colnames(st_7)<-c("Ensembl gene ID","Gene","Tissue","r2","F_stat","nsnp")
st_7<-distinct(st_7)

st_6<-merge(st_6,st_7,by=c("Ensembl gene ID","Tissue","Gene"),all.x = TRUE)

# Supplementary table 9 ---------------------------------------------------

st_8<-fread("3_MR/MR_results_druggable_genome.csv")
st_8<-st_8 %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
st_8$Gene<-substr(st_8$Gene,1,15)

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(st_8$Gene), mart=ensembl)
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)


st_8<-merge(st_8,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)

st_8$OR<-exp(st_8$b)
st_8<-st_8[order(st_8$pval),]

st_8<-dplyr::select(st_8,Gene,hgnc_symbol,Tissue,subtype,OR,pval)

st_8$OR<-round(st_8$OR,digits=2)
st_8$pval<-format(st_8$pval,digits=3)


colnames(st_8)<-c("Ensembl gene ID","Gene","Tissue","Subtype","Odds ratio","P-value")
st_8<-distinct(st_8)

st_9<-fread("3_MR/f_stat_and_r2_druggable.csv")
nsnps<-fread("3_MR/MR_results_druggable_genome.csv")
nsnps<-nsnps %>% tidyr::separate(exposure,c("Gene","Tissue"),sep=";")
nsnps<-dplyr::select(nsnps,Gene,Tissue,nsnp)
st_9<-merge(st_9,nsnps,by=c("Gene","Tissue"))
st_9$r2<-round(st_9$r2,digits=2)
st_9$F_stat<-round(st_9$F_stat,digits=2)
st_9$Gene<-substr(st_9$Gene,1,15)
st_9<-distinct(st_9)
min(st_9$F_stat)
max(st_9$F_stat)
median(st_9$F_stat)

st_9<-merge(st_9,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)

st_9<-dplyr::select(st_9,Gene,hgnc_symbol,Tissue,r2,F_stat,nsnp)
colnames(st_9)<-c("Ensembl gene ID","Gene","Tissue","r2","F_stat","nsnp")

st_9<-distinct(st_9)

st_8<-merge(st_8,st_9,by=c("Ensembl gene ID","Tissue","Gene"),all.x = TRUE)


# Supplementary table 10 --------------------------------------------------

st_10<-fread("2_Colocalization/Coloc_results_abf.csv")
st_10<-st_10 %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(st_10$Gene), mart=ensembl)
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)
st_10<-distinct(st_10)

st_10<-merge(st_10,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
st_10$hgnc_symbol[st_10$Gene=="ENSG00000213753"]<-"CENPBD2P"

st_10<-st_10[rev(order(st_10$PP.H4.abf)),]
st_10$PP.H0.abf<-format(round(st_10$PP.H0.abf,digits=3),nsmall=3)
st_10$PP.H1.abf<-format(round(st_10$PP.H1.abf,digits=3),nsmall=3)
st_10$PP.H2.abf<-format(round(st_10$PP.H2.abf,digits=3),nsmall=3)
st_10$PP.H3.abf<-format(round(st_10$PP.H3.abf,digits=3),nsmall=3)
st_10$PP.H4.abf<-format(round(st_10$PP.H4.abf,digits=3),nsmall=3)

st_10<-dplyr::select(st_10,Gene,hgnc_symbol,Tissue,subtype,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps)
st_10<-distinct(st_10)
colnames(st_10)<-c("Ensembl gene ID","Gene","Tissue","Subtype","H0","H1","H2","H3","H4","Number of SNPs")
st_10<-distinct(st_10)

# Supplementary table 11 --------------------------------------------------

st_11<-fread("2_Colocalization/Coloc_results_abf_splice.csv")
st_11<-st_11 %>% tidyr::separate(id,c("Chromosome","Start","End","Tissue"),sep=":")

st_11<-st_11 %>%
  mutate(Chromosome = str_remove_all(Chromosome, "chr"))

events<-dplyr::select(st_1,`Ensembl gene ID`,Gene,`Splice event`,`Splice event type`,`Splice event chromosome`,`Splice event start position`,`Splice event end position`)
colnames(events)<-c("Ensembl gene ID","Gene","Splice event","Splice event type","Chromosome","Start","End")

st_11<-merge(st_11,events,by=c("Chromosome","Start","End"),all.x=T)

st_11<-st_11[rev(order(st_11$PP.H4.abf)),]
st_11$PP.H0.abf<-format(round(st_11$PP.H0.abf,digits=3),nsmall=3)
st_11$PP.H1.abf<-format(round(st_11$PP.H1.abf,digits=3),nsmall=3)
st_11$PP.H2.abf<-format(round(st_11$PP.H2.abf,digits=3),nsmall=3)
st_11$PP.H3.abf<-format(round(st_11$PP.H3.abf,digits=3),nsmall=3)
st_11$PP.H4.abf<-format(round(st_11$PP.H4.abf,digits=3),nsmall=3)

st_11<-dplyr::select(st_11,`Ensembl gene ID`,Gene,`Splice event`,`Splice event type`,Chromosome,Start,End,Tissue,subtype,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps)
st_6$Gene[st_6$`Ensembl gene ID`=="ENSG00000213753"]<-"CENPBD2P"

st_11<-distinct(st_11)
colnames(st_11)<-c("Ensembl gene ID","Gene","Splice event","Splice event type","Splice event chromosome","Splice event start position","Splice event end position","Tissue","Subtype","H0","H1","H2","H3","H4","Number of SNPs")
st_11<-distinct(st_11)

# Supplementary table 12 --------------------------------------------------

st_12<-fread("2_Colocalization/Coloc_results_abf_druggable.csv")
st_12<-st_12 %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
st_12$Gene<-substr(st_12$Gene,1,15)
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(st_12$Gene), mart=ensembl)
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)


st_12<-merge(st_12,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
st_12<-st_12[rev(order(st_12$PP.H4.abf)),]
st_12$PP.H0.abf<-format(round(st_12$PP.H0.abf,digits=3),nsmall=3)
st_12$PP.H1.abf<-format(round(st_12$PP.H1.abf,digits=3),nsmall=3)
st_12$PP.H2.abf<-format(round(st_12$PP.H2.abf,digits=3),nsmall=3)
st_12$PP.H3.abf<-format(round(st_12$PP.H3.abf,digits=3),nsmall=3)
st_12$PP.H4.abf<-format(round(st_12$PP.H4.abf,digits=3),nsmall=3)
st_12<-distinct(st_12)
st_12<-dplyr::select(st_12,Gene,hgnc_symbol,Tissue,subtype,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps)
colnames(st_12)<-c("Ensembl gene ID","Gene","Tissue","Subtype","H0","H1","H2","H3","H4","Number of SNPs")
st_12<-distinct(st_12)

# Supplementary table 13 --------------------------------------------------

st_13<-readxl::read_excel("4_Drug target prioritization/CMAP/Results/Overall CMAP TWAS results.xlsx")
st_13$input<-"Overall CRC"
colnames(st_13)[which(colnames(st_13)=="pert_itome")]<-"pert_itime"
st_13<-distinct(st_13)


# Supplementary table 14 --------------------------------------------------

st_14a<-read.csv("spliceai.csv")
colnames(st_14a)<-c("Chromosome","Position","Splice event ID","Reference allele","Alternative allele","del1","del2","del3","Gene","DS_AG","DS_AL","DS_DL","DP_AG","DP_DG","Reference transcript","Strand","Cryptic_Acceptor_activation	Cryptic_Donor_activation	Any_splicing_aberration	bp_5prime	bp_3prime	Partial_intron_retention	Partial_exon_deletion	Partial_exon_start	Partial_exon_end	Partial_frameshift	Partial_intron_retention_aaseq	Partial_exon_deletion_aaseq	Gained_exon_size	Pseudoexon_activation	Pseudoexon_start	Pseudoexon_end	Pseudoexon_frameshift	Pseudoexon_intron	Pseudoexon_activation_aaseq	Exon_skipping	Lost_exons	Exon_skipping_frameshift	Exon_skipping_aaseq	Retained_intron_size	Intron_retention	Retained_intron	Intron_retention_frameshift	Intron_retention_aaseq")

# Supplementary table 15 --------------------------------------------------

st_14<-read.csv("CMap results.csv")
st_14<-distinct(st_14)
colnames(st_14)<-c("Gene","Colorectal cancer essentiality")



# Supplementary table 16 --------------------------------------------------

st_18a<-fread("6_Risk factors/Risk_factor_colocalization_results.csv")
st_18b<-fread("6_Risk factors/Risk_factor_druggable_colocalization_results.csv")
st_18<-rbind(st_18a,st_18b)

st_18<-st_18 %>% tidyr::separate(id,c("Gene","Tissue"),sep=";")
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(st_18$Gene), mart=ensembl)
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
genes<-distinct(genes)


st_18<-merge(st_18,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
st_18<-st_18[rev(order(st_18$PP.H4.abf)),]
st_18$PP.H0.abf<-format(round(st_18$PP.H0.abf,digits=3),nsmall=3)
st_18$PP.H1.abf<-format(round(st_18$PP.H1.abf,digits=3),nsmall=3)
st_18$PP.H2.abf<-format(round(st_18$PP.H2.abf,digits=3),nsmall=3)
st_18$PP.H3.abf<-format(round(st_18$PP.H3.abf,digits=3),nsmall=3)
st_18$PP.H4.abf<-format(round(st_18$PP.H4.abf,digits=3),nsmall=3)

st_18<-dplyr::select(st_18,risk_factor,Gene,hgnc_symbol,Tissue,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps)
colnames(st_18)<-c("Risk Factor","Ensembl gene ID","Gene","Tissue","H0","H1","H2","H3","H4","Number of SNPs")
st_18<-distinct(st_18)

# Supplementary table 17 --------------------------------------------------
st_19<-fread("3_MR/Exposure_data.csv")
st_19<-st_19 %>% tidyr::separate(Exposure,c("Gene","Tissue"),sep=";")
st_19<-st_19 %>% tidyr::separate(Gene,c("temp","Gene"),sep=",")
st_19$`Splice event`<-st_19$temp
st_19$`Splice event`[is.na(st_19$Gene)]<-NA
st_19$Gene[is.na(st_19$Gene)]<-st_19$temp
st_19<-dplyr::select(st_19,Gene,`Splice event`,Tissue,SNP,Effect_allele,Other_allele,EAF,Beta,SE,`P-value`)

st_19<-merge(st_19,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
st_19$hgnc_symbol[st_19$Gene=="ENSG00000213753"]<-"CENPBD2P"

colnames(st_19)[colnames(st_19)=="Gene"]<-"Ensembl gene ID"
colnames(st_19)[colnames(st_19)=="hgnc_symbol"]<-"Gene"

st_19<-dplyr::select(st_19,`Ensembl gene ID`,Gene,`Splice event`,Tissue,SNP,Effect_allele,Other_allele,EAF,Beta,SE,`P-value`)
st_19<-distinct(st_19)



# Labelling missing gene names --------------------------------------------

#missing<-unique(c(st_1$`Ensembl gene ID`[is.na(st_1$Gene)],st_2$`Ensembl gene ID`[is.na(st_2$Gene)],st_3$`Ensembl gene ID`[is.na(st_3$Gene)],st_5$gene_name[is.na(st_5$gene)],st_6$`Ensembl gene ID`[is.na(st_6$Gene)],st_8$`Ensembl gene ID`[is.na(st_8$Gene)],st_10$`Ensembl gene ID`[is.na(st_10$Gene)],st_11$`Ensembl gene ID`[is.na(st_11$Gene)],st_12$`Ensembl gene ID`[is.na(st_12$Gene)],st_14$`Ensembl gene ID`[is.na(st_14$Gene)],st_15$`Ensembl gene ID`[is.na(st_15$Gene)],st_16$`Ensembl gene ID`[is.na(st_16$Gene)],st_18$`Ensembl gene ID`[is.na(st_18$Gene)],st_19$`Ensembl gene ID`[is.na(st_19$Gene)]))

#write.table(missing,"missing.txt",row.names = FALSE,quote=FALSE)

missing<-fread("missing_with_gene_names.txt",header = FALSE)

missing_still<-missing$V1[missing$V2==""]

missing$V2[missing$V1=="ENSG00000277427"]<-"XXbac-BPG154L12.5"
missing$V2[missing$V1=="ENSG00000116957"]<-"TBCE"
missing$V2[missing$V1=="ENSG00000182109"]<-"RP11-69E11.4"
missing$V2[missing$V1=="ENSG00000228798"]<-"AP000473.5"
missing$V2[missing$V1=="ENSG00000242687"]<-"AC004893.11"
missing$V2[missing$V1=="ENSG00000255864"]<-"RP11-444D3.1"
missing$V2[missing$V1=="ENSG00000224441"]<-"AC068831.3"
missing$V2[missing$V1=="ENSG00000187695"]<-"RP11-723O4.6"
missing$V2[missing$V1=="ENSG00000240800"]<-"ATP8A2P1"
missing$V2[missing$V1=="ENSG00000272053"]<-"RP11-367G6.3"
missing$V2[missing$V1=="ENSG00000284413"]<-"BTBD8"
missing$V2[missing$V1=="ENSG00000240050"]<-"RP1-93H18.1"
missing$V2[missing$V1=="ENSG00000225302"]<-"RP11-539I5.1"
missing$V2[missing$V1=="ENSG00000231335"]<-"AC107072.2"
missing$V2[missing$V1=="ENSG00000228512"]<-"RP11-281A20.2"
missing$V2[missing$V1=="ENSG00000262222"]<-"RP11-876N24.4"
missing$V2[missing$V1=="ENSG00000256164"]<-"CCND2-AS1"

colnames(missing)<-colnames(genes)

still_missing_again<-genes[genes$ensembl_gene_id %in% setdiff(genes$ensembl_gene_id[genes$hgnc_symbol==""],missing$ensembl_gene_id),]

still_missing_again$hgnc_symbol[still_missing_again$ensembl_gene_id=="ENSG00000261888"]<-"Lnc-METRNL-1"
still_missing_again$hgnc_symbol[still_missing_again$ensembl_gene_id=="ENSG00000262003"]<-"Lnc-TIMM22-4"
still_missing_again$hgnc_symbol[still_missing_again$ensembl_gene_id=="ENSG00000272334"]<-"Lnc-EPM2AIP1-3"
still_missing_again$hgnc_symbol[still_missing_again$ensembl_gene_id=="ENSG00000273619"]<-"Lnc-CABLES2-7"
still_missing_again$hgnc_symbol[still_missing_again$ensembl_gene_id=="ENSG00000275437"]<-"Lnc-RBBP8NL-2"

missing<-rbind(missing,still_missing_again)

#still_missing_again_again<-unique(c(st_1$`Ensembl gene ID`[st_1$Gene==""],st_2$`Ensembl gene ID`[st_2$Gene==""],st_3$`Ensembl gene ID`[st_3$Gene==""],st_6$`Ensembl gene ID`[st_6$Gene==""],st_7$`Ensembl gene ID`[st_7$Gene==""],st_8$`Ensembl gene ID`[st_8$Gene==""],st_9$`Ensembl gene ID`[st_9$Gene==""],st_10$`Ensembl gene ID`[st_10$Gene==""],st_11$`Ensembl gene ID`[st_11$Gene==""],st_12$`Ensembl gene ID`[st_12$Gene==""],st_13$`Ensembl gene ID`[st_13$Gene==""],st_14$`Ensembl gene ID`[st_14$Gene==""],st_15$`Ensembl gene ID`[st_15$Gene==""],st_16$`Ensembl gene ID`[st_16$Gene==""],st_17$`Ensembl gene ID`[st_17$Gene==""],st_18$`Ensembl gene ID`[st_18$Gene==""],st_19$`Ensembl gene ID`[st_19$Gene==""]))

#write.table(still_missing_again_again,"still_missing.txt",row.names = FALSE,quote=FALSE)

missing_2<-fread("still_missing_with_gene_names.txt",header = FALSE)
colnames(missing_2)<-colnames(missing)
missing<-rbind(missing,missing_2)


write.table(missing,"missing_with_gene_names_all.txt",row.names = FALSE,quote=FALSE)




# Build supplementary table Excel file ------------------------------------

library(openxlsx)

#xl_lst <- list('Supplementary table 1'=st_1a,'Supplementary table 2' =tb_2, 'Supplementary table 3' = st_1, 'Supplementary table 4' = st_2, 'Supplementary table 5' = st_3, 'Supplementary table 7' = st_5, 'Supplementary table 8' = st_6, 'Supplementary table 9' = st_8, 'Supplementary table 10' = st_10, 'Supplementary table 11' = st_11, 'Supplementary table 12' = st_12, 'Supplementary table 13' = st_13,'Supplementary table 14' = st_14a, 'Supplementary table 15' = st_14, 'Supplementary table 16' = st_18, 'Supplementary table 17' = st_19)


xl_lst <- list('Supplementary table 1' = st_1, 'Supplementary table 2' = st_2, 'Supplementary table 3' = st_3, 'Supplementary table 4' = st_6, 'Supplementary table 5' = st_8, 'Supplementary table 6' = st_10, 'Supplementary table 7' = st_11, 'Supplementary table 8' = st_12, 'Supplementary table 9' = st_14a, 'Supplementary table 10' = st_14, 'Supplementary table 11' = st_13,'Supplementary table 12' = st_18, 'Supplementary table 13'=st_1a,'Supplementary table 14' = st_19)


write.xlsx(xl_lst, file = "Supplementary_tables.xlsx")




# Files for Zenodo --------------------------------------------------------


write.xlsx(st_4, file = "All_sQTLs_S-PrediXcan.xlsx")
write.xlsx(st_5, file = "All_eQTLs_S-PrediXcan.xlsx")



# New supplementary tables after review -----------------------------------

#splicing MR results
# Supplementary table 7 ---------------------------------------------------
# Supplementary table 7 ---------------------------------------------------
spl_1<-fread("data/CRC_TWAS/MR_results_sqtls.csv")
spl_1<-spl_1 %>% tidyr::separate(exposure,c("Splice event","Tissue"),sep=";")

genes<-fread("~/Downloads/temp.csv")
genes<-dplyr::select(genes,`Ensembl gene ID`,`Splice event`)

spl_1<-merge(spl_1,genes,by="Splice event",all.x=T)

spl_1$OR<-exp(spl_1$b)

spl_1<-dplyr::select(spl_1,`Ensembl gene ID`,`Splice event`,Tissue,subtype,OR,pval,nsnp)

spl_1$OR<-round(spl_1$OR,digits=2)
spl_1$pval<-format(spl_1$pval,digits=3)

spl_1<-spl_1[order(spl_1$pval),]

colnames(spl_1)<-c("Ensembl gene ID","Splice event","Tissue","Subtype","Odds ratio","P-value","nsnp")
spl_1<-distinct(spl_1)

spl_2<-fread("data/CRC_TWAS/f_stat_and_r2_sqtls.csv")

spl_2$`Splice event`<-spl_2$Gene

spl_2<-merge(spl_1,spl_2,by=c("Splice event","Tissue"),all.x=TRUE)
spl_2$r2<-round(spl_2$r2,digits=2)
spl_2$F_stat<-round(spl_2$F_stat,digits=2)
spl_2<-distinct(spl_2)
min(spl_2$F_stat)
max(spl_2$F_stat)
median(spl_2$F_stat)

spl_2<-dplyr::select(spl_2,`Ensembl gene ID`,`Splice event`,Tissue,Subtype,`Odds ratio`,`P-value`,nsnp,r2,F_stat)
spl_2<-distinct(spl_2)

fwrite(spl_2,"~/Downloads/results.csv")

#Work out figure 8
no_out<-fread("data/CRC_TWAS/no_outcome_druggable.csv",col.names = c("Gene","Tissue"))
no_harm<-fread("data/CRC_TWAS/no_harmonise_druggable.csv",col.names = c("Gene","Tissue"))

no_harm<-no_harm[!no_harm$Gene %in% no_out$Gene,]

mr<-fread("data/CRC_TWAS/MR_results_druggable_genome.csv")




