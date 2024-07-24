#module load R/4.2.0

library(data.table)
library(dplyr)
library(tidyr)

# SMultiXcan results ------------------------------------------------------

all<-list.files()
all<-all[!grepl("splicing",all)]

colon<-all[grepl("colon",all)]
distal<-all[grepl("distal",all)]
female<-all[grepl("female",all)]
male<-all[grepl("male",all)]
male<-male[!male %in% female]
overall<-all[grepl("overall",all)]
proximal<-all[grepl("proximal",all)]
rectal<-all[grepl("rectal",all)]

colon_res<-fread(colon)
colon_res$subtype<-"Colon"
distal_res<-fread(distal)
distal_res$subtype<-"Distal"
female_res<-fread(female)
female_res$subtype<-"Female"
male_res<-fread(male)
male_res$subtype<-"Male"
overall_res<-fread(overall)
overall_res$subtype<-"Overall"
proximal_res<-fread(proximal)
proximal_res$subtype<-"Proximal"
rectal_res<-fread(rectal)
rectal_res$subtype<-"Rectal"
all_res<-rbind(colon_res,distal_res,female_res,male_res,overall_res,proximal_res,rectal_res)

strong_p<-0.05/(7*length(unique(all_res$gene)))
weak_p<-0.05/length(unique(all_res$gene))

strong_all<-all_res[all_res$pvalue<strong_p,]
weak_all<-all_res[all_res$pvalue<weak_p,]

strong_all<-strong_all[order(strong_all,pvalue),]
weak_all<-weak_all[order(weak_all,pvalue),]

strong_genes<-unique(strong_all$gene_name)
weak_genes<-unique(weak_all$gene_name)

fwrite(strong_all,"Results/SMultiXcan_Expression_All_Strong.csv")
fwrite(weak_all,"Results/SMultiXcan_Expression_All_Weak.csv")
fwrite(all_res,"Results/SMultiXcan_Expression_All.csv")


# JTI ---------------------------------------------------------------------
all<-list.files()

all_res<-data.frame()
for (a in all){
  x<-fread(a)
  x$file<-paste(a)
  all_res<-rbind(all_res,x)
}

all_res$tissue<-NA
all_res$tissue[grepl("Adipose_Subcutaneous",all_res$file)]<-"Adipose_Subcutaneous"
all_res$tissue[grepl("Adipose_Visceral_Omentum",all_res$file)]<-"Adipose_Visceral_Omentum"
all_res$tissue[grepl("Cells_EBV-transformed_lymphocytes",all_res$file)]<-"Cells_EBV-transformed_lymphocytes"
all_res$tissue[grepl("Colon_Sigmoid",all_res$file)]<-"Colon_Sigmoid"
all_res$tissue[grepl("Colon_Transverse",all_res$file)]<-"Colon_Transverse"
all_res$tissue[grepl("Whole_Blood",all_res$file)]<-"Whole_Blood"

all_res$subtype<-NA
all_res$subtype[grepl("colon",all_res$file)]<-"colon"
all_res$subtype[grepl("distal",all_res$file)]<-"distal"
all_res$subtype[grepl("male",all_res$file)]<-"male"
all_res$subtype[grepl("female",all_res$file)]<-"female"
all_res$subtype[grepl("overall",all_res$file)]<-"overall"
all_res$subtype[grepl("proximal",all_res$file)]<-"proximal"
all_res$subtype[grepl("rectal",all_res$file)]<-"rectal"

all_res<-all_res[!all_res$subtype=="",]

strong_p<-0.05/(7*6*length(unique(all_res$gene)))
weak_p<-0.05/length(unique(all_res$gene))

strong_all<-all_res[all_res$pvalue<strong_p,]
weak_all<-all_res[all_res$pvalue<weak_p,]

strong_all<-strong_all[order(strong_all,pvalue),]
weak_all<-weak_all[order(weak_all,pvalue),]

strong_genes<-unique(strong_all$gene_name)
weak_genes<-unique(weak_all$gene_name)

fwrite(strong_all,"Results/JTI_All_Strong.csv")
fwrite(weak_all,"Results/JTI_All_Weak.csv")
fwrite(all_res,"Results/JTI_All.csv")


# Splicing ----------------------------------------------------------------

all<-list.files()
all<-all[grepl("splicing",all)]

colon<-all[grepl("colon",all)]
distal<-all[grepl("distal",all)]
female<-all[grepl("female",all)]
male<-all[grepl("male",all)]
male<-male[!male %in% female]
overall<-all[grepl("overall",all)]
proximal<-all[grepl("proximal",all)]
rectal<-all[grepl("rectal",all)]

colon_res<-fread(colon)
colon_res$subtype<-"Colon"
distal_res<-fread(distal)
distal_res$subtype<-"Distal"
female_res<-fread(female)
female_res$subtype<-"Female"
male_res<-fread(male)
male_res$subtype<-"Male"
overall_res<-fread(overall)
overall_res$subtype<-"Overall"
proximal_res<-fread(proximal)
proximal_res$subtype<-"Proximal"
rectal_res<-fread(rectal)
rectal_res$subtype<-"Rectal"
all_res<-rbind(colon_res,distal_res,female_res,male_res,overall_res,proximal_res,rectal_res)


splice_strong<-all_res[all_res$pvalue<(0.05/(7*length(unique(all_res$gene)))),]
  
  
#Annotate with ENSG
genes<-all_res %>% tidyr::separate(gene, c("indel", "chr", "start", "end"))
genes<-genes[!duplicated(genes),]
genes$chr<-paste("chr",genes$chr,sep="")

as<-fread("GTEx/GTEx_Analysis_v8_sQTL_groups/Adipose_Subcutaneous.leafcutter.phenotype_groups.txt",header = F)
av<-fread("GTEx/GTEx_Analysis_v8_sQTL_groups/Adipose_Visceral_Omentum.leafcutter.phenotype_groups.txt",header = F)
ly<-fread("GTEx/GTEx_Analysis_v8_sQTL_groups/Cells_Cultured_fibroblasts.leafcutter.phenotype_groups.txt",header = F)
cs<-fread("GTEx/GTEx_Analysis_v8_sQTL_groups/Colon_Sigmoid.leafcutter.phenotype_groups.txt",header = F)
ct<-fread("GTEx/GTEx_Analysis_v8_sQTL_groups/Colon_Transverse.leafcutter.phenotype_groups.txt",header = F)
wb<-fread("GTEx/GTEx_Analysis_v8_sQTL_groups/Whole_Blood.leafcutter.phenotype_groups.txt",header = F)
all<-rbind(as,av,ly,cs,ct,wb)
all<-all %>% tidyr::separate(V1, c("chr", "start", "end","id","gene"),sep=":")
all$gene<-substr(all$gene,1,15)
all<-select(all,chr,start,end,gene)
all<-all[!duplicated(all),]

all_res<-merge(genes,all,all.x=T)
all_res<-all_res[!is.na(all_res$pvalue),]

strong_p<-0.05/(7*length(unique(all_res$gene)))
weak_p<-0.05/length(unique(all_res$gene))
all_res$pvalue<-as.numeric(all_res$pvalue)

strong_all<-all_res[all_res$pvalue<strong_p,]
weak_all<-all_res[all_res$pvalue<weak_p,]

#Annotate the ones which are missing
strong_all$gene[strong_all$gene_name=="intron_11_111356433_111357445"]<-"ENSG00000110777"
strong_all$gene[strong_all$gene_name=="intron_11_61802462_61802801"]<-"ENSG00000149485"
strong_all$gene[strong_all$gene_name=="intron_12_50130694_50130863"]<-"ENSG00000139624"
strong_all$gene[strong_all$gene_name=="intron_9_89418118_89424476"]<-"ENSG00000187764"


fwrite(strong_all,"Results/SMultiXcan_Splicing_Expression_All_Strong.csv")
fwrite(weak_all,"Results/SMultiXcan_Splicing_Expression_All_Weak.csv")
fwrite(all_res,"Results/SMultiXcan_Splicing_Expression_All.csv")

# Crossovers --------------------------------------------------------------

jti<-fread("Results/JTI_All_Strong.csv")
exp<-fread("Results/SMultiXcan_Expression_All_Strong.csv")
spl<-fread("Results/SMultiXcan_Splicing_Expression_All_Strong.csv")
jti<-fread("Results/1_TWAS/JTI_All_Strong.csv")
exp<-fread("Results/1_TWAS/SMultiXcan_Expression_All_Strong.csv")
spl<-fread("Results/1_TWAS/SMultiXcan_Splicing_Expression_All_Strong.csv")
exp$gene<-substr(exp$gene,1,15)
spl<-spl[!spl$gene=="",]

tot<-c(unique(exp$gene),unique(spl$gene),unique(jti$gene))
length(unique(tot))

length(unique(jti$gene))
length(unique(exp$gene))
length(unique(spl$gene))
spl$id<-paste(spl$gene,spl$gene_name,sep="_")
length(unique(spl$id))

all3<-intersect(intersect(unique(exp$gene),unique(spl$gene)),unique(jti$gene))
length(unique(all3))

expspl<-intersect(unique(exp$gene),unique(spl$gene))
expspl<-setdiff(unique(expspl),unique(jti$gene))
length(unique(expspl))

expjti<-intersect(unique(exp$gene),unique(jti$gene))
expjti<-setdiff(unique(expjti),unique(spl$gene))
length(unique(expjti))

jtispl<-intersect(unique(spl$gene),unique(jti$gene))
jtispl<-setdiff(unique(jtispl),unique(exp$gene))
length(unique(jtispl))

length(unique(c(expspl,expjti,jtispl)))

length(unique(c(jti$gene,exp$gene,spl$gene)))
length(unique(c(jti$gene,exp$gene)))

jtionly<-setdiff(jti$gene,c(exp$gene,spl$gene))
length(unique(jtionly))

exponly<-setdiff(exp$gene,c(jti$gene,spl$gene))
length(unique(exponly))

splonly<-setdiff(spl$gene,c(exp$gene,jti$gene))
length(unique(splonly))

length(unique(c(jtionly,exponly,splonly)))

both_genes<-intersect(unique(exp$gene_name),unique(jti$gene_name))
all_genes<-unique(c(jti$gene_name,exp$gene_name))


twas_exp_multi_up<-exp[exp$gene %in% expjti & exp$z_mean>0]
twas_exp_multi_down<-exp[exp$gene %in% expjti & exp$z_mean<0]
twas_exp_jti_up<-jti[jti$gene %in% expjti & jti$effect_size>0]
twas_exp_jti_down<-jti[jti$gene %in% expjti & jti$effect_size<0]

no1<-unique(twas_exp_multi_up$gene %in% twas_exp_jti_down$gene)
no2<-unique(twas_exp_jti_down$gene %in% twas_exp_multi_up$gene)
no3<-unique(twas_exp_multi_down$gene %in% twas_exp_jti_up$gene)
no4<-unique(twas_exp_jti_up$gene %in% twas_exp_multi_down$gene)

