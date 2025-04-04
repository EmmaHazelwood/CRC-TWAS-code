library(dplyr)
library(ggforestplot)
library(ggplot2)
library(data.table)
library(qqman)
library(tidyr)

firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Supplementary figure 1 --------------------------------------------------

fg_1<-fread("1_TWAS/JTI_All.csv")
a1 <-fg_1[fg_1$gene_name=="AAMP" & fg_1$tissue=="Adipose_Subcutaneous",]
a1$panel<-"A. AAMP & Adipose Subcutaneous"
b1 <-fg_1[fg_1$gene_name=="AAMP" & fg_1$tissue=="Adipose_Visceral_Omentum",]
b1$panel<-"B. AAMP & Adipose Visceral Omentum"
c1 <-fg_1[fg_1$gene_name=="AAMP" & fg_1$tissue=="Colon_Sigmoid",]
c1$panel<-"C. AAMP & Colon Sigmoid"
d1 <-fg_1[fg_1$gene_name=="COLCA1" & fg_1$tissue=="Colon_Transverse",]
d1$panel<-"D. COLCA1 & Colon Transverse"
e1 <-fg_1[fg_1$gene_name=="EPM2AIP1" & fg_1$tissue=="Adipose_Visceral_Omentum",]
e1$panel<-"E. EPMA2IP1 & Adipose Visceral Omentum"
f1 <-fg_1[fg_1$gene_name=="EPM2AIP1" & fg_1$tissue=="Whole_Blood",]
f1$panel<-"F. EPMA2IP1 & Whole Blood"
g1 <-fg_1[fg_1$gene_name=="LAMC1" & fg_1$tissue=="Whole_Blood",]
g1$panel<-"G. LAMC1 & Whole Blood"
h1 <-fg_1[fg_1$gene_name=="MLH1" & fg_1$tissue=="Adipose_Subcutaneous",]
h1$panel<-"H. MLH1 & Adipose Subcutaneous"
i1 <-fg_1[fg_1$gene_name=="MLH1" & fg_1$tissue=="Adipose_Visceral_Omentum",]
i1$panel<-"I. MLH1 & Adipose Visceral Omentum"
j1 <-fg_1[fg_1$gene_name=="MLH1" & fg_1$tissue=="Whole_Blood",]
j1$panel<-"J. MLH1 & Whole Blood"
k1 <-fg_1[fg_1$gene_name=="MLH1" & fg_1$tissue=="Cells_EBV-transformed_lymphocytes",]
k1$panel<-"K. MLH1 & lymphocytes"
l1 <-fg_1[fg_1$gene_name=="AC011816.2" & fg_1$tissue=="Adipose_Subcutaneous",]
l1$panel<-"L. RP11-129K12.1 & Adipose Subcutaneous"
m1 <-fg_1[fg_1$gene_name=="AC011816.2" & fg_1$tissue=="Colon_Transverse",]
m1$panel<-"M. RP11-129K12.1 & Colon Transverse"
n1 <-fg_1[fg_1$gene_name=="CCM2" & fg_1$tissue=="Whole_Blood",]
n1$panel<-"N. CCM2 & Whole Blood"

data<-rbind(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1)
data$se<-data$effect_size/data$zscore
data$subtype<-firstup(data$subtype)
data$subtype<-factor(data$subtype,levels=c("Overall","Colon","Distal","Proximal","Rectal","Female","Male"))
data<-data[order(data$subtype),]


p1<- ggforestplot::forestplot(
  df = data,
  name=subtype,
  estimate = effect_size,
  se=se,
  pvalue = pvalue,
  psignif = 6.01*10^(-8),
  xlab = "Effect size (95% CI)",
  logodds=FALSE
)+
  facet_wrap(~ panel, scales = "free_y", ncol = 4)
p1

#Output source data
fwrite(data,"Source data figure 6.txt")

ggplot2::ggsave(filename="Results/Supplementary figure 1.png", plot=ggplot2::last_plot(),width = 0.85*500, height = 0.85*400, units = "mm",bg="white",dpi=1000)


# Supplementary figure 2 --------------------------------------------------

fg_2b<-fread("1_TWAS/SMultiXcan_Expression_All.csv")
fg_2c<-fread("1_TWAS/SMultiXcan_Splicing_Expression_All.csv")
fg_2b$gene<-substr(fg_2b$gene,1,15)
fg_2b$gene_name[fg_2b$gene=="ENSG00000272334"]<-"AC011816.2"

a2 <-fg_2b[fg_2b$gene_name=="ABCC2",]
a2$panel<-"A. ABCC2 expression"
a2$analysis<-"expression"
b2 <-fg_2b[fg_2b$gene_name=="MLH1",]
b2$panel<-"B. MLH1 expression"
b2$analysis<-"expression"
c2 <-fg_2b[fg_2b$gene_name=="AC011816.2",]
c2$panel<-"C. RP11-129K12.1 expression"
c2$analysis<-"expression"
d2 <-fg_2c[fg_2c$gene=="ENSG00000163466",]
d2$panel<-"D. ARPC2 splicing"
d2$analysis<-"splicing"
d2<-d2[rev(order(d2$pvalue)),]
d2<-d2[d2$gene_name==d2$gene_name[1],]
e2 <-fg_2c[fg_2c$gene=="ENSG00000076650",]
e2<-e2[order(e2$pvalue),]
e2<-e2[e2$gene_name==e2$gene_name[1],]
e2$panel<-"E. GPATCH1 splicing"
e2$analysis<-"splicing"

data<-rbind(a2,b2,c2,d2,e2,fill=TRUE)
data$subtype<-factor(data$subtype,levels=c("Overall","Colon","Distal","Proximal","Rectal","Female","Male"))
data<-data[order(data$subtype),]
data$z_sd[is.na(data$z_sd)]<-0
data$pval_fill<-1
data$pval_fill[data$analysis=="expression" & data$pvalue<3.91*10^(-7)]<-0
data$pval_fill[data$analysis=="splicing" & data$pvalue<5.49*10^(-7)]<-0

p2<- ggforestplot::forestplot(
  df = data,
  name=subtype,
  estimate = z_mean,
  se=z_sd,
  pvalue = pval_fill,
  psignif = 0.05,
  xlab = "Mean Z-score (standard deviation)",
  logodds=FALSE
)+
  facet_wrap(~ panel, scales = "free_y", ncol = 3)
p2

#Output source data
fwrite(data,"Source data figure 7.txt")

ggplot2::ggsave(filename="Results/Supplementary figure 2.png", plot=ggplot2::last_plot(),width = 0.7*500, height = 0.5*400, units = "mm",bg="white",dpi=1000)

