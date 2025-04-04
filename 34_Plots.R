library(dplyr)
library(ggforestplot)
library(ggplot2)
library(data.table)
library(qqman)
library(biomaRt)
library(tidyr)



# Gene names --------------------------------------------------------------

tb_2<-fread("Consistent_all_3.csv")
hgnc_names<-distinct(dplyr::select(tb_2,gene,hgnc_symbol))


# Figures 3 and 4----------------------------------------------------------------

fg_1a<-fread("1_TWAS/JTI_All.csv")
fg_1b<-fread("1_TWAS/SMultiXcan_Expression_All.csv")
fg_1c<-fread("1_TWAS/SMultiXcan_Splicing_Expression_All.csv")

#Limit to ensembl ID and p-value and subtype
fg_1a<-dplyr::select(fg_1a,gene,pvalue,subtype)
fg_1b<-dplyr::select(fg_1b,gene,pvalue,subtype)
fg_1b$gene<-substr(fg_1b$gene,1,15)
fg_1c<-dplyr::select(fg_1c,gene,pvalue,subtype)
fg_1<-rbind(fg_1a,fg_1b,fg_1c)
all_rbind<-fg_1




fg_1$subtype<-tolower(fg_1$subtype)

#annotate with TSS (chromosome and position)
list<-c(unique(fg_1$gene))
options(biomart.update = FALSE)
options(httr.use.curl = FALSE)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","chromosome_name","transcription_start_site","strand","transcript_mane_select"), values=list, mart=ensembl)

fg_1<-merge(fg_1,genes,by.x="gene",by.y="ensembl_gene_id",all.x=T,allow.cartesian = T)

fg_1$transcription_start_site<-as.numeric(fg_1$transcription_start_site)
fg_1<-fg_1[!is.na(fg_1$transcription_start_site),]
fg_1$pvalue<-as.numeric(fg_1$pvalue)
fg_1<-fg_1[!is.na(fg_1$pvalue),]

#Group by subtype
overall<-fg_1[fg_1$subtype=="overall",]
colon<-fg_1[fg_1$subtype=="colon",]
distal<-fg_1[fg_1$subtype=="distal",]
female<-fg_1[fg_1$subtype=="female",]
male<-fg_1[fg_1$subtype=="male",]
proximal<-fg_1[fg_1$subtype=="proximal",]
rectal<-fg_1[fg_1$subtype=="rectal",]

fg_1<-fg_1[order(fg_1$pvalue),]
fg_1<-fg_1[match(unique(fg_1$gene), fg_1$gene),]

fg_1<-fg_1[rev(order(fg_1$transcript_mane_select)),]
fg_1<-fg_1[match(unique(fg_1$gene), fg_1$gene),]



#Get list of genes to annotate
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(tb_2$gene), mart=ensembl)
missing<-fread("missing_with_gene_names_all.txt")
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]
hgnc<-distinct(genes)
colnames(hgnc)<-c("gene","hgnc_symbol")
hgnca<-hgnc
list1_all<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Paper/Table 1 r output.csv")
list2_all<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Paper/Table 2 r output.csv")

list<-unique(c(list1_all$`Ensembl ID`,list2_all$`Ensembl ID`))
list<-c(list,c("ENSG00000093167","ENSG00000137310"))

hgnc<-rbind(hgnc,c("ENSG00000102554","KLF5"))
hgnc<-rbind(hgnc,c("ENSG00000099326","MZF1"))
hgnc<-rbind(hgnc,c("ENSG00000093167","LRRFIP2"))
hgnc$hgnc_symbol[hgnc$hgnc_symbol=="AC011816.2"]<-"RP11-129K12.1"
hgnc$hgnc_symbol[hgnc$hgnc_symbol=="AL121832.2"]<-"RPS21-DT"
hgnc<-distinct(hgnc)

fg_1<-merge(fg_1,hgnc,by="gene",all.x=T)
fg_1<-distinct(fg_1)



#Organise by p-value
overall<-overall[order(overall$pvalue),]
colon<-colon[order(colon$pvalue),]
distal<-distal[order(distal$pvalue),]
female<-female[order(female$pvalue),]
male<-male[order(male$pvalue),]
proximal<-proximal[order(proximal$pvalue),]
rectal<-rectal[order(rectal$pvalue),]

#Just keep one with lowest p-value
overall<-overall[match(unique(overall$gene), overall$gene),]
colon<-colon[match(unique(colon$gene), colon$gene),]
distal<-distal[match(unique(distal$gene), distal$gene),]
female<-female[match(unique(female$gene), female$gene),]
male<-male[match(unique(male$gene), male$gene),]
proximal<-proximal[match(unique(proximal$gene), proximal$gene),]
rectal<-rectal[match(unique(rectal$gene), rectal$gene),]



#Make manhattan plot
library(ggrepel)

# Prepare the dataset
don <- fg_1 %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(fg_1, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))

annotate<-don[don$is_annotate=="yes",]
names<-annotate[,8]
names<-names[order(names)]

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

options(ggrepel.max.overlaps = Inf)

# Make the plot
all_fg1<-ggplot(don, aes(x=transcription_start_sitecum, y=-log10(pvalue))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#DCDCDC", "#82D8F8"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="#D11D22", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,force=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylab("-log10 P-value")+
  scale_y_continuous(expand = c(0,0.5))
all_fg1

#Output source data
fwrite(don,"~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Paper/Source data figures 8 and 9.txt")


# Make for poster
all_fg1_poster<-ggplot(don, aes(x=transcription_start_sitecum, y=-log10(pvalue))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#4A4B9F", "#00BFD7"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="#B41E3B", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=3,force=50) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylab("-log10 P-value")+
  scale_y_continuous(expand = c(0,0.5))
all_fg1_poster

annotatelist<-unique(don$gene[don$is_annotate=="yes"])
annotatelistdf<-don[don$is_annotate=="yes",]
missing<-setdiff(list,annotatelist)

#overall
list<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Overall"],list2_all$`Ensembl ID`[list2_all$Subtype=="Overall"]))
list<-c(list,"ENSG00000137310")
hgnc<-hgnc[hgnc$gene %in% list,]

overall<-merge(overall,hgnc,by="gene",all.x=T)
overall<-distinct(overall)

# Prepare the dataset
don <- overall %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(overall, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Make the plot
overall_fg1<-ggplot(don, aes(x=transcription_start_sitecum, y=-log10(pvalue))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#DCDCDC", "#82D8F8"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="#D11D22", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  )+
  ylab("-log10 P-value")+
  scale_y_continuous(expand = c(0,0.9))
overall_fg1

#rectal and colon
list1<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Rectal"],list2_all$`Ensembl ID`[list2_all$Subtype=="Rectal"]))

hgnc1<-hgnca[hgnca$gene %in% list1,]


list2<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Colon"],list2_all$`Ensembl ID`[list2_all$Subtype=="Colon"]))
hgnc2<-hgnca[hgnca$gene %in% list2,]


hgnc<-rbind(hgnc1,hgnc2)
list<-c(list1,list2)

colon<-merge(colon,hgnc,by="gene",all.x=T)
colon<-distinct(colon)

rectal<-merge(rectal,hgnc,by="gene",all.x=T)
rectal<-distinct(rectal)

# Prepare the dataset
don_rectal <- rectal %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(rectal, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))%>%
  mutate( is_highlight1=ifelse(gene %in% list1, "yes", "no"))%>%
  mutate( is_highlight2=ifelse(gene %in% list2, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Prepare the dataset
don_colon <- colon %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(colon, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))%>%
  mutate( is_highlight1=ifelse(gene %in% list1, "yes", "no"))%>%
  mutate( is_highlight2=ifelse(gene %in% list2, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Make the plot
don_colon$log10p<--log10(don_colon$pvalue)
don_rectal$log10p<--1*(-log10(don_rectal$pvalue))

don_colon$gene_subtype<-paste(don_colon$gene,"_colon",sep="")
don_rectal$gene_subtype<-paste(don_rectal$gene,"_rectal",sep="")

don<-rbind(don_colon,don_rectal)

list1<-paste(list1,"_rectal",sep="")
list2<-paste(list2,"_colon",sep="")

don <- don %>%
  mutate( is_highlight1=ifelse(gene_subtype %in% list1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(gene_subtype %in% list2, "yes", "no")) %>%
  mutate( is_highlight3=ifelse(gene %in% list, "yes", "no")) 

colon_rectal_fg1<-ggplot(don, aes(x=transcription_start_sitecum, y=log10p)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#DCDCDC", "#82D8F8"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=0)+
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight3=="yes"), color="black", size=2) +
  geom_point(data=subset(don, is_highlight1=="yes"), color="#D11D22", size=2) +
  geom_point(data=subset(don, is_highlight2=="yes"), color="#D11D22", size=2) +
  
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes" & subtype=="colon"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,ylim=c(0,NA),force=13) +
  geom_label_repel( data=subset(don, is_annotate=="yes" & subtype=="rectal"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,ylim=c(NA,0),force=1.5) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  )+
  ylab("-log10 P-value")+
  scale_y_continuous(expand = c(0,0.5),breaks=c(-20,-10,0,10,20),labels=c("20","10","0","10","20")) 
colon_rectal_fg1

#proximal and distal
#proximal

list1<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Proximal"],list2_all$`Ensembl ID`[list2_all$Subtype=="Proximal"]))
list1<-c(list1,"ENSG00000093167")
hgnc1<-hgnca[hgnca$gene %in% list1,]
hgnc1<-rbind(hgnc1,c("ENSG00000093167","LRRFIP2"))


list2<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Distal"],list2_all$`Ensembl ID`[list2_all$Subtype=="Distal"]))
hgnc2<-hgnca[hgnca$gene %in% list2,]


hgnc<-rbind(hgnc1,hgnc2)
list<-c(list1,list2)

distal<-merge(distal,hgnc,by="gene",all.x=T)
distal<-distinct(distal)

proximal<-merge(proximal,hgnc,by="gene",all.x=T)
proximal<-distinct(proximal)

# Prepare the dataset
don_proximal <- proximal %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(proximal, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))%>%
  mutate( is_highlight1=ifelse(gene %in% list1, "yes", "no"))%>%
  mutate( is_highlight2=ifelse(gene %in% list2, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Prepare the dataset
don_distal <- distal %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(distal, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))%>%
  mutate( is_highlight1=ifelse(gene %in% list1, "yes", "no"))%>%
  mutate( is_highlight2=ifelse(gene %in% list2, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Make the plot
don_distal$log10p<--log10(don_distal$pvalue)
don_proximal$log10p<--1*(-log10(don_proximal$pvalue))

don_distal$gene_subtype<-paste(don_distal$gene,"_distal",sep="")
don_proximal$gene_subtype<-paste(don_proximal$gene,"_proximal",sep="")

don<-rbind(don_distal,don_proximal)

list1<-paste(list1,"_proximal",sep="")
list2<-paste(list2,"_distal",sep="")

don <- don %>%
  mutate( is_highlight1=ifelse(gene_subtype %in% list1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(gene_subtype %in% list2, "yes", "no")) %>%
  mutate( is_highlight3=ifelse(gene %in% list, "yes", "no")) 

distal_proximal_fg1<-ggplot(don, aes(x=transcription_start_sitecum, y=log10p)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#DCDCDC", "#82D8F8"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=0)+
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight3=="yes"), color="black", size=2) +
  geom_point(data=subset(don, is_highlight1=="yes"), color="#D11D22", size=2) +
  geom_point(data=subset(don, is_highlight2=="yes"), color="#D11D22", size=2) +
  
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes" & subtype=="distal"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,ylim=c(0,NA),force=1.5) +
  geom_label_repel( data=subset(don, is_annotate=="yes" & subtype=="proximal"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,ylim=c(NA,0),force=1.5) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  )+
  ylab("-log10 P-value")+
  scale_y_continuous(expand = c(0,0.5),breaks=c(10,0,-10,-20),labels=c("10","0","10","20")) 
distal_proximal_fg1

#male and female
#male

list1<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Male"],list2_all$`Ensembl ID`[list2_all$Subtype=="Male"]))
hgnc1<-hgnca[hgnca$gene %in% list1,]


list2<-unique(c(list1_all$`Ensembl ID`[list1_all$Subtype=="Female"],list2_all$`Ensembl ID`[list2_all$Subtype=="Female"]))
hgnc2<-hgnca[hgnca$gene %in% list2,]

hgnc<-rbind(hgnc1,hgnc2)
list<-c(list1,list2)

female<-merge(female,hgnc,by="gene",all.x=T)
female<-distinct(female)

male<-merge(male,hgnc,by="gene",all.x=T)
male<-distinct(male)

# Prepare the dataset
don_male <- male %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(male, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))%>%
  mutate( is_highlight1=ifelse(gene %in% list1, "yes", "no"))%>%
  mutate( is_highlight2=ifelse(gene %in% list2, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Prepare the dataset
don_female <- female %>% 
  
  # Compute chromosome size
  group_by(chromosome_name) %>% 
  summarise(chr_len=max(transcription_start_site)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(female, ., by=c("chromosome_name"="chromosome_name")) %>%
  
  # Add a cumulative position of each gene
  arrange(chromosome_name, transcription_start_site) %>%
  mutate( transcription_start_sitecum=transcription_start_site+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(gene %in% list, "yes", "no")) %>%
  mutate( is_annotate=ifelse(gene %in% list, "yes", "no"))%>%
  mutate( is_highlight1=ifelse(gene %in% list1, "yes", "no"))%>%
  mutate( is_highlight2=ifelse(gene %in% list2, "yes", "no"))

# Prepare X axis
axisdf <- don %>% group_by(chromosome_name) %>% summarize(center=( max(transcription_start_sitecum) + min(transcription_start_sitecum) ) / 2 )

# Make the plot
don_female$log10p<--log10(don_female$pvalue)
don_male$log10p<--1*(-log10(don_male$pvalue))

don_female$gene_subtype<-paste(don_female$gene,"_female",sep="")
don_male$gene_subtype<-paste(don_male$gene,"_male",sep="")

don<-rbind(don_female,don_male)

list1<-paste(list1,"_male",sep="")
list2<-paste(list2,"_female",sep="")

don <- don %>%
  mutate( is_highlight1=ifelse(gene_subtype %in% list1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(gene_subtype %in% list2, "yes", "no")) %>%
  mutate( is_highlight3=ifelse(gene %in% list, "yes", "no")) 

female_male_fg1<-ggplot(don, aes(x=transcription_start_sitecum, y=log10p)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#DCDCDC", "#82D8F8"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=0)+
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight3=="yes"), color="black", size=2) +
  geom_point(data=subset(don, is_highlight1=="yes"), color="#D11D22", size=2) +
  geom_point(data=subset(don, is_highlight2=="yes"), color="#D11D22", size=2) +
  
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes" & subtype=="female"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,ylim=c(0,NA),force=1.5) +
  geom_label_repel( data=subset(don, is_annotate=="yes" & subtype=="male"), aes(label=hgnc_symbol),min.segment.length = unit(0, 'lines'), size=2,ylim=c(NA,0),force=1.5) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  )+
  ylab("-log10 P-value")+
  scale_y_continuous(expand = c(0,0.5),breaks=c(10,0,-10,-20),labels=c("10","0","10","20")) 
female_male_fg1


library(ggpubr)
fg_1<-ggarrange(overall_fg1,NULL,colon_rectal_fg1,NULL,distal_proximal_fg1, NULL,female_male_fg1,labels=c("Overall","","Colon and rectal","","Distal and proximal","","Female and male"),font.label = list(size = 10),ncol=1,heights=c(1,0.05,1.1,0.05,1,0.05,1))
fg_1

ggplot2::ggsave(filename="Results/TWAS_Results_subtype.png", plot=ggplot2::last_plot(),width = 280, height = 520, units = "mm",bg="white",dpi=1000)

all_fg1
ggplot2::ggsave(filename="Results/TWAS_Results_all_combined.png", plot=ggplot2::last_plot(),width = 0.8*280, height = 0.8*230, units = "mm",bg="white",dpi=1000)

#Wide for presentations
ggplot2::ggsave(filename="Results/TWAS_Results_all_combined_wide.png", plot=ggplot2::last_plot(),width = 0.6*400, height = 0.6*0.8*230, units = "mm",bg="white",dpi=300,limitsize=FALSE)
#Wider for posters
all_fg1_poster
ggplot2::ggsave(filename="Results/TWAS_Results_all_combined_wider.png", plot=ggplot2::last_plot(),width = 0.7*400, height = 0.4*0.8*280, units = "mm",bg="white",dpi=300,limitsize=FALSE)

# Figure 5 ----------------------------------------------------------------
firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#MR results
fg_4<-fread("3_MR/MR_results.csv")
fg_4<- fg_4 %>% separate(exposure,c("Gene","Tissue","Subtype"),sep=";")

#Annotate with gene names
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=fg_4$Gene, mart=ensembl)
missing<-fread("missing_with_gene_names_all.txt")
genes<-rbind(genes,missing)
genes<-genes[!genes$hgnc_symbol=="",]

fg_4<-merge(fg_4,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
fg_4<-distinct(fg_4)

#fg_4$hgnc_symbol[fg_4$Gene=="ENSG00000261888"]<-"Novel transcript, Lnc-METRNL-1"
#fg_4$hgnc_symbol[fg_4$Gene=="ENSG00000262003"]<-"LOC101927727"
#fg_4$hgnc_symbol[fg_4$Gene=="ENSG00000272334"]<-"Novel transcript, Lnc-EPM2AIP1-3"
#fg_4$hgnc_symbol[fg_4$Gene=="ENSG00000272368"]<-"Novel Transcript, Antisense To CERS5"
#fg_4$hgnc_symbol[fg_4$Gene=="ENSG00000273619"]<-"Novel Transcript, Antisense To RPS21"
#fg_4$hgnc_symbol[fg_4$Gene=="ENSG00000275437"]<-"Novel Transcript, Sense Intronic To CABLES2"

fg_4$Subtype<-firstup(fg_4$Subtype)

fg_4$Tissue<-gsub("_"," ",fg_4$Tissue)

fg_4$Subtype<-factor(fg_4$Subtype,levels=c("Overall","Colon","Distal","Proximal","Rectal","Female","Male"))

fg_4<-fg_4[order(fg_4$b),]
fg_4$hgnc_symbol<-factor(fg_4$hgnc_symbol,levels=unique(fg_4$hgnc_symbol))

#Facet by subtype, row by gene, colour by tissue

p<- ggforestplot::forestplot(
  df = fg_4,
  name=hgnc_symbol,
  estimate = b,
  pvalue = pval,
  se=se,
  psignif = 3.924647e-05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  colour = Tissue
)+
  ggforce::facet_col(
    facets = ~Subtype,
    scales = "free_y",
    space = "free"
  ) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#82D8F8","#FFB800","#5B9366","#8CE19C")))
p

#Output source data
fwrite(fg_4,"~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Paper/Source data figure 3.txt")

ggplot2::ggsave(filename="Results/MR_results_figure_4.png", plot=ggplot2::last_plot(),width = 280, height = 400, units = "mm",bg="white",dpi=1000)




# Figure 6 ----------------------------------------------------------------

#Druggable MR results
fg_5<-fread("3_MR/MR_results_druggable_genome.csv")
fg_5<- fg_5 %>% separate(exposure,c("Gene","Tissue"),sep=";")
fg_5$Gene<-substr(fg_5$Gene,1,15)

#Annotate with gene names
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=fg_5$Gene, mart=ensembl)

fg_5<-merge(fg_5,genes,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
fg_5<-distinct(fg_5)
fg_5$hgnc_symbol[fg_5$hgnc_symbol==""]<-fg_5$Gene[fg_5$hgnc_symbol==""]

druggable_list<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Final scripts and results/Results/Druggable_list.csv")
fg_5<-fg_5[fg_5$Gene %in% druggable_list$ENSG00000117620,]


fg_5$Tissue<-gsub("_"," ",fg_5$Tissue)

fg_5$subtype<-factor(fg_5$subtype,levels=c("Overall","Colon","Distal","Proximal","Rectal","Female","Male"))

fg_5<-fg_5[fg_5$Gene %in%
             fg_5$Gene[(fg_5$pval<(0.05/380))],]

#Facet by subtype, row by gene, colour by tissue

p<- ggforestplot::forestplot(
  df = fg_5,
  name=hgnc_symbol,
  estimate = b,
  pvalue = pval,
  se=se,
  psignif = (0.05/380),
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  colour = Tissue
)+
  ggforce::facet_col(
    facets = ~subtype,
    scales = "free_y",
    space = "free"
  ) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#82D8F8","#5B9366","#8CE19C")))
p

#Output source data
fwrite(fg_5,"~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Paper/Source data figure 5.txt")

ggplot2::ggsave(filename="Results/Druggable_MR_results_figure_5.png", plot=ggplot2::last_plot(),width = 280, height = 260, units = "mm",bg="white",dpi=1000)





# Presentations -----------------------------------------------------------


res<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Final scripts and results/Results/Consistent_all_3.csv")
all_res<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Final scripts and results/Results/all_3_together.csv")

res$OR_MR<-exp(res$b)
res$LCI<-exp(res$b-1.96*res$se)
res$UCI<-exp(res$b+1.96*res$se)
res$CI_MR<-paste(res$LCI," to ",res$UCI,sep="")
res<-dplyr::select(res,gene,hgnc_symbol,t_i_best,subtype,z_mean,pvalue,analysis,PP.H4.abf,OR_MR,LCI,UCI,CI_MR,pval,b,se)
colnames(res)<-c("Ensembl","Gene","Tissue","Subtype","Z_TWAS","P_TWAS","TWAS","H4_coloc","OR_MR","CI_MR","LCI","UCI","P_MR","Beta_MR","SE_MR")

res$id<-paste(res$Ensembl,res$Tissue,res$Subtype,sep=";")
res<-res[match(unique(res$id), res$id),]

res$Analysis_type<-"MR"
res$Analysis_type[is.na(res$OR_MR) & res$TWAS!="splice"]<-"TWAS"
res$P_MR[is.na(res$OR_MR) & res$TWAS!="splice"]<-res$P_TWAS[is.na(res$OR_MR) & res$TWAS!="splice"]
res$SE_MR[is.na(res$OR_MR) & res$TWAS!="splice"]<-0
res$OR_MR[is.na(res$OR_MR) & res$TWAS!="splice"]<-res$Z_TWAS[is.na(res$OR_MR) & res$TWAS!="splice"]

all_res$OR_MR<-exp(all_res$b)
all_res$LCI<-exp(all_res$b-1.96*all_res$se)
all_res$UCI<-exp(all_res$b+1.96*all_res$se)
all_res$CI_MR<-paste(all_res$LCI," to ",all_res$UCI,sep="")
all_res<-dplyr::select(all_res,gene,hgnc_symbol,t_i_best,subtype,z_mean,pvalue,analysis,PP.H4.abf,OR_MR,LCI,UCI,CI_MR,pval,b,se)
colnames(all_res)<-c("Ensembl","Gene","Tissue","Subtype","Z_TWAS","P_TWAS","TWAS","H4_coloc","OR_MR","CI_MR","LCI","UCI","P_MR","Beta_MR","SE_MR")

all_res$id<-paste(all_res$Ensembl,all_res$Tissue,all_res$Subtype,sep=";")
all_res<-all_res[match(unique(all_res$id), all_res$id),]

all_res$Analysis_type<-"MR"
all_res$Analysis_type[is.na(all_res$OR_MR) & all_res$TWAS!="splice"]<-"TWAS"
all_res$P_MR[is.na(all_res$OR_MR) & all_res$TWAS!="splice"]<-all_res$P_TWAS[is.na(all_res$OR_MR) & all_res$TWAS!="splice"]
all_res$SE_MR[is.na(all_res$OR_MR) & all_res$TWAS!="splice"]<-0
all_res$OR_MR[is.na(all_res$OR_MR) & all_res$TWAS!="splice"]<-all_res$Z_TWAS[is.na(all_res$OR_MR) & all_res$TWAS!="splice"]

all_res$new_id<-paste(all_res$Ensembl,all_res$Tissue,sep=";")

overall<-res[res$Subtype=="overall",]
overall<-overall[order(overall$OR_MR),]

subtypes<-res[res$Subtype=="colon" | res$Subtype=="distal" |res$Subtype=="rectal" |res$Subtype=="proximal",]
subtypes<-subtypes[order(subtypes$OR_MR),]
subtypes$new_id<-paste(subtypes$Ensembl,subtypes$Tissue,sep=";")
list<-c(paste(subtypes$new_id,"colon",sep=";"),paste(subtypes$new_id,"distal",sep=";"),paste(subtypes$new_id,"proximal",sep=";"),paste(subtypes$new_id,"rectal",sep=";"))
subtypes<-all_res[all_res$id %in% list,]


sex<-res[res$Subtype=="male" | res$Subtype=="female",]
sex<-sex[order(sex$OR_MR),]
sex$new_id<-paste(sex$Ensembl,sex$Tissue,sep=";")
list<-c(paste(sex$new_id,"male",sep=";"),paste(sex$new_id,"female",sep=";"))
sex<-all_res[all_res$id %in% list,]

p<- ggforestplot::forestplot(
  df = res,
  name=Gene,
  estimate = OR_MR,
  pvalue = P_MR,
  se=SE_MR,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI) or Z mean",
  logodds=FALSE,
  colour = Tissue,
  shape = Analysis_type
)+
  ggforce::facet_col(
    facets = ~Subtype,
    scales = "free_y",
    space = "free"
  )
p


ggplot2::ggsave(filename="Results/All_results.png", plot=ggplot2::last_plot(),width = 280, height = 200, units = "mm",bg="white",dpi=1000)

p<- ggforestplot::forestplot(
  df = overall,
  name=Gene,
  estimate = OR_MR,
  pvalue = P_MR,
  se=SE_MR,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI) or Z mean",
  logodds=FALSE,
  colour = Tissue,
  shape = Analysis_type
)
p


ggplot2::ggsave(filename="Results/Overall_results.png", plot=ggplot2::last_plot(),width = 280, height = 200, units = "mm",bg="white",dpi=1000)

p<- ggforestplot::forestplot(
  df = subtypes,
  name=Gene,
  estimate = OR_MR,
  pvalue = P_MR,
  se=SE_MR,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI) or Z mean",
  logodds=FALSE,
  colour = Tissue,
  shape = Analysis_type
)+
  ggforce::facet_col(
    facets = ~Subtype,
    scales = "free_y",
    space = "free"
  )
p


ggplot2::ggsave(filename="Results/Subtypes_results.png", plot=ggplot2::last_plot(),width = 280, height = 200, units = "mm",bg="white",dpi=1000)

p<- ggforestplot::forestplot(
  df = sex,
  name=Gene,
  estimate = OR_MR,
  pvalue = P_MR,
  se=SE_MR,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI) or Z mean",
  logodds=FALSE,
  colour = Tissue,
  shape = Analysis_type
)+
  ggforce::facet_col(
    facets = ~Subtype,
    scales = "free_y",
    space = "free"
  )
p


ggplot2::ggsave(filename="Results/Sex_results.png", plot=ggplot2::last_plot(),width = 280, height = 200, units = "mm",bg="white",dpi=1000)


# Druggable ---------------------------------------------------------------
res<-fread("Results/Druggable_genes_strong.csv")
res$p_one_sided<-res$`P value`/2
res$z<-qnorm(res$p_one_sided)
res$se<-log(res$OR)/res$z

library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=unique(res$`Ensembl ID`), mart=ensembl)
res<-merge(res,genes,by.x="Ensembl ID",by.y="ensembl_gene_id",all.x=T)

p<- ggforestplot::forestplot(
  df = res,
  name=hgnc_symbol,
  estimate = OR,
  pvalue = `P value`,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=FALSE,
  colour = Tissue,
  se=se
)+
  ggforce::facet_col(
    facets = ~Subtype,
    scales = "free_y",
    space = "free"
  )
p
ggplot2::ggsave(filename="Results/Druggable_MR_results.png", plot=ggplot2::last_plot(),width = 280, height = 200, units = "mm",bg="white",dpi=1000)

#CCND2
res<-fread("~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Final scripts and results/Results/3_MR/MR_results_druggable_genome.csv")
ccnd2<-res[grep("ENSG00000118971",res$exposure),]

p<- ggforestplot::forestplot(
  df = ccnd2,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05/(length(unique(res$exposure)))*7,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  colour=subtype,
  se=se
)
p
ggplot2::ggsave(filename="Results/Druggable_MR_results_CCND2.png", plot=ggplot2::last_plot(),width = 280, height = 200, units = "mm",bg="white",dpi=1000)



# Splicing MR -------------------------------------------------------------

library(readxl)
library(ggplot2)
library(dplyr)

fg_4<-read_excel("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/Year 3/QIMR/Paper/Submitted/Review/Supplementary_tables.xlsx",sheet=6)
fg_4<-dplyr::select(fg_4,`Ensembl gene ID`,`Gene`,`Splice event`)
res<-data.table::fread("~/Downloads/MR_results_sqtls.csv")
res <- tidyr::separate(data = res, 
                       col = exposure,
                       into = c("Splice event", "Tissue"),  # Name your new columns
                       sep = ";")
res<-merge(res,fg_4,by="Splice event")

res$Gene_event<-paste(res$Gene,res$`Splice event`,sep="; ")

res$subtype<-firstup(res$subtype)
res$subtype<-factor(res$subtype,levels=c("Overall","Distal","Colon" ,"Proximal","Rectal","Male"))
res<-res[order(res$Gene_event),]
res$Gene_event<-factor(res$Gene_event,levels=unique(res$Gene_event))
res$Tissue<-gsub("_"," ",res$Tissue)

p<- ggforestplot::forestplot(
  df = res,
  name=`Gene_event`,
  estimate = b,
  pvalue = pval,
  se=se,
  psignif = (0.05/42),
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  colour = Tissue
)+
  ggforce::facet_col(
    facets = ~subtype,
    scales = "free_y",
    space = "free"
  ) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#82D8F8","#FFB800","#5B9366","#8CE19C")))
p

#Output source data
fwrite(res,"~/OneDrive - University of Bristol/Documents/Year 3/QIMR/Paper/Source data figure 4.txt")

ggplot2::ggsave(filename="Results/sQTL_MR_results_figure.png", plot=ggplot2::last_plot(),width = 0.8*280, height = 0.8*310, units = "mm",bg="white",dpi=1000)
