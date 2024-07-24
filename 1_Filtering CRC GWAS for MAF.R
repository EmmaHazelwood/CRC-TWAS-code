#module load R/4.2.0

library(data.table)
library(dplyr)

files<-list.files()

for (a in files){
paste(a)
gwas<-fread(a)
gwas<-gwas[gwas$Freq1>0.1,]
gwas<-gwas[gwas$Freq1<0.9,]
fwrite(gwas,paste("Formatted_GECCO/formatted_",a,sep=""))
}
