sink("Rscript.txt")
library(data.table)
library(dplyr)
#Download from OpenGWAS
#wget https://objectstorage.us-ashburn-1.oraclecloud.com/n/idrvm4tkz2a8/b/OpenGWAS/o/ieu-b/ieu-b-40/ieu-b-40.vcf.gz
#wget https://objectstorage.us-ashburn-1.oraclecloud.com/n/idrvm4tkz2a8/b/OpenGWAS/o/ieu-b/ieu-b-4877/ieu-b-4877.vcf.gz
#wget https://objectstorage.us-ashburn-1.oraclecloud.com/n/idrvm4tkz2a8/b/OpenGWAS/o/ieu-a/ieu-a-974/ieu-a-974.vcf.gz
#wget https://objectstorage.us-ashburn-1.oraclecloud.com/n/idrvm4tkz2a8/b/OpenGWAS/o/ieu-a/ieu-a-785/ieu-a-785.vcf.gz
setwd("CRC_risk_factors")

com<-fread("fat-distn.giant.ukbb.meta-analysis.whr.combined.tbl")
com<-com %>% tidyr::separate(SNP, c("SNP", "effect_allele","other"),sep=":")
com<-select(com,SNP,Allele1,Allele2,`P-value`,Freq1,Effect)
fwrite(com,"fat-distn.giant.ukbb.meta-analysis.whr.combined_formatted.txt",sep="\t",na="NA",quote=FALSE)

fem<-fread("fat-distn.giant.ukbb.meta-analysis.whr.females.tbl")
fem<-fem %>% tidyr::separate(SNP, c("SNP", "effect_allele","other"),sep=":")
fem<-select(fem,SNP,Allele1,Allele2,`P-value`,Freq1,Effect)
fwrite(fem,"fat-distn.giant.ukbb.meta-analysis.whr.female_formatted.txt",sep="\t",na="NA",quote=FALSE)

mal<-fread("fat-distn.giant.ukbb.meta-analysis.whr.males.tbl")
mal<-mal %>% tidyr::separate(SNP, c("SNP", "effect_allele","other"),sep=":")
mal<-select(mal,SNP,Allele1,Allele2,`P-value`,Freq1,Effect)
fwrite(mal,"fat-distn.giant.ukbb.meta-analysis.whr.male_formatted.txt",sep="\t",na="NA",quote=FALSE)


crc<-fread("Formatted_GECCO/formatted_overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_overall_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)

crc<-fread("Formatted_GECCO/formatted_colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_colon_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)

crc<-fread("Formatted_GECCO/formatted_distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_distal_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)

crc<-fread("Formatted_GECCO/formatted_female_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_female_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)

crc<-fread("Formatted_GECCO/formatted_male_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_male_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)

crc<-fread("Formatted_GECCO/formatted_proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_proximal_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)

crc<-fread("Formatted_GECCO/formatted_rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
crc<-select(crc,SNP,Allele1,Allele2,P.value,Freq1,Effect)
fwrite(crc,"Formatted_GECCO/formatted_rectal_CRC_GWAS_noUKBio_summary_stats_annotated_space.txt",sep="\t",na="NA",quote=FALSE)


library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(genetics.binaRies)
library(data.table)

set_plink()
set_bcftools()

vcf <- readVcf("ieu-a-785.vcf")
gwas <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, type="exposure")
fwrite(gwas,"ieu-a-785.csv",sep="\t",na="NA",quote=FALSE)

vcf <- readVcf("ieu-b-40.vcf")
gwas <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, type="exposure")
fwrite(gwas,"ieu-b-40.csv",sep="\t",na="NA",quote=FALSE)

vcf <- readVcf("ieu-a-974.vcf")
gwas <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, type="exposure")
fwrite(gwas,"ieu-a-974.csv",sep="\t",na="NA",quote=FALSE)

vcf <- readVcf("ieu-b-4877.vcf")
gwas <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, type="exposure")
fwrite(gwas,"ieu-b-4877.csv",sep="\t",na="NA",quote=FALSE)



sink()