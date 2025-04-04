sink("coloc_abf_risk_factors_CRC.txt")
library(coloc)
library(dplyr)
library(gwasglue)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(plinkbinr)


# Get analyses needed -----------------------------------------------------
list<-fread("data/CRC_TWAS/Consistent_all_3.csv")

# Read in exposure data ---------------------------------------------------
bmi<-fread("data/CRC_risk_factors/ieu-b-40.csv")
whr<-fread("data/CRC_risk_factors/fat-distn.giant.ukbb.meta-analysis.whr.combined_formatted.txt")
ac<-fread("data/CRC_risk_factors/DRINKS_PER_WEEK_GWAS.txt")
tu<-fread("data/CRC_risk_factors/ieu-b-4877.csv")
bmi_f<-fread("data/CRC_risk_factors/ieu-a-974.csv")
whr_f<-fread("data/CRC_risk_factors/fat-distn.giant.ukbb.meta-analysis.whr.female_formatted.txt")
bmi_m<-fread("data/CRC_risk_factors/ieu-a-785.csv")
whr_m<-fread("data/CRC_risk_factors/fat-distn.giant.ukbb.meta-analysis.whr.male_formatted.txt")
