library(data.table)
library(dplyr)

overall<-fread("Genetic correlation overall.csv")
colon<-fread("Genetic correlation colon.csv")
distal<-fread("Genetic correlation distal.csv")
female<-fread("Genetic correlation female.csv")
male<-fread("Genetic correlation male.csv")
proximal<-fread("Genetic correlation proximal.csv")
rectal<-fread("Genetic correlation rectal.csv")

all<-rbind(overall,colon,distal,female,male,proximal,rectal)
all_weak<-all[p<0.05,]
all_strong<-all[p<(0.05/4),]
