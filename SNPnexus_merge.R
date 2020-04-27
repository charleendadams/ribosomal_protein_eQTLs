###################################################################################################
#### Merge the SNPnexus annotation with the RP data for qval<0.05 and P<5x10-6 
#### Author: Charleen D. Adams
###################################################################################################
setwd('C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/')
#devtools::install_github('MRCIEU/TwoSampleMR')
library(sqldf)
library(CMplot)
library(oncofunco)
library(qqman)
library(dplyr)
library(tidyr)
library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(SNPlocs.Hsapiens.dbSNP.20120608)
library(data.table)
library(DataCombine)
library(Cairo)
library(devtools)
library(RColorBrewer)
library(MRInstruments)
library(TwoSampleMR)
library(RadialMR)
library(gplots)
library(tidyverse)
library(packcircles)
library(ggplot2)
library(viridis)
library(igraph)

load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/V8_master_tissues_n.Rdata")

#### Subset to the RP genes
RP_master=master_v8_n[which(master_v8_n$gene_name %in% myvars),]
RP_master$gene_name <- factor(RP_master$gene_name)
#write.csv(RP_master, 'RP_v8.csv')

#### Get the significant RPs by qval<0.05 and nominal P<5x10-6
RP_sig=RP_master[which(RP_master$qval<=0.05),]
stringent_RP_sig_beta=RP_sig[which(RP_sig$pval<0.000005),]

#### Get the list of unique eQTLs (set above is counts all tissues)
snp=unique(stringent_RP_sig_beta$SNP)
length(snp)
snp2=na.omit(snp)
length(snp2)
dim(stringent_RP_sig_beta)

#### The set of unique eQTLs was fed into SNPnexus to obtain annotation for MAF and genomic consequences
#### (coding/non-coding/UTR and deleteriousness, etc..)

#### Read in the MAF annotation
my_data <- read.table("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/gen_coords_v8_unique.txt", 
  header=TRUE, sep='\t',stringsAsFactors = FALSE)

dim(my_data)
985/996 #99%
str(my_data)
head(my_data)
my_data <- na.omit(my_data)

#### Create a dummy variable for whether the minor allele is the reference allele (1=yes)
my_data<-within(my_data,match <- ifelse (Minor.Allele == REF.Allele,1,0))
colnames(my_data)
table(my_data$match)
my_data$match=as.character(my_data$match)

#### Create the EAF variable by reasoning that if the minor allele is the references allele, then
#### 1) the EAF is 1-MAF for our data 
#### 2) otherwise, the EAF is the MAF

my_data$Minor.Allele.Global.Frequency=as.numeric(my_data$Minor.Allele.Global.Frequency)
my_data$EAF <- my_data$Minor.Allele.Global.Frequency
my_data$EAF <- ifelse(my_data$match=="1", 1-my_data$Minor.Allele.Global.Frequency, 
                  my_data$Minor.Allele.Global.Frequency)

summary(my_data$EAF)
hist(my_data$EAF)
write.csv(my_data, 'match.csv')

#### merge the MAF into the 'stringent_RP_sig_beta' file
my_data$dbSNP=as.character(my_data$dbSNP)

merged=sqldf("select * 
         from stringent_RP_sig_beta left join my_data 
         on stringent_RP_sig_beta.SNP = my_data.dbSNP", 
      row.names = TRUE)
head(merged)

merged <- merged[order(merged$SNP),]
hist(merged$EAF)
summary(merged$EAF)
merged$EAF=merged$EAF
merged_EAF=merged
head(merged_EAF$SNP)

#### One SNP has two rsids: change to the rsid that has MAF info: rs12012747
merged_EAF$SNP=as.character(merged_EAF$SNP)
merged_EAF$SNP=ifelse(merged_EAF$SNP=='rs4909,rs12012747', 'rs12012747', merged_EAF$SNP)
merged_EAF$match=ifelse(merged_EAF$SNP=='rs12012747', 0, merged_EAF$match)
merged_EAF$EAF=ifelse(merged_EAF$SNP=='rs12012747', 0.331126, merged_EAF$EAF)
table(merged_EAF$EAF)

### Rename the 'eaf' variable to reflect that it is the MAF in GTex (not the actual EAF)
colnames(merged_EAF)[31] <- "maf_in_GTex"
write.csv(merged_EAF, 'merged_master_RP_SNPnexus_EAF.csv')

#### Create a copy of my_data to use in the Markdown document to show what I did to my_data
#my_data2 <- read.table("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/gen_coords_v8_unique.txt", 
#                      header=TRUE, sep='\t',stringsAsFactors = FALSE)

save(myvars, myvars_all_mito, myvars_mito_39S, myvars_40S_60S, myvars_40S, myvars_60S, ao, master_v8_n, RP_master,
     stringent_RP_sig_beta, RP_sig, my_data, my_data2, merged_EAF,
     file="./V8_master_EAF.Rdata")
