##############################################################################################
#### GTex samples annotation files
#### Contains a function to read files from a url: read.url
#### Contains Canonical Correlation Analysis CCA
#### Author: Charleen D. Adams

###############################################################################################
##### V8 subject and sample annotation files 
###############################################################################################
#https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
#https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx

setwd('C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/')
library(CCA)
library(dplyr)
library(tidyr)
library(EnsDb.Hsapiens.v79) 
library(biomaRt)
library(SNPlocs.Hsapiens.dbSNP.20120608)
library(data.table)
library(Cairo)
library(devtools)
library(MRInstruments) 
library(TwoSampleMR)
library(RadialMR)
library(data.table)
library(tidyverse)

read.url <- function(url, ...){
  tmpFile <- tempfile()
  download.file(url, destfile = tmpFile, method = "curl")
  url.data <- fread(tmpFile, ...)
  return(url.data)
}

samples_v8=read.url("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
write.csv(samples_v8, 'samples_v8.csv')
head(samples_v8$SAMPID, n=50)
head(samples_v8)
#breast=samples_v8[which(samples_v8$SMTSD=="Breast - Mammary Tissue"),]

#### The 'SMTSD' variable stands for 'Tissue Site Detail field' contains the counts of samples per tissue
head(samples_v8$SMTSD)
freq=table(samples_v8$SMTSD)
write.csv(freq, 'v8_tissue_frequencies.csv')

#### Annotations
subject_phenotypes=read.csv("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.csv")
write.csv(subject_phenotypes, 'subject_phenotypes.csv')
head(subject_phenotypes)

#### Separate SAMPID into subject, tissue, and aliquot

df <- data.frame(x = samples_v8$SAMPID)
samples_v83=df %>% separate(x, c("GTex","Donor_ID", "tissue_site_ID", "SM","aliquot_ID"))
dim(samples_v83)
head(samples_v83)
colnames(samples_v83)

#### cbind the separated id data with the samples_v8 data
samples_v8_final <- cbind(samples_v8, samples_v83)
colnames(samples_v8_final)
head(samples_v8_final)
table(samples_v8_final$Donor_ID)
table(samples_v8_final$tissue_site_ID)
donor_tissue=table(samples_v8_final$Donor_ID, samples_v8_final$SMTSD)
write.csv(donor_tissue, 'donor_tissue.csv')

#### Determine how many tissue types each donor gave
donor_tissue2=read.csv("donor_tissue2.csv")
str(donor_tissue2)
rownames(donor_tissue2)=donor_tissue2$Donor
donor_tissue2$Donor <- NULL
head(donor_tissue) 
cols=colnames(donor_tissue2)
donor_tissue2$new_col <- apply(donor_tissue2[,cols],1,sum)
head(donor_tissue2)
donor_tissue2 <- donor_tissue2[order(donor_tissue2$new_col),]
summary(donor_tissue2$new_col)

#### There is a potential problem with the correlation matrix of tissue by RP gene: 
#### The correlations could be confounded by some hidden structure related to the which 
#### tissues donor gave

#### To try to sort through this, I will a canonical correlation matrix approach, which is similar 
#### to PCA
#### Obtain a canonical correlation matrix

#### First I need the two matrices aligned in terms of their columns
#### The columns will be the tissue names
#### Since I already removed some tissues from the tissue by gene data (tissues without any RP eQTLs), 
#### I will remove those tissues for the the tissue by donor matrix

#### To do so, I need to drop the rows in samples_v8_final that correspond to unused tissues
#### Table samples_v8_final$SMTSD to recall how the tissues are named
tissues=table(samples_v8_final$SMTSD)
tissues_to_remove=as.data.frame(tissues)
write.csv(tissues_to_remove, 'tissues_to_remove.csv')

#### Open in Excel and identify which of the tissues to drop (compare with the Frequency data saved previously)

my_drop=c('Bladder', "Cells - Leukemia cell line (CML)", "Cervix - Ectocervix",
          "Cervix - Endocervix", "Fallopian Tube", "Kidney - Medulla")

samples_v8_final[!my_drop]
samples_v8_final_clean <- samples_v8_final[!which(samples_v8_final$SMTSD %in% my_drop),]
dim(samples_v8_final)
dim(samples_v8_final_clean)
samples_v8_final_clean_tab=table(samples_v8_final_clean$SMTSD)
write.csv(samples_v8_final_clean_tab, 'samples_v8_final_clean_tab.csv')
#### Compare again in Excel -- (#All good)

donor_tissue_clean=table(samples_v8_final_clean$SMTSD, samples_v8_final_clean$Donor_ID)
donor_tissue_clean_mat<- cor(t(donor_tissue_clean))
head(donor_tissue_clean_mat)



