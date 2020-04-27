##############################################################################
#### Example of the pre-formatting done for each tissue for GTex v8
#### GTEX new release - Breast_Mammary_Tissue
#### Helpful code: gsub, paste, and use of EnsDb.Hsapiens.v79 to fetch annotations
##############################################################################
setwd('C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/breast')
#install.packages(c("curl", "httr"))
devtools::install_github('MRCIEU/TwoSampleMR')

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

### Load RP eQTLs
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/ribosomal_protein_lists.Rdata")
#ao <- available_outcomes()

Breast_Mammary_Tissue <- read.table('C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/breast/Breast_Mammary_Tissue.v8.egenes.txt', sep="\t", header=TRUE)
colnames(Breast_Mammary_Tissue)[1] <- "ensemblsIDS"
ensemblsIDS=Breast_Mammary_Tissue$ensemblsIDS
head(Breast_Mammary_Tissue)

#### Remove the version numbers at the end of the ensembl names
Breast_Mammary_Tissue$ensembl_gene_id <- gsub('\\..+$', '', Breast_Mammary_Tissue$ensemblsIDS)
head(Breast_Mammary_Tissue$ensembl_gene_id)
head(Breast_Mammary_Tissue)
colnames(Breast_Mammary_Tissue)

###########################################################################

#### Fetch the gene names: convert from ensembl.gene to gene.symbol
ensembl.genes=as.character(Breast_Mammary_Tissue$ensembl_gene_id)
head(ensembl.genes)
str(ensembl.genes)

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, 
  keytype = "GENEID", columns = c("SYMBOL","GENEID"))
head(geneIDs1)
colnames(Breast_Mammary_Tissue)[34]="GENEID"

#### Merge the gene names with the Breast_Mammary_Tissue tissue data 
Breast_Mammary_Tissue <- merge(x = geneIDs1, 
  y = Breast_Mammary_Tissue, 
  by.x="GENEID",
  by.y="GENEID")
head(Breast_Mammary_Tissue)

colnames(Breast_Mammary_Tissue)[18] <- "other_allele"
colnames(Breast_Mammary_Tissue)[19] <- "effect_allele"
colnames(Breast_Mammary_Tissue)[21] <- "SNP"
colnames(Breast_Mammary_Tissue)[24] <- "eaf"
colnames(Breast_Mammary_Tissue)[27] <- "beta"
colnames(Breast_Mammary_Tissue)[28] <- "se"
colnames(Breast_Mammary_Tissue)[26] <- "pval" #check

#### Separate out the variant ID to obtain the position to use for getting rsids
#### also gets the alleles

df <- data.frame(x = Breast_Mammary_Tissue$variant_id)
Breast_Mammary_Tissue3=df %>% separate(x, c("chromsome", "position", "other_allele", "effect_allele", "version"))
dim(Breast_Mammary_Tissue3)
head(Breast_Mammary_Tissue3)
Breast_Mammary_Tissue$position=Breast_Mammary_Tissue3$position
Breast_Mammary_Tissue3$start=Breast_Mammary_Tissue3$position
Breast_Mammary_Tissue3$stop=Breast_Mammary_Tissue3$position

#### merge the separated id data with the Breast_Mammary_Tissue tissue data
Breast_Mammary_Tissue <- merge(x = Breast_Mammary_Tissue3, 
  y = Breast_Mammary_Tissue, 
  by.x="position",
  by.y="position")
head(Breast_Mammary_Tissue)

#### Paste chromosome number and colon onto position
Breast_Mammary_Tissue$snp_position=(paste("1:", Breast_Mammary_Tissue$position, sep=""))
head(Breast_Mammary_Tissue$snp_position)

#### Rename other variables for MR-Base
colnames(Breast_Mammary_Tissue)[3] <- "other_allele"
colnames(Breast_Mammary_Tissue)[4] <- "effect_allele"
#colnames(Breast_Mammary_Tissue)[9] <- "gene_name"

head(Breast_Mammary_Tissue$variant_id)
myvars

save(myvars, myvars_all_mito, myvars_mito_39S, myvars_40S_60S, myvars_40S, myvars_60S, Breast_Mammary_Tissue,
     file="./V8_master_Breast_Mammary_Tissue_RP_lists.Rdata")