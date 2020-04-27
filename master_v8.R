#######################################################################################
#### Merge the v8 GTEx files
#### Author: Charleen D. Adams
#######################################################################################
setwd('C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8')
#devtools::install_github('MRCIEU/TwoSampleMR')

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

ao <- available_outcomes()
#### Adipose Subcutaneous
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/adipose_sub/V8_master_adi_subcut_RP_lists.Rdata")
adi_subcut
adi_subcut$tissue="Adipose_Subcutaneous"
adi_subcut$n=763
colnames(adi_subcut)
#### Adipose Visceral Omentum
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/adipose_vo/V8_master_adipose_vo_RP_lists.Rdata")
adipose_vo
adipose_vo$tissue="Adipose_Visceral_Omentum"
adipose_vo$n=564
#### Adrenal
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/adrenal/V8_master_Adrenal_Gland_RP_lists.Rdata")
Adrenal_Gland
Adrenal_Gland$tissue="Adrenal_Gland"
Adrenal_Gland$n=275
#### Amygdala
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/amygdala/V8_master_amygdala_RP_lists.Rdata")
amygdala
amygdala$tissue="Amygdala"
amygdala$n=177
#### Artery -- aorta
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/artery_aorta/V8_master_artery_aorta_RP_lists.Rdata")
artery_aorta
artery_aorta$tissue="Artery_Aorta"
artery_aorta$n=450
#### Artery -- tibial
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/Artery_Tibial/V8_master_Artery_Tibial_RP_lists.Rdata")
Artery_Tibial
Artery_Tibial$tissue='Artery_Tibial'
Artery_Tibial$n=770
#### Artery_Coronary
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/artery_coronary/V8_master_Artery_Coronary_RP_lists.Rdata")
Artery_Coronary
Artery_Coronary$tissue='Artery_Coronary'
Artery_Coronary$n=253
#### Brain_Anterior_cingulate_cortex_BA24
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/brain_anterior_cin_cortex/V8_master_Brain_Anterior_cingulate_cortex_BA24_RP_lists.Rdata")
Brain_Anterior_cingulate_cortex_BA24
Brain_Anterior_cingulate_cortex_BA24$tissue='Brain_Anterior_cingulate_cortex_BA24'
Brain_Anterior_cingulate_cortex_BA24$n=213
#### Brain_Caudate_basal_ganglia
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/caudate_BG/V8_master_Brain_Caudate_basal_ganglia_RP_lists.Rdata")
Brain_Caudate_basal_ganglia
Brain_Caudate_basal_ganglia$tissue='Brain_Caudate_basal_ganglia'
Brain_Caudate_basal_ganglia$n=291
#### Brain_Cerebellar_Hemisphere
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/cerebellar_hem/V8_master_Brain_Cerebellar_Hemisphere_RP_lists.Rdata")
Brain_Cerebellar_Hemisphere
Brain_Cerebellar_Hemisphere$tissue='Brain_Cerebellar_Hemisphere'
Brain_Cerebellar_Hemisphere$n=263
#### Brain_Cerebellum
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/cerebellum/V8_master_Brain_Cerebellum_RP_lists.Rdata")
Brain_Cerebellum
Brain_Cerebellum$tissue='Brain_Cerebellum'
Brain_Cerebellum$n=298
#### Brain_Cortex
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/cortex/V8_master_Brain_Cortex_RP_lists.Rdata")
Brain_Cortex
Brain_Cortex$tissue='Brain_Cortex'
Brain_Cortex$n=325
#### Brain_Frontal_Cortex_BA9
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/prefrontal/V8_master_Brain_Frontal_Cortex_BA9_RP_lists.Rdata")
Brain_Frontal_Cortex_BA9
Brain_Frontal_Cortex_BA9$tissue='Brain_Frontal_Cortex_BA9'
Brain_Frontal_Cortex_BA9$n=425
#### Brain_Hippocampus
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/hippocampus/V8_master_Brain_Hippocampus_RP_lists.Rdata")
Brain_Hippocampus
Brain_Hippocampus$tissue='Brain_Hippocampus'
Brain_Hippocampus$n=243
#### Brain_Hypothalamus
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/hypothalamus/V8_master_Brain_Hypothalamus_RP_lists.Rdata")
Brain_Hypothalamus
Brain_Hypothalamus$tissue='Brain_Hypothalamus'
Brain_Hypothalamus$n=236
#### Brain_Nucleus_accumbens_basal_ganglia
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/nucleus_accum_BG/V8_master_Brain_Nucleus_accumbens_basal_ganglia_RP_lists.Rdata")
Brain_Nucleus_accumbens_basal_ganglia
Brain_Nucleus_accumbens_basal_ganglia$tissue='Brain_Nucleus_accumbens_basal_ganglia'
Brain_Nucleus_accumbens_basal_ganglia$n=277
#### Brain_Putamen_basal_ganglia
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/putamen_BG/V8_master_Brain_Putamen_basal_ganglia_RP_lists.Rdata")
Brain_Putamen_basal_ganglia
Brain_Putamen_basal_ganglia$tissue='Brain_Putamen_basal_ganglia'
Brain_Putamen_basal_ganglia$n=232
#### Brain_Spinal_cord_cervical
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/spinal_cord_cervical/V8_master_Brain_Spinal_cord_cervical_RP_lists.Rdata")
Brain_Spinal_cord_cervical
Brain_Spinal_cord_cervical$tissue='Brain_Spinal_cord_cervical'
Brain_Spinal_cord_cervical$n=182
#### Brain_Substantia_nigra
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/substantia_nigra/V8_master_Brain_Substantia_nigra_RP_lists.Rdata")
Brain_Substantia_nigra
Brain_Substantia_nigra$tissue='Brain_Substantia_nigra'
Brain_Substantia_nigra$n=164
#### Breast_Mammary_Tissue
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/breast/V8_master_Breast_Mammary_Tissue_RP_lists.Rdata")
Breast_Mammary_Tissue
Breast_Mammary_Tissue$tissue='Breast_Mammary_Tissue'
Breast_Mammary_Tissue$n=480
#### Cells_Cultured_fibroblasts
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/cells_fibroblasts/V8_master_Cells_Cultured_fibroblasts_RP_lists.Rdata")
Cells_Cultured_fibroblasts
Cells_Cultured_fibroblasts$tissue='Cells_Cultured_fibroblasts'
Cells_Cultured_fibroblasts$n=527
#### Cells_EBV_transformed_lymphocytes
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/cells_trans_lymph/V8_master_Cells_EBV_transformed_lymphocytes_RP_lists.Rdata")
Cells_EBV_transformed_lymphocytes
Cells_EBV_transformed_lymphocytes$tissue='Cells_EBV_transformed_lymphocytes'
Cells_EBV_transformed_lymphocytes$n=192
#### Colon_Sigmoid
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/colon_sigmoid/V8_master_Colon_Sigmoid_RP_lists.Rdata")
Colon_Sigmoid
Colon_Sigmoid$tissue='Colon_Sigmoid'
Colon_Sigmoid$n=389
#### Colon_Transverse
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/colon_transverse/V8_master_Colon_Transverse_RP_lists.Rdata")
Colon_Transverse
Colon_Transverse$tissue='Colon_Transverse'
Colon_Transverse$n=432
#### Esophagus_Gastroesophageal_Junction
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/esophagus_gastr_junction/V8_master_Esophagus_Gastroesophageal_Junction_RP_lists.Rdata")
Esophagus_Gastroesophageal_Junction
Esophagus_Gastroesophageal_Junction$tissue='Esophagus_Gastroesophageal_Junction'
Esophagus_Gastroesophageal_Junction$n=401
#### Esophagus_Mucosa
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/esophagus_mucosa/V8_master_Esophagus_Mucosa_RP_lists.Rdata")
Esophagus_Mucosa
Esophagus_Mucosa$tissue='Esophagus_Mucosa'
Esophagus_Mucosa$n=622
#### Esophagus_Muscularis
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/esophagus_muscularis/V8_master_Esophagus_Muscularis_RP_lists.Rdata")
Esophagus_Muscularis
Esophagus_Muscularis$tissue='Esophagus_Muscularis'
Esophagus_Muscularis$n=559
#### Heart_Atrial_Appendage
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/heart_atrial_app/V8_master_Heart_Atrial_Appendage_RP_lists.Rdata")
Heart_Atrial_Appendage
Heart_Atrial_Appendage$tissue='Heart_Atrial_Appendage'
Heart_Atrial_Appendage$n=452
#### Heart_Left_Ventricle
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/heart_left_ventricle/V8_master_Heart_Left_Ventricle_RP_lists.Rdata")
Heart_Left_Ventricle
Heart_Left_Ventricle$tissue='Heart_Left_Ventricle'
Heart_Left_Ventricle$n=689
#### Kidney_Cortex
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/kidney/V8_master_Kidney_Cortex_RP_lists.Rdata")
Kidney_Cortex
Kidney_Cortex$tissue='Kidney_Cortex'
Kidney_Cortex$n=100
#### Liver
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/liver/V8_master_Liver_RP_lists.Rdata")
Liver
Liver$tissue='Liver'
Liver$n=251
#### Lung 
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/lung/V8_master_Lung_RP_lists.Rdata")
Lung 
Lung$tissue='Lung'
Lung$n=867
#### Minor_Salivary_Gland
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/minor_salivary_gland/V8_master_Minor_Salivary_Gland_RP_lists.Rdata")
Minor_Salivary_Gland
Minor_Salivary_Gland$tissue='Minor_Salivary_Gland'
Minor_Salivary_Gland$n=181
#### Muscle_Skeletal
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/muscle_skeletal/V8_master_Muscle_Skeletal_RP_lists.Rdata")
Muscle_Skeletal
Muscle_Skeletal$tissue='Muscle_Skeletal'
Muscle_Skeletal$n=1132
#### Nerve_Tibial
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/nerve_tibial/V8_master_Nerve_Tibial_RP_lists.Rdata")
Nerve_Tibial
Nerve_Tibial$tissue='Nerve_Tibial'
Nerve_Tibial$n=722
#### Ovary
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/ovary/V8_master_Ovary_RP_lists.Rdata")
Ovary
Ovary$tissue='Ovary'
Ovary$n=195
#### Pancreas
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/pancreas/V8_master_Pancreas_RP_lists.Rdata")
Pancreas
Pancreas$tissue='Pancreas'
Pancreas$n=360
#### Pituitary
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/pituitary/V8_master_Pituitary_RP_lists.Rdata")
Pituitary
Pituitary$tissue='Pituitary'
Pituitary$n=301
#### Prostate
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/prostate/V8_master_Prostate_RP_lists.Rdata")
Prostate
Prostate$tissue='Prostate'
Prostate$n=262
#### Skin_Not_Sun_Exposed_Suprapubic
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/skin_suprapub/V8_master_Skin_Not_Sun_Exposed_Suprapubic_RP_lists.Rdata")
Skin_Not_Sun_Exposed_Suprapubic
Skin_Not_Sun_Exposed_Suprapubic$tissue='Skin_Not_Sun_Exposed_Suprapubic'
Skin_Not_Sun_Exposed_Suprapubic$n=638
#### Skin_Sun_Exposed_Lower_leg
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/skin_leg/V8_master_Skin_Sun_Exposed_Lower_leg_RP_lists.Rdata")
Skin_Sun_Exposed_Lower_leg
Skin_Sun_Exposed_Lower_leg$tissue='Skin_Sun_Exposed_Lower_leg'
Skin_Sun_Exposed_Lower_leg$n=849
#### Small_Intestine_Terminal_Ileum
Small_Intestine_Terminal_Ileum
Small_Intestine_Terminal_Ileum$tissue='Small_Intestine_Terminal_Ileum'
Small_Intestine_Terminal_Ileum$n=193
str(Small_Intestine_Terminal_Ileum)
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/small_intestine_tl/V8_master_Small_Intestine_Terminal_Ileum_RP_lists.Rdata")
#### Spleen
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/spleen/V8_master_Spleen_RP_lists.Rdata")
Spleen
Spleen$tissue='Spleen'
Spleen$n=260
#### Stomach
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/stomach/V8_master_Stomach_RP_lists.Rdata")
Stomach
Stomach$tissue='Stomach'
Stomach$n=381
#### Testis
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/testis/V8_master_Testis_RP_lists.Rdata")
Testis
Testis$tissue='Testis'
Testis$n=406
#### Thyroid
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/thyroid/V8_master_Thyroid_RP_lists.Rdata")
Thyroid
Thyroid$tissue='Thyroid'
Thyroid$n=812
#### Uterus
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/uterus/V8_master_Uterus_RP_lists.Rdata")
Uterus
Uterus$tissue='Uterus'
Uterus$n=166
#### Vagina
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/vagina/V8_master_Vagina_RP_lists.Rdata")
Vagina
Vagina$tissue='Vagina'
Vagina$n=173
#### Whole_Blood
load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/blood/V8_master_Whole_Blood_RP_lists.Rdata")
Whole_Blood
Whole_Blood$tissue='Whole_Blood'
Whole_Blood$n=3288
#### rbind

master_v8_n=rbind(adi_subcut, adipose_vo, Adrenal_Gland, amygdala, artery_aorta, Artery_Tibial, Artery_Coronary, Brain_Caudate_basal_ganglia,
                Brain_Anterior_cingulate_cortex_BA24,Brain_Cerebellar_Hemisphere,Brain_Cerebellum,Brain_Cortex, Brain_Frontal_Cortex_BA9,
                Brain_Hippocampus,Brain_Hypothalamus, Brain_Nucleus_accumbens_basal_ganglia,Brain_Putamen_basal_ganglia,Brain_Spinal_cord_cervical,
                Brain_Substantia_nigra,Breast_Mammary_Tissue,Cells_Cultured_fibroblasts,Cells_EBV_transformed_lymphocytes, Colon_Sigmoid, 
                Colon_Transverse,Esophagus_Gastroesophageal_Junction,Esophagus_Mucosa,Esophagus_Muscularis,Heart_Atrial_Appendage,
                Heart_Left_Ventricle,Kidney_Cortex,Liver,Lung,Minor_Salivary_Gland,Muscle_Skeletal,Nerve_Tibial,Ovary,Pancreas,Pituitary,
                Prostate,Skin_Not_Sun_Exposed_Suprapubic,Skin_Sun_Exposed_Lower_leg,Small_Intestine_Terminal_Ileum,Spleen,Stomach,Testis,
                Thyroid,Uterus,Vagina,Whole_Blood)
save(myvars, myvars_all_mito, myvars_mito_39S, myvars_40S_60S, myvars_40S, myvars_60S, master_v8_n,ao,
     file="./V8_master_tissues_n.Rdata")
write.csv(master_v8_n, 'v8_master_tissues_n.csv')

###########################################################################################
##### Get just the RP data 
###########################################################################################
setwd('C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8')

load("C:/Users/charl/Dropbox/Harvard/ribosomal_proteins/Main_Tables_Paper/v8/V8_master_tissues_n.Rdata")

#### Subset to the RP genes
RP_master=master_v8_n[which(master_v8_n$gene_name %in% myvars),]
RP_master$gene_name <- factor(RP_master$gene_name)
table(RP_master$gene_name)
#write.csv(RP_master, 'RP_v8.csv')

#### Look at the data
RP_sig=RP_master[which(RP_master$qval<=0.05),]
summary(RP_sig$pval)
#write.csv(RP_sig, 'RP_sig.csv')

stringent_RP_sig_beta=RP_sig[which(RP_sig$pval<0.000005),]
#write.csv(stringent_RP_sig_beta, 'stringent_RP_sig_beta.csv')

tissue_gene=table(stringent_RP_sig_beta$tissue, stringent_RP_sig_beta$gene_name)
write.csv(tissue_gene, 'tissue_gene.csv')
#NOTE: there are more than 1 eQTL for some genes

#### Heatmap
tissue_gene_mat <- cor(t(tissue_gene))
head(tissue_gene_mat)

data <- as.matrix(tissue_gene_mat)
write.csv(data, 'matrix_tissue_correlation_by_RPs_across_genes.csv')
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

heatmap.2(data, 
          col=coul, 
          main="Between-tissue comparison of ribosomal protein eQTLs by gene (P<0.000005)",
          distfun=function(x) dist(x, method="euclidean"), 
          hclustfun=function(x) hclust(x, method="ward.D2"))

dev.off()
#### This code increases the margins for the labels

par(mar=c(7,4,4,2)+0.1) 
png(filename='test.png', width=1200, height=1200)

heatmap.2(data, 
          col=coul, 
          main="Between-tissue comparison of RP eQTLs (P<0.000005)",
          distfun=function(x) dist(x, method="euclidean"), 
          hclustfun=function(x) hclust(x, method="ward.D2"), cexRow=1,cexCol=1,margins=c(20,20),trace="none")

graphics.off()
#####

tissue_gene2=table(stringent_RP_sig_beta$tissue, stringent_RP_sig_beta$gene_name)

# Library
#library(igraph)

tissue_gene_mat <- cor(t(tissue_gene2))
rownames(tissue_gene_mat)

# Keep only high correlations
tissue_gene_mat[tissue_gene_mat<0.60] <- 0

# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix(tissue_gene_mat, weighted=T, 
                                       mode="undirected", diag=F)

# Basic chart
plot(network)

# title and legend
text(0,1.25,"Tissue networks by shared ribosomal protein eQTLs",col="black", cex=2.5)

#par(mfrow=c(2,3)) #Code to have more than one plot in the plot space
dev.off()







