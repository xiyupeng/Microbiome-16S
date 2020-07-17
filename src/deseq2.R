### this script is going to test differential expressed genus and family under different treatment

### require phyloseq object of 
#### data_0
#### taxa_abundace_g_i taxa_abundace_f_i
#### taxa_abundace_g_c taxa_abundace_f_c
#### taxa_abundace_g_f taxa_abundace_f_f

## set the address of working space 
library(phyloseq)
library(DESeq2)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")
source("/Users/xiyupeng/Documents/Research/metagenomics/src/functions.R")

#----------------------------------------------------------------
##Differential Abundance for Microbiome Data

### Analysis across all sample_types
#Use DESeq2 for the comparison metagenomics analysis

# if we only consider colonic and ileal samples 
#data0<-data_ci

# consider the genus level
levels(sample_data(data0)$Trt)<-c("U","U1","W","X","Y","Z")
data0<-tax_glom(data0,"Genus")    # camparison analysis on genus level
data0<-tax_glom(data0,"Family")   # comparison analysis on family level
data0<-tax_glom(data0,"Phylum")   # comparison analysis on phylum level 

# remove OTUs that are only seen in under 10% total samples.
prev0 <-apply(X = otu_table(data0),
              MARGIN = ifelse(taxa_are_rows(data0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
data0 <- prune_taxa((prev0 > 0.1 * nsamples(data0) ), data0)
data0

des_trt = phyloseq_to_deseq2(data0,design = ~NULL)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(des_trt), 1, gm_mean)
des_trt = estimateSizeFactors(des_trt, geoMeans = geoMeans)
design.group = paste0(des_trt$Trt, des_trt$Sample_type)

# if include feces 
design.group[59:71]<-"BasalFeces"

des_trt$group <-factor(design.group)
design(des_trt) <- ~ 0+group
#pre_filtering ?
#keep <-rowSums(counts(des_trt)) >100
#des_trt<-des_trt[keep,]
diagdds = DESeq(des_trt,fitType = "parametric",test= "Wald")
plotDispEsts(diagdds)
# access the normalized data
counts_norm<-counts(diagdds,normalized=TRUE)
resultsNames(diagdds)

#-----------------------------------------------------------
# Analysis based on the samples of Ileal digesta
alpha<-1

#genus level
max_abun<-apply(taxa_abundace_g_i[,1:6],1,max)
topg_i<-taxa_abundace_g_i$Genus[which(max_abun>1)]   #>1%
#table_g_i_1<-subset(taxa_abundace_g_i, Genus %in% topg_i)
#write.csv(table_g_i_1,file="abundance_Ileal_genus.csv",row.names = TRUE)
# compare "groupU1Ileal.digesta"  VS  "groupUIleal.digesta" 

# family level 
max_abun<-apply(taxa_abundace_f_i[,1:6],1,max)
topf_i<-taxa_abundace_f_i$Family[which(max_abun>1)]   #>1%

# consider feces, ileal, colonic at the same time 
contrast<-rep(0,17)
contrast[3]<-1
contrast[5]<--1

# only consider ileal and colonic 
#contrast<-rep(0,12)
#contrast[2]<-1
#contrast[4]<--1

# genus 
U1_U_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
# family 
U1_U_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

U1_U_Ileal.ind<-rep("U1_vs_U",nrow(U1_U_Ileal))
U1_U_Ileal.mean<-group_mean(diagdds,rownames(U1_U_Ileal),counts_norm,
                            "U1Ileal digesta","UIleal digesta",NULL)
U1_U_Ileal<-cbind(U1_U_Ileal,U1_U_Ileal.ind,U1_U_Ileal.mean)
colnames(U1_U_Ileal)[(ncol(U1_U_Ileal)-2):ncol(U1_U_Ileal)]<-c("comparison",
                                                               "group1_mean","group2_mean")

# compare "groupWIleal.digesta" VS "groupU1Ileal.digesta"
contrast<-rep(0,17)
contrast[3]<--1
contrast[8]<-1

#contrast<-rep(0,12)
#contrast[2]<--1
#contrast[6]<-1

#genus
W_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
#family
W_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

W_U1_Ileal.ind<-rep("W_vs_U1",nrow(W_U1_Ileal))
W_U1_Ileal.mean<-group_mean(diagdds,rownames(W_U1_Ileal),counts_norm,
                            "WIleal digesta","U1Ileal digesta",NULL)
W_U1_Ileal<-cbind(W_U1_Ileal,W_U1_Ileal.ind,W_U1_Ileal.mean)
colnames(W_U1_Ileal)[(ncol(W_U1_Ileal)-2):ncol(W_U1_Ileal)]<-c("comparison",
                                                               "group1_mean","group2_mean")

# compare "groupXIleal.digesta" VS "groupU1Ileal.digesta"
contrast<-rep(0,17)
contrast[3]<--1
contrast[11]<-1

#contrast<-rep(0,12)
#contrast[2]<--1
#contrast[8]<-1

#genus
X_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
#family
X_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

X_U1_Ileal.ind<-rep("X_vs_U1",nrow(X_U1_Ileal))
X_U1_Ileal.mean<-group_mean(diagdds,rownames(X_U1_Ileal),counts_norm,
                            "XIleal digesta","U1Ileal digesta",NULL)
X_U1_Ileal<-cbind(X_U1_Ileal,X_U1_Ileal.ind,X_U1_Ileal.mean)
colnames(X_U1_Ileal)[(ncol(X_U1_Ileal)-2):ncol(X_U1_Ileal)]<-c("comparison",
                                                               "group1_mean","group2_mean")

# compare "groupYIleal.digesta" VS "groupU1Ileal.digesta"
contrast<-rep(0,17)
contrast[3]<--1
contrast[14]<-1

#contrast<-rep(0,12)
#contrast[2]<--1
#contrast[10]<-1
# genus
Y_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
#family
Y_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

Y_U1_Ileal.ind<-rep("Y_vs_U1",nrow(Y_U1_Ileal))
Y_U1_Ileal.mean<-group_mean(diagdds,rownames(Y_U1_Ileal),counts_norm,
                            "YIleal digesta","U1Ileal digesta",NULL)
Y_U1_Ileal<-cbind(Y_U1_Ileal,Y_U1_Ileal.ind,Y_U1_Ileal.mean)
colnames(Y_U1_Ileal)[(ncol(Y_U1_Ileal)-2):ncol(Y_U1_Ileal)]<-c("comparison",
                                                               "group1_mean","group2_mean")

# compare "groupZIleal.digesta" VS "groupU1Ileal.digesta"
contrast<-rep(0,17)
contrast[3]<--1
contrast[17]<-1

#contrast<-rep(0,12)
#contrast[2]<--1
#contrast[12]<-1
#genus
Z_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
#family
Z_U1_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

Z_U1_Ileal.ind<-rep("Z_vs_U1",nrow(Z_U1_Ileal))
Z_U1_Ileal.mean<-group_mean(diagdds,rownames(Z_U1_Ileal),counts_norm,
                            "ZIleal digesta","U1Ileal digesta",NULL)
Z_U1_Ileal<-cbind(Z_U1_Ileal,Z_U1_Ileal.ind,Z_U1_Ileal.mean)
colnames(Z_U1_Ileal)[(ncol(Z_U1_Ileal)-2):ncol(Z_U1_Ileal)]<-c("comparison",
                                                               "group1_mean","group2_mean")
# compare "groupXIleal.digesta" VS "groupWIleal.digesta"
contrast<-rep(0,17)
contrast[8]<--1
contrast[11]<-1

#contrast<-rep(0,12)
#contrast[6]<--1
#contrast[8]<-1
# genus
X_W_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
# family
X_W_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

X_W_Ileal.ind<-rep("X_vs_W",nrow(X_W_Ileal))
X_W_Ileal.mean<-group_mean(diagdds,rownames(X_W_Ileal),counts_norm,"XIleal digesta",
                           "WIleal digesta",NULL)
X_W_Ileal<-cbind(X_W_Ileal,X_W_Ileal.ind,X_W_Ileal.mean)
colnames(X_W_Ileal)[(ncol(X_W_Ileal)-2):ncol(X_W_Ileal)]<-c("comparison",
                                                            "group1_mean","group2_mean")

# compare "groupZIleal.digesta" VS "groupYIleal.digesta"
contrast<-rep(0,17)
contrast[14]<--1
contrast[17]<-1

#contrast<-rep(0,12)
#contrast[10]<--1
#contrast[12]<-1
#genus
Z_Y_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topg_i,"genus")
#family
Z_Y_Ileal<-contra_ana(diagdds,contrast,alpha,data0,topf_i,"family")

Z_Y_Ileal.ind<-rep("Z_vs_Y",nrow(Z_Y_Ileal))
Z_Y_Ileal.mean<-group_mean(diagdds,rownames(Z_Y_Ileal),counts_norm,
                           "ZIleal digesta","YIleal digesta",NULL)
Z_Y_Ileal<-cbind(Z_Y_Ileal,Z_Y_Ileal.ind,Z_Y_Ileal.mean)
colnames(Z_Y_Ileal)[(ncol(Z_Y_Ileal)-2):ncol(Z_Y_Ileal)]<-c("comparison",
                                                            "group1_mean","group2_mean")

Ileal.res<-rbind(U1_U_Ileal,W_U1_Ileal,X_U1_Ileal,
                 Y_U1_Ileal,Z_U1_Ileal,X_W_Ileal,Z_Y_Ileal)
#Ileal.res<-cbind(rownames(Ileal.res),Ileal.res)
Ileal.res<-Ileal.res[,1:(ncol(Ileal.res)-2)]
write.csv(as.data.frame(Ileal.res),file="Ileal.csv",row.names = TRUE)

#----------------------------------------------------------
# Analysis based on the samples of feces

alpha <- 1

# genus level 
max_abun<-apply(taxa_abundace_g_f[,1:5],1,max)
topg_f<-taxa_abundace_g_f$Genus[which(max_abun>1)]
topg_f<-topg_f[topg_f!="g__"]
#table_g_f_1<-subset(taxa_abundace_g_f, Genus %in% topg_f)
#write.csv(table_g_f_1,file="abundance_feces_genus.csv",row.names = TRUE)

# family level 
max_abun<-apply(taxa_abundace_f_f[,1:5],1,max)
topf_f<-taxa_abundance_f_f$Family[which(max_abun>1)]
topf_f<-topf_f[topf_f!="f__"]


# compare "groupWFeces" VS "groupU1Feces" & "groupUFeces"
contrast<-rep(0,17)
contrast[1]<--1
contrast[7]<-1
# genus
SF_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family
SF_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

SF_Basal_feces<-cbind(SF_Basal_feces,rep("SF_vs_Basal",nrow(SF_Basal_feces)))
colnames(SF_Basal_feces)[ncol(SF_Basal_feces)]<-"comparison"

# compare "groupXFeces" VS "groupU1Feces" & "groupUFeces"
contrast<-rep(0,17)
contrast[1]<--1
contrast[10]<-1
# genus
SFE_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family
SFE_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

SFE_Basal_feces<-cbind(SFE_Basal_feces,rep("SFE_vs_Basal",nrow(SFE_Basal_feces)))
colnames(SFE_Basal_feces)[ncol(SFE_Basal_feces)]<-"comparison"

# compare "groupYFeces" VS "groupU1Feces" & "groupUFeces"
contrast<-rep(0,17)
contrast[1]<--1
contrast[13]<-1
# genus
IF_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family
IF_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

IF_Basal_feces<-cbind(IF_Basal_feces,rep("IF_vs_Basal",nrow(IF_Basal_feces)))
colnames(IF_Basal_feces)[ncol(IF_Basal_feces)]<-"comparison"

# compare "groupZFeces" VS "groupU1Feces" & "groupUFeces"
contrast<-rep(0,17)
contrast[1]<--1
contrast[16]<-1
# genus
IFE_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family
IFE_Basal_feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

IFE_Basal_feces<-cbind(IFE_Basal_feces,rep("IFE_vs_Basal",nrow(IFE_Basal_feces)))
colnames(IFE_Basal_feces)[ncol(IFE_Basal_feces)]<-"comparison"

# compare "groupXFeces" VS "groupWFeces"
contrast<-rep(0,17)
contrast[7]<--1
contrast[10]<-1
# genus 
SFE_SF_feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family
SFE_SF_feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

SFE_SF_feces<-cbind(SFE_SF_feces,rep("SFE_vs_SF",nrow(SFE_SF_feces)))
colnames(SFE_SF_feces)[ncol(SFE_SF_feces)]<-"comparison"

# compare "groupZFeces" VS "groupYFeces"
contrast<-rep(0,17)
contrast[13]<--1
contrast[16]<-1
# genus
IFE_IF_feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family 
IFE_IF_feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

IFE_IF_feces<-cbind(IFE_IF_feces,rep("IFE_vs_IF",nrow(IFE_IF_feces)))
colnames(IFE_IF_feces)[ncol(IFE_IF_feces)]<-"comparison"

feces.res<-rbind(SF_Basal_feces,SFE_Basal_feces,IF_Basal_feces,
                 IFE_Basal_feces,SFE_SF_feces,IFE_IF_feces)
#feces.res<-cbind(rownames(feces.res),feces.res)
write.csv(as.data.frame(feces.res),file="feces2.csv",row.names = TRUE)

#--------------------------------------------------------
# Analysis based on the samples of Colonic

# genus level
max_abun<-apply(taxa_abundace_g_c[,1:6],1,max)
topg_c<-taxa_abundace_g_c$Genus[which(max_abun>1)]   
topg_c<-topg_c[topg_c!="g__"]
#table_g_c_1<-subset(taxa_abundace_g_c, Genus %in% topg_c)
#write.csv(table_g_c_1,file="abundance_Colonic_genus.csv",row.names = TRUE)

# family level 
max_abun<-apply(taxa_abundance_f_c[,1:6],1,max)
topf_c<-taxa_abundance_f_c$Family[which(max_abun>1)]
topf_c<-topf_c[topf_c!="f__"]

# compare "groupU1Colonic.digesta"  VS  "groupUColonic.digesta" 
contrast<-rep(0,17)
contrast[2]<-1
contrast[4]<--1

#contrast<-rep(0,12)
#contrast[1]<-1
#contrast[3]<--1

# genus
U1_U_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family 
U1_U_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topf_c,"family")

U1_U_Colonic<-cbind(U1_U_Colonic,rep("U1_vs_U",nrow(U1_U_Colonic)))
colnames(U1_U_Colonic)[ncol(U1_U_Colonic)]<-"comparison"

# compare "groupWColonic.digesta" VS "groupU1Colonic.digesta"
contrast<-rep(0,17)
contrast[2]<--1
contrast[6]<-1

#contrast<-rep(0,12)
#contrast[1]<--1
#contrast[5]<-1

# genus
W_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family
W_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topf_c,"family")


W_U1_Colonic<-cbind(W_U1_Colonic,rep("W_vs_U1",nrow(W_U1_Colonic)))
colnames(W_U1_Colonic)[ncol(W_U1_Colonic)]<-"comparison"

# compare "groupXColonic.digesta" VS "groupU1Colonic.digesta"
contrast<-rep(0,17)
contrast[2]<--1
contrast[9]<-1

#contrast<-rep(0,12)
#contrast[1]<--1
#contrast[7]<-1

# genus
X_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family
X_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topf_c,"family")


X_U1_Colonic<-cbind(X_U1_Colonic,rep("X_vs_U1",nrow(X_U1_Colonic)))
colnames(X_U1_Colonic)[ncol(X_U1_Colonic)]<-"comparison"

# compare "groupYColonic.digesta" VS "groupU1Colonic.digesta"
contrast<-rep(0,17)
contrast[2]<--1
contrast[12]<-1

#contrast<-rep(0,12)
#contrast[1]<--1
#contrast[9]<-1

# genus
Y_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family
Y_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topf_c,"family")

Y_U1_Colonic<-cbind(Y_U1_Colonic,rep("Y_vs_U1",nrow(Y_U1_Colonic)))
colnames(Y_U1_Colonic)[ncol(Y_U1_Colonic)]<-"comparison"

# compare "groupZColonic.digesta" VS "groupU1Colonic.digesta"
contrast<-rep(0,17)
contrast[2]<--1
contrast[15]<-1

#contrast<-rep(0,12)
#contrast[1]<--1
#contrast[11]<-1

# genus
Z_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family
Z_U1_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topf_c,"family")


Z_U1_Colonic<-cbind(Z_U1_Colonic,rep("Z_vs_U1",nrow(Z_U1_Colonic)))
colnames(Z_U1_Colonic)[ncol(Z_U1_Colonic)]<-"comparison"

# compare "groupXColonic.digesta" VS "groupWColonic.digesta"
contrast<-rep(0,17)
contrast[6]<--1
contrast[9]<-1

#contrast<-rep(0,12)
#contrast[5]<--1
#contrast[7]<-1

# genus
X_W_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family
X_W_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topf_c,"family")


X_W_Colonic<-cbind(X_W_Colonic,rep("X_vs_W",nrow(X_W_Colonic)))
colnames(X_W_Colonic)[ncol(X_W_Colonic)]<-"comparison"

# compare "groupZColonic.digesta" VS "groupYColonic.digesta"
contrast<-rep(0,17)
contrast[12]<--1
contrast[15]<-1

#contrast<-rep(0,12)
#contrast[9]<--1
#contrast[11]<-1
# genus
Z_Y_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"genus")
# family
Z_Y_Colonic<-contra_ana(diagdds,contrast,alpha,data0,topg_c,"family")

Z_Y_Colonic<-cbind(Z_Y_Colonic,rep("Z_vs_Y",nrow(Z_Y_Colonic)))
colnames(Z_Y_Colonic)[ncol(Z_Y_Colonic)]<-"comparison"


Colonic.res<-rbind(U1_U_Colonic,W_U1_Colonic,X_U1_Colonic,
                   Y_U1_Colonic,Z_U1_Colonic,X_W_Colonic,Z_Y_Colonic)
#Colonic.res<-cbind(rownames(Colonic.res),Colonic.res)
write.csv(as.data.frame(Colonic.res),file="Colonic.csv",row.names = TRUE)

#----------------------------------------------------------
# Comparing NC and PC group before challenge (only use feces sample)

alpha <- 1

# genus level 
max_abun<-apply(taxa_abundace_g_f[,1:5],1,max)
topg_f<-taxa_abundace_g_f$Genus[which(max_abun>1)]
topg_f<-topg_f[topg_f!="g__"]
#table_g_f_1<-subset(taxa_abundace_g_f, Genus %in% topg_f)
#write.csv(table_g_f_1,file="abundance_feces_genus.csv",row.names = TRUE)

# family level 
max_abun<-apply(taxa_abundace_f_f[,1:5],1,max)
topf_f<-taxa_abundace_f_f$Family[which(max_abun>1)]
topf_f<-topf_f[topf_f!="f__"]


# compare "groupU1Feces" VS "groupUFeces"
contrast<-rep(0,18)
contrast[2]<-1
contrast[5]<--1

# genus
U1_U_Feces<-contra_ana(diagdds,contrast,alpha,data0,topg_f,"genus")
# family
U1_U_Feces<-contra_ana(diagdds,contrast,alpha,data0,topf_f,"family")

U1_U_Feces<-cbind(U1_U_Feces,rep("U1_vs_U",nrow(U1_U_Feces)))
colnames(U1_U_Feces)[ncol(U1_U_Feces)]<-"comparison"

# Genes
write.csv(as.data.frame(U1_U_Feces),file="Feces_U1vsU_Genus.csv",row.names = TRUE)
# family 
write.csv(as.data.frame(U1_U_Feces),file="Feces_U1vsU_Famliy.csv",row.names = TRUE)
