# install packages 

library(phyloseq)

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggsci)

source("/Users/xiyupeng/Documents/Research/metagenomics/src/functions.R")

theme_set(theme_bw())

setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME")
data<-import_biom("qiimegoodformat.biom.txt",treefilename= "rep.tre.txt",
                  refseqfilename = "rep.fasta.txt")
sample.design<-import_qiime_sample_data(mapfilename = "qiime.mapping.txt")

y<-tax_table(data)
colnames(y)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
sample.design$Trt
levels(sample.design$Trt)<-c("NC","PC","SF-","SF+","IF-","IF+")
sample.design$Diet
levels(sample.design$Diet)<-c("Basal","IF-","IF+","SF-","SF+")
#diets.name<-c("Basal","SF","SFE","IF","IFE")
#trts.name<-c("NC","PC","SBP","SBP+E","DDGS","DDGS+E")
data0<-phyloseq(otu_table(data),tax_table(y),sample_data(sample.design),
                refseq(data),phy_tree(data))
data0 <- prune_taxa(taxa_sums(data0) > 0, data0)


#-------------------------------------------------------------------------
#Distribution of number of reads in samples
readsumsdf = data.frame(nreads = sort(taxa_sums(data0), TRUE), sorted = 1:ntaxa(data0), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(data0), 
                                                        TRUE), sorted = 1:nsamples(data0), type = "SampleID"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#-------------------------------------------------------------------------
# Filtering due to experimental performance

#Here, we may delete the data of pig 25 and 26, 
#based on the bad performance in the experiment. 

#data0<-subset_samples(data0, Pen_ID!=25)
#data0<-subset_samples(data0, Pen_ID!=26)
#data0 <- prune_taxa(taxa_sums(data0) > 0, data0)
#data0

#-----------------------------------------------------------------------
# relative abundance

#separate the OTU table for three different samples
data_c<-subset_samples(data0, Sample_type=="Colonic digesta")
data_c <- prune_taxa(taxa_sums(data_c) > 0, data_c)
data_f<-subset_samples(data0, Sample_type=="Feces")
data_f <- prune_taxa(taxa_sums(data_f) > 0, data_f)
data_i<-subset_samples(data0, Sample_type=="Ileal digesta")
data_i <- prune_taxa(taxa_sums(data_i) > 0, data_i)
data_c
data_f
data_i
data_i.otu <-rownames(otu_table(data_i))
data_f.otu <-rownames(otu_table(data_f))
data_c.otu <-rownames(otu_table(data_c))

#----------------------------------------------------------------

# if we only consider sample type ileal and colonic 
samples<-c("Colonic digesta", "Ileal digesta")
data_ori<-data0      # keep a copy of data0
data_ci<-subset_samples(data0, Sample_type %in% samples)
data_ci <- prune_taxa(taxa_sums(data_ci) > 0, data_ci)
data_ci


########  codes below just try different package and method (NOT USED !!!!!!!!!!!!!!!!!!!!!!!!!)
#--------------------------------------------------------------------------------------------------------
# permutation test at phylum level 
# Not used here 

#data0<-data_ci

# consider the genus level
levels(sample_data(data0)$Trt)<-c("U","U1","W","X","Y","Z")
#data0<-tax_glom(data0,"Genus")    # camparison analysis on genus level
#data0<-tax_glom(data0,"Family")   # comparison analysis on family level
data0<-tax_glom(data0,"Phylum")   # comparison analysis on phylum level 

# remove phylums that are only seen in under 10% total samples.
prev0 <-apply(X = otu_table(data0),
              MARGIN = ifelse(taxa_are_rows(data0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
data0 <- prune_taxa((prev0 > 0.1 * nsamples(data0) ), data0)
data0

##------------------------------------------------------------------------------
# microbiomeSeq

library(microbiomeSeq)
library(ggplot2)
library(adespatial)

# alpha_diversity
data5<-data_i   # use phyloseq object
#data5<-data_c

physeq<-phyloseq(otu_table(t(as.matrix(otu_table(data5))),taxa_are_rows = FALSE),tax_table(data5),sample_data(data5),refseq(data5),phy_tree(data5))

data_phylum <- taxa_level(physeq, "Phylum")

p <- plot_anova_diversity(physeq, method = c("richness", "simpson", "shannon"), 
                          grouping_column = "Trt", pValueCutoff = 0.05)
print(p)

# beta_diversity contribution
data5<-data_c_2   # After filtering
physeq<-phyloseq(otu_table(t(as.matrix(otu_table(data5))),taxa_are_rows = FALSE),tax_table(data5),sample_data(data5),refseq(data5),phy_tree(data5))
physeq <- taxa_level(physeq, "Genus")
physeq <- normalise_data(physeq, norm.method = "relative")
p <- plot_taxa(physeq, grouping_column = "Trt", method = "hellinger", number.taxa = 25, 
               filename = NULL)
print(p)  # 11 for ileal, 25 for colonic 

##-----------------------------------------------------------------
# have bugs, not work
ord.res <- ordination(physeq, distance = "bray", method = "PCoA", grouping_column = "Trt", 
                      pvalue.cutoff = 0.05)

p <- plot_ordination(ord.res, method = "PCoA", pvalue.cutoff = 0.05, show.pvalues = T, 
                     num.signi.groups = NULL)
print(p)

###---------------------------------------------------------------------
# try plot heatmap 

##ileal
data_i_map <- tax_glom(data_i, "Family")
data_i_map<-transform_sample_counts(data_i_map,function(x) 100 * x/sum(x))
GP2 = prune_taxa(names(sort(taxa_sums(data_i_map), TRUE)[1:30]), data_i_map)
p<-plot_heatmap(GP2, "PCoA", "unifrac", "Trt", "Family",weighted = TRUE,sample.order="Trt")

png(filename = "heatmap_ileal_family.png",
    width = 3000, height = 1200, res = 300)

p

dev.off()

###colonic
data_c_map <- tax_glom(data_c, "Family")
data_c_map<-transform_sample_counts(data_c_map,function(x) 100 * x/sum(x))
GP2 = prune_taxa(names(sort(taxa_sums(data_c_map), TRUE)[1:30]), data_c_map)
p<-plot_heatmap(GP2, "PCoA", "unifrac", "Trt", "Family",weighted = TRUE,sample.order="Trt", trans = log_trans(4))

png(filename = "heatmap_colonic_family.png",
    width = 3000, height = 1200, res = 300)

p

dev.off()

###----------------------------------------------------------------------
### mixOmics

library(mixOmics)

### import environment factor
Ileal_genes <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/Ileal_genes.csv")
colon_vfa <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/colon_vfa.csv")

rownames(Ileal_genes)<-Ileal_genes[,1]
Ileal_genes<-Ileal_genes[,1:17]
rownames(colon_vfa)<-colon_vfa[,1]
colon_vfa<-colon_vfa[,c(1,10:15)]

###ileal

### prepare sample data 
sample.data<-sample_data(data_i)[,c("Pen_ID","Trt")]
sample.data$names<-rownames(sample.data)
sample.data<-as.data.frame(cbind(sample.data@.Data[[1]],sample.data@.Data[[3]],as.character(sample.data@.Data[[2]])))
vfa.sample<-merge(colon_vfa,sample.data,by.x = "Pen",by.y = "V1")
rownames(vfa.sample)<-vfa.sample$V2
vfa.sample<-vfa.sample[order(vfa.sample$V2),]

### prepare otu data
#data_i.rela<-transform_sample_counts(data_i,function(x) 100 * x/sum(x))
#ileal.otu<-t(otu_table(data_i.rela))
ileal.otu<-as.matrix(otu_table(data_i))

###colonic

####---------------------------------------------------------------------
## try plot network

# ileal 
ig = make_network(data_i, type = "samples", distance = "bray", max.dist = 0.3)
plot_network(ig, data_i, color = "Trt",line_weight = 0.4, 
             label = NULL)


####other
##-------------------------------------------------------------------
### ileal 

data_i_map <- data_i
data_i_map<-transform_sample_counts(data_i_map,function(x) 100 * x/sum(x))
GP2 = prune_taxa(names(sort(taxa_sums(data_i_map), TRUE)[1:100]), data_i_map)
GP.wUF.ord = ordinate(GP2, "PCoA", "unifrac", weighted = TRUE) 
p2 = plot_ordination(GP2, GP.wUF.ord, type = "samples", color = "Trt")
p2 + geom_point(size = 5) 

