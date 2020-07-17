### this script is going to test differential expressed phylum under different treatment

### require phyloseq object of 
#### data_c
#### data_i
#### data_f


## set the address of working space 
library(phyloseq)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")

##------------------------------------------------------------------
## differential analysis in ileal 

data5<-data_i  
physeq<-phyloseq(otu_table(t(as.matrix(otu_table(data5))),taxa_are_rows = FALSE),tax_table(data5),sample_data(data5),refseq(data5),phy_tree(data5))
physeq <- tax_glom(physeq, "Family")  # Phylum, Genus
# data transformation
#physeq_log_re <- normalise_data(physeq, norm.method = "log-relative")
physeq_log_re<- transform_sample_counts(physeq,function(x) 100 * (x)/sum(x))

# filtering 
prev0 <-apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
physeq <- prune_taxa((prev0 > 0.1 * nsamples(physeq) ), physeq_log_re)

# We just think about just calculate abundant phylum (max > 0.5%)
max_abun<-apply(taxa_abundace_p_i[,1:6],1,max)
topp_i<-taxa_abundace_p_i$Phylum[which(max_abun>0.5)]
#topp<-c("p__Proteobacteria","p__Actinobacteria","p__Bacteroidetes","p__Firmicutes")
#physeq<-prune_taxa(Phylum %in% topp, physeq)
physeq = subset_taxa(physeq, Phylum %in% topp_i)
physeq = subset_taxa(physeq, Family %in% topf_i)
physeq = subset_taxa(physeq, Genus %in% topg_i)

# compare NC and PC 
samples<-c("NC", "PC")
physeq_NC_PC<-subset_samples(physeq, Trt %in% samples)

result_NC_PC<-kw_test(physeq_NC_PC,"BH")

result_NC_PC<-data.frame(cbind(result_NC_PC,rep("NC vs PC",nrow(result_NC_PC)))) 
colnames(result_NC_PC)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_NC_PC)<-c("Family","p.value"," p.adj","comparison")
result_NC_PC
topp_i
result_NC_PC$Phylum <-  topp_i
topf_i
result_NC_PC$Family <- topf_i

# compare PC and SF-
samples<-c("PC", "SF-")
physeq_PC_W<-subset_samples(physeq, Trt %in% samples)

result_PC_W<-kw_test(physeq_PC_W,"BH")

result_PC_W<-data.frame(cbind(result_PC_W,rep("PC vs SF-",nrow(result_NC_PC))))
colnames(result_PC_W)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_PC_W)<-c("Family","p.value"," p.adj","comparison")

result_PC_W
topp_i
result_PC_W$Phylum <-  topp_i
topf_i
result_PC_W$Family <-  topf_i

# compare PC and SF+
samples<-c("PC", "SF+")
physeq_PC_X<-subset_samples(physeq, Trt %in% samples)

result_PC_X<-kw_test(physeq_PC_X,"BH")

result_PC_X<-data.frame(cbind(result_PC_X,rep("PC vs SF+",nrow(result_NC_PC))))
colnames(result_PC_X)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_PC_X)<-c("Family","p.value"," p.adj","comparison")
result_PC_X
topp_i
result_PC_X$Phylum <-  topp_i
topf_i
result_PC_X$Family <- topf_i

# compare PC and IF-
samples<-c("PC", "IF-")
physeq_PC_Y<-subset_samples(physeq, Trt %in% samples)

result_PC_Y<-kw_test(physeq_PC_Y,"BH")

result_PC_Y<-data.frame(cbind(result_PC_Y,rep("PC vs IF-",nrow(result_NC_PC))))
colnames(result_PC_Y)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_PC_Y)<-c("Family","p.value"," p.adj","comparison")
result_PC_Y
topp_i
result_PC_Y$Phylum <-  topp_i
topf_i
result_PC_Y$Family <- topf_i

# compare PC and IF+
samples<-c("PC", "IF+")
physeq_PC_Z<-subset_samples(physeq, Trt %in% samples)

result_PC_Z<-kw_test(physeq_PC_Z,"BH")

result_PC_Z<-data.frame(cbind(result_PC_Z,rep("PC vs IF+",nrow(result_NC_PC))))
colnames(result_PC_Z)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_PC_Z)<-c("Family","p.value"," p.adj","comparison")
result_PC_Z
topp_i
result_PC_Z$Phylum <-  topp_i
topf_i
result_PC_Z$Family <- topf_i


# compare SF- and SF+
samples<-c("SF-", "SF+")
physeq_W_X<-subset_samples(physeq, Trt %in% samples)

result_W_X<-kw_test(physeq_W_X,"BH")

result_W_X<-data.frame(cbind(result_W_X,rep("SF- vs SF+",nrow(result_NC_PC))))
colnames(result_W_X)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_W_X)<-c("Family","p.value"," p.adj","comparison")
result_W_X
topp_i
result_W_X$Phylum <-  topp_i
topf_i
result_W_X$Family <-topf_i

# compare IF- and IF+
samples<-c("IF-", "IF+")
physeq_Y_Z<-subset_samples(physeq, Trt %in% samples)

result_Y_Z<-kw_test(physeq_Y_Z,"BH")

result_Y_Z<-data.frame(cbind(result_Y_Z,rep("IF- vs IF+",nrow(result_NC_PC))))
colnames(result_Y_Z)<-c("Phylum","p.value"," p.adj","comparison")
colnames(result_Y_Z)<-c("Family","p.value"," p.adj","comparison")
result_Y_Z
topp_i
result_Y_Z$Phylum <-  topp_i
topf_i
result_Y_Z$Family <-topf_i

result<-rbind(result_NC_PC,result_PC_W,result_PC_X,result_PC_Y,result_PC_Z,result_W_X,result_Y_Z)

write.csv(as.data.frame(result),file="ileal_family_ks_test.csv",row.names = FALSE)

##------------------------------------------------------------------
## differential analysis in colonic

data5<-data_c  
physeq<-phyloseq(otu_table(t(as.matrix(otu_table(data5))),taxa_are_rows = FALSE),tax_table(data5),sample_data(data5),refseq(data5),phy_tree(data5))
physeq <- taxa_level(physeq, "Phylum")
# data transformation
physeq_log_re <- normalise_data(physeq, norm.method = "log-relative")

# filtering unprevalence phylum
prev0 <-apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
physeq <- prune_taxa((prev0 > 0.1 * nsamples(physeq) ), physeq_log_re)

# We just think about just calculate abundant phylum (max > 0.5%)
topp<-c("p__Proteobacteria","p__Actinobacteria","p__Bacteroidetes","p__Firmicutes","p__Tenericutes","p__Spirochaetae",
        "p__Verrucomicrobia","p__Fusobacteria","p__Saccharibacteria")
physeq<-prune_taxa(colnames(otu_table(physeq)) %in% topp,physeq)

# compare NC and PC 
samples<-c("NC", "PC")
physeq_NC_PC<-subset_samples(physeq, Trt %in% samples)

result_NC_PC<-kw_test(physeq_NC_PC,"BH")

result_NC_PC<-data.frame(cbind(result_NC_PC,rep("NC vs PC",9)))
colnames(result_NC_PC)<-c("Phylum","p.value"," p.adj","comparison")
result_NC_PC

# compare PC and SF-
samples<-c("PC", "SF-")
physeq_PC_W<-subset_samples(physeq, Trt %in% samples)

result_PC_W<-kw_test(physeq_PC_W,"BH")

result_PC_W<-data.frame(cbind(result_PC_W,rep("PC vs SF-",9)))
colnames(result_PC_W)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_W

# compare PC and SF+
samples<-c("PC", "SF+")
physeq_PC_X<-subset_samples(physeq, Trt %in% samples)

result_PC_X<-kw_test(physeq_PC_X,"BH")

result_PC_X<-data.frame(cbind(result_PC_X,rep("PC vs SF+",9)))
colnames(result_PC_X)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_X


# compare PC and IF-
samples<-c("PC", "IF-")
physeq_PC_Y<-subset_samples(physeq, Trt %in% samples)

result_PC_Y<-kw_test(physeq_PC_Y,"BH")

result_PC_Y<-data.frame(cbind(result_PC_Y,rep("PC vs IF-",9)))
colnames(result_PC_Y)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_Y

# compare PC and IF+
samples<-c("PC", "IF+")
physeq_PC_Z<-subset_samples(physeq, Trt %in% samples)

result_PC_Z<-kw_test(physeq_PC_Z,"BH")

result_PC_Z<-data.frame(cbind(result_PC_Z,rep("PC vs IF+",9)))
colnames(result_PC_Z)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_Z

# compare SF- and SF+
samples<-c("SF-", "SF+")
physeq_W_X<-subset_samples(physeq, Trt %in% samples)

result_W_X<-kw_test(physeq_W_X,"BH")

result_W_X<-data.frame(cbind(result_W_X,rep("SF- vs SF+",9)))
colnames(result_W_X)<-c("Phylum","p.value"," p.adj","comparison")
result_W_X

# compare IF- and IF+
samples<-c("IF-", "IF+")
physeq_Y_Z<-subset_samples(physeq, Trt %in% samples)

result_Y_Z<-kw_test(physeq_Y_Z,"BH")

result_Y_Z<-data.frame(cbind(result_Y_Z,rep("IF- vs IF+",9)))
colnames(result_Y_Z)<-c("Phylum","p.value"," p.adj","comparison")
result_Y_Z

result<-rbind(result_NC_PC,result_PC_W,result_PC_X,result_PC_Y,result_PC_Z,result_W_X,result_Y_Z)

write.csv(as.data.frame(result),file="colonic_phylum_ks_test.csv",row.names = FALSE)

###------------------------------------
## feces
library(microbiomeSeq)

data5<-data_f 
physeq<-phyloseq(otu_table(t(as.matrix(otu_table(data5))),taxa_are_rows = FALSE),tax_table(data5),sample_data(data5),refseq(data5),phy_tree(data5))
physeq <- taxa_level(physeq, "Phylum")
# data transformation
physeq_log_re <- normalise_data(physeq, norm.method = "log-relative")

# filtering unprevalence phylum
prev0 <-apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
physeq <- prune_taxa((prev0 > 0.1 * nsamples(physeq) ), physeq_log_re)

# We just think about just calculate abundant phylum (max > 0.5%)
topp<-c("p__Proteobacteria","p__Actinobacteria","p__Bacteroidetes","p__Firmicutes","p__Tenericutes","p__Spirochaetae",
        "p__Verrucomicrobia","p__Planctomycetes","p__Synergistetes")
physeq<-prune_taxa(colnames(otu_table(physeq)) %in% topp,physeq)

# Compare NC and PC before challenge
samples<-c("NC", "PC")
physeq_PC_NC<-subset_samples(physeq, Trt %in% samples)

result_PC_NC<-kw_test(physeq_PC_NC,"BH")

result_PC_NC<-data.frame(cbind(result_PC_NC,rep("PC vs NC",9)))
colnames(result_PC_NC)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_NC

# phylum
write.csv(as.data.frame(result_PC_NC),file="Feces_U1vsU_Phylum.csv",row.names = TRUE)

# compare Basal and SF-
samples<-c("Basal", "SF-")
physeq_PC_W<-subset_samples(physeq, Diet %in% samples)

result_PC_W<-kw_test(physeq_PC_W,"BH")

result_PC_W<-data.frame(cbind(result_PC_W,rep("Basal vs SF-",9)))
colnames(result_PC_W)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_W

# compare Basal and SF+
samples<-c("Basal", "SF+")
physeq_PC_X<-subset_samples(physeq, Diet %in% samples)

result_PC_X<-kw_test(physeq_PC_X,"BH")

result_PC_X<-data.frame(cbind(result_PC_X,rep("Basal vs SF+",9)))
colnames(result_PC_X)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_X


# compare Basal and IF-
samples<-c("Basal", "IF-")
physeq_PC_Y<-subset_samples(physeq, Diet %in% samples)

result_PC_Y<-kw_test(physeq_PC_Y,"BH")

result_PC_Y<-data.frame(cbind(result_PC_Y,rep("Basal vs IF-",9)))
colnames(result_PC_Y)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_Y

# compare Basal and IF+
samples<-c("Basal", "IF+")
physeq_PC_Z<-subset_samples(physeq, Diet %in% samples)

result_PC_Z<-kw_test(physeq_PC_Z,"BH")

result_PC_Z<-data.frame(cbind(result_PC_Z,rep("Basal vs IF+",9)))
colnames(result_PC_Z)<-c("Phylum","p.value"," p.adj","comparison")
result_PC_Z

# compare SF- and SF+
samples<-c("SF-", "SF+")
physeq_W_X<-subset_samples(physeq, Diet %in% samples)

result_W_X<-kw_test(physeq_W_X,"BH")

result_W_X<-data.frame(cbind(result_W_X,rep("SF- vs SF+",9)))
colnames(result_W_X)<-c("Phylum","p.value"," p.adj","comparison")
result_W_X

# compare IF- and IF+
samples<-c("IF-", "IF+")
physeq_Y_Z<-subset_samples(physeq, Diet %in% samples)

result_Y_Z<-kw_test(physeq_Y_Z,"BH")

result_Y_Z<-data.frame(cbind(result_Y_Z,rep("IF- vs IF+",9)))
colnames(result_Y_Z)<-c("Phylum","p.value"," p.adj","comparison")
result_Y_Z

result<-rbind(result_PC_W,result_PC_X,result_PC_Y,result_PC_Z,result_W_X,result_Y_Z)

write.csv(as.data.frame(result),file="feces_phylum_ks_test.csv",row.names = FALSE)

