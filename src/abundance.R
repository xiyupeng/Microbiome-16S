### After we have upload the data, this script is going to 
### draw relative abundance table of all groups
### require phyloseq object of 
#### data_i
#### data_c
#### data_f

## set the address of working space 
library(phyloseq)
library(RColorBrewer)
library(ggsci)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")

#-----------------------------------------------------------------------------
# Feces
data_f_m1<-merge_samples(data_f,"Diet")
diets.name<-c("Basal","SF-","SF+","IF-","IF+")

# Phylum
# phylum >0.5% in feces
data_f_m2<-tax_glom(data_f_m1,"Phylum")
data_f_ra<-transform_sample_counts(data_f_m2,function(x){100*x/sum(x)})

taxa_abundace_p_f<- as.data.frame(cbind(t(otu_table(data_f_ra)),tax_table(data_f_ra)))
max_abun<-apply(taxa_abundace_p_f[,1:5],1,max)
#topp_f<-taxa_abundace_p_f$Phylum[which(max_abun>0.5)]  #>0.5%
topp_f<-taxa_abundace_p_f$Phylum[which(max_abun>0.5)] 
topp_f<-topp_f[topp_f!="p__"]
table_p_f<-subset(taxa_abundace_p_f,max_abun>0.5)
data_p_f<-prune_taxa(rownames(subset(table_p_f, Phylum %in% topp_f)), data_f_ra)
#data_p_f<-prune_taxa(rownames(table_p_f), data_f_ra)

# plot relative abundance 
plot_bar(data_p_f,fill="Phylum")+xlab("Diet")+ylab("Percentage of Sequences")+
  scale_fill_igv()+ylim(0,100)#+scale_x_discrete(limits= diets.name)
#subset(taxa_abundace_p_f, Phylum %in% topp_f)
table_p_f

##output the relative abundance table
#write.csv(subset(table_p_f, Phylum %in% topp_f),file = "feces_phylum_relative_abundance.csv",row.names = FALSE)

# family
# family >1% and have been identified in feces
data_f_m2_2<-tax_glom(data_f_m1,"Family")
data_f_ra_2<-transform_sample_counts(data_f_m2_2,function(x){100*x/sum(x)})
taxa_abundace_f_f<- data.frame(cbind(t(otu_table(data_f_ra_2)),tax_table(data_f_ra_2)))
max_abun<-apply(taxa_abundace_f_f[,1:5],1,max)
topf_f<-taxa_abundace_f_f$Family[which(max_abun>1)]  
topf_f<-topf_f[topf_f!="f__"]   #still not show unidentified family
table_f_f<-subset(taxa_abundace_f_f,max_abun>1)
# sort table  450* 600 
table_f_f<-table_f_f[order(table_f_f$Order),]
table_f_f<-table_f_f[order(table_f_f$Class),]
table_f_f<-table_f_f[order(table_f_f$Phylum),]

data_f_f<-prune_taxa(rownames(subset(table_f_f, Family %in% topf_f)), data_f_ra_2)

# plot 
plot_bar(data_f_f,fill="Family")+xlab("Diet")+ylab("Percentage of Sequences")+
  ylim(0,100)+scale_x_discrete(limits= diets.name)+scale_fill_igv()
# taxonomy table
#table_f_f
#write.csv(subset(table_f_f, Family %in% topf_f),file = "feces_family_relative_abundance.csv",row.names = FALSE)

# genus >1% and have been identified in feces
data_f_m2_3<-tax_glom(data_f_m1,"Genus")
data_f_ra_3<-transform_sample_counts(data_f_m2_3,function(x){100*x/sum(x)})
taxa_abundace_g_f<- data.frame(cbind(t(otu_table(data_f_ra_3)),tax_table(data_f_ra_3)))
max_abun<-apply(taxa_abundace_g_f[,1:5],1,max)
topg_f<-taxa_abundace_g_f$Genus[which(max_abun>1)]   #>1%
topg_f<-topg_f[topg_f!="g__"]   #still not show unidentified family
table_g_f<-subset(taxa_abundace_g_f,max_abun>1)
# sort table 
table_g_f<-table_g_f[order(table_g_f$Family),]
table_g_f<-table_g_f[order(table_g_f$Order),]
table_g_f<-table_g_f[order(table_g_f$Class),]
table_g_f<-table_g_f[order(table_g_f$Phylum),]

#data_g_f<-prune_taxa(rownames(table_g_f), data_f_ra_3)
#table_g_f<-subset(taxa_abundace_g_f, Genus %in% topg_f)
data_g_f<-prune_taxa(rownames(subset(table_g_f, Genus %in% topg_f)), data_f_ra_3)

# plot 800*600, 600*450
plot_bar(data_g_f,fill="Genus")+xlab("Diet")+ylab("Percentage of Sequences")+ 
  ylim(0,100)+scale_fill_igv()+scale_x_discrete(limits= diets.name)
# taxonomy table
#table_g_f
#write.csv(subset(table_g_f, Genus %in% topg_f),file="feces_genus_relative_abundance.csv",row.names = FALSE)

#---------------------------------------------------------------
#Ileal 
data_i_m1<-merge_samples(data_i,"Trt")
trts.name<-c("NC","PC","SF-","SF+","IF-","IF+")

# phylum 
# phylum >0.5% in Ileal
data_i_m2<-tax_glom(data_i_m1,"Phylum")
data_i_ra<-transform_sample_counts(data_i_m2,function(x){100*x/sum(x)})
taxa_abundace_p_i<- data.frame(cbind(t(otu_table(data_i_ra)),tax_table(data_i_ra)))
max_abun<-apply(taxa_abundace_p_i[,1:6],1,max)
topp_i<-taxa_abundace_p_i$Phylum[which(max_abun>0.5)]  #>0.5%
table_p_i<-subset(taxa_abundace_p_i,max_abun>0.5)
data_p_i<-prune_taxa(rownames(table_p_i), data_i_ra)
plot_bar(data_p_i,fill="Phylum")+xlab("Treatment")+ylab("Percentage of Sequences")+
  ylim(0,100)+scale_fill_igv()+scale_x_discrete(limits= trts.name)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust=0))
#subset(taxa_abundace_p_i, Phylum %in% topp_i)   ### output 600 * 450
#table_p_i
#write.csv(table_p_i,file = "ileal_phylum_relative_abundance.csv",row.names = FALSE)

# family 
# family >2% and have been identified in Ileal
data_i_m2_2<-tax_glom(data_i_m1,"Family")
data_i_ra_2<-transform_sample_counts(data_i_m2_2,function(x){100*x/sum(x)})
taxa_abundace_f_i<- data.frame(cbind(t(otu_table(data_i_ra_2)),tax_table(data_i_ra_2)))
max_abun<-apply(taxa_abundace_f_i[,1:6],1,max)
topf_i<-taxa_abundace_f_i$Family[which(max_abun>1)]   #>1%
topf_i<-topf_i[topf_i!="f__"]   #still not show unidentified family
table_f_i<-subset(taxa_abundace_f_i,max_abun>1)
# sort table 
table_f_i<-table_f_i[order(table_f_i$Order),]
table_f_i<-table_f_i[order(table_f_i$Class),]
table_f_i<-table_f_i[order(table_f_i$Phylum),]

data_f_i<-prune_taxa(rownames(table_f_i), data_i_ra_2)
#table_f_i<-subset(taxa_abundace_f_i, Family %in% topf_i)
plot_bar(data_f_i,fill="Family")+xlab("Treatment")+ylab("Percentage of Sequences")+
  ylim(0,100)+scale_fill_igv()+scale_x_discrete(limits= trts.name)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust=0))
# taxonomy table
#table_f_i
#write.csv(table_f_i,file = "ileal_family_relative_abundance.csv",row.names = FALSE)

# genus
# genus >2% and have been identified in Ileal
data_i_m2_3<-tax_glom(data_i_m1,"Genus")
data_i_ra_3<-transform_sample_counts(data_i_m2_3,function(x){100*x/sum(x)})
taxa_abundace_g_i<- data.frame(cbind(t(otu_table(data_i_ra_3)),tax_table(data_i_ra_3)))
max_abun<-apply(taxa_abundace_g_i[,1:6],1,max)
topg_i<-taxa_abundace_g_i$Genus[which(max_abun>1)]   #>1%
topg_i<-topg_i[topg_i!="g__"]   #still not show unidentified family
table_g_i<-subset(taxa_abundace_g_i,max_abun>1)
# sort table 
table_g_i<-table_g_i[order(table_g_i$Family),]
table_g_i<-table_g_i[order(table_g_i$Order),]
table_g_i<-table_g_i[order(table_g_i$Class),]
table_g_i<-table_g_i[order(table_g_i$Phylum),]

data_g_i<-prune_taxa(rownames(table_g_i), data_i_ra_3)
#table_g_i<-subset(taxa_abundace_g_i, Genus %in% topg_i)
plot_bar(data_g_i,fill="Genus")+xlab("Treatment")+ylab("Percentage of Sequences")+ylim(0,100)+
  scale_fill_igv()+scale_x_discrete(limits= trts.name)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust=0))
# taxonomy table
#table_g_i
#write.csv(table_g_i,file = "ileal_genus_relative_abundance.csv",row.names = FALSE)


#-------------------------------------------------------------
# Colonic 

# Colonic
data_c_m1<-merge_samples(data_c,"Trt")
trts.name<-c("NC","PC","SF-","SF+","IF-","IF+")

# Phylum
# Phylum >0.5% in Colonic
data_c_m2<-tax_glom(data_c_m1,"Phylum")
data_c_ra<-transform_sample_counts(data_c_m2,function(x){100*x/sum(x)})
taxa_abundace_p_c<- data.frame(cbind(t(otu_table(data_c_ra)),tax_table(data_c_ra)))
max_abun<-apply(taxa_abundace_p_c[,1:6],1,max)
#topp_c<-taxa_abundace_p_c$Phylum[which(max_abun>0.5)]  #>0.5%
table_p_c<-subset(taxa_abundace_p_c,max_abun>0.5)
data_p_c<-prune_taxa(rownames(table_p_c), data_c_ra)
plot_bar(data_p_c,fill="Phylum")+xlab("Treatment")+ylab("Percentage of Sequences")+
  ylim(0,100)+scale_fill_igv()+scale_x_discrete(limits= trts.name)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust=0))
#subset(taxa_abundace_p_c, Phylum %in% topp_c)
#table_p_c
#write.csv(table_p_c,file = "colonic_phylum_relative_abundance.csv",row.names = FALSE)

# Family
# Famlily >1% and have been identified in Colonic
data_c_m2_2<-tax_glom(data_c_m1,"Family")
data_c_ra_2<-transform_sample_counts(data_c_m2_2,function(x){100*x/sum(x)})
taxa_abundance_f_c<- data.frame(cbind(t(otu_table(data_c_ra_2)),
                                      tax_table(data_c_ra_2)))
max_abun<-apply(taxa_abundance_f_c[,1:6],1,max)
topf_c<-taxa_abundance_f_c$Family[which(max_abun>1)]   #>1%
topf_c<-topf_c[topf_c!="f__"]   #still not show unidentified family
#data_f_c<-prune_taxa(rownames(subset(taxa_abundace_f_c, which(max_abun>2))), data_c_ra_2)
table_f_c<-subset(taxa_abundance_f_c,max_abun>1)
# sort table 
table_f_c<-table_f_c[order(table_f_c$Order),]
table_f_c<-table_f_c[order(table_f_c$Class),]
table_f_c<-table_f_c[order(table_f_c$Phylum),]

data_f_c<-prune_taxa(rownames(subset(table_f_c, Family %in% topf_c)), data_c_ra_2)
#table_f_c<-subset(taxa_abundace_f_c, Family %in% topf_c)
plot_bar(data_f_c,fill="Family")+xlab("Treatment")+ylab("Percentage of Sequences")+ylim(0,100)+
  scale_fill_igv()+scale_x_discrete(limits= trts.name)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust=0))
# taxonomy table
# subset(table_f_c, Family %in% topf_c)
#write.csv(subset(table_f_c, Family %in% topf_c),file = "colonic_family_relative_abundance.csv",row.names = FALSE)

# Genus
# genus >2% and have been identified in Colonic
data_c_m2_3<-tax_glom(data_c_m1,"Genus")
data_c_ra_3<-transform_sample_counts(data_c_m2_3,function(x){100*x/sum(x)})
taxa_abundace_g_c<- data.frame(cbind(t(otu_table(data_c_ra_3)),
                                     tax_table(data_c_ra_3)))
max_abun<-apply(taxa_abundace_g_c[,1:6],1,max)
topg_c<-taxa_abundace_g_c$Genus[which(max_abun>1)]   #>1%
topg_c<-topg_c[topg_c!="g__"]   #still not show unidentified family
#data_g_c<-prune_taxa(rownames(subset(taxa_abundace_g_c, Genus %in% topg_c)), data_c_ra_3)
#table_g_c<-subset(taxa_abundace_g_c, Genus %in% topg_c)
table_g_c<-subset(taxa_abundace_g_c, max_abun>1)
table_g_c<-subset(table_g_c, Genus %in% topg_c)
# sort table 
table_g_c<-table_g_c[order(table_g_c$Family),]
table_g_c<-table_g_c[order(table_g_c$Order),]
table_g_c<-table_g_c[order(table_g_c$Class),]
table_g_c<-table_g_c[order(table_g_c$Phylum),]

data_g_c<-prune_taxa(rownames(subset(table_g_c, Genus %in% topg_c)), data_c_ra_3)
plot_bar(data_g_c,fill="Genus")+xlab("Treatment")+ylab("Percentage of Sequences")+
  ylim(0,100)+scale_fill_igv()+scale_x_discrete(limits= trts.name)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust=0))
# taxonomy table
#table_g_c
#subset(table_g_c, Genus %in% topg_c)
#write.csv(subset(table_g_c, Genus %in% topg_c),file = "colonic_genus_relative_abundance.csv",row.names = FALSE)

###### ---------------------------------------------------
## heatmap of relative abundance at diffferent level in ileal 


### Family
data_i_map <- tax_glom(data_i, "Family")
data_i_map<-transform_sample_counts(data_i_map,function(x) 100 * x/sum(x))
GP2 = subset_taxa(data_i_map, Family %in% topf_i)
p<-plot_heatmap(GP2, "PCoA", "unifrac", "Trt", "Family",weighted = TRUE,sample.order="Trt")

png(filename = "heatmap_ileal_family.png",
    width = 3000, height = 1200, res = 300)

p

dev.off()


### Genus
data_i_map <- tax_glom(data_i, "Genus")
data_i_map<-transform_sample_counts(data_i_map,function(x) 100 * x/sum(x))
GP2 = subset_taxa(data_i_map, Genus %in% topg_i)
p<-plot_heatmap(GP2, "PCoA", "unifrac", "Trt", "Genus",weighted = TRUE,sample.order="Trt")

png(filename = "heatmap_ileal_genus.png",
    width = 3000, height = 1200, res = 300)

p

dev.off()

### output relative abundance of ileal taxa at phylum, family and genus level


### Family
data_i_map <- tax_glom(data_i, "Family")
data_i_map<-transform_sample_counts(data_i_map,function(x) 100 * x/sum(x))
GP2 = subset_taxa(data_i_map, Family %in% topf_i)
psmelt(GP2)->ileal_famlily_rela_abun

### Genus
data_i_map <- tax_glom(data_i, "Genus")
data_i_map<-transform_sample_counts(data_i_map,function(x) 100 * x/sum(x))
GP2 = subset_taxa(data_i_map, Genus %in% topg_i)
psmelt(GP2)->ileal_genus_rela_abun




