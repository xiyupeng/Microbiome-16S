### this script is going to draw graph of beta diversity 
### and do some statistics.(ANOSIM AND PERMANOVA)

### require phyloseq object of 
#### data_0

## set the address of working space 
library(phyloseq)
library(vegan)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")


#-----------------------------------------------------------------------
# Prevalence filtering

# if we only consider colonic and ileal samples 
#data0<-data_ci

##########################################
prev0 <-apply(X = otu_table(data0),
              MARGIN = ifelse(taxa_are_rows(data0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevalenceThreshold <- 0.05 * nsamples(data0)  # threshold as 5% of total samples
prevalenceThreshold
# filter Phylum that have OTUs smaller than 6
prevdf <- data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(data0),
                     tax_table(data0))
keepPhyla <- table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
#top5phyla = sort(keepPhyla, TRUE)[1:5]
prevdf1 <- subset(prevdf, Phylum %in% names(keepPhyla))
# filter taxa that that were rarely seen in the samples.
data1 <- prune_taxa((prev0 > prevalenceThreshold), data0)
data1
# filter taxa with rare phylum
data2 <- subset_taxa(data1, Phylum %in% names(keepPhyla))
data2
ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 1, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") + facet_wrap(~Phylum)

#-----------------------------------------------------------------
# Beta Diversity 

library(vegan)
set.seed(123)

# based on weighted unifrac 
#data_unifrac = rarefy_even_depth(data2)
data_unifrac<-transform_sample_counts(data2, function(x) 1E6 * x/sum(x))
# change the label of data_unifrac
names(sample_data(data_unifrac))[9]<-"Sample"
names(sample_data(data_unifrac))[6]<-"Treatment"
ordu_w <- ordinate(data_unifrac,"PCoA","unifrac",weighted=TRUE)
#plot_ordination(data_unifrac,ordu_w,color = "Sample",shape="Treatment",label = "Treatment")+xlab("PCoA1 [44.4%]")+ylab("PCoA2 [19%]")
plot_ordination(data_unifrac,ordu_w,color = "Sample",shape="Treatment")+xlab("PCoA1 [43.9%]")+ylab("PCoA2 [14.7%]")#+facet_wrap(~Treatment)#

#### new ------------------------------------------------
plot_ordination(data_unifrac,ordu_w,color = "Sample",shape="Treatment")+xlab("PCoA1 [43.9%]")+ylab("PCoA2 [14.7%]")+facet_wrap(~Treatment)+stat_ellipse(aes(fill=Sample), alpha=.2,type='t', geom="polygon")
# size 800 * 500
#plot_ordination(data_unifrac,ordu_w,color = "Trt",shape="Sample_type")
### 
ordu_w <- ordinate(data_unifrac,"PCoA","unifrac",weighted=TRUE)
plot_ordination(data_unifrac,ordu_w,color = "Sample",shape="Gender")+xlab("PCoA1 [43.9%]")+ylab("PCoA2 [14.7%]")+facet_wrap(~Gender)+stat_ellipse(aes(fill=Sample), alpha=.2,type='t', geom="polygon")

# anosim test for weighted unifrac distance of microbio comunity
distance_w = distance(data_unifrac, "wunifrac")
# Sample_Type
sample_group = get_variable(data_unifrac, "Sample")
sample_ano_w = anosim(distance_w, sample_group)
sample_ano_w
#pairwise comparison
pairwise.adonis(distance_w,sample_group,'BH')

## sex 
gender_group = get_variable(data_unifrac, "Gender")
gender_ano_w = anosim(distance_w, gender_group)
gender_ano_w
#pairwise comparison
pairwise.adonis(distance_w,gender_group,'BH')


# Treatment
trt_group = get_variable(data_unifrac, "Treatment")
trt_ano_w = anosim(distance_w,trt_group)
trt_ano_w
# pairwise comparison
pairwise.adonis(distance_w,trt_group,'BH')

# combine treatment and sample_type
group <- factor(paste0(sample_group, "_", trt_group))
# pairwise comparison (only show beta diversity in the Ileal)
pairwise.adonis(distance_w,group,'BH')

# based on unweighted unifrac distance
ordu_u = ordinate(data_unifrac,"PCoA","unifrac",weighted=FALSE)
plot_ordination(data_unifrac,ordu_u,color = "Sample",shape="Treatment")+xlab("PCoA1 [39.8%]")+ylab("PCoA2 [4.3%]")#+facet_wrap(~Treatment)
plot_ordination(data_unifrac,ordu_u,color = "Sample",shape="Treatment")+xlab("PCoA1 [39.8%]")+ylab("PCoA2 [4.3%]")+facet_wrap(~Treatment)+stat_ellipse(aes(fill=Sample), alpha=.2,type='t', geom="polygon")
## sex 
plot_ordination(data_unifrac,ordu_u,color = "Sample",shape="Gender")+xlab("PCoA1 [39.8%]")+ylab("PCoA2 [4.3%]")+facet_wrap(~Gender)+stat_ellipse(aes(fill=Sample), alpha=.2,type='t', geom="polygon")


# size 800 * 500
#plot_ordination(data_unifrac,ordu_u,color = "Trt",shape="Sample_type")

# anosim test for unweighted unifrac distance of microbio comunity
distance_uw = distance(data_unifrac, "unifrac")
# Sample_Type
sample_ano_uw = anosim(distance_uw, sample_group)
sample_ano_uw
# pairwise comparison
pairwise.adonis(distance_uw,sample_group,'BH')

gender_group = get_variable(data_unifrac, "Gender")
gender_ano_uw = anosim(distance_uw, gender_group)
gender_ano_uw
#pairwise comparison
pairwise.adonis(distance_w,gender_group,'BH')

# Treatment
trt_ano_uw = anosim(distance_uw, trt_group)
trt_ano_uw
# pairwise comparison
pairwise.adonis(distance_uw,trt_group,'BH')

# combine Sample_type and Treatment
# pairwise comparison (only show beta diversity in the ileal)
pairwise.adonis(distance_uw,group,'BH')

# try to do PCoA separately

#-----------------------------------------------------------------------------------
# Feces ####################

#data_f_unifrac= rarefy_even_depth(data_f)
data_f_2<-subset_samples(data2, Sample_type=="Feces")
data_f_2 <- prune_taxa(taxa_sums(data_f_2) > 0, data_f_2)

# weighted unifrac 
data_f_unifrac<-transform_sample_counts(data_f_2,function(x) 1E6 * x/sum(x))
ordu_f_w = ordinate(data_f_unifrac,"PCoA","unifrac",weighted=TRUE)
plot_ordination(data_f_unifrac,ordu_f_w,color = "Diet",label = "Diet")+xlab("PCoA1 [21.6%]")+ylab("PCoA2 [15.2%]")#+facet_wrap(~Diet)

# anosim test for weighted unifrac distance of microbio comunity in feces
distance_w_feces = distance(data_f_unifrac, "wunifrac")
diet_group_feces = get_variable(data_f_unifrac, "Diet")
diet_ano_w_feces = anosim(distance_w_feces,diet_group_feces)
diet_ano_w_feces
# pairwise comparison 
pairwise.adonis(distance_w_feces,diet_group_feces,'BH')
### gender 
gender_group_feces = get_variable(data_f_unifrac, "Gender")
distance_w_feces = distance(data_f_unifrac, "wunifrac")
pairwise.adonis(distance_w_feces,gender_group_feces,'BH')

# unweighted unifrac 
ordu_f_uw = ordinate(data_f_unifrac,"PCoA","unifrac",weighted=FALSE)
plot_ordination(data_f_unifrac,ordu_f_uw,color = "Diet",label = "Diet")+xlab("PCoA1 [13.5%]")+ylab("PCoA2 [8.5%]")#+facet_wrap(~Diet)

# anosim test for weighted unifrac distance of microbio comunity in feces
distance_uw_feces = distance(data_f_unifrac, "wunifrac")
diet_group_feces = get_variable(data_f_unifrac, "Diet")
diet_ano_uw_feces = anosim(distance_uw_feces,diet_group_feces)
diet_ano_uw_feces
# pairwise comparison 
pairwise.adonis(distance_uw_feces,diet_group_feces,'BH')
### gender 
gender_group_feces = get_variable(data_f_unifrac, "Gender")
distance_uw_feces = distance(data_f_unifrac, "unifrac")
pairwise.adonis(distance_uw_feces,gender_group_feces,'BH')



#----------------------------------------------------------------------------------
# Ileal ######################

data_i_2<-subset_samples(data2, Sample_type=="Ileal digesta")
data_i_2<-prune_taxa(taxa_sums(data_i_2) > 0, data_i_2)
#data_i_unifrac= rarefy_even_depth(data_i_2)

# weighted unifrac
data_i_unifrac<-transform_sample_counts(data_i_2,function(x) 1E6 * x/sum(x))
names(sample_data(data_i_unifrac))[6]<-"Treatment"

ordu_i_w = ordinate(data_i_unifrac,"PCoA","unifrac",weighted=TRUE)
plot_ordination(data_i_unifrac,ordu_i_w,color = "Treatment",label = "Treatment")+xlab("PCoA1 [57.2%]")+ylab("PCoA2 [12.3%]")
####new ---------------------------
plot_ordination(data_i_unifrac,ordu_i_w,color = "Treatment",label = "Treatment")+xlab("PCoA1 [57.2%]")+ylab("PCoA2 [12.3%]")+stat_ellipse(aes(fill=Treatment), alpha=0.2,type='t', geom="polygon")
#---------------------------------
#plot_ordination(data_i_unifrac,ordu_i_w,shape = "Treatment")+facet_wrap(~Treatment)+xlab("PCoA1 [57.2%]")+ylab("PCoA2 [12.3%]")

# anosim test 
distance_w_ileal = distance(data_i_unifrac, "wunifrac")
trt_group_ileal = get_variable(data_i_unifrac, "Treatment")
trt_ano_w_ileal = anosim(distance_w_ileal,trt_group_ileal)
trt_ano_w_ileal
# pairwise comparison 
pairwise.adonis(distance_w_ileal,trt_group_ileal,'BH')
### gender 
gender_group_ileal = get_variable(data_i_unifrac, "Gender")
distance_w_ileal = distance(data_i_unifrac, "wunifrac")
pairwise.adonis(distance_w_ileal,gender_group_ileal,'BH')


# Use Bray-Curtis dissimilarity
#distance_bc_ileal = distance(data_i_unifrac,"bray")
#trt_group_ileal = get_variable(data_i_unifrac, "Trt")
#trt_ano_bc_ileal = anosim(distance_bc_ileal,trt_group_ileal)
#trt_ano_bc_ileal
# pairwise comparison 
#pairwise.adonis(distance_bc_ileal,trt_group_ileal,'BH')

# unweighted unifrac
ordu_i_uw = ordinate(data_i_unifrac,"PCoA","unifrac",weighted=FALSE)
plot_ordination(data_i_unifrac,ordu_i_uw,color = "Treatment",label = "Treatment")+xlab("PCoA1 [13.6%]")+ylab("PCoA2 [9.1%]")
plot_ordination(data_i_unifrac,ordu_i_uw,shape = "Treatment")+facet_wrap(~Treatment)+xlab("PCoA1 [13.6%]")+ylab("PCoA2 [9.1%]")

# anosim test
distance_uw_ileal = distance(data_i_unifrac, "unifrac")
trt_ano_uw_ileal = anosim(distance_uw_ileal,trt_group_ileal)
trt_ano_uw_ileal
# pairwise comparison 
pairwise.adonis(distance_uw_ileal,trt_group_ileal,'BH')
### gender 
gender_group_ileal = get_variable(data_i_unifrac, "Gender")
distance_uw_ileal = distance(data_i_unifrac, "unifrac")
pairwise.adonis(distance_uw_ileal,gender_group_ileal,'BH')


#-----------------------------------------------------------------------
# Colonic
#data_c_unifrac= rarefy_even_depth(data_c)
data_c_2<-subset_samples(data2, Sample_type=="Colonic digesta")
data_c_2 <- prune_taxa(taxa_sums(data_c_2) > 0, data_c_2)

data_c_unifrac<-transform_sample_counts(data_c_2,function(x) 1E6 * x/sum(x))
names(sample_data(data_c_unifrac))[6]<-"Treatment"

ordu_c_w = ordinate(data_c_unifrac,"PCoA","unifrac",weighted=TRUE)
plot_ordination(data_c_unifrac,ordu_c_w,color = "Treatment",label = "Treatment")+xlab("PCoA1 [18%]")+ylab("PCoA2 [15.1%]")
plot_ordination(data_c_unifrac,ordu_c_w,shape = "Treatment")+facet_wrap(~Treatment)+xlab("PCoA1 [18%]")+ylab("PCoA2 [15.1%]")

# anosim test for weighted unifrac distance of microbio comunity in colonic
distance_w_colonic = distance(data_c_unifrac, "wunifrac")
trt_group_colonic = get_variable(data_c_unifrac, "Treatment")
trt_ano_w_colonic = anosim(distance_w_colonic,trt_group_colonic)
trt_ano_w_colonic
# pairwise comparison
pairwise.adonis(distance_w_colonic,trt_group_colonic,'BH')
### gender 
gender_group_colonic = get_variable(data_c_unifrac, "Gender")
distance_w_colonic = distance(data_c_unifrac, "wunifrac")
pairwise.adonis(distance_w_colonic,gender_group_colonic,'BH')

# unweighted unifrac 
ordu_c_uw = ordinate(data_c_unifrac,"PCoA","unifrac",weighted=FALSE)
plot_ordination(data_c_unifrac,ordu_c_uw,color = "Treatment",label = "Treatment")+xlab("PCoA1 [16.1%]")+ylab("PCoA2 [8.8%]")
plot_ordination(data_c_unifrac,ordu_c_uw,shape = "Treatment")+facet_wrap(~Treatment)+xlab("PCoA1 [16.1%]")+ylab("PCoA2 [8.8%]")

# anosim test for weighted unifrac distance of microbio comunity in colonic
distance_uw_colonic = distance(data_c_unifrac, "unifrac")
trt_group_colonic = get_variable(data_c_unifrac, "Treatment")
trt_ano_uw_colonic = anosim(distance_uw_colonic,trt_group_colonic)
trt_ano_uw_colonic
# pairwise comparison
pairwise.adonis(distance_uw_colonic,trt_group_colonic,'BH')
### gender 
gender_group_colonic = get_variable(data_c_unifrac, "Gender")
distance_uw_colonic = distance(data_c_unifrac, "unifrac")
pairwise.adonis(distance_uw_colonic,gender_group_colonic,'BH')

