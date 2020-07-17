### this script is going to draw graph of alpha diversity 
### and do some statistics.

### require phyloseq object of 
#### data_i
#### data_c
#### data_f

## set the address of working space
library(phyloseq)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")

# Alpha diversity

# Feces

### new way  550 * 350
diets.name<-c("Basal","SF-","SF+","IF-","IF+")
plot_richness(data_f,x = "Diet",measures=c("Chao1", "Shannon","InvSimpson"))+geom_boxplot()+xlab("Diet")+scale_x_discrete(limits= diets.name)
## difference between sex
plot_richness(data_f,x = "Gender",measures=c("Chao1", "Shannon","InvSimpson"))+geom_boxplot()+xlab("Gender")


### old way 
richness.feces<-estimate_richness(data_f, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
diet_index_f<-sample_data(data_f)$Diet
richness.feces<-cbind(richness.feces,diet_index_f)

plot_f_c<-ggplot(richness.feces,aes(diet_index_f,Chao1))+geom_boxplot()+xlab("Diets")+
  ylim(0,1000)+labs(title ="Feces_Chao1")+theme(axis.text.x = element_text(angle = 45))

plot_f_s<-ggplot(richness.feces,aes(diet_index_f,Shannon))+geom_boxplot()+xlab("Diets")+
  ylim(0,4.5)+labs(title="Feces_Shannon")+theme(axis.text.x = element_text(angle = 45))

plot_f_i<-ggplot(richness.feces,aes(diet_index_f,InvSimpson))+geom_boxplot()+xlab("Diets")+
  ylim(0,40)+labs(title="Feces_InvSimpson")+theme(axis.text.x = element_text(angle = 45))

# Ileal 

### new way  550 * 350
plot_richness(data_i,x = "Trt",measures=c("Chao1", "Shannon","InvSimpson"))+geom_boxplot()+xlab("Treatment")
plot_richness(data_i,x = "Gender",measures=c("Chao1", "Shannon","InvSimpson"))+geom_boxplot()+xlab("Gender")


#### old way 
richness.ileal<-estimate_richness(data_i, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
trt_index_i<-sample_data(data_i)$Trt
richness.ileal<-cbind(richness.ileal,trt_index_i)
plot_i_c<-ggplot(richness.ileal,aes(trt_index_i,Chao1))+geom_boxplot()+xlab("Treatment")+
  ylim(0,1000)+labs(title="Ileal")+theme(axis.text.x = element_text(angle = 45))
plot_i_s<-ggplot(richness.ileal,aes(trt_index_i,Shannon))+geom_boxplot()+xlab("Treatment")+
  ylim(0,4.5)+labs(title="Ileal")+theme(axis.text.x = element_text(angle = 45))
plot_i_i<-ggplot(richness.ileal,aes(trt_index_i,InvSimpson))+geom_boxplot()+xlab("Treatment")+
  ylim(0,40)+labs(title="Ileal")+theme(axis.text.x = element_text(angle = 45))

# Colonic

#### new way  550 * 350
plot_richness(data_c,x = "Trt",measures=c("Chao1", "Shannon","InvSimpson"))+geom_boxplot()+xlab("Treatment")
plot_richness(data_c,x = "Gender",measures=c("Chao1", "Shannon","InvSimpson"))+geom_boxplot()+xlab("Gender")


### old way 
richness.colonic<-estimate_richness(data_c, measures=c("Observed",
                                                       "InvSimpson", "Shannon", "Chao1"))
trt_index_c<-sample_data(data_c)$Trt
richness.colonic<-cbind(richness.colonic,trt_index_c)
plot_c_c<-ggplot(richness.colonic,aes(trt_index_c,Chao1))+geom_boxplot()+xlab("Treatment")+
  ylim(0,1000)+labs(title="Colonic")+theme(axis.text.x = element_text(angle = 45))
plot_c_s<-ggplot(richness.colonic,aes(trt_index_c,Shannon))+geom_boxplot()+xlab("Treatment")+
  ylim(0,4.5)+labs(title="Colonic")+theme(axis.text.x = element_text(angle = 45))
plot_c_i<-ggplot(richness.colonic,aes(trt_index_c,InvSimpson))+geom_boxplot()+xlab("Treatment")+
  ylim(0,40)+labs(title="Colonic")+theme(axis.text.x = element_text(angle = 45))
#ggplot(richness.colonic,aes(diet_index_c,InvSimpson))+geom_boxplot()+xlab("Treatment")

#grid.arrange(nrow=3,plot_f_c,plot_i_c,plot_c_c,plot_f_s,plot_i_s,plot_c_s,plot_f_i,plot_i_i,plot_c_i)
grid.arrange(nrow=3,plot_i_c,plot_c_c,plot_i_s,plot_c_s,plot_i_i,plot_c_i)
# Chao1  // current output 600*400 eps graph
grid.arrange(nrow=1,plot_i_c,plot_c_c)
# Shannon
grid.arrange(nrow=1,plot_i_s,plot_c_s)
# InvSimpson
grid.arrange(nrow=1,plot_i_i,plot_c_i)


# statistics analysis for the alpha diversity
# between sample types
richness.all<-estimate_richness(data0,measures = c("Observed", "InvSimpson", "Shannon", "Chao1"))
# Chao 1
pairwise.t.test(richness.all$Chao1,sample_data(data0)$Sample_type,"BH")
pairwise.t.test(richness.all$Chao1,sample_data(data0)$Gender,"BH") # 0.49

# Shannon
pairwise.t.test(richness.all$Shannon,sample_data(data0)$Sample_type,"BH")
pairwise.t.test(richness.all$Shannon,sample_data(data0)$Gender,"BH") # 0.67

# InvSimpson
pairwise.t.test(richness.all$InvSimpson,sample_data(data0)$Sample_type,"BH")
pairwise.t.test(richness.all$InvSimpson,sample_data(data0)$Gender,"BH") # 0.25


# Between diets in feces
# Chao 1
diet_index_f<-sample_data(data_f)$Diet
lm(richness.feces$Chao1 ~ diet_index_f)->r
summary(r)
pairwise.t.test(richness.feces$Chao1,diet_index_f,"BH")
pairwise.t.test(richness.feces$Chao1,sample_data(data_f)$Gender,"BH")  # 0.11
# Shannon
lm(richness.feces$Shannon ~ diet_index_f)->r
summary(r)
pairwise.t.test(richness.feces$Shannon,diet_index_f,"BH")
pairwise.t.test(richness.feces$Shannon,sample_data(data_f)$Gender,"BH")  # 0.081
# InvSimpson
lm(richness.feces$InvSimpson ~ diet_index_f)->r
summary(r)
pairwise.t.test(richness.feces$InvSimpson,diet_index_f,"BH")
pairwise.t.test(richness.feces$InvSimpson,sample_data(data_f)$Gender,"BH")  # 0.18


# Between treatment in ileal
# Chao 1
lm(richness.ileal$Chao1 ~ trt_index_i)->r
summary(r)
pairwise.t.test(richness.ileal$Chao1,trt_index_i,"BH")
pairwise.t.test(richness.ileal$Chao1,sample_data(data_i)$Gender,"BH")  # 0.4
# Shannon
lm(richness.ileal$Shannon ~ trt_index_i)->r
summary(r)
pairwise.t.test(richness.ileal$Shannon,trt_index_i,"BH")
pairwise.t.test(richness.ileal$Shannon,sample_data(data_i)$Gender,"BH")  # 0.67
# InvSimpson
lm(richness.ileal$InvSimpson ~ trt_index_i)->r
summary(r)
pairwise.t.test(richness.ileal$InvSimpson,trt_index_i,"BH")
pairwise.t.test(richness.ileal$InvSimpson,sample_data(data_i)$Gender,"BH")  # 0.86


# Between treatment in Colonic
# Chao 1
lm(richness.colonic$Chao1 ~ trt_index_c)->r
summary(r)
pairwise.t.test(richness.colonic$Chao1,trt_index_c,"BH")
pairwise.t.test(richness.colonic$Chao1,sample_data(data_c)$Gender,"BH")  # 0.46
# Shannon
lm(richness.colonic$Shannon ~ trt_index_c)->r
summary(r)
pairwise.t.test(richness.colonic$Shannon,trt_index_c,"BH")
pairwise.t.test(richness.colonic$Shannon,sample_data(data_c)$Gender,"BH")  # 0.6
# InvSimpson
lm(richness.colonic$InvSimpson ~ trt_index_c)->r
summary(r)
pairwise.t.test(richness.colonic$InvSimpson,trt_index_c,"BH")
pairwise.t.test(richness.colonic$InvSimpson,sample_data(data_c)$Gender,"BH")  # 0.42

