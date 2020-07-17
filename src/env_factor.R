### this script tries to link VFA and genus

### require phyloseq object of 
#### data_i_unifrac
#### data_c_unifrac
#### data_i
#### data_c

## set the address of working space 
library(phyloseq)
library(vegan)
library(mixOmics)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")

### import environment factor
Ileal_genes <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/Ileal_genes.csv")
colon_vfa <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/colon_vfa.csv")


### for ileal 
###------------------------------------------------
## mixOmics

rownames(Ileal_genes)<-Ileal_genes[,1]
Ileal_genes<-Ileal_genes[,1:17]
rownames(colon_vfa)<-colon_vfa[,1]
colon_vfa<-colon_vfa[,c(1,10:15)]

### prepare sample data 
sample.data<-sample_data(data_i)[,c("Pen_ID","Trt")]
sample.data$names<-rownames(sample.data)
sample.data<-as.data.frame(cbind(sample.data@.Data[[1]],sample.data@.Data[[3]],as.character(sample.data@.Data[[2]])))

##### dataset combination vfa and sample 
vfa.sample<-merge(colon_vfa,sample.data,by.x = "Pen",by.y = "V1")
rownames(vfa.sample)<-vfa.sample$V2
vfa.sample<-vfa.sample[order(vfa.sample$V2),]
vfa.sample.data<-vfa.sample[,2:7]


###ileal

### prepare otu data
#data_i.rela<-transform_sample_counts(data_i,function(x) 100 * x/sum(x))
#ileal.otu<-t(otu_table(data_i.rela))
data_i.cca<-tax_glom(data_i_unifrac, "Genus")
#data_i.cca<-data_i_unifrac
ileal.otu<-t(as.matrix(otu_table(data_i.cca)@.Data))
ileal.data<-data.frame(cbind(otu_table(data_i.cca),
                             tax_table(data_i.cca)))
names<-names(sort(taxa_sums(data_i.cca), TRUE))
#rownames(ileal.data)<-ileal.data$Genus
ileal.data.sub30<-ileal.data[names,]


#### CCA and rCCA
X<-as.matrix(ileal.otu[,names]) # (only plot top 30 genus)
#colnames(X)<-ileal.data.sub30$Genus
Y<-as.matrix(vfa.sample.data)
dim(X)
dim(Y)
head(cbind(rownames(X), rownames(Y)))

png(filename = "cor_ileal_top30_genus_vs_colon_vfa.png",
    width = 2500, height = 2500, res = 300)
colnames(X)<-ileal.data.sub30$Genus
imgCor(X[,1:30],Y)

dev.off()

imgCor(X[,1:30],Y,"separate")

### shrinkage
colnames(X)<-rownames(ileal.data.sub30)
il.vfa.shrink<-rcc(X, Y, ncomp = 3, method = 'shrinkage')
plot(il.vfa.shrink, scree.type = "barplot")

### tune cv

grid1 <- seq(0.05, 0.2, length = 5) 
grid2 <- seq(0.05, 0.2, length = 5)

cv <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo")

#### rCCA
il.vfa.rcc <- rcc(X,Y, ncomp = 3,  lambda1 = cv$opt.lambda1, 
                  lambda2 = cv$opt.lambda2)

plotIndiv(il.vfa.rcc, comp = 1:2, ind.names = vfa.sample$V3,
          group = vfa.sample$V3, rep.space = "XY-variate",
          legend = TRUE, title = 'rCCA XY-space')

plotIndiv(il.vfa.rcc, comp = 1:2,  ind.names = vfa.sample$V3,
          group = vfa.sample$V3, 
          legend = TRUE, title = 'each subspace')

plotVar(il.vfa.rcc, comp = 1:2, cutoff = 0.5, var.names = c(TRUE, TRUE),
        cex = c(4, 4), title = 'rCCA comp 1 - 2')

png(filename = "cca_ileal_genus_vs_colon_vfa.png",
    width = 3000, height = 4000, res = 300)

cim(il.vfa.rcc, comp = 1:3, margins = c(5, 6))

dev.off()

###------------------------------------------------------
## cca or rda by vegan  
library(vegan)

#X<-log(X+1)
decorana(X)

sp0<-cca(as.data.frame(X) ~ ., as.data.frame(Y))
sp0

plot(sp0)

# rda plot constrained on vfa
plot_ordination(data_i_unifrac, sp0, type="samples",color="Treatment", label="Treatment")

# rda plot constrianed on treatment
ileal.rda<-ordinate(data_i_unifrac ~ Treatment, "RDA")
plot_ordination(data_i_unifrac, ileal.rda, type="samples",color="Treatment", label="Treatment")

###---------------------------------------------------------------------------
## for colon

####----------------------------------------------------------------
### colon

### import environment factor
#Ileal_genes <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/Ileal_genes.csv")
#colon_vfa <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/colon_vfa.csv")

#rownames(Ileal_genes)<-Ileal_genes[,1]
#Ileal_genes<-Ileal_genes[,1:17]
#rownames(colon_vfa)<-colon_vfa[,1]
#colon_vfa<-colon_vfa[,c(1,10:15)]

### prepare sample data 
sample.data<-sample_data(data_c)[,c("Pen_ID","Trt")]
sample.data$names<-rownames(sample.data)
sample.data<-as.data.frame(cbind(sample.data@.Data[[1]],sample.data@.Data[[3]],as.character(sample.data@.Data[[2]])))

##### dataset combination vfa and sample 
vfa.sample<-merge(colon_vfa,sample.data,by.x = "Pen",by.y = "V1")
rownames(vfa.sample)<-vfa.sample$V2
vfa.sample<-vfa.sample[order(vfa.sample$V2),]
vfa.sample.data<-vfa.sample[,2:7]

### prepare otu data
#data_i.rela<-transform_sample_counts(data_i,function(x) 100 * x/sum(x))
#ileal.otu<-t(otu_table(data_i.rela))
data_c.cca<-tax_glom(data_c_unifrac, "Genus")
#data_c.cca<-transform_sample_counts(data_c.cca,function(x) log(x+1))
#data_i.cca<-data_i_unifrac
colon.otu<-t(as.matrix(otu_table(data_c.cca)@.Data))
colon.data<-data.frame(cbind(otu_table(data_c.cca),
                             tax_table(data_c.cca)))
names<-names(sort(taxa_sums(data_c.cca), TRUE))
#rownames(ileal.data)<-ileal.data$Genus
colon.data.sub30<-colon.data[names,]


#### CCA and rCCA
X<-as.matrix(colon.otu[,names]) # (only consider top 30 genus)
colnames(X)<-colon.data.sub30$Genus
Y<-as.matrix(vfa.sample.data)
dim(X)
dim(Y)
head(cbind(rownames(X), rownames(Y)))

png(filename = "cor_colon_top30_genus_vs_colon_vfa.png",
    width = 2500, height = 2500, res = 300)
colnames(X)<-colon.data.sub30$Genus
imgCor(X[,1:30],Y)

dev.off()

colnames(X)<-colon.data.sub30$Genus
imgCor(X[,1:30],Y,"separate")

### shrinkage
colnames(X)<-rownames(colon.data.sub30)
co.vfa.shrink<-rcc(X, Y, ncomp = 3, method = 'shrinkage')
plot(co.vfa.shrink, scree.type = "barplot")

### tune cv

grid1 <- seq(0.05, 0.2, length = 5) 
grid2 <- seq(0.05, 0.2, length = 5)

cv <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo")

#### rCCA
co.vfa.rcc <- rcc(X,Y, ncomp = 3,  lambda1 = cv$opt.lambda1, 
                  lambda2 = cv$opt.lambda2)

plotIndiv(co.vfa.rcc, comp = 1:2, ind.names = vfa.sample$V3,
          group = vfa.sample$V3, rep.space = "XY-variate",
          legend = TRUE, title = 'rCCA XY-space')

plotIndiv(co.vfa.rcc, comp = 1:2,  ind.names = vfa.sample$V3,
          group = vfa.sample$V3, 
          legend = TRUE, title = 'each subspace')

plotVar(co.vfa.rcc, comp = 1:2, cutoff = 0.5, var.names = c(TRUE, TRUE),
        cex = c(4, 4), title = 'rCCA comp 1 - 2')

png(filename = "cca_colon_genus_vs_colon_vfa.png",
    width = 3000, height = 4000, res = 300)

cim(co.vfa.rcc, comp = 1:3, margins = c(5, 6),xlab = "VFA",ylab = "Genus")

dev.off()


###------------------------------------------------------
## cca or rda by vegan  
library(vegan)

#X<-log(X+1)
decorana(X)

sp0<-rda(as.data.frame(X) ~ ., as.data.frame(Y))
sp0

plot(sp0)


### rda plot constrained on vfa 
plot_ordination(data_c_unifrac, sp0, type="samples",color="Treatment", label="Treatment")

### rda on trt
colon.rda<-ordinate(data_c_unifrac ~ Treatment, "RDA")
plot_ordination(data_c_unifrac, colon.rda, type="samples",color="Treatment", label="Treatment")


plot(sp0,choices=c(1,2),display=c("wa","bp"),type="points",xlim=c(-4,1.5),scaling=2)
points(sp0,disp="sites",pch=21,col="red",bg="red",cex=1)
text(sp0,"sites",pos=3,axis.bp=TRUE)



