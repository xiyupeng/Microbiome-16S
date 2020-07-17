### this new script tries to link VFA and genus

### require phyloseq object of 
#data_c
#topg_c   top genus list in colon

## set the address of working space 
library(phyloseq)
library(vegan)
library(corrplot)
setwd("/Users/xiyupeng/Documents/Research/metagenomics/QIIME/working_space")

### import environment factor vfa in colon
colon_vfa <- read.csv("~/Documents/Research/metagenomics/VFA and gene expression data/colon_vfa.csv")
rownames(colon_vfa)<-colon_vfa[,1]
colon_vfa<-colon_vfa[,c(1,10:15)]

### prepare sample data 
sample.data<-sample_data(data_c)[,c("Pen_ID","Trt")]
sample.data$names<-rownames(sample.data)
sample.data<-as.data.frame(cbind(sample.data@.Data[[1]],sample.data@.Data[[3]],as.character(sample.data@.Data[[2]])))

##### dataset combination vfa and sample 
vfa.sample<-merge(colon_vfa,sample.data,by.x = "Pen",by.y = "V1")
rownames(vfa.sample)<-vfa.sample$V2
vfa.sample<-vfa.sample[order(vfa.sample$V2),]
vfa.sample.data<-vfa.sample[,2:7]

#### prepare genus table (select max abundance > 1%)
topg_c<-topg_c[topg_c!="g__uncultured"]
data_c_unifrac<-transform_sample_counts(data_c,function(x) 1E3 * x/sum(x))
data_c.cca<-tax_glom(data_c_unifrac, "Genus")
colon.otu<-t(as.matrix(otu_table(data_c.cca)@.Data))
colon.data<-data.frame(cbind(otu_table(data_c.cca),
                             tax_table(data_c.cca)))
colon.data.select<-colon.otu[,colon.data$Genus %in% topg_c]

##X and Y matrix
X<-as.matrix(colon.data.select)
colnames(X)<-colon.data$Genus[colon.data$Genus %in% topg_c]
colnames(X)<-colon.data.sub30$Genus
Y<-as.matrix(vfa.sample.data)
dim(X)
dim(Y)
head(cbind(rownames(X), rownames(Y)))

### plot correlation map 
res <- cor(X,Y)
col<- colorRampPalette(c("blue", "white", "red"))(10)
#heatmap(x = res, col =col,margins = c(8,8),scale = "column",xlab = "VFA")
corrplot(t(res),col = col)

write.csv(res,file = "cc_colon_genus_vs_vfa.csv")

res.p <- t(cor.mtest(cbind(X,Y),conf.level = .95)$p[36:41,1:35])
colnames(res.p)<-colnames(res)
rownames(res.p)<-rownames(res)

write.csv(res.p,file = "cc_pvalue_colon_genus_vs_vfa.csv")

## p value adjust 
res.adjp<-p.adjust(array(res.p), method = "BH")
dim(res.adjp)<-dim(res.p)
colnames(res.adjp)<-colnames(res)
rownames(res.adjp)<-rownames(res)

write.csv(res.adjp,file = "cc_adjp_colon_genus_vs_vfa.csv")

### RDA analysis----------------------------------------------------------------
library(vegan)

#X<-log(X+1)
decorana(X)

sp0<-rda(as.data.frame(X) ~ ., as.data.frame(Y))
sp0

## plot1
plot(sp0,type = "points")

## plot2
plot(sp0,choices=c(1,2),display=c("wa","bp","sites"),type="points",xlim=c(-4,1.5),scaling=2)
points(sp0,disp="sites",pch=21,col="red",bg="red",cex=1)

##plot3
plot(sp0, type="n", scaling="sites")
text(sp0, dis="cn",scaling = "sites")
points(sp0, pch=21, col="red", bg="yellow", cex=0.8,scaling = "sites")
text(sp0, "species", col="blue", cex=0.8,scaling = "sites")


## plot4 
### data prepare
new<-sp0$CCA
new
samples<-data.frame(sample=row.names(new$u),RDA1=new$u[,1],RDA2=new$u[,2])
samples
species<-data.frame(spece=row.names(new$v),RDA1=new$v[,1],RDA2=new$v[,2])
species
envi<-data.frame(en=row.names(new$biplot),RDA1=new$biplot[,1],RDA2=new$biplot[,2])
envi
line_x =c(0,envi[1,2],0,envi[2,2],0,envi[3,2],0,envi[4,2],0,envi[5,2],0,envi[6,2])
line_x
line_y =c(0,envi[1,3],0,envi[2,3],0,envi[3,3],0,envi[4,3],0,envi[5,3],0,envi[6,3])
line_y
line_g =c("Acetic","Propionic","Isobutyric","Butyric","Isovaleric","Valeric")
line_g
line_data =data.frame(x=line_x,y=line_y,group=line_g)

line_data

# plot
library(ggplot2)
ggplot(data=samples,aes(RDA1,RDA2)) +geom_point(size=1,shape = 3) +
  geom_point(data=species,aes(col=spece)) +
  geom_text(data=envi,aes(label=en),color="blue") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  theme_bw() + theme(panel.grid=element_blank())
ggsave("RDA2.PDF")

