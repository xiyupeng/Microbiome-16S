pairwise.adonis <- function(x,factors, p.adjust.m ='BH')
{
  co = combn(levels(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(as.matrix(x)[factors %in% c(co[1,elem],co[2,elem]),
                             factors %in% c(co[1,elem],co[2,elem])] ~ 
                  factors[factors %in% c(co[1,elem],co[2,elem])]);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

contra_ana<-function(dds,contrast,alpha,ori_data,topg,level){
  res<-results(dds, contrast =contrast)
  res <- res[order(res$padj, na.last=NA), ]
  sigtab <- res
  sigtab <- res[(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ori_data)[rownames(sigtab), ],"matrix"))
  # only consider genus when abundance > 2% # use topg_i (caculated from previous analysis )
  if(level == "genus")
    return(subset(sigtab, Genus %in% topg))
  if(level == "family")
    return(subset(sigtab, Family %in% topg))
  #return(sigtab)
}

group_mean<-function(diagdds,otulist,counts_norm,group1,group2,group_other){
  
  group1.count<-counts_norm[otulist,which(colData(diagdds)$group == group1)]
  group1.mean<-apply(group1.count,1,mean)
  if (length(group_other)!=0){
    group2.count<-counts_norm[otulist,which(colData(diagdds)$group == group2)]
    group_other.count<-counts_norm[otulist,which(colData(diagdds)$group == group_other)]
    group2.count<-cbind(group2,group_other)
  }else
    group2.count<-counts_norm[otulist,which(colData(diagdds)$group == group2)]
  group2.mean<-apply(group2.count,1,mean)
  abundance<-as.data.frame(cbind(group1.mean,group2.mean))
  return(abundance)
}


kw_test<-function(physeq,padj_method){
  
  label<-sample_data(physeq)$Trt
  
  n<-ncol(otu_table(physeq))

  p.value<-NULL
  
  for (i in 1:n){
    abun_data<-otu_table(physeq)[,i]
    p.value[i]<-kruskal.test(abun_data ~ label)$p.value
  }
  
  p.adjusted = p.adjust(p.value,method= padj_method)
  
  result<-data.frame(colnames(otu_table(physeq)),p.value,p.adjusted)
  
  return(result)
}