---
  title: "CLIC_Analysis"
author: "PTDolan"
date: "8/21/2018" #Happy Birthday, Mom.
output: pdf_document
---
  ```{r setup, include=FALSE}
library(data.table)
knitr::opts_chunk$set(echo = TRUE)
```
## R Markdown
This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>

    Installing GEO functions from bioconductor.
```{r GEO, echo=FALSE}
# source("https://bioconductor.org/biocLite.R")
# biocLite("GEOquery",ask = F)
```

Add GEO tools to library
```{r libGEO}
library(GEOquery)
library(data.table)
library(ggplot2)
library(viridis)
```


Get GEO DB entry. 
```{r}
retrieveGEO<-function(ID){
  retrieved<- GEOquery::getGEOSuppFiles(GEO="GSE97987") 
  
  return(retrieved)
}

generateTable<-function(data,COLUMNS,TABLE){
  N<-ncol(TABLE)
  DT<-data.table(TABLE)
  
  Z<-scale(DT[,3:N])
  Filter<-apply(MARGIN = 1,scale(DT[,3:N]),FUN = function(X){any(abs(X)>1)})
  
  input<-na.omit(DT[Filter,])
  #smacof::mds(dist(input[,3:6]))
  return(input)
}

plotPopStructure<-function(input,COLUMNS,N){
  #Dsample<-dist(t(input[,3:(N)]),'manhattan')
  Dsample<-as.dist(1-cor(input[,3:N]))
  scaleDS<-cmdscale(Dsample)
  popStructure<-data.frame(scaleDS,COLUMNS)
  
  
  for (el in 3:(ncol(popStructure))){
    ggsave(width=10,height=6,paste("GDS",ID,"_",colnames(popStructure)[el],"_plot.pdf",sep = ""),ggplot(popStructure)+geom_point(aes(X1,X2,color=popStructure[,el])))
  }
  
  PC<-princomp(cor=T,na.omit(input[,3:N]))
  Sinput<-cbind(input,PC$scores)
  sM<-melt(Sinput,id.vars = c(1:2,(N+1):ncol(Sinput)))
  MsM<-melt(sM,id.vars = c(1,2,(ncol(sM)-1):(ncol(sM))),variable.name = "Comp")
  
  ggplot(MsM)+geom_bar(aes(reorder(IDENTIFIER,value.1),value.1),stat='identity')+facet_grid(~Comp)
  
  PCinput<-data.frame(scaleDS,scale(PC$loadings[]))
  PCinput<-melt(PCinput,measure.vars = 3:(ncol(PCinput)))
  
  ggsave(width=10,height=6,paste("GDS",ID,"_PCplot.pdf",sep = ""),
         ggplot(PCinput)+geom_point(aes(X1,X2,color=value))+facet_wrap(~variable)+scale_color_viridis())
  
  for (comp in 1:ncol(PC$scores)){
    pdf(paste(ID,"Comp_",comp,"_heatmap.pdf",sep=""),width = 5,height=8)
    gplots::heatmap.2(col = viridis(100),labCol=COLUMNS[,2],cexRow=0.2,cexCol = 1,labRow=input$IDENTIFIER, scale='row', trace='none',
                      as.matrix(input[which(abs(scale(PC$scores[,comp]))>2),][,-c(1:2)]))
    dev.off()  
  }
  
  Scores<-melt(data.frame(input[,1:2],scores=scale(PC$scores)))
  write.csv(Scores,file = paste("GDS",ID,"_PCscores.csv",sep = ""))
  
  ggsave(paste("GDS",ID,"_PCscores.pdf",sep = ""),ggplot(Scores[abs(Scores$value)>2,])+geom_bar(aes(IDENTIFIER,value),stat='identity')+facet_wrap(~variable))
}

plotClusters<-function(input,COLUMNS,ID,N){
  mInput<-melt(input,measure.vars = 3:(N))
  mDF<-merge(mInput,COLUMNS,by.x = "variable",by.y = "sample")
  
  ggsave(paste("GDS",ID,"_clusterplot.pdf",sep = ""),ggplot(mDF)+
           geom_hline(yintercept = c(0,2,-2))+
           #geom_line(aes(variable,value,color=cluster,group=ID_REF),alpha=0.2)+
           stat_summary(geom ="path",fun.y = mean,aes(variable,value,color=cluster,group=IDENTIFIER))+
           facet_wrap(~cluster)+scale_y_log10())
  
  ggsave(paste("GDS",ID,"_clusters.pdf",sep = ""),ggplot(mDF)+
           geom_hline(yintercept = c(0,2,-2))+
           geom_boxplot(aes(get(colnames(mDF[,6])),value,color=cluster),alpha=0.2)+
           facet_wrap(~cluster)+scale_y_log10()+theme(axis.text.x = element_text(angle=90)))
  
}

collectData<-function(input,allData,COLUMNS,ID){
  NR<-nrow(allData)
  input$dataset<-paste("GDS",ID,sep = "")
  N<-ncol(input)
  mDATA<-melt(input,id.vars = c(1:2,(N-1):N))
  mDATA<-merge(mDATA,COLUMNS,by.x='variable',by.y='sample')
  allData<-rbind(mDATA,allData,fill=TRUE)
  return(allData)
}

clusterGenes<-function(input){
  print("computing correlation")
  #print(head(input))
  N=ncol(input)
  Dgenes<- as.dist(1-(cor(t(na.omit(input[,3:N])))))
  HC<- hclust(Dgenes,method = "ward.D2")
  #plot(HC)
  cutHC<-cutree(HC,k = ceiling(sqrt(N)))
  input[,cluster:=as.factor(cutHC)]
  return(input)
}

PHNanalysis<-function(input,ID,COLUMNS){
  #PHN<-input[toupper(input$IDENTIFIER)%in%toupper(MGI[GO_ID%in%c("GO:0034620","GO:0006457")]$Symbol)]
  PHN<-input[toupper(input$IDENTIFIER)%in%toupper(chaperome$Symbol)]
  centeredFC<-apply(PHN[,-c("IDENTIFIER","ID_REF","cluster")],MARGIN = 1,FUN = function(r){ scale(r)})
  pdf(paste(ID,"heatmapPHN.pdf",sep=""))
  #print(PHN)
  gplots::heatmap.2(margins = c(10,5),trace='none',scale="none",col=viridis(200),cexRow = 0.3,labRow=PHN$IDENTIFIER, labCol=paste(COLUMNS$description),cexCol = .5,colRow = viridis(n = length(levels(PHN$cluster)),option = "E")[as.integer(PHN$cluster)],x = as.matrix(centeredFC))
  legend(fill = brewer.pal("Set1",n = 4),legend = c(1:4),x = 0.925,y=.98 ,cex = 0.7)
  dev.off()
  return(PHN)
}
