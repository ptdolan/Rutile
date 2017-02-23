#pairwiseScatter

DF<-data.frame(randu)
N<-paste("name",1:length(DF),sep="-")
colnames(DF)<-N

##
##pairPlot
##
X<-DF
pairPlot<-function(X=DF,mask=NULL, ggstuff=NULL,cor=TRUE){
  X[mask]<-NA
  stacks<-data.frame()
  names<-colnames(X)
  for (i in 1:(ncol(X)-1)){
    for (j in (i+1):ncol(X)){
      pair=na.omit(X[,c(i,j)])
      block=data.frame(Xval=pair[,1],Yval=pair[,2],Xvar=names[i],Yvar=names[j],corr=corr(pair))
      stacks<-rbind(stacks,block)
    }
  }
  stacks$Xvar<-factor(stacks$Xvar)
  stacks$Yvar<-factor(stacks$Yvar)
  G<-ggplot(stacks)+geom_point(size=0.4,aes(Xval,Yval,color=corr))+facet_grid(Xvar~Yvar)+scale_color_gradient(high = "red",low="darkgrey")
  
  plot(G+ggstuff)
  
  return(stacks)
}

mask <- DF<0.05
pairPlot(DF,mask = mask)
