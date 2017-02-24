#R_useful
DF<-data.frame(iris[,1:4])
X<-DF
#pairwiseScatter


# #
# pairPlot: pairwise correlation plots for ggplotters.
# usage: pairPlot(dataframe with columns)
# #

pairPlot<-function(X,mask=NA, ggstuff=NULL,cor=TRUE){
  if(!is.na(mask)){
    X[mask]<-NA
  }
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
  G<-ggplot(stacks)+geom_point(size=0.4,aes(Xval,Yval,color=corr))+facet_grid(Xvar~Yvar)+scale_color_gradient2(high = "red",low = 'blue',mid="darkgrey")
  
  plot(G+ggstuff)
  return(stacks)
}

pairPlot(DF)
