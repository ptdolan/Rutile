
# # R_useful.R
# Collection of useful functions in R.
# Edited: 2-24-17

# #
# pairPlot: pairwise correlation plots for ggplotters.
# usage: pairPlot(X, mask, cor=TRUE)
# X: dataframe with columns for each set of observations. 
# mask: Optional Boolean DF of same dimensions as X used to mask observations based on specific criteria
#
# Will generate pairwise scatter of all columns with Correlation coef.
# #
require(boot)
pairPlot<-function(X,mask=NA, ggstuff=NULL,cor=TRUE,al=0.7){
  if(!is.na(mask)){
    X[mask]<-NA
  }
  stacks<-data.frame()
  names<-colnames(X)
  for (i in 1:(ncol(X)-1)){
    for (j in (i+1):ncol(X)){
      #print(i)
      #print(j)
      pair=na.omit(X[,c(i,j)])
      block=data.frame(Yval=pair[,1],Xval=pair[,2],Xvar=names[i],Yvar=names[j],corr=corr(X[,c(i,j)]))
      stacks<-rbind(stacks,block)
    }
  }
  stacks$Xvar<-factor(stacks$Xvar)
  stacks$Yvar<-factor(stacks$Yvar)
  if(cor==TRUE){
    G<-ggplot(stacks)+geom_point(size=0.4,alpha = al,aes(Xval,Yval,color=corr))+facet_grid(Xvar~Yvar)+scale_color_gradient2(high = "red",low = 'blue',mid="darkgrey")
  }
  else{    
    G<-ggplot(stacks)+geom_point(size=0.4,alpha = al,aes(Xval,Yval))+facet_grid(Xvar~Yvar)
  }
  plot(G+ggstuff)
  return(stacks)
}

colGen<-function(z){#https://stackoverflow.com/questions/39117827/colorful-plot-using-persp
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)
  return(z.facet.range)
}
