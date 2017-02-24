#import from R_useful on GitHub
source("https://raw.githubusercontent.com/ptdolan/Rutile/master/R_useful.R")

#make a numeric(!) DF

DF<-data.frame(iris[,1:4])
X<-DF

#deploy pairwise scatter
pairPlot(DF)
