# Rutile
Useful r scripts

<h1>Rutile: R utilities resource</h1>

<h3>pairPlot(X,mask=NA,ggstuff=NULL,cor=TRUE)</h3>

>X: Numeric dataframe where each column is the observations for a given set. (e.g. FC values for treatment ~ gene)

>mask: Boolean DF of same dimensions as X. e.g. P-values for FC values < 0.05  mask = Pmatrix<0.05

>ggstuff: a single layer to ggplot2. e.g.: "geom_vline(0.05)"

![alt text](https://github.com/ptdolan/Rutile/raw/master/pairPlot_IrisExample.png "pairPlot() Iris Example")




~~~source("https://raw.githubusercontent.com/ptdolan/Rutile/master/R_useful.R)

