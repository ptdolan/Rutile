# # # # Basic Circle Plot
Circos<-ggplot()+coord_polar()+theme_void()+
  scale_x_continuous(breaks = function(limits){signif(seq(limits[1],limits[2],length.out = 1),digits = 4)})+
  ylim(0,NA)+
  ylab("Annotation Track")+
  xlab("Position")+
  scale_color_brewer(palette = "Paired")
