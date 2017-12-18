# ########
# #
# # RCircles in GGplot - V4: now with for loops across files and pbapply for TE and EVE assignments.
# # Patrick Dolan
# #
# #######

#Requirements
library(ggplot2)
library(data.table)
library(limma)
library(reshape2)
library(pbapply)
library(viridis)

# # # # # # Functions # # # # # # # # #
stagger<-function(input,scaler){
  flip=1
  staggerVec<-NULL
  for(x in 1:length(input)){
    staggerVec<-append(staggerVec, (scaler*flip) )
    flip= -flip
  }
  return(staggerVec)}
posVecs<-function(X){X['position5']:X['position3']}

piRNAmap<-function(piRNAs){
  piRNAs<-data.table(piRNAs,stringsAsFactors = T)
  sortedIndex<-allContigs[allContigs$ID%in%c("Aag2_Contigs"),]
  piRNAs[,contigID:=as.factor(strsplit2(as.character(contig),split = "|",fixed=T)[,1])]
  piRNAs[,adjStart:=(position5 + allContigs$startpositions[as.character(allContigs$CONTIG)==as.character(contigID)][1]),contigID]
  filteredpiRNAs<-piRNAs
  filteredpiRNAs[,l:=(nchar(as.character(sequence)))]
  filteredpiRNAs[,direc:=(-2*(as.integer(filteredpiRNAs$strand)-1.5))]
  filteredpiRNAs[,position3:=((direc*(l-1))+position5)]
  return(filteredpiRNAs)
}

CZ<-function(N){
  CZplot<-ggplot()+theme_minimal()+ylab("Relative TE Density")+xlab("Position")+scale_fill_brewer(palette = "Paired")+
    geom_segment(data=allContigs[allContigs$ID==lev,][N,],aes(x=startpositions,xend=stoppositions,y = stagger/4-.2,yend=stagger/4-.2),lwd=2)+
    geom_area(stat="bin",position = "fill",binwidth=25000,data=allTEs[allTEs$ID==lev&allTEs$contigID%in%allContigs$CONTIG[N],], aes(adjStart,fill=subclass),lwd=.25,color="black",alpha=1)
  ggsave(paste(inDir,"NewFigures/DensityFill",N,"_ContigTE.pdf",sep=""),CZplot,width = 8,height =3,dpi = 600)
}

RNAtally<-function(piRNAs){
  tally<-table(piRNAs$adjStart)
  piTally=data.table(tally,positions=as.integer(rownames(tally)))
  piTally[,Freq:=N/sum(piTally$N)]
  piTally[,logFreq:=log10(Freq)]
  piTally[,logCount:=log10(N)]
  piTally[,adjFreq:=(N/max(piTally$N))]
  piTally[,EVEtarget:=positions%in%allEVEpositions]
  return(piTally)
}

contigZoom<-function(X,limits=c(sortedIndex$startpositions[X],sortedIndex$stoppositions[X]),args=""){### Uses local contig position!!!
  selectPosTE<-plusTEs[as.character(plusTEs$contigID)==as.character(sortedIndex$CONTIG[X]),]
  selectNegTE<-minusTEs[as.character(minusTEs$contigID)==as.character(sortedIndex$CONTIG[X]),]
  levelN<-length(levels(selectNegTE$family))
  levels<-1:levelN/(.25*levelN)
  Zoom<-ggplot()+scale_color_brewer(palette="Paired")+
    coord_cartesian(xlim=limits,ylim=c(-4,12))+
    ggtitle(args)+
    geom_hline(aes(yintercept=c(7+levels,levels)),alpha=0.3,lwd=0.1)+
    geom_segment(data=sortedIndex[X,],aes(x=startpositions,xend=stoppositions,y = 12,yend=12),lwd=1)+
    
    geom_segment(data=plusEVEs,aes(x=adjStart,xend=adjEnd,y = 5.75,yend=5.75),alpha=0.5,lwd=0.5,arrow = arrow(ends = 'last',length = unit(.1,"cm"),type = "closed"))+
    geom_segment(data=minusEVEs,aes(x=adjStart,xend=adjEnd,y = 5.25,yend=5.25),alpha=0.5,lwd=0.5,arrow = arrow(ends = 'first',length = unit(.1,"cm"),type = "closed"))+
    
    geom_linerange(data=pluspiTally[pluspiTally$positions>=sortedIndex$startpositions[X]&pluspiTally$positions<=sortedIndex$stoppositions[X],],aes(x=positions,ymin = -2,ymax=-2+(logCount)/4),lwd=.3)+
    geom_linerange(data=minuspiTally[minuspiTally$positions>=sortedIndex$startpositions[X]&minuspiTally$positions<=sortedIndex$stoppositions[X],],aes(x=positions,ymax = -2,ymin=-2-(logCount)/4),lwd=.3)+
    
    geom_hline(yintercept = c(5.25,5.75,-2,-1,-3),lwd=.1)+
    
    geom_segment(data=selectPosTE,aes(x=adjStart, xend=adjEnd, y=7+as.integer(family)/(.25*levelN),color=subclass,yend=7+as.integer(family)/(.25*levelN)),arrow = arrow(length = unit(0.05,"cm"),type = "closed"))+
    geom_segment(data=selectNegTE,aes(x=adjStart, xend=adjEnd,y=as.integer(family)/(.25*levelN),color=subclass,yend=as.integer(family)/(.25*levelN)),arrow = arrow(ends = 'first',length = unit(0.05,"cm"),type = "closed"))+
    xlab(paste("Contig", sortedIndex[X,1]))+
    ylab("")+theme_minimal()+theme(axis.text.y=element_blank())
  
  ggsave(filename = paste("Zoom_Contig",X,limits[1],limits[2],"_new.pdf",sep="_"),Zoom+theme_bw(),width=6,height=3  )
}

strandZoom<-function(X,limits=c(sortedIndex$startpositions[X],sortedIndex$stoppositions[X]),args="")### Uses local contig position!!!
{
  selectPosTE<-plusTEs[as.character(plusTEs$contigID)==as.character(sortedIndex$CONTIG[X]),]
  selectNegTE<-minusTEs[as.character(minusTEs$contigID)==as.character(sortedIndex$CONTIG[X]),]
  selectTFsites<-superMerge[as.character(superMerge$Contig)%in%as.character(sortedIndex$CONTIG[X]),]
  levelN<-length(levels(selectNegTE$TE_ID))
  levels<-1:levelN/(.25*levelN)
  Zoom<-ggplot()+
    scale_color_brewer(palette = "Spectral")+
    scale_fill_brewer(palette = "Dark2")+
    coord_cartesian(xlim=limits,ylim=c(-5,9))+
    ggtitle(args)+
    ylab("log10(piRNA)")+
    xlab(paste("Contig", sortedIndex[X,1]))+
    geom_segment(data=selectPosTE ,aes(x=adjStart, xend=adjEnd, y=6,yend=6),alpha=0.8,lwd=1.0,arrow = arrow(length = unit(0.15,"cm"),type = "closed"))+
    geom_segment(data=selectNegTE,aes(x=adjStart, xend=adjEnd,y=6,yend=6),alpha=0.8,lwd=1.0,arrow = arrow(ends = 'first',length = unit(0.15,"cm"),type = "closed"))+
    geom_segment(data=plusEVEs,aes(x=adjStart,xend=adjEnd,y = 6,yend=6,color=family),lwd=1.0,arrow = arrow(ends = 'last',length = unit(.1,"cm"),type = "closed"))+
    geom_segment(data=minusEVEs,aes(x=adjStart,xend=adjEnd,y = 6,yend=6,color=family),lwd=1.0,arrow = arrow(ends = 'first',length = unit(.1,"cm"),type = "closed"))+
    
    geom_linerange(data=pluspiTally[pluspiTally$positions>=sortedIndex$startpositions[X]&pluspiTally$positions<=sortedIndex$stoppositions[X],],aes(x=positions,ymin = 0,ymax=0+(logCount)),lwd=.3)+
    geom_linerange(data=minuspiTally[minuspiTally$positions>=sortedIndex$startpositions[X]&minuspiTally$positions<=sortedIndex$stoppositions[X],],aes(x=positions,ymax = 0,ymin=0-(logCount)),lwd=.3)+
    geom_point(data=selectTFsites,aes(adjStart,-5+(3-as.integer(TFdir))/2,pch=TFdir,color=TF),cex=1)+
    geom_text(hjust = 0,angle=50,cex = 2,data=plusEVEs,aes(x=(adjStart+adjEnd)/2,y = 6,label=species,color=family),nudge_y = .25)+
    geom_text(hjust = 0,angle=50,cex = 2,data=minusEVEs,aes(x=(adjStart+adjEnd)/2,y = 6,label=species,color=family),nudge_y = .25)+
    theme_minimal()
  ggsave(filename = paste("StrandZoom",X,limits[1],limits[2],"_new.pdf",sep="_"),Zoom+theme_bw(),width=6,height=3)
}

CZplot<-function(N){
  G<-ggplot()+theme_minimal()+ylab("Relative TE Density")+xlab("Position")+scale_fill_brewer(palette = "Paired")+
    #geom_segment(data=allContigs[allContigs$ID==lev,][N,],aes(x=startpositions,xend=stoppositions,y = stagger/4-.2,yend=stagger/4-.2),lwd=2)+
    geom_area(stat="bin",position = "fill",binwidth=25000,data=allTEs[allTEs$contigID%in%allContigs$CONTIG[N],], aes(adjStart,fill=subclass),lwd=.25,color="black",alpha=1)+ylim(0,1)
  
  ggsave(paste(inDir,"NewFigures/fill",N,"_ContigTEDensity_new.pdf",sep=""),G,width = 6,height=2)
}

# # # # # # Directory information # # #
inDir<-"~/Research/EVEsAndpiRNA/FrozenData_4-27-17/"
setwd(inDir)
dir.create(paste(inDir,"NewFigures/"))

# 
# load("Contigs.RData")
# load("TEs_Aag2.RData")
# load("~/Research/EVEsAndpiRNA/FrozenData_4-27-17/allpiRNAs_new.RData")
# load("filtered_piRNAs.Rdata")
# 


# # # # # # Import CONTIG info # # # # #

ColNames<-c("CONTIG","LENGTH","I1","I2","I3")
allContigs<-NULL
for(filename in list.files(recursive = T,path = inDir,pattern = ".fa.fai",full.names = T)){
  print(filename)
  contigIndex<-fread(stringsAsFactors = T,filename,header = F,sep="\t")
  colnames(contigIndex)<-ColNames
  sortedIndex<-contigIndex[order(contigIndex$LENGTH,decreasing = T),]
  sortedIndex$startpositions<- c(1,cumsum(as.numeric(sortedIndex$LENGTH))[1:length(sortedIndex$LENGTH)-1])
  sortedIndex$stoppositions<- cumsum(as.numeric(sortedIndex$LENGTH))+1
  sortedIndex$stagger<- stagger(sortedIndex$startpositions,.2)
  print(sortedIndex$ID<-as.factor(strsplit2(strsplit2(filename,"/")[length(strsplit2(filename,"/"))],".fa")[1]))
  allContigs<-rbind(allContigs,sortedIndex)
}
#save(allContigs,file = "Contigs.RData")

# # # # # # Import TE Data # # # # #
allTEs<-NULL
for(filename in list.files(recursive = T,path = inDir,pattern = "TEs.bed",full.names = T)[1]){
  print(filename)
  TEs_All<-fread(stringsAsFactors = T,filename,header = F,sep="\t")
  colnames(TEs_All)<-c("contigID","start","end","info","score","strand","TE_ID")
  splitIDs<-strsplit2(TEs_All$TE_ID,"/")
  TEs_All$subclass<-as.factor(splitIDs[,1])
  TEs_All$family<-as.factor(splitIDs[,2])

  filteredTEs<-TEs_All
  filteredTEs$ID<-as.factor(strsplit2(strsplit2(filename,"/")[length(strsplit2(filename,"/"))],"_TEs.bed")[1])
  
  setContigs<- allContigs[as.character(allContigs$ID)==as.character(filteredTEs$ID[1]),]
  print(paste("mapping",filteredTEs$ID[1],"TEs to reordered",setContigs[1,9]))
  
  filteredTEs[,adjStart:=(start+setContigs$startpositions[as.character(setContigs$CONTIG)==as.character(contigID)][1]),contigID]
  filteredTEs$adjEnd<-filteredTEs$adjStart+(filteredTEs$end-filteredTEs$start)
  filteredTEs$ID<-as.factor(strsplit2(strsplit2(filename,"/")[length(strsplit2(filename,"/"))],"_TEs.bed")[1])
  allTEs<-rbind(allTEs, filteredTEs)
}

#save(allTEs , file = "TEs_Aag2.RData")
allTEs$length<-abs(allTEs$end-allTEs$start)

# # # # Import EVE info # # # # # #
allEVEsOut<-NULL
#for(filename in list.files(recursive = T,path = inDir,pattern = "Aag2_Contigs_EVEs_sorted.bed_withTaxonomy.txt",full.names = T)){
for(filename in "/Users/ptdolan/Research/EVEsAndpiRNA/FrozenData_4-27-17/Aag_Contigs_EVEs_NEW_cut_sorted.bed_withTaxonomy.txt"){
  print(filename)
  EVEs_Set<-read.delim(filename,header = T,sep="\t")
  
  filteredEVEs<-EVEs_Set
  filteredEVEs$adjStart<- -2
  
  filteredEVEs$ID<-"Aag2_Contigs"
  setContigs<- allContigs[as.character(allContigs$ID)==as.character(filteredEVEs$ID[1]),]
  print(paste("mapping",filteredEVEs$ID[1],"to",setContigs$ID[1]))
  for( row in 1:nrow(setContigs)){
    sI<-setContigs[row,]
    #print(head(sI))
    hits<-which(as.character(filteredEVEs$ContigEVE) == as.character(sI[[1]]))
    if(length(hits)>0){
      contigPos<-as.integer(sI[[6]])                           #start of EVE on Contig
      EVEPos<-filteredEVEs$EVEstart[hits]
      adjustedPos<-EVEPos+contigPos
      filteredEVEs$adjStart[hits]<-adjustedPos
    }
  }
  filteredEVEs$adjEnd<-filteredEVEs$adjStart+(filteredEVEs$EVEend-filteredEVEs$EVEstart)
  allEVEsOut<-rbind(allEVEsOut, filteredEVEs)
}

filtEVEs<- allEVEsOut[ allEVEsOut$adjStart!=-2,]
allEVEpositions<-NULL
for(eve in 1:nrow(filtEVEs)){
  allEVEpositions<-c(allEVEpositions,filtEVEs$adjStart[eve]:filtEVEs$adjEnd[eve])
}
allEVEpositions<-unique(allEVEpositions)
length(allEVEpositions)

# # # # Read piRNAs # # #
piRNAsPB<-fread("~/Research/EVEsAndpiRNA/FrozenData_4-27-17/I234_dsFluc-B_Aag2-PB_Frozen_v1.csv",stringsAsFactors = T,header=F)
#write.csv(piRNAsPB, "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/piRNA_PB.map",sep = "\t",row.names = F,col.names=F,quote = F)
fwrite(piRNAsPB, "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/piRNA_PB.table")
totN<-sum(allContigs$LENGTH[allContigs$ID=="Aag2_Contigs"])
colnames(piRNAsPB)<-c("sequence", "strand", "contig", "position5", "complement")
PBpiDT<-data.table(piRNAsPB,Platform=factor("PacBio"),totalN=totN)
# 
piRNAsLP<-fread("~/Research/EVEsAndpiRNA/FrozenData_4-27-17/I234_dsFluc-B_Aag-LP.map_piRNA.csv",header=F)
colnames(piRNAsLP)<-c("sequence", "strand", "contig", "position5", "complement")
LPpiDT<-data.table(piRNAsLP,Platform=factor("Sanger"),totalN=1309407401)
# 
allpis<-rbind(PBpiDT,LPpiDT)

##
#
# piRNAs<-fread("~/Dropbox/Aag2_Genome/MT/dsFluc-B_S2_c_Aag2_m1v1.map",stringsAsFactors = T)
# piRNAs[,ID:=strsplit2(V1,"-")[,1],]
# piRNAs[,counts:=strsplit2(V1,"-")[,2],]
# colnames(piRNAs)[1:5]<-c("num","strand","contig","position","seq")
# 
#
##

#save(allpis,file = "allpiRNAs_new.RData")
fwrite(allpis,file = "allpiRNAs_new.table")
filteredpiRNAs<-piRNAmap(piRNAsPB) # UNCOMMENT: if you want to look at only PB mapped reads.
#filteredpiRNAs<-piRNAmap(allpis) #UNCOMMENT:if you want to look at all mapped reads.
#save(filteredpiRNAs,file = "filtered_piRNAs.Rdata")
#catpositions<-pbapply(filteredpiRNAs,MARGIN = 1,FUN = posVecs)
#allpositions<-unlist(catpositions)
#piPositions<-unique(allpositions)


# Aag2TEs<-allTEs[allTEs$ID=="Aag2_Contigs",]
# allTEpositions<-NULL
# allpositions<-apply(Aag2TEs,1,function(X){
#   allTEpositions<-c(allTEpositions,X[2]:X[3])
# })
# 
# for(te in 1:nrow(Aag2TEs)){
#   print(te)
#   allTEpositions<-unique(c(allTEpositions,Aag2TEs$adjStart[eve]:Aag2TEs$adjEnd[te]))
# }

pluspiRNAs<-filteredpiRNAs[filteredpiRNAs$strand=="+",]
minuspiRNAs<-filteredpiRNAs[filteredpiRNAs$strand=="-",]

pluspiTally<-RNAtally(pluspiRNAs)
minuspiTally<-RNAtally(minuspiRNAs)

piClusters<-read.delim(file = "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/proTRAC_piRNAs.map_2016y11m1d22h53m3s/Aag2_piClusters.bed",header = F)

colnames(piClusters)<-c("contig","start",'stop','directionality','strand','length','pipositions')
head(piClusters)
for(piC in 1:nrow(piClusters)){
  piClusters$adjStart[piC]<-piClusters$start[piC]+allContigs$startpositions[as.character(allContigs$CONTIG)==as.character(piClusters$contig[piC])]                          #start of EVE on Conti
}
piClusters$adjEnd<-piClusters$adjStart+piClusters$length

# # # # #
# PLOTTING
# # # # #

# # # # Basic Circle Plot
Circos<-ggplot()+coord_polar()+theme_void()+scale_x_continuous(breaks = function(limits){signif(seq(limits[1],limits[2],length.out = 1),digits = 4)})+ylim(0,NA)+ylab("Annotation Track")+xlab("Genome Position")+scale_color_brewer(palette = "Paired")

# # # # Fig 1 - Alignment


contigSize<-ggplot(allContigs)+theme_bw()+scale_fill_brewer(palette = "Set1")+geom_histogram(aes(LENGTH,fill=ID,weight=LENGTH),position='dodge')+scale_x_log10()+ylab("Bases")+xlab("Contig Length")
ggsave(contigSize,filename = "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/NewFigures/contigSize_new.pdf",width = 5,height = 2)

#plot TEs per contig 
out<-dcast(allTEs,ID+contigID~.) #count
TEperContig<-ggplot(out)+theme_bw()+geom_density(aes(.,fill=ID,weight=.))+scale_x_log10()+xlab("TEs per Contig")
ggsave(TEperContig,filename = "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/NewFigures/TEperContig_new.pdf",width = 5,height = 2)

out<-dcast(allTEs,ID+contigID~., value.var="length",fun.aggregate = sum) #count
TEperContig<-ggplot(out)+theme_bw()+scale_fill_brewer(palette = "Set1")+geom_density(aes(.,fill=ID,weight=.))+scale_x_log10()+xlab("TE bases per Contig")
ggsave(TEperContig,filename = "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/NewFigures/TEbpperContig_new.pdf",width = 5,height = 2)


# # # # Fig 2 - TE content

lev="Aag2_Contigs"

allTEs$subclass<-factor(allTEs$subclass,ordered = T, c("DNA", "LINE", "LTR", "Unknown","MITEs","UD","Helitrons","SINE","Penelope","RC"))
allTEs$family<-with(allTEs,reorder(family,as.integer(subclass),mean))

tallytable<-dcast(allTEs,contigID+ID~.)

for(lev in "Aag2_Contigs"){
  print(paste("Plotting",lev,"..."))
  contigs<-allContigs[allContigs$ID==lev,][1:10,]
  TEtrack=3
  contigTrack=9.5
  
  tePlot<-function(contigs){Cplot<-Circos+
    geom_segment(data=contigs,aes(x=startpositions,xend=stoppositions,y = stagger+contigTrack,yend=stagger+contigTrack),lwd=6)+
    geom_point(size=I(0.1),
               data=allTEs[allTEs$ID==lev&allTEs$adjStart!=-2 & allTEs$contigID%in%contigs$CONTIG,],
               aes(x=adjStart,y=TEtrack+as.integer(family)/(.18*length(levels(allTEs$family[allTEs$ID==lev]))),
                   color=subclass),alpha=0.5,lwd=0.1)
  return(Cplot)
  }
  
  Cplot<-tePlot(contigs)
  ggsave(paste(inDir,"NewFigures/",lev,"_top10TEdotsAll_new.pdf",sep=""),Cplot,width = 10,height = 8,dpi = 600)
  ggsave(paste(inDir,"NewFigures/",lev,"_top10TEdotsAll.tiff",sep=""),Cplot,width = 10,height = 8,dpi = 600)
  
  C100plot<-tePlot(allContigs[allContigs$ID==lev,][1:100,])
  ggsave(paste(inDir,"NewFigures/",lev,"_top100TEdotsAll_new.pdf",sep=""),C100plot,width = 10,height = 8,dpi = 1200)
  
  Call_plot<-tePlot(allContigs[allContigs$ID==lev,])
  ggsave(paste(inDir,"NewFigures/",lev,"_TEdotsAll.tiff",sep=""),Call_plot,width = 10,height = 8,dpi = 1200)
  
  #
  for(i in 1:20){
    N=i
    CZ(N)
  }
  piechart<-ggplot( data=allTEs[allTEs$ID==lev&allTEs$adjStart!=-2 & allTEs$contigID%in%contigs$CONTIG,],aes(x=1,fill=subclass),color='black',lwd=1 )+scale_fill_brewer(palette = "Paired")+
    coord_polar(theta = "y")+geom_bar(width=1)+theme_void()+xlim( -1,2)
  
  barchart<-ggplot( data=allTEs[allTEs$ID==lev&allTEs$adjStart!=-2 & allTEs$contigID%in%contigs$CONTIG,],aes(x=1,fill=subclass),color='black',lwd=1 )+scale_fill_brewer(palette = "Paired")+geom_bar(width=1)+theme_void()
  
  ggsave(paste(inDir,"NewFigures/",lev,"_pieChart_new.pdf",sep=""),piechart,width = 8,height =3,dpi = 600)
  print("Done.")
  
}

# # # # Fig 3 - Global piRNA and EVEs
ggsave("ReadsMappedByPlatform.pdf",ggplot(allpis)+geom_bar(aes(Platform,fill=Platform))+ylab("piRNA Reads Mapped"),width=3,height=4)


# Fig3PiandEVEs<-Circos+scale_color_brewer(palette="Spectral")+
#   ylim(0, NA)+
#   #geom_segment(data=allContigs[allContigs$ID=="Aag2_Contigs",][c(9,24,372),],aes(x=startpositions,xend=stoppositions,y = stagger/2+2.5,yend=stagger/2+2.5),lwd=1)+
#   #geom_segment(data=allContigs[allContigs$ID=="Aag2_Contigs",],aes(x=startpositions,xend=stoppositions,y = stagger/2+2.5,yend=stagger/2+2.5),lwd=1,alpha=0.5)+
#   #geom_rect(data=piClusters,aes(ymin = 5.5,ymax = 8.5, xmin = adjStart, xmax = adjEnd),fill= "gold",color="gold",alpha=0.5,lwd=0.1)+
#   geom_hline(yintercept = 3+(2*(1:nlevels(filtEVEs$family))/nlevels(filtEVEs$family)),lwd=1.5,alpha=0.25)+
#   geom_linerange(data=filtEVEs,aes(x=adjStart,
#                                    ymin=3-1/length(levels(filtEVEs$family))+(2*as.integer(family)/length(levels(filtEVEs$family))),
#                                    ymax=3+1/length(levels(filtEVEs$family))+
#                                      (2*as.integer(family)/length(levels(filtEVEs$family))),
#                                    color=family),lwd=.7,alpha=0.8)+
#   geom_linerange(data=pluspiTally,aes(x=positions,ymin = 7,ymax=7+(logCount/4)),lwd=.2,alpha=0.4)+ #scaled to 4log10!!
#   geom_linerange(data=minuspiTally,aes(x=positions,ymax= 7,ymin=7-(logCount/4)),lwd=.2,alpha=0.4)+
#   geom_linerange(data=pluspiTally[pluspiTally$EVEtarget==T,],aes(x=positions,ymin = 7,ymax=7+(logCount)/4),lwd=.2,color='brown2')+ #scaled to 4log10!!
#   geom_linerange(data=minuspiTally[minuspiTally$EVEtarget==T,],aes(x=positions,ymax= 7,ymin=7-(logCount)/4),lwd=.2,color='brown2')+
#   geom_hline(yintercept = c(6,6.5,7,7.5,8),lwd=.2,alpha=0.2)+
#   geom_hline(yintercept = c(7),lwd=.5,color="grey")
# 
# ggsave("NewFigures/Fig3_piRNAsandEVEs2_new.pdf",Fig3PiandEVEs,width = 6,height=4.5,dpi = 800)

#plot rank of sites, activity and EVE status!

#Fig3 - Zoom to areas of interest
sortedIndex<-allContigs[allContigs$ID=="Aag2_Contigs",]
plusTEs<- allTEs[allTEs$ID=="Aag2_Contigs"&allTEs$strand=="+",]
minusTEs<-allTEs[allTEs$ID=="Aag2_Contigs"&allTEs$strand=="-",]

X<-10 #Contig Number (as rank)
limits=c(2.1E6,2.5E6) #NT limits

filtEVEs$family<-as.character(filtEVEs$family)
filtEVEs$family[filtEVEs$family==""]<-"Unannotated"
filtEVEs$family<-factor(filtEVEs$family)
filtEVEs$family<-with(filtEVEs,reorder(family,family,length))
EVEs<-filtEVEs[filtEVEs$ID=="Aag2_Contigs",]
plusEVEs<-EVEs[EVEs$EVEstrand=="+",]
minusEVEs<-EVEs[EVEs$EVEstrand=="-",]

piechartCount<-ggplot( data=filtEVEs,aes(x=1,fill=family) )+
  xlim(-1,1.6)+
  scale_fill_brewer(palette = "Spectral")+
  coord_polar(theta = "y")+geom_bar(stat='count',position = 'stack',aes(x=1,fill=family),color='black',lwd=.25)+theme_void()

ggsave(paste(inDir,"NewFigures/EVE_barChart_new.pdf",sep=""),barchart,width = 4,height =3,dpi = 600)
ggsave(paste(inDir,"NewFigures/EVE_pieChartCount_new.pdf",sep=""),piechartCount,width = 4,height =3,dpi = 600)

ggplot(filtEVEs)+theme_bw()+geom_histogram(aes(EVEend-EVEstart,fill=family),position = 'stack',bins=15,color="black")+scale_x_log10()+scale_fill_brewer(palette = "Spectral")+geom_vline(xintercept = median(filtEVEs$EVEend-filtEVEs$EVEstart))+xlab("Length")
ggplot(filtEVEs)+theme_bw()+geom_histogram(aes(EVEend-EVEstart,fill=family,weight=EVEend-EVEstart),position = 'stack',bins=15,color="black")+scale_x_log10()+scale_fill_brewer(palette = "Spectral")+geom_vline(xintercept = median(filtEVEs$EVEend-filtEVEs$EVEstart))+xlab("Length")+ylab("Number of bps")
median(filtEVEs$EVEend-filtEVEs$EVEstart)
sum(filtEVEs$EVEend-filtEVEs$EVEstart)


dir.create("EVEZooms")
setwd("EVEZooms")
for(EVEcontig in 1:nrow(filteredEVEs)){
  contigZoom(which(as.character(allContigs$CONTIG)==filteredEVEs$ContigEVE[EVEcontig]))
  strandZoom(which(as.character(allContigs$CONTIG)==filteredEVEs$ContigEVE[EVEcontig]),c(filteredEVEs$adjStart[EVEcontig]-1E4,filteredEVEs$adjEnd[EVEcontig]+1E4))
}

#for(cont in 1:3000){try(contigZoom(cont))} #NT limits
contigZoom(10) #NT limits
contigZoom(10,c(4.9E7,4.95E7))
strandZoom(10,c(4.9E7,4.95E7))

contigZoom(25)
contigZoom(25,c(1.02E8,1.03E8))
strandZoom(25,c(102460000,102525000))

contigZoom(373)
contigZoom(373,c(8.61E8,8.615E8))
strandZoom(373,c(8.61E8,8.6115E8))


contigZoom(459,c(9.754E8,9.757E8))
strandZoom(459,c(9.754E8,9.757E8))

contigZoom(459,c(9.754E8,9.757E8)) #whole + strand region of piCluster
strandZoom(459,c(9.754E8,9.757E8))

strandZoom(459,c(9.75508E8,9.75515E8))

contigZoom(378,c(870255E3,870275E3))


CZplot(10)
CZplot(25)
CZplot(372)
CZplot(933)
lev="Aag2_Contigs"
for(cluster in 1:nrow(piClusters)){
  print(cluster)
  hits<-which(as.character(allContigs$CONTIG)==as.character(piClusters$contig[cluster]))
  contigPos<-setContigs$startpositions[hits]
  CZ(hits)
  contigZoom(hits)
  strandZoom(hits,c(piClusters$adjStart[cluster],piClusters$adjEnd[cluster]),args=paste(piClusters$contig[cluster],piClusters$start[cluster],"-",piClusters$stop[cluster]))
}



as.integer(superMerge$TFdir)
CZ(906)
contigZoom()
CZ(10)
strandZoom(3)

genePred<-read.delim(comment.char = ">",file = "~/Research/EVEsAndpiRNA/FrozenData_4-27-17/Genomes/snapAnnotation.ann",header=F)
colnames(genePred)<-c("ElementType","start","end","Strand", "Score", "5overhang", "3overhang", "Frame", "identifier")


superMerge<-superMerge[superMerge$TF%in%c("SOX9","Gata4")]
superMerge[,TFdir:=factor(ifelse(seq%in%c("AGATAAC","AGATAAG","AACAATAA","AACAATAG","AACAATGA","AACAATGG"),"For","Rev")),]



#####KIMURA

dat<-fread("~/Research/EVEsAndpiRNA/FrozenData_4-27-17/Kimura/Aag2_TEsClosestToEVEs_overlapOrNearest_withKimura.txt",stringsAsFactors = T)
dat$family[(is.na(dat$family))]
dat$family<-with(dat,reorder(family,family,length))
dat$TEclass<-factor(dat$TEclass,ordered = T, c("DNA", "LINE", "LTR", "Unknown","MITEs","UD","Helitrons","SINE","Penelope","RC"))
dat$TEfamily<-with(dat,reorder(TEfamily,as.integer(TEclass),mean))
ggplot(dat)+geom_point(aes(family,KimuraDivergence,color=TEclass))+theme(axis.text.x = element_text(angle=90))

ggplot(dat[dat$EVEstrand==dat$TEstrand,])+
  geom_dotplot(aes(KimuraDivergence,fill=family),stackgroups = T,binwidth = 1,method = 'histodot')+
  scale_fill_brewer(palette = "Spectral")+theme_minimal()

ggplot(dat[dat$EVEstrand==dat$TEstrand,])+
  geom_dotplot(aes(KimuraDivergence,fill=TEclass),stackdir = 'down',stackgroups = T,binwidth = 1,method = 'histodot')+
  scale_fill_brewer(palette = "Paired")+theme_minimal()

plot(dat$TEclass+dat$family)
dcast(dat[dat$EVEstrand==dat$TEstrand,],formula = EVEdescription+TEclass+TEfamily~.,value.var = "KimuraDivergence")

