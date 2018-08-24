# ---
# title: "CodonAlign.R"
# author: "PTDolan"
# date: "8/23/2018"
# description: Codon-align CDS sequences in fasta format. 
# usage: 
# ---


# UnComment for Installation

# source("https://bioconductor.org/biocLite.R")
# biocLite("limma",ask = F)
# biocLite("muscle",ask = F)
# biocLite("Biostrings",ask = F)
# biocLite("seqinr",ask = F)
# biocLite("seqinr",ask = F)
# install.packages("lettercase")
# install.packages("data.table")

library(data.table)
library(limma)
library(muscle)
library(Biostrings)
library(seqinr)
library(lettercase)

args = commandArgs(trailingOnly=TRUE)#Takes sequence from command line

ifelse(length(args)>0,
       CDS<-Biostrings::readDNAStringSet(args), #CDS file from NCBI (fasta format)
       "No input CDS fasta.")

#CDS<-Biostrings::readDNAStringSet("/Users/ptdolan/Downloads/sequence.txt") #CDS file from NCBI (fasta format)
TLN<- Biostrings::translate(CDS,if.fuzzy.codon = "solve") #translate your CDS.

aligned<-muscle(TLN)

AlignedSet<-AAStringSet(aligned)
outputNT<-data.frame()

for(i in 1:length(AlignedSet)){
  print(i/length(AlignedSet),digits = 2)
  aaAligned <- (AlignedSet[i] %>% toString() %>% strsplit2(split=""))
  codonTemp <- (CDS[i] %>% toString() %>% strsplit2(split=""))
  aligner <- NULL
  pos=1
  for (j in aaAligned){
    if (j != "-"){
      #print(j)
      aligner<-append(aligner,codonTemp[((pos*3)-2):(pos*3)])
      pos=pos+1
    }
    else{
      aligner<-append(aligner,"---")
    }
  }
  codonAlignedSeq<- aligner %>% as.character() %>% paste(collapse="")
  codonAligned<-data.frame(seq=codonAlignedSeq,name=names(AlignedSet[i]))
  outputNT<-rbind.data.frame(outputNT,codonAligned)
}

DMSA<-DNAStringSet(outputNT$seq)



#NCBI Name Parsing

ns1<-limma::strsplit2(outputNT$name,split = "lcl\\|| ")[,2]
ns2<-(strsplit2(outputNT$name,split = "locus_tag=")[,2] %>%   strsplit2(split = "\\]|\\["))[,1]

names(DMSA)<-ns2
writeXStringSet(format = 'fasta',filepath = ".fa",x = DMSA)



