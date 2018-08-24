# ---
# title: "CodonAlign.R"
# author: "PTDolan"
# date: "8/23/2018"
# description: Codon-align CDS sequences in fasta format. 
# usage: "RScript CodonAlign.R <myFasta.fa>"
# ---

# UnComment for Installation

# source("https://bioconductor.org/biocLite.R")
# biocLite("limma",ask = F)
# biocLite("muscle",ask = F)
# biocLite("Biostrings",ask = F)
# biocLite("seqinr",ask = F)
# install.packages("lettercase")
# install.packages("data.table")

library(data.table,quietly = T)
library(limma,quietly = T)
library(muscle,quietly = T)
library(Biostrings,quietly = T)
library(seqinr,quietly = T)
library(lettercase,quietly = T)

args = commandArgs(trailingOnly=TRUE) #Takes sequence from command line

infiles=list.files(pattern = ".fa|.fasta")

ifelse(length(args)>0,
       infiles<-args, #CDS files from NCBI (fasta format)
       "No input CDS fasta. Searching in Working Directory...")

for (infile in infiles){
  print(infile)
  
  CDS<-Biostrings::readDNAStringSet(infile) #CDS file from NCBI (fasta format)
  TLN<- Biostrings::translate(CDS,if.fuzzy.codon = "solve") #translate your CDS.
  
  aligned<-muscle(TLN)
  
  AlignedSet<-AAStringSet(aligned)
  outputNT<-data.frame()
  
  #Main Function
  for(i in 1:length(AlignedSet)){
    flush.console()
    cat(".")
    #print(i/length(AlignedSet),digits = 2)
    
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
  outfile<-paste(strsplit2(infile,"\\.")[1],"_codAln.fa",sep = "")
  writeXStringSet(format = 'fasta',filepath = paste(strsplit2(infile,"\\.")[1],"_codAln.fa",sep = ""),x = DMSA)
  print(paste("Wrote ",length(DMSA)," sequences to ",outfile,sep = ""))
}

print("Done.")



