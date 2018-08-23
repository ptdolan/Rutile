library(data.table)
library(limma)
library(muscle)
library(Biostrings)
library(seqinr)
library(lettercase)

CDS<-Biostrings::readDNAStringSet("/Users/ptdolan/Downloads/sequence.txt") #CDS file from NCBI (fasta format)
TLN<- Biostrings::translate(CDS,if.fuzzy.codon = "solve") #translate your CDS.

aligned<-muscle(TLN)

AlignedSet<-AAStringSet(aligned)
outputNT<-data.frame()
for(i in 1:length(AlignedSet)){
  print(i)
  aaAligned <- (AlignedSet[i] %>% toString() %>% strsplit2(split=""))[[1]]
  codonTemp <- (CDS[i] %>% toString() %>% strsplit2(split=""))[[1]]
  aligner <- NULL
  pos=1
  for (j in aaAligned){
    if (j != "-"){
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

#name parsing

ns1<-limma::strsplit2(outputNT$name,split = "lcl\\|| ")[,2]
ns2<-(strsplit2(outputNT$name,split = "locus_tag=")[,2] %>%   strsplit2(split = "\\]|\\["))[,1]

write.fasta(DMSA,file.out = "codonAlign.fa",names = ns2)



