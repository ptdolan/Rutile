#Annot data
library(seqinr)
library(limma)
library(data.table)
Coord<-fread(input = "PV_VaccineCohortData/MappingStatstics_M4a/Coordinates.txt",stringsAsFactors = T)
nOPV2="~/PV_VaccineCohortData/References/S2Cre5-DM.fa"#809:7430
Sabin2="~/PV_VaccineCohortData/References/Sabin2.fa"#
Mahoney="~/PV_VaccineCohortData/References/Mahoney_PV1.fasta.txt"#
Dengue2="~/PV_VaccineCohortData/References/Dengue16681.fasta"#97	10272

codons=as.character(c('TTT', 'TTC', 'TTA', 'TTG',
                      'TCT', 'TCC', 'TCA', 'TCG',
                      'TAT', 'TAC', 'TAA', 'TAG',
                      'TGT', 'TGC', 'TGA', 'TGG',
                      'CTT', 'CTC', 'CTA', 'CTG',
                      'CCT', 'CCC', 'CCA', 'CCG',
                      'CAT', 'CAC', 'CAA', 'CAG',
                      'CGT', 'CGC', 'CGA', 'CGG',
                      'ATT', 'ATC', 'ATA', 'ATG',
                      'ACT', 'ACC', 'ACA', 'ACG',
                      'AAT', 'AAC', 'AAA', 'AAG',
                      'AGT', 'AGC', 'AGA', 'AGG',
                      'GTT', 'GTC', 'GTA', 'GTG',
                      'GCT', 'GCC', 'GCA', 'GCG',
                      'GAT', 'GAC', 'GAA', 'GAG',
                      'GGT', 'GGC', 'GGA', 'GGG', NA))

residues = c( "F","F","L","L",
              "S","S","S","S",
              "Y","Y","X","X",
              "C","C","X","W",
              "L","L","L","L",
              "P","P","P","P",
              "H","H","Q","Q",
              "R","R","R","R",
              "I","I","I","M",
              "T","T","T","T",
              "N","N","K","K",
              "S","S","R","R",
              "V","V","V","V",
              "A","A","A","A",
              "D","D","E","E",
              "G","G","G","G", "U")


translateAndMutate<-function(fasta,startN=1,endN=0,regionTable=NA){
  inFasta<-read.fasta(fasta)
  LEN<-length(inFasta[[1]])
  if(endN==0){endN=LEN}
  faTrans<-getTrans.SeqFastadna(inFasta[[1]][startN:endN])
  mutInfo<-data.table(pos=1:length(inFasta[[1]]),wt=inFasta[[1]])
  mutInfo[pos%in%startN:endN,resPos:=floor((pos-startN+3)/3),]#Residue Positions
  mutInfo[pos%in%startN:endN,codon:=toupper(paste(wt,collapse = "")),resPos]#Wt Codons
  mutInfo[pos%in%startN:endN,wtRes:=residues[which(codons==codon)],resPos]#Wt Residue Positions
  
  mutTable<-data.table(pos=rep(1:length(inFasta[[1]]),each=4),mut=rep(c("a","c","g","t"),length(inFasta[[1]]))) #define all SNVs
  
  outputTable<-merge.data.table(mutInfo,mutTable,by = "pos")  ##merge mutation table with wt table
  
  outputTable[pos%in%startN:endN,mutCodon:=toupper(paste(collapse = "",replace(strsplit2(codon[1],split = ""),((pos-startN+3)%%3)+1,mut))),by=c("pos","mut")]
  outputTable[pos%in%startN:endN,mutRes:=residues[which(codons==mutCodon)],by=mutCodon]
  outputTable[pos%in%startN:endN,mutType:=ifelse(wtRes==mutRes,"Syn","NonSyn"),by=pos]
  outputTable[pos%in%startN:endN,mutType:=ifelse(wt==mut,"WT",mutType),by=pos]
  
  if(!is.na(regionTable)){
    regionTable[,outputTable[pos>=Position,reg:=Feature],by=Position]
    #annotate residues? need sequence of each virus. 
  }
  return(outputTable)
}

nOPV2out<-translateAndMutate(nOPV2,809,7429,Coord[virus=="nOPV2"])
Sabin2out<-translateAndMutate(Sabin2,748,7368,Coord[virus=="Sabin2"])
Mahoneyout<-translateAndMutate(Mahoney,746,7372,Coord[virus=="Mahoney"])
DENV2out<-translateAndMutate(Dengue2,97,10273)

