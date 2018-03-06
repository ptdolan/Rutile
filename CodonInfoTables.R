subTypeTable<-read.delim(file = "~/Research/CirSeq/ResSubTable_Pechmann2014.txt" , head = T)
rownames(subTypeTable)<-colnames(subTypeTable)
states=as.array(c('TTT', 'TTC', 'TTA', 'TTG',
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

alphaVec<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X","U")

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

res3 = c( "PHE","PHE","LEU","LEU",
          "SER","SER","SER","SER",
          "TYR","TYR","XXX","XXX",
          "CYS","CYS","XXX","TRP",
          "LEU","LEU","LEU","LEU",
          "PRO","PRO","PRO","PRO",
          "HIS","HIS","GLN","GLN",
          "ARG","ARG","ARG","ARG",
          "ILE","ILE","ILE","MET",
          "THR","THR","THR","THR",
          "ASN","ASN","LYS","LYS",
          "SER","SER","ARG","ARG",
          "VAL","VAL","VAL","VAL",
          "ALA","ALA","ALA","ALA",
          "ASP","ASP","GLU","GLU",
          "GLY","GLY","GLY","GLY", "UTR")

resCharacter=c('acidic',
               'acyclic',
               'aliphatic',
               'aromatic',
               'basic',
               'buried',
               'charged',
               'cyclic',
               'hydrophobic',
               'large',
               'medium',
               'negative',
               'neutral',
               'polar',
               'positive',
               'small',
               'surface',
               'stop')

resList=c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X")

resAnnotMatrix=##Binary Vector -- could be mproved to be continuous data, see K. Chou's stuff.
matrix(  
  #A  R  N  D  C  E  Q  G  H  I  L  K  M  F  P  S  T  W  Y  V  * U
  c(0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,# acidic
  1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0,# acyclic
  1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,# aliphatic
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,# aromatic
  0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,# basic
  1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0,# buried
  0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,# charged
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,# cyclic
  1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0,# hydrophobic
  0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0,# large
  0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,# medium
  0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,# negative
  1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,# neutral
  0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,# polar
  0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,# positive
  1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,# small
  0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0,# surface
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),# stop
  byrow=T,ncol=22,dimnames = list(resCharacter,c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X","U")))


