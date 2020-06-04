#Set the variables
wd <- "C:/Users/reole/Desktop/Research Summer 2020/Tipp/W1_contigs_tipp/markers"
#Creates a phylogenetic tree
running <- function(){
  library("ctv") 
  library(ape)
  #library(phytools)
  library(ggplot2)
  library(phangorn)
  library(seqinr)
  library(msa)
  library(tidyverse)
  library(ggtree)  
  library(ggimage)
  library(treeio)
}
running()
setwd(wd)

#Phylogenetics for the individual markers
file <- "tipp_dnaG_placement.json"
align <- "tipp_dnaG_alignment.fasta"
classification <- "tipp_dnaG_classification.txt"
tree <- read.tree(file)
alignment <- read.alignment(align, format = "fasta")
class <- read.table(classification, sep =",", header = FALSE)

msaplot(p=ggtree(tree), fasta = align, offset = 14, window=c(200, 300)) + geom_tiplab(size = 3, align = TRUE)

