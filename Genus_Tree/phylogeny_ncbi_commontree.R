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

#Set the variables
wd <- "C:/Users/reole/Desktop/Research Summer 2020/Tipp/Genus_Tree"
file <- "R4_contigs_2_tree.phy" #Tree created using NCBI common taxonomy database
annotation <- "C:/Users/reole/Desktop/Research Summer 2020/Tipp/R4_contigs_tipp/abundance.genus.csv"
sample_name <- "R4_contigs_2"
pic_name <- "R4_contigs_2.png"


#Create gene annotation file, can take first column for NCBI
setwd(wd)
gene <- read.table(annotation, header = TRUE, sep = "\t")
gene <- subset(gene, taxa != "unclassified" & abundance >= 0.002 & taxa != "Candidatus Liberibacter")
colnames(gene)[2] <- "abundance"
breaks <- quantile(gene$abundance)
breaks[1] <- breaks[1]-0.0001
breaks <- unique(breaks)  
labels <- vector()
for(i in 1:length(breaks)-1){
  labels[i] <- paste(breaks[i],"-",breaks[i+1])
}
gene$category <- cut(gene$abundance, breaks = breaks, labels = labels)
write.csv(gene[1],paste(paste(sample_name, "genus", sep="_"),".csv",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
png((paste("density_", sample_name,".png",sep="")), width = 200, height = 200, units='mm', res = 200)
plot(density(gene$abundance), main=(paste("Density Plot", sample_name))) #visualize density of abundance information
dev.off()
#AFTER MAKING NCBI NEED TO FIND AND REPLACE NAMES WITH: (\n\w+)\s(\w+\s) then sub for \1_\2

#Reads in alignment and trees
z <- readLines(file)
z <- gsub("'(\\w+)\\s(\\w+)'","\\1_\\2",z)
writeLines(z, file)
tree <- read.tree(file)
gene$taxa <- gsub("\\s","_",gene$taxa)
gene1<- merge(gene, tree[["tip.label"]], by.x = "taxa", by.y = "y", all.y = TRUE, all.x = FALSE)

ggtree(tree) %<+% gene1 + geom_tiplab(offset = 0.3, size=2.5, align=FALSE, linesize=.5) +
  geom_tippoint(aes(size = category, color = category)) + 
  theme(legend.position = "right") +  xlim(0, 25) +
  scale_size_discrete(range = c(1,4))

png(pic_name, width = 465, height = 275, units='mm', res = 300)
dev.off()


#If you want to compare results from tsv vs results from tipp
tsv_report <- "C:/Users/reole/Desktop/Research Summer 2020/centrifuge_report.tsv"
tsv <- read.table(tsv_report, header = TRUE, sep = "\t")
tsv <- subset(tsv, taxRank == "species")
tsv <- tsv[,c(1,7)]
colnames(tsv) <- c("taxa", "tsv_abundance")
tsv <- subset(tsv, tsv_abundance >= 0.0000001)
new_df <- merge(W4_contig, gene)



