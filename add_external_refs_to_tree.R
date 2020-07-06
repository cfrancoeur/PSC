##code to create Figure S13

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(vegan)
library(tidyverse)
library(dplyr)
library(ape)
library(ggpubr)
library(DECIPHER)
library(phangorn)



###############add external references ############# based off https://github.com/joey711/phyloseq/issues/1150

##PROTEOBACTERIA
#create a phyloseq object with an empty tree slot
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1

#take subset of ASVs of interest
Prot = subset_taxa(ps1, Phylum %in% c("Proteobacteria")) #taxa of interest

#read in reference file/create DNA string set
rs <- readDNAStringSet(file = "reference.txt",
                       format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#make sample data and OTU table for reference sequences
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Prot)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#add taxonomic levels with random strings for 6 ranks
randos <- ids::random_id(624,4) #104 taxa x 6 taxonomic ranks = 624
wholetax <- replicate(6,sample(randos, 104, rep=TRUE))
tax_df <- data.frame(wholetax) ## 6 columns for 88 references
rownames(tax_df) <- names(rs)
names(tax_df)[names(tax_df) == 'X1'] <- 'Kingdom'
names(tax_df)[names(tax_df) == 'X2'] <- 'Phylum'
names(tax_df)[names(tax_df) == 'X3'] <- 'Class'
names(tax_df)[names(tax_df) == 'X4'] <- 'Order'
names(tax_df)[names(tax_df) == 'X5'] <- 'Family'
names(tax_df)[names(tax_df) == 'X6'] <- 'Genus'
tax_df$Genus <- names(rs)

ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

#tree calculation
ps.merged <- merge_phyloseq(Prot, ps.ref) #number of taxa in ps.merged should be the sum of the two phyloseq objects
pseqs <- refseq(ps.merged)
palignment <- AlignSeqs(pseqs, anchor = NA, normPower=0) #normPower added to account for differences in 16S length
psalignment <- StaggerAlignment(palignment) #added this step as recommended in "The Art of Multiple Sequence Alignment in R by Erick S. Wright"
pphang.align <- phyDat(as(psalignment, "matrix"), type="DNA")
pdm <- dist.ml(pphang.align)
ptreeNJ <- NJ(pdm) # Note, tip order != sequence order
pfit = pml(ptreeNJ, data = pphang.align)
pfitGTR <- update(pfit, k = 4, inv = 0.2)
pfitGTR <- optim.pml(pfitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

#calcalate bootstraps (from "Estimating phylogenetic trees with phangorn" by Klaus P. Schliep)

bs <-bootstrap.pml(pfitGTR, bs = 250, optNni=TRUE, control=pml.control(trace=0))
bs_tree <- plotBS(midpoint(pfitGTR$tree), bs, p =50, type ='p')

detach("package:phangorn", unload = TRUE)

# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(pfitGTR$tree))
plot_tree(ps.merged, label.tips='taxa_names', ladderize='left') ##no bootstraps

##export newick and edit in FigTree
tree1 = phy_tree(ps.merged)
ape::write.tree(tree1, "./tree2.asv.newick")

##bootstrap tree export for FigTree
tree2 = phy_tree(bs_tree)
ape::write.tree(tree2, "./bs_burkholderia_tree.asv.newick")

######BY GENUS#########
##same code as above with different levels of filtering

###Burkholderia
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1
Burk2 = subset_taxa(ps1, Genus %in% c("Burkholderia-Caballeronia-Paraburkholderia")) #taxa of interest

#read in reference file/create DNA string set
rs <- readDNAStringSet(file = "burkholderia_reference_2.txt",
                       format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#make sample data and OTU table for reference sequences
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Burk2)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#add taxonomic levels with random strings for 6 ranks, 6 x # of reference seq
randos <- ids::random_id(132,4)
wholetax <- replicate(6,sample(randos, 22, rep=TRUE))
tax_df <- data.frame(wholetax) ## 6 columns for X references
rownames(tax_df) <- names(rs)
names(tax_df)[names(tax_df) == 'X1'] <- 'Kingdom'
names(tax_df)[names(tax_df) == 'X2'] <- 'Phylum'
names(tax_df)[names(tax_df) == 'X3'] <- 'Class'
names(tax_df)[names(tax_df) == 'X4'] <- 'Order'
names(tax_df)[names(tax_df) == 'X5'] <- 'Family'
names(tax_df)[names(tax_df) == 'X6'] <- 'Genus'
tax_df$Genus <- names(rs)
ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

#tree calculation
library(phangorn)
ps.merged <- merge_phyloseq(Burk2, ps.ref) #number of taxa in ps.merged should be the sum of the two phyloseq objects
pseqs <- refseq(ps.merged)
palignment <- AlignSeqs(pseqs, anchor = NA, normPower=0)
psalignment <- StaggerAlignment(palignment)
pphang.align <- phyDat(as(psalignment, "matrix"), type="DNA")
pdm <- dist.ml(pphang.align)
ptreeNJ <- NJ(pdm) # Note, tip order != sequence order
pfit = pml(ptreeNJ, data = pphang.align)
pfitGTR <- update(pfit, k = 4, inv = 0.2)
pfitGTR <- optim.pml(pfitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

bs <-bootstrap.pml(pfitGTR, bs =250, optNni=TRUE, control=pml.control(trace=0))
bs_tree <- plotBS(midpoint(pfitGTR$tree), bs, p =50, type ='p')

detach("package:phangorn", unload = TRUE)
# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(pfitGTR$tree))
plot_tree(ps.merged, label.tips='taxa_names', ladderize='left')

##export newick and edit in FigTree
tree1 = phy_tree(ps.merged)
ape::write.tree(tree1, "./burkholderia_tree.asv.newick")

tree2 = phy_tree(bs_tree)
ape::write.tree(tree2, "./bs_burkholderia_tree.asv.newick")


###Pantoea
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1
Pantoea = subset_taxa(ps1, Genus %in% c("Pantoea")) #taxa of interest

#read in reference file/create DNA string set
rs <- readDNAStringSet(file = "pantoea_reference.txt",
                       format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#make sample data and OTU table for reference sequences
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Pantoea)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#add taxonomic levels with random strings for 6 ranks, 6 x # of reference seq
randos <- ids::random_id(66,4)
wholetax <- replicate(6,sample(randos, 11, rep=TRUE))
tax_df <- data.frame(wholetax) ## 6 columns for X references
rownames(tax_df) <- names(rs)
names(tax_df)[names(tax_df) == 'X1'] <- 'Kingdom'
names(tax_df)[names(tax_df) == 'X2'] <- 'Phylum'
names(tax_df)[names(tax_df) == 'X3'] <- 'Class'
names(tax_df)[names(tax_df) == 'X4'] <- 'Order'
names(tax_df)[names(tax_df) == 'X5'] <- 'Family'
names(tax_df)[names(tax_df) == 'X6'] <- 'Genus'
tax_df$Genus <- names(rs)
ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

#tree calculation
library(phangorn)
ps.merged <- merge_phyloseq(Pantoea, ps.ref) #number of taxa in ps.merged should be the sum of the two phyloseq objects
pseqs <- refseq(ps.merged)
palignment <- AlignSeqs(pseqs, anchor = NA, normPower=0)
psalignment <- StaggerAlignment(palignment)
pphang.align <- phyDat(as(psalignment, "matrix"), type="DNA")
pdm <- dist.ml(pphang.align)
ptreeNJ <- NJ(pdm) # Note, tip order != sequence order
pfit = pml(ptreeNJ, data = pphang.align)
pfitGTR <- update(pfit, k = 4, inv = 0.2)
pfitGTR <- optim.pml(pfitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

bs <-bootstrap.pml(pfitGTR, bs = 250, optNni=TRUE, control=pml.control(trace=0))
bs_tree <- plotBS(midpoint(pfitGTR$tree), bs, p =50, type ='p')

detach("package:phangorn", unload = TRUE)
# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(pfitGTR$tree))

plot_tree(ps.merged, label.tips='taxa_names', ladderize='left')

##export newick and edit in FigTree
tree1 = phy_tree(ps.merged)
ape::write.tree(tree1, "./pantoe_tree.asv.newick")

tree2 = phy_tree(bs_tree)
ape::write.tree(tree2, "./bs_pantoea_tree.asv.newick")

###Pseudomonas
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1
Pseudomonas = subset_taxa(ps1, Genus %in% c("Pseudomonas")) #taxa of interest

#read in reference file/create DNA string set
rs <- readDNAStringSet(file = "pseudomonas_reference.txt",
                       format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#make sample data and OTU table for reference sequences
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Pseudomonas)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#add taxonomic levels with random strings for 6 ranks, 6 x # of reference seq
randos <- ids::random_id(30,4)
wholetax <- replicate(6,sample(randos, 5, rep=TRUE))
tax_df <- data.frame(wholetax) ## 6 columns for X references
rownames(tax_df) <- names(rs)
names(tax_df)[names(tax_df) == 'X1'] <- 'Kingdom'
names(tax_df)[names(tax_df) == 'X2'] <- 'Phylum'
names(tax_df)[names(tax_df) == 'X3'] <- 'Class'
names(tax_df)[names(tax_df) == 'X4'] <- 'Order'
names(tax_df)[names(tax_df) == 'X5'] <- 'Family'
names(tax_df)[names(tax_df) == 'X6'] <- 'Genus'
tax_df$Genus <- names(rs)
ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

#tree calculation
library(phangorn)
ps.merged <- merge_phyloseq(Pseudomonas, ps.ref) #number of taxa in ps.merged should be the sum of the two phyloseq objects
pseqs <- refseq(ps.merged)
palignment <- AlignSeqs(pseqs, anchor = NA, normPower=0)
psalignment <- StaggerAlignment(palignment)
pphang.align <- phyDat(as(psalignment, "matrix"), type="DNA")
pdm <- dist.ml(pphang.align)
ptreeNJ <- NJ(pdm) # Note, tip order != sequence order
pfit = pml(ptreeNJ, data = pphang.align)
pfitGTR <- update(pfit, k = 4, inv = 0.2)
pfitGTR <- optim.pml(pfitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

bs <-bootstrap.pml(pfitGTR, bs = 250, optNni=TRUE, control=pml.control(trace=0))
bs_tree <- plotBS(midpoint(pfitGTR$tree), bs, p =50, type ='p')

detach("package:phangorn", unload = TRUE)
# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(pfitGTR$tree))
plot_tree(ps.merged, label.tips='taxa_names', ladderize='left')

##export newick and edit in FigTree
tree1 = phy_tree(ps.merged)
ape::write.tree(tree1, "./pseudomonas_tree.asv.newick")

tree2 = phy_tree(bs_tree)
ape::write.tree(tree2, "./bs_pseudomonas_tree.asv.newick")

###Klebsiella
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1
Klebsiella = subset_taxa(ps1, Genus %in% c("Klebsiella")) #taxa of interest

#read in reference file/create DNA string set
rs <- readDNAStringSet(file = "klebsiella_reference.txt",
                       format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#make sample data and OTU table for reference sequences
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Klebsiella)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#add taxonomic levels with random strings for 6 ranks, 6 x # of reference seq
randos <- ids::random_id(42,4)
wholetax <- replicate(6,sample(randos, 7, rep=TRUE))
tax_df <- data.frame(wholetax) ## 6 columns for X references
rownames(tax_df) <- names(rs)
names(tax_df)[names(tax_df) == 'X1'] <- 'Kingdom'
names(tax_df)[names(tax_df) == 'X2'] <- 'Phylum'
names(tax_df)[names(tax_df) == 'X3'] <- 'Class'
names(tax_df)[names(tax_df) == 'X4'] <- 'Order'
names(tax_df)[names(tax_df) == 'X5'] <- 'Family'
names(tax_df)[names(tax_df) == 'X6'] <- 'Genus'
tax_df$Genus <- names(rs)
ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

#tree calculation
library(phangorn)
ps.merged <- merge_phyloseq(Klebsiella, ps.ref) #number of taxa in ps.merged should be the sum of the two phyloseq objects
pseqs <- refseq(ps.merged)
palignment <- AlignSeqs(pseqs, anchor = NA, normPower=0)
psalignment <- StaggerAlignment(palignment)
pphang.align <- phyDat(as(psalignment, "matrix"), type="DNA")
pdm <- dist.ml(pphang.align)
ptreeNJ <- NJ(pdm) # Note, tip order != sequence order
pfit = pml(ptreeNJ, data = pphang.align)
pfitGTR <- update(pfit, k = 4, inv = 0.2)
pfitGTR <- optim.pml(pfitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

#calculate bootstraps
bs <-bootstrap.pml(pfitGTR, bs = 250, optNni=TRUE, control=pml.control(trace=0))
bs_tree <- plotBS(midpoint(pfitGTR$tree), bs, p =50, type ='p')
detach("package:phangorn", unload = TRUE)

# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(pfitGTR$tree))
plot_tree(ps.merged, label.tips='taxa_names', ladderize='left')

##export newick and edit in FigTree
tree1 = phy_tree(ps.merged)
ape::write.tree(tree1, "./klebsiella_tree.asv.newick")

tree2 = phy_tree(bs_tree)
ape::write.tree(tree2, "./bs_klebsiella_tree.asv.newick")

###Enterobacter
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1
Ent = subset_taxa(ps1, Genus %in% c("Enterobacter", "Erwinia")) #taxa of interest

#read in reference file/create DNA string set
rs <- readDNAStringSet(file = "enterobacter_reference.txt",
                       format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#make sample data and OTU table for reference sequences
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Ent)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#add taxonomic levels with random strings for 6 ranks, 6 x # of reference seq
randos <- ids::random_id(64,4)
wholetax <- replicate(6,sample(randos, 8, rep=TRUE))
tax_df <- data.frame(wholetax) ## 6 columns for X references
rownames(tax_df) <- names(rs)
names(tax_df)[names(tax_df) == 'X1'] <- 'Kingdom'
names(tax_df)[names(tax_df) == 'X2'] <- 'Phylum'
names(tax_df)[names(tax_df) == 'X3'] <- 'Class'
names(tax_df)[names(tax_df) == 'X4'] <- 'Order'
names(tax_df)[names(tax_df) == 'X5'] <- 'Family'
names(tax_df)[names(tax_df) == 'X6'] <- 'Genus'
tax_df$Genus <- names(rs)
#tax_df <- tax_df %>% separate(wholetax, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

#tree calculation
library(phangorn)
ps.merged <- merge_phyloseq(Ent, ps.ref) #number of taxa in ps.merged should be the sum of the two phyloseq objects
pseqs <- refseq(ps.merged)
palignment <- AlignSeqs(pseqs, anchor = NA, normPower=0)
psalignment <- StaggerAlignment(palignment)
pphang.align <- phyDat(as(psalignment, "matrix"), type="DNA")
pdm <- dist.ml(pphang.align)
ptreeNJ <- NJ(pdm) # Note, tip order != sequence order
pfit = pml(ptreeNJ, data = pphang.align)
pfitGTR <- update(pfit, k = 4, inv = 0.2)
pfitGTR <- optim.pml(pfitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

bs <-bootstrap.pml(pfitGTR, bs = 250, optNni=TRUE, control=pml.control(trace=0))
bs_tree <- plotBS(midpoint(pfitGTR$tree), bs, p =50, type ='p')

detach("package:phangorn", unload = TRUE)
# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(pfitGTR$tree))
plot_tree(ps.merged, label.tips='taxa_names', ladderize='left')

##export newick and edit in FigTree
tree1 = phy_tree(ps.merged)
ape::write.tree(tree1, "./enterobacter_tree.asv.newick")

tree2 = phy_tree(bs_tree)
ape::write.tree(tree2, "./bs_enterobacter_tree.asv.newick")

