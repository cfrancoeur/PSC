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

#Figure S5A, based on code from https://benjjneb.github.io/dada2/tutorial.html

top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
ps.top50 <- sort(ps.top50)
p<-plot_bar(ps.top50, x="compound_dose", fill="Phylum") + 
  scale_x_discrete(limits=c("Control5", "Control25", "ap5","ap25", "lina5", "lina25","noants")) 
p
pd <- p$data %>% as_tibble %>%
  mutate(Phylum = as.character(Phylum)) %>%
  replace_na(list(Phylum = "Phylum"))
genus_abun <- pd %>%
  group_by(Phylum) %>%
  summarize(Abundance = sum(Abundance)) %>%
  arrange(Abundance)
genus_levels <- genus_abun$Phylum
pd0 <- pd %>%
  mutate(Phylum = factor(Phylum, genus_levels))
colony.labels <- c('RM120223-02', 'CF180405-02', 'CF180406-01', 'HH180403-03', 'CR14')
names(colony.labels) <- c('Dora', 'Joan', "Leslie", 'ML','Regina')
ggplot(pd0, aes(x = compound_dose, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~colony, scales="free_x", labeller=labeller(colony=colony.labels)) + scale_x_discrete(limits=c("Control5", "Control25", "ap5","ap25", "lina5", "lina25","noants"))


##Figure S5B, based on code from https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html###
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Genus"){
  p1f = subset_taxa(physeq, Genus %in% c("Pantoea", "Klebsiella", "Pseudomonas","Serratia", "Burkholderia-Caballeronia-Paraburkholderia", "Enterobacter","Erwinia","Acinetobacter","Weissella"))
  #p1f= subset_taxa(physeq,Phylum %in% c("Proteobacteria"))
  p2f=subset_taxa(p1f, Genus %in% c("Pantoea", "Klebsiella", "Pseudomonas", "Serratia","Burkholderia-Caballeronia-Paraburkholderia", "Enterobacter", "Erwinia", "Bacillus"))
  #ps2 = prune_taxa(keepTaxa, ps)
  mphyseq = psmelt(p2f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  colony.labels <- c('RM120223-02', 'CF180405-02', 'CF180406-01', 'HH180403-03', 'CR14')
  names(colony.labels) <- c('Dora', 'Joan', "Leslie", 'ML','Regina')
  ggplot(data = mphyseq, mapping = aes_string(x = "Genus",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.2,
               position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, labeller=labeller(colony=colony.labels)) + scale_y_log10() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size =12), strip.text.x = element_text(size = 12)) + ylab("Abundance (log scale)")
  #+ theme(legend.position="none")
}
ps3ra = transform_sample_counts(ps, function(x){x / sum(x)})
plotGenus = plot_abundance(ps3ra,"", Color='Genus', Facet='colony') #+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
plotGenus
ggsave("select_genus.pdf", plot=plotGenus)
