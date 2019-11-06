##test stuff from https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html ##
rank_names(ps)

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#phylum
prevdf1 = subset(prevdf,Genus %in% get_taxa_unique(ps, "Genus"))
prevdf_proteo=subset(prevdf1, Phylum == "Proteobacteria")
prevdf2=subset(prevdf, Phylum%in% get_taxa_unique(ps, "Phylum"))
ggplot(prevdf_proteo, aes(TotalAbundance, Prevalence / nsamples(ps),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")

#filter out low prevalence taxa##
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold 
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)


length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  #p1f = subset_taxa(physeq, Genus %in% c("Acinetobacter",""))
  #ps2 = prune_taxa(keepTaxa, ps)
  mphyseq = psmelt(ps)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "compound_dose",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}



ps3ra = transform_sample_counts(ps, function(x){x / sum(x)})
plotnotrans=plot_abundance(ps, "", Facet="Family")
plotFamily = plot_abundance(ps3ra,"", Facet="Family") + scale_x_discrete(limits=c("Control5", "Control25", "ap5","ap25", "lina5", "lina25","noants")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plotGenus = plot_abundance(ps3ra,"", Facet="Genus") + scale_x_discrete(limits=c("Control5", "Control25", "ap5","ap25", "lina5", "lina25","noants")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plotPhylum = plot_abundance(ps3ra,"", Facet="Phylum") + scale_x_discrete(limits=c("Control5", "Control25", "ap5","ap25", "lina5", "lina25","noants")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

grid.arrange(nrow = 2,  plotnotrans, plotFamily)
plotnotrans
plotFamily
plotGenus
plotPhylum


##by colony - Figure 5B and 5C code###
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  #p1f = subset_taxa(physeq, Genus %in% c("Pantoea", "Klebsiella", "Pseudomonas","Serratia", "Burkholderia-Caballeronia-Paraburkholderia", "Enterobacter","Erwinia","Acinetobacter","Weissella"))
  p1f= subset_taxa(physeq,Phylum %in% c("Proteobacteria"))
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
plotGenus = plot_abundance(ps3ra,"", Color='Genus', Facet='colony') #+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

plotGenus
ggsave("select_genus.pdf", plot=plotGenus)
plotPhylum = plot_abundance(ps3ra,"", Facet="colony", Color="Genus") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plotPhylum


# Open a pdf file
pdf("rplot.pdf", width=10, height=6) 
# 2. Create a plot
plot_abundance(ps3ra,"", Color='Genus', Facet='colony') 
# Close the pdf file
dev.off() 

