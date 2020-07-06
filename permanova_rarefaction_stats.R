##continuation of dada2 pipeline

##get OTU table
OTU1=as(otu_table(ps),"matrix")
if(taxa_are_rows(ps)){OTU1<-t(OTU1)}
OTUdf=as.data.frame(OTU1)

#permanova using vegan
BC.dist=vegdist(OTUdf, distance="bray")
adonis(BC.dist ~ compound_dose, data = samdf_manual, permutations = 1000) 
adonis(BC.dist ~ lina, data = samdf_manual, permutations = 1000) 
adonis(BC.dist ~ ap, data = samdf_manual, permutations = 1000)
adonis(BC.dist ~ control, data = samdf_manual, permutations = 1000)
adonis(BC.dist ~ colony*compound_dose, data = samdf_manual, permutations = 1000)

#permanova using phyloseq distance
psbray <- phyloseq::distance(ps, method="bray")
adonis(psbray ~ compound, data=samdf)

##rarefaction - vegan
rarecurve(OTUdf,step=100, lwd=2, ylab="ASVs", label=T)
