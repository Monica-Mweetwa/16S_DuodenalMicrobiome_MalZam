#title: "Assigning Taxonomy - BEED and SEEM"
#author: Monica
#Last updated: 01/01/2025
#*Adapted from https://benjjneb.github.io/dada2/tutorial.html*
  

library(dada2); packageVersion("dada2")
setwd("~/Dropbox/BEECH_16s_Analysis/Meta-Analysis")


#SEEM
seqtab.agg_seem <- readRDS("Data/seqtab.agg_seem.RDS")
taxa.new <- assignTaxonomy(seqs = seqtab.agg_seem, refFasta = "~/Dropbox/BEECH_16s_Analysis/Data/Tax/silva_nr99_v138.1_train_set.fa.gz", tryRC = T, multithread=F)
taxa.new2 <- addSpecies(taxtab = taxa.new, refFasta =  "~/Dropbox/BEECH_16s_Analysis/Data/Tax/silva_species_assignment_v138.1.fa.gz", tryRC = T)

seqtab.nochim.seem <- seqtab.agg_seem
taxa.seem <- taxa.new2
save(seqtab.nochim.seem, file="~/Dropbox/BEECH_16s_Analysis/Meta-Analysis/Data/ExtDat_SEEM_CountTab.RData") #Count table
save(taxa.seem, file="~/Dropbox/BEECH_16s_Analysis/Meta-Analysis/Data/EXtDat_SEEM_TaxTab.RData") #Taxonomy table


#BEED
seqtab.agg_beed <- readRDS("Data/seqtab.agg_beed.RDS")
taxa.new <- assignTaxonomy(seqs = seqtab.agg_beed, refFasta = "~/Dropbox/BEECH_16s_Analysis/Data/Tax/silva_nr99_v138.1_train_set.fa.gz", tryRC = T, multithread=F)
taxa.new2 <- addSpecies(taxtab = taxa.new, refFasta =  "~/Dropbox/BEECH_16s_Analysis/Data/Tax/silva_species_assignment_v138.1.fa.gz", tryRC = T)

seqtab.nochim.bang <- seqtab.agg_beed
taxa.bang <- taxa.new2
save(seqtab.nochim.bang, file="~/Dropbox/BEECH_16s_Analysis/Meta-Analysis/Data/ExtDat_Bang_CountTab.RData") #Count table
save(taxa.bang, file="~/Dropbox/BEECH_16s_Analysis/Meta-Analysis/Data/EXtDat_Bang_TaxTab.RData") #Taxonomy table

