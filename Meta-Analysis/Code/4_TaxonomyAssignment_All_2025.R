#title: "Assigning Taxonomy to all samples used in meta-analysis"
#author: Monica
#Last updated: 01/01/2025
#*Adapted from https://benjjneb.github.io/dada2/tutorial.html*
  
rm(list=ls())
library(dada2); packageVersion("dada2")
setwd("~/Dropbox/BEECH_16s_Analysis/Meta-Analysis") #Mac
#setwd("C:/Users/Monica/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Windows

#Importing abundance tables 
seqtab.agg_seem <- readRDS("Data/RData/seqtab.agg_seem.RDS")
seqtab.agg_beed <- readRDS("Data/RData/seqtab.agg_beed.RDS")

#merge BEED and SEEM and remove chimeras. This step was already done for the AFRIBIOTA dataset
st.bs <- mergeSequenceTables(seqtab.agg_seem, seqtab.agg_beed) # dim 619 7181
# Check orientations (looking for all-TACG variants, not CCTG variants)
table(substr(colnames(st.bs), start = 0, stop = 4))
#AACA AACG CACC CACG CTTC GACA GACG GGCG TACA TACC TACG TACT TATG TGCG TTCG 
#2   65    1   36    1    7   23    1  257    1 7565    2    3    2    1 

# Remove chimera
st.bs.int <- apply(st.bs, MARGIN = c(1, 2), function(x) as.integer(x))
seqtab.nochim <- removeBimeraDenovo(st.bs.int, method="consensus", multithread=TRUE, verbose = TRUE)

#merge with AFRIBIOTA and Zambia dataset
load("Data/RData/ExtDat_Afr_CountTab.RData") #seqtab.nochim.afr
load("Data/RData/seqtab.agg_Zam.RDS") #asv_final
seqtab.nochim_all <- mergeSequenceTables(seqtab.nochim, seqtab.nochim.afr, as.matrix(asv_final))
save(seqtab.nochim_all, file="Data/RData/All_CountTab2025.RData") #Count table

#Assign Taxonomy
load("Data/RData/All_CountTab2025.RData") #seqtab.nochim.afr
taxa.new <- assignTaxonomy(seqs = seqtab.nochim_all, refFasta = "Data/Tax/silva_nr99_v138.1_train_set.fa", tryRC = T, multithread=F)
taxa.ext_all <- addSpecies(taxtab = taxa.new, refFasta =  "Data/Tax/silva_species_assignment_v138.1.fa", tryRC = T)
save(taxa.ext_all, file="~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis/Data/RData/All_TaxTab2025.RData") #Taxonomy table

#BEED - Assign using Greegenes database
seqtab.agg_beed <- readRDS("Data/RData/seqtab.agg_beed.RDS")
st.beed <- apply(seqtab.agg_beed, MARGIN = c(1, 2), function(x) as.integer(x))
seqtab.nochim.beed <- removeBimeraDenovo(st.beed, method="consensus", multithread=TRUE, verbose = TRUE)
save(seqtab.nochim.beed, file="Data/RData/ExtDat_BEED_CountTab2025_Greengenes.RData") #Count table
taxa.beed_silva <- assignTaxonomy(seqs = seqtab.nochim.beed, refFasta = "Data/Tax/gg_13_8_train_set_97.fa.gz", tryRC = T, multithread=F)
save(taxa.beed_silva, file="Data/RData/EXtDat_BEED_TaxTab_Greengenes.RData") #Taxonomy table

#DONE!