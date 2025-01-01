#!/usr/bin/env -S Rscript --vanilla

# Load libraries

library("Rcpp")
library("dada2")
library("stringi")
library("optparse")

sI <- sessionInfo()
print(sI, locale=F)

# Set random seed for reproducibility purposes
set.seed(100)

# Global variables for file naming
datestring <- format(Sys.time(), "%y%m%d_%H%M")
#experiment <- gsub("^.*/", "", getwd()) #get from working directory

# Local functions

## Reverse complement
rc <- function(seq) {
  return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", seq)))
}

## Pre-check merging
merge_check_single <- function(dadaF, derepF, dadaR, derepR) {
  if (is(dadaF, "dada")) 
    dadaF <- list(dadaF)
  if (is(derepF, "derep")) 
    derepF <- list(derepF)
  if (is(dadaR, "dada")) 
    dadaR <- list(dadaR)
  if (is(derepR, "derep")) 
    derepR <- list(derepR)
  if (!(dada2:::is.list.of(dadaF, "dada") && dada2:::is.list.of(derepF, "derep") && 
        dada2:::is.list.of(dadaR, "dada") && dada2:::is.list.of(derepR, "derep"))) {
    stop("This function requires dada-class and derep-class input arguments.")
  }
  nrecs <- c(length(dadaF), length(derepF), length(dadaR), 
             length(derepR))
  if (length(unique(nrecs)) > 1) 
    stop("The dadaF/derepF/dadaR/derepR arguments must be the same length.")
  rval <- list()
  
  for (i in seq_along(dadaF)) {
    mapF <- derepF[[i]]$map
    mapR <- derepR[[i]]$map
    if (!(is.integer(mapF) && is.integer(mapR))) 
      stop("Incorrect format of $map in derep-class arguments.")
    if (!(length(mapF) == length(mapR) && max(mapF) == length(dadaF[[i]]$map) && 
          max(mapR) == length(dadaR[[i]]$map))) {
      stop("Non-corresponding maps and dada-outputs.")
    }
    rF <- dadaF[[i]]$map[mapF]
    rR <- dadaR[[i]]$map[mapR]
    pairdf <- data.frame(sequence = "", abundance = 0, forward = rF, 
                         reverse = rR)
    ups <- unique(pairdf)
    if (nrow(ups[!is.na(ups$forward) & !is.na(ups$reverse),]) > 1) {
      return(TRUE)
      #      print("TRUE")
    } else {
      return(FALSE)
      #      print("FALSE")
    }
  }
}


# Arguments ---------------------------------------------------------------

option_list = list(
  make_option(c("-t", "--study"), type="character", default="", 
              help="addition to file name for tagging purposes [default= %default]", dest="study")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Set filtered file location
filt_path <- paste(getwd(), "bbmap_beed/filtered", sep = "/")

# Identify filtered R1 and R2 files (in case file pairs didn't pass filter)
filtFs <- sort(list.files(filt_path, pattern=".*R1_(515F|806R).matched.trimmed.paired.filtered.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern=".*R2_(515F|806R).matched.trimmed.paired.filtered.fastq.gz", full.names = TRUE))

sample_names2 <- intersect(unique(sapply(strsplit(basename(filtFs), "_R1"), `[`, 1)), unique(sapply(strsplit(basename(filtRs), "_R2"), `[`, 1)))

# Learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

# Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

save.image(file = file.path(getwd(), "bbmap_beed/dada2_beed.Rdata"))

# Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

save.image(file = file.path(getwd(), "bbmap_beed/dada2_beed.Rdata"))

# Merging
seqtab.agg <- matrix(0, dimnames = list("empty", "empty"))

for (i in 1:length(sample_names2)) {
  sample <- sample_names2[i]
  print(paste(sample, " ", i, "/", length(sample_names2), sep = ""))
  # The following logical structure for seqtab 1 and seqtab 2 allows instances where only one orientation exists
  if (!is.null(unlist(dadaFs[paste0(sample, "_R1_515F.matched.trimmed.paired.filtered.fastq.gz")]))) {
    if (merge_check_single(dadaFs[paste0(sample, "_R1_515F.matched.trimmed.paired.filtered.fastq.gz")],
                           derepFs[paste0(sample, "_R1_515F.matched.trimmed.paired.filtered.fastq.gz")],
                           dadaRs[paste0(sample, "_R2_806R.matched.trimmed.paired.filtered.fastq.gz")],
                           derepRs[paste0(sample, "_R2_806R.matched.trimmed.paired.filtered.fastq.gz")])) {
      mergers.1 <- mergePairs(dadaFs[paste0(sample, "_R1_515F.matched.trimmed.paired.filtered.fastq.gz")],
                              derepFs[paste0(sample, "_R1_515F.matched.trimmed.paired.filtered.fastq.gz")],
                              dadaRs[paste0(sample, "_R2_806R.matched.trimmed.paired.filtered.fastq.gz")],
                              derepRs[paste0(sample, "_R2_806R.matched.trimmed.paired.filtered.fastq.gz")],
                              verbose = TRUE)
      seqtab.1 <- makeSequenceTable(mergers.1)
    } else {
      seqtab.1 <- matrix(0, dimnames = list("empty", "empty"))
    }
  } else {
    seqtab.1 <- matrix(0, dimnames = list("empty", "empty"))
  }
  if (!is.null(unlist(dadaFs[paste0(sample, "_R1_806R.matched.trimmed.paired.filtered.fastq.gz")]))) {
    if (merge_check_single(dadaFs[paste0(sample, "_R1_806R.matched.trimmed.paired.filtered.fastq.gz")],
                           derepFs[paste0(sample, "_R1_806R.matched.trimmed.paired.filtered.fastq.gz")],
                           dadaRs[paste0(sample, "_R2_515F.matched.trimmed.paired.filtered.fastq.gz")],
                           derepRs[paste0(sample, "_R2_515F.matched.trimmed.paired.filtered.fastq.gz")])) {
      mergers.2 <- mergePairs(dadaFs[paste0(sample, "_R1_806R.matched.trimmed.paired.filtered.fastq.gz")],
                              derepFs[paste0(sample, "_R1_806R.matched.trimmed.paired.filtered.fastq.gz")],
                              dadaRs[paste0(sample, "_R2_515F.matched.trimmed.paired.filtered.fastq.gz")],
                              derepRs[paste0(sample, "_R2_515F.matched.trimmed.paired.filtered.fastq.gz")],
                              verbose = TRUE)
      seqtab.2 <- makeSequenceTable(mergers.2)
      colnames(seqtab.2) <- rc(colnames(seqtab.2))
    } else {
      seqtab.2 <- matrix(0, dimnames = list("empty", "empty"))
    }
  } else {
    seqtab.2 <- matrix(0, dimnames = list("empty", "empty"))
  }
  seqtab <- t(rowsum(t(cbind(seqtab.1, seqtab.2)), colnames(cbind(seqtab.1, seqtab.2))))
  rownames(seqtab) <- sample
  seqtab.agg <- mergeSequenceTables(seqtab.agg, seqtab)
}

seqtab.agg_beed <- seqtab.agg[rownames(seqtab.agg) != "empty", colnames(seqtab.agg) != "empty"]
#table(nchar(colnames(seqtab)))

# Set seqtab storage location
dir.create(paste(getwd(), "bbmap_beed/seqtabs", sep = "/"))

# Save complete R session in case user would like to inspect intermediate outputs
save.image(file = file.path(getwd(), "bbmap_beed/seqtabs/dada2_beed.Rdata"))

# Save ASV table to serialized file
saveRDS(seqtab.agg_beed, file = paste(getwd(),"bbmap_beed/seqtabs/seqtab.agg_beed.RDS", sep = "/"))
