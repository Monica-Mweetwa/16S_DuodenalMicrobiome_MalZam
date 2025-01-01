#!/usr/bin/env -S Rscript --vanilla

# Load libraries

library("Rcpp")
library("dada2")
library("ggplot2")
library("stringi")
library("optparse")
library("dplyr")

.libPaths()

sI <- sessionInfo()
print(sI, locale=F)

# Set random seed for reproducibility purposes
set.seed(100)

# Global variables for file naming
datestring <- format(Sys.time(), "%y%m%d_%H%M") # Example: 230314_1418
experiment <- gsub("^.*/", "", getwd()) #get from working directory

# Arguments ---------------------------------------------------------------

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, help="sample name", dest="sample"),
  make_option("--truncLenR1", type="integer", default=200, help="Trim length for R1 sequences", dest="truncLenR1"),
  make_option("--truncLenR2", type="integer", default=150, help="Trim length for R2 sequences", dest="truncLenR2"),
  make_option("--maxEER1", type="integer", default=2, help="After truncation, R1 sequences with more than this number of expected errors will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))", dest="maxEER1"),
  make_option("--maxEER2", type="integer", default=2, help="After truncation, R2 sequences with more than this number of expected errors will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))", dest="maxEER2")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (sample file)\n", call.=FALSE)
}

# TODO: Path definition below should be passed through as an option when dada2_filter.R is run (using the ${PAIRED_DIR} value 
#       defined in config.sh), rather than being hard-coded
# Path for properly-paired, demultiplexed data
path <- paste(getwd(), "bbmap_seem/paired", sep = "/")

# identify R1 and R2 files (NOTE: the "starts with" carrot (^) is required to make sure the right sample is selected
print("Finding R1 and R2 files for processing...")
fnFs <- sort(list.files(path, pattern=paste("^", opt$sample, "_R1_(515F|806R).matched.trimmed.paired.fastq.gz", sep = ""), full.names = TRUE))
fnFs <- fnFs[file.size(fnFs) > 100] # removes files with < 70 bytes file size. These files were originally empty in fastq format but increase to 69bytes after compression.
fnRs <- sort(list.files(path, pattern=paste("^", opt$sample, "_R2_(515F|806R).matched.trimmed.paired.fastq.gz", sep = ""), full.names = TRUE))
fnRs <- fnRs[file.size(fnRs) > 100] # removes files with < 70 bytes file size. These files were originally empty in fastq format but increase to 69bytes after compression.

if (length(fnFs) == 0 || length(fnRs) == 0) {
  print(paste(opt$sample, "one or more input files contain no data after preprocessing.", sep = "    "))
  quit()
}

# Extract sample names
# To-Do: Talk to Matt about why this is necessary; will there ever be situations where the extracted sample name will not simply match
#        the sample name passed via the options sent to .R script? Did this have something to do with the commented out make_option()
#        variation in which file names might have additions to allow for study designations?
sample_names <- intersect(unique(sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)), unique(sapply(strsplit(basename(fnRs), "_R2"), `[`, 1)))

# Set filtered file location
filt_path <- paste(getwd(), "bbmap_seem/filtered", sep = "/")

# Set up data structure to capture filtering performance
#out <- data.frame(reads.in = numeric(), reads.out = numeric())
out <- data.frame(reads.in = numeric(),
                  reads.out = numeric(),
                  reads.removed = numeric(),
                  pct.filter.removed = numeric(),
                  pct.filter.retained = numeric())
                  
print("Performing filtering...")
# Will there ever be more than one sample name stored in sample_names? Is a loop needed here?
for (i in 1:length(sample_names)) {
  sample <- sample_names[i]
  print(paste("Processing sample ", i, "/", length(sample_names), sep = ""))
  print("Processing PS1 orientation...")
  if (file.exists(file.path(path, paste0(sample, "_R1_515F.matched.trimmed.paired.fastq.gz"))) &
      (file.size(file.path(path, paste0(sample, "_R1_515F.matched.trimmed.paired.fastq.gz"))) > 0 )) {
    out.1 <- filterAndTrim(file.path(path, paste0(sample, "_R1_515F.matched.trimmed.paired.fastq.gz")),
                           file.path(filt_path, paste0(sample, "_R1_515F.matched.trimmed.paired.filtered.fastq.gz")),
                           file.path(path, paste0(sample, "_R2_806R.matched.trimmed.paired.fastq.gz")),
                           file.path(filt_path, paste0(sample, "_R2_806R.matched.trimmed.paired.filtered.fastq.gz")),
                           truncLen = c(opt$truncLenR1, opt$truncLenR2),
                           maxN = 0,
                           maxEE = c(opt$maxEER1, opt$maxEER2),
                           truncQ = 2,
                           rm.phix=TRUE,
                           compress=FALSE,
                           multithread=TRUE)
    row.names(out.1) <- paste(sample, "R1_515F_R2_806R", sep = "_")
  } else {
    print("No data in PS1 orientation...")
    out.1 <- matrix(data = c(NA, NA), nrow = 1, ncol = 2, dimnames = list(paste(sample, "R1_515F_R2_806R", sep = "_"), c("reads.in", "reads.out")))
  }
  out.1.df <- as.data.frame(out.1)
  out.1.df <- out.1.df %>% mutate(reads.removed = reads.in - reads.out, 
                                  pct.filter.removed = reads.removed/reads.in, 
                                  pct.filter.retained = reads.out/reads.in)
  
  print("Processing PS2 orientation...")
  if (file.exists(file.path(path, paste0(sample, "_R1_806R.matched.trimmed.paired.fastq.gz"))) &
      (file.size(file.path(path, paste0(sample, "_R1_806R.matched.trimmed.paired.fastq.gz"))) > 0 )) {
    out.2 <- filterAndTrim(file.path(path, paste0(sample, "_R1_806R.matched.trimmed.paired.fastq.gz")), 
                           file.path(filt_path, paste0(sample, "_R1_806R.matched.trimmed.paired.filtered.fastq.gz")),
                           file.path(path, paste0(sample, "_R2_515F.matched.trimmed.paired.fastq.gz")),
                           file.path(filt_path, paste0(sample, "_R2_515F.matched.trimmed.paired.filtered.fastq.gz")),
                           truncLen = c(opt$truncLenR1, opt$truncLenR2),
                           maxN = 0, 
                           maxEE = c(opt$maxEER1, opt$maxEER2), 
                           truncQ = 2, 
                           rm.phix=TRUE,
                           compress=FALSE, 
                           multithread=TRUE)  
    row.names(out.2) <- paste(sample, "R1_806R_R2_515F", sep = "_")
  } else {
    print("No data in PS2 orientation...")
    out.2 <- matrix(data = c(NA, NA), nrow = 1, ncol = 2, dimnames = list(paste(sample, "R1_806R_R2_515F", sep = "_"), c("reads.in", "reads.out")))
  }
  
  # Create df from matrix for easier addition of new columns
  out.2.df <- as.data.frame(out.2)
  out.2.df <- out.2.df %>% mutate(reads.removed = reads.in - reads.out, 
                                  pct.filter.removed = reads.removed/reads.in, 
                                  pct.filter.retained = reads.out/reads.in)

  out <- rbind(out, out.1.df, out.2.df)
}

write.table(out, file.path(filt_path, paste(datestring, experiment, opt$sample, "filtering_stats.txt", sep = "_")), quote = F, sep = "\t", col.names=NA)
