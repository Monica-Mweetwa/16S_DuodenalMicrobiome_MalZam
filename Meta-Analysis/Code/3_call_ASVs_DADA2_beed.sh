#!/bin/bash

################################################
# 3. Perform ASV calling using DADA2
################################################

# Description:
#   Using the DADA2 R package, merge and perform ASV calling on the read pairs that have 
#   survived preprocessing (primer matching/removal, quality-based filtering), and create an output
#   file containing the ASV table (encoded as an R object) 

export DADA2_LOGS_BASENAME='3_call_ASVs_DADA2'
export LOGS_DIR=${WORK_DIR}/bbmap_beed/logs

STEP_LOGS_DIR=${LOGS_DIR}/${DADA2_LOGS_BASENAME}
if [ ! -d ${STEP_LOGS_DIR} ]; then mkdir -p ${STEP_LOGS_DIR}; fi

# Set the path for a locally installed R version, if needed
# Update the path to point to your locally installed R version or use system R
export PATH="/usr/local/bin:$PATH"

# Ensure the correct R version is loaded
R_VERSION="4.3.2"
if ! Rscript --version | grep -q $R_VERSION; then
    echo "Warning: Installed R version does not match required version ($R_VERSION)."
fi

#execute calling ASVs

./3_call_ASVs_DADA2_beed.R

echo "ASV calling done"

