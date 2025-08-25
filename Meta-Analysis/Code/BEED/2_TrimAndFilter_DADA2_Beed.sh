#!/bin/bash

#############################################
# Perform DADA2-appropriate read filtering
#############################################

# Description:
# This script trims and filters R1 and R2 sequences using the DADA2 R package. It calculates 
# the number of tasks based on the mapping file and directly executes the R script with 
# the specified parameters for local execution.

# Load configuration variables
WORK_DIR=$( pwd )
export WORK_DIR
export LOGS_DIR=${WORK_DIR}/bbmap_beed/logs
export FILTER_LOGS_BASENAME='2_trim_and_filter_DADA2'
export FILTER_TRUNCLEN_R1=200
export FILTER_TRUNCLEN_R2=150
export FILTER_MAXEE_R1=2
export FILTER_MAXEE_R2=2
export PAIRED_DIR=${WORK_DIR}/bbmap_beed/paired       # This is not yet customizable; do not change
export FILTERED_DIR=${WORK_DIR}/bbmap_beed/filtered   # This is not yet customizable (or used); do not change
export MAPPING_FILE_PATH=${WORK_DIR}/bbmap_beed/mapping_file_beed.txt

# Count non-header lines in the mapping file
JOB_COUNT=$(grep "^[^#]" ${MAPPING_FILE_PATH} | wc -l | awk '{print $1}')

# Set up logs directory
STEP_LOGS_DIR=${LOGS_DIR}/${FILTER_LOGS_BASENAME}
if [ ! -d "${STEP_LOGS_DIR}" ]; then
    mkdir -p "${STEP_LOGS_DIR}"
fi

set -e

# Set the path for a locally installed R version, if needed
# Update the path to point to your locally installed R version or use system R
export PATH="/usr/local/bin:$PATH"

# Ensure the correct R version is loaded
R_VERSION="4.3.2"
if ! Rscript --version | grep -q $R_VERSION; then
    echo "Warning: Installed R version does not match required version ($R_VERSION)."
fi

# Loop through each sample and process
for TASK_ID in $(seq 1 ${JOB_COUNT}); do
    SAMPLEID=$(grep "^[^#]" ${MAPPING_FILE_PATH} | sed -n ${TASK_ID}p | cut -f 1)
    echo "Processing Sample: ${SAMPLEID}"

    # Call the R script for trimming and filtering
    ./2_trim_and_filter_DADA2_beed.R \
        -s "${SAMPLEID}" \
        --truncLenR1 ${FILTER_TRUNCLEN_R1} \
        --truncLenR2 ${FILTER_TRUNCLEN_R2} \
        --maxEER1 ${FILTER_MAXEE_R1} \
        --maxEER2 ${FILTER_MAXEE_R2} \
        > "${STEP_LOGS_DIR}/${FILTER_LOGS_BASENAME}_${TASK_ID}.txt" 2>&1

    if [ $? -eq 0 ]; then
        echo "Sample ${SAMPLEID} processed successfully."
    else
        echo "Error processing Sample ${SAMPLEID}. Check logs for details."
    fi
done

echo "All samples processed."