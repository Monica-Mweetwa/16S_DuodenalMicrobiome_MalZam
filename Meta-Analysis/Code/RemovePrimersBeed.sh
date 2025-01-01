#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

# Load bbmap locally
export PATH="./bbmap@38.63:$PATH" 

# General (project-level)
#########################
WORK_DIR=$( pwd )  # Set the current directory as the working directory
export WORK_DIR

# Define directories
export DEMUX_DIR=${WORK_DIR}/bbmap_beed/demux
export P515F_DIR=${WORK_DIR}/bbmap_beed/515F
export P806R_DIR=${WORK_DIR}/bbmap_beed/806R
export PAIRED_DIR=${WORK_DIR}/bbmap_beed/paired       # Directory for paired reads
export MAPPING_FILE_PATH=${WORK_DIR}/bbmap_beed/mapping_file_beed.txt  # Ensure correct mapping file path

# Create necessary directories if they don't exist
mkdir -p "${DEMUX_DIR}" "${P515F_DIR}" "${P806R_DIR}" "${PAIRED_DIR}"

# 1_remove_primer_seqs_BBTools
##############################
export PREPROCESSING_CPU='2'          # Adjust CPU cores for preprocessing
export PREPROCESSING_MEM='1G'         # Adjust memory for preprocessing
export PREPROCESSING_LOGS_BASENAME='1_remove_primer_seqs_BBTools'

# Define variables manually
SAMPLE_LIST=${1}  # Input file containing sample information
THREADS=4         # Number of threads to use

if [[ ! -f $SAMPLE_LIST ]]; then
    echo "Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

# Loop through each line in the sample list
while IFS=$'\t' read -r SAMPLEID MATCHID; do
    echo "Processing Sample: $SAMPLEID"

    # Read 1 Orientation 1 (515F)
    echo -e "\n\nProcessing R1-PS1 orientation (515F_R1)"
    bbduk.sh \
            in="${DEMUX_DIR}/${SAMPLEID}"_R1.fastq.gz \
            restrictleft=30 \
            k=19 \
            hdist=1 \
            literal=GTGYCAGCMGCCGCGGTAA \
            copyundefined=t \
            outm="${P515F_DIR}/${SAMPLEID}"_R1_515F.matched.fastq \
            out="${P515F_DIR}/${SAMPLEID}"_R1_515F.nomatch.fastq \
            stats="${P515F_DIR}/${SAMPLEID}"_R1_515F.matched.sort_stats \
            -Xmx2g \
            ow=t \;

    bbduk.sh \
            in="${P515F_DIR}/${SAMPLEID}"_R1_515F.matched.fastq \
            restrictleft=30 \
            k=19 \
            ktrim=l \
            hdist=1 \
            literal=GTGYCAGCMGCCGCGGTAA \
            copyundefined=t \
            out="${P515F_DIR}/${SAMPLEID}"_R1_515F.matched.trimmed.fastq \
            stats="${P515F_DIR}/${SAMPLEID}"_R1_515F.matched.trim_stats \
            -Xmx2g \
            ow=t \;

    # Read 1 Orientation 2 (806R)
    echo -e "\n\nProcessing R1-PS2 orientation (806R_R1)"
    bbduk.sh \
            in="${DEMUX_DIR}/${SAMPLEID}"_R1.fastq.gz \
            restrictleft=30 \
            k=20 \
            hdist=1 \
            literal=GGACTACNVGGGTWTCTAAT \
            copyundefined=t \
            outm="${P806R_DIR}/${SAMPLEID}"_R1_806R.matched.fastq \
            out="${P806R_DIR}/${SAMPLEID}"_R1_806R.nomatch.fastq \
            stats="${P806R_DIR}/${SAMPLEID}"_R1_806R.matched.sort_stats \
            -Xmx2g \
            ow=t \;

    bbduk.sh \
            in="${P806R_DIR}/${SAMPLEID}"_R1_806R.matched.fastq \
            restrictleft=30 \
            k=20 \
            ktrim=l \
            hdist=1 \
            literal=GGACTACNVGGGTWTCTAAT \
            copyundefined=t \
            out="${P806R_DIR}/${SAMPLEID}"_R1_806R.matched.trimmed.fastq \
            stats="${P806R_DIR}/${SAMPLEID}"_R1_806R.matched.trim_stats \
            -Xmx2g \
            ow=t \;

    # Read 2 Orientation 2 (515F)
    echo -e "\n\nProcessing R2-PS2 orientation (515F_R2)"
    bbduk.sh \
            in="${DEMUX_DIR}/${SAMPLEID}"_R2.fastq.gz \
            restrictleft=30 \
            k=19 \
            hdist=1 \
            literal=GTGYCAGCMGCCGCGGTAA \
            copyundefined=t \
            outm="${P515F_DIR}/${SAMPLEID}"_R2_515F.matched.fastq \
            out="${P515F_DIR}/${SAMPLEID}"_R2_515F.nomatch.fastq \
            stats="${P515F_DIR}/${SAMPLEID}"_R2_515F.matched.sort_stats \
            -Xmx2g \
            ow=t \;

    bbduk.sh \
            in="${P515F_DIR}/${SAMPLEID}"_R2_515F.matched.fastq \
            restrictleft=30 \
            k=19 \
            ktrim=l \
            hdist=1 \
            literal=GTGYCAGCMGCCGCGGTAA \
            copyundefined=t \
            out="${P515F_DIR}/${SAMPLEID}"_R2_515F.matched.trimmed.fastq \
            stats="${P515F_DIR}/${SAMPLEID}"_R2_515F.matched.trimmed.trim_stats \
            -Xmx2g \
            ow=t \;

    # Read 2 Orientation 1 (806R)
    echo -e "\n\nProcessing R2-PS1 orientation (806R_R2)"
    bbduk.sh \
            in="${DEMUX_DIR}/${SAMPLEID}"_R2.fastq.gz \
            restrictleft=30 \
            k=20 \
            hdist=1 \
            literal=GGACTACNVGGGTWTCTAAT \
            copyundefined=t \
            outm="${P806R_DIR}/${SAMPLEID}"_R2_806R.matched.fastq \
            out="${P806R_DIR}/${SAMPLEID}"_R2_806R.nomatch.fastq \
            stats="${P806R_DIR}/${SAMPLEID}"_R2_806R.matched.sort_stats \
            -Xmx2g \
            ow=t \;

    bbduk.sh \
            in="${P806R_DIR}/${SAMPLEID}"_R2_806R.matched.fastq \
            restrictleft=30 \
            k=20 \
            ktrim=l \
            hdist=1 \
            literal=GGACTACNVGGGTWTCTAAT \
            copyundefined=t \
            out="${P806R_DIR}/${SAMPLEID}"_R2_806R.matched.trimmed.fastq \
            stats="${P806R_DIR}/${SAMPLEID}"_R2_806R.matched.trimmed.trim_stats \
            -Xmx2g \
            ow=t \;

    # Re-pairing both sets of reads after primer-finding
    # Orientation 1
    echo "Checking read pairing for PS1 orientation"
    repair.sh \
            in="${P515F_DIR}/${SAMPLEID}"_R1_515F.matched.trimmed.fastq \
            in2="${P806R_DIR}/${SAMPLEID}"_R2_806R.matched.trimmed.fastq \
            out="${PAIRED_DIR}/${SAMPLEID}"_R1_515F.matched.trimmed.paired.fastq \
            out2="${PAIRED_DIR}/${SAMPLEID}"_R2_806R.matched.trimmed.paired.fastq \
            -Xmx2g \
            ow=t \;

    # Orientation 2
    echo "Checking read pairing for PS2 orientation"
    repair.sh \
            in="${P806R_DIR}/${SAMPLEID}"_R1_806R.matched.trimmed.fastq \
            in2="${P515F_DIR}/${SAMPLEID}"_R2_515F.matched.trimmed.fastq \
            out="${PAIRED_DIR}/${SAMPLEID}"_R1_806R.matched.trimmed.paired.fastq \
            out2="${PAIRED_DIR}/${SAMPLEID}"_R2_515F.matched.trimmed.paired.fastq \
            -Xmx2g \
            ow=t \;

done < <(grep "^[^#]" ${SAMPLE_LIST} | cut -f 1,3)

echo "All samples processed."

#compress fastq files to save space
find bbmap_beed  -type f -name "*.fastq" -exec gzip {} +

echo "All fastq files compressed"