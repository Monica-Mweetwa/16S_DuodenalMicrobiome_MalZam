# 16S_DuodenalMicrobiome_MalZam
This repository contains the code used to analyse 16S RNA sequencing data from duodenal samples of malnourished children in Zambia

## Summary
This repository holds the code for 2 analyses:

* A retrospective cross-sectional study using samples collected during the conduct of previously reported studies on hospitalized children with SAM at the university teaching hospital (Lusaka, Zambia) and stunted children from a high-density community in Lusaka known to have with a high risk of malnutrition and enteropathy.
* A meta-analysis of all publicly available duodenal microbiome 16s data on undernourished children retrieved from SRA Archive.

## Description of Datasets:
### Zambian Datasets:
1. BEECH - Duodenal samples were collected from 53 children with linear growth stunting living in an urban slum in Lusaka, Zambia and used for 16s paired amplicon amplicon sequencing at the Gordon Lab in Washington Uiversity
2. EEBiomarker - Duodenal samples were collected from 24 children with SAM admitted at UTH in Lusaka, Zambia with no clinical explanation for chronic diarrhoea and used for 16s paired amplicon amplicon sequencing at the Gordon Lab in Washington Uiversity
### Additional datasets for Meta-Analysis:
3. BEED	Duodenal samples were collected from 36 children with linear growth stunting, living in an urban slum in Dhaka, Bangladesh, who had biopsy-confirmed EED.
   Project ID PRJEB32184	(16s paired amplicon sequenced at the Gordon Lab in Washington Uiversity)
5. SEEM	Duodenal samples were collected from 43 children with linear growth wasting, living in Pakistan 3-6 months of age.
   Project ID - PRJEB75181	(16s paired amplicon sequenced at the Gordon Lab in Washington Uiversity)
7. Afribiota - Duodenal samples were collected from 46 stunted children aged 2-5 years living in Bangui, Central African Republic and in Antananarivo, Madagascar.
   Project ID - PRJEB27868	(single-ended amplicons sequenced by Microbiome insights)

## SAM vs Stunting Analysis
62% of children with SAM were stunted but stunted children were not severely malnourished i.e severe wasting (WLZ < -3). The ASV counts and taxonomic identification from fastQ files was done by the WashU team. The code deposited here, details the analysis done to relate this data to metdata (clincal features for these children). A guidance document providing deatiled description of this analysis and data used to implement this is provided: **Stunting_SAM/16S_DuodenalMicrobiome_MalZam_GuidanceDocument_Apr.docx**  

## Meta-Analysis of Duodenal Microbiome of stunted children.
Publicly available fastq files (16S) of the duodenal microbiome were retrieved from SRA archive and compared with the Zambian dataset to understand if geography had an impact of the compositon and diversity. This analysis was restried to children with LAZ < -2. 
A guidance document providing deatiled description of this analysis and data used to implement this is provided: **Meta-Analysis/16S_DuodenalMicrobiome_Meta-Analysis_GuidanceDocument_Apr.docx**  
The scripts for filtering primers, and DADA2 implementaton are only modified to run on a local compture instead of HPC. All other parameters were not modified and is synonymous to the original scripts: https://gitlab.com/Gordon_Lab/amplicon_sequencing_pipeline
