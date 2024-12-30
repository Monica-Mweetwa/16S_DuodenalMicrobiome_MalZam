# 16S_DuodenalMicrobiome_MalZam
This repository contains the code used to analyse 16S RNA sequencing data from duodenal samples of malnourished children in Zambia

## Summary
This repository holds the code for 2 analyses:
•	A retrospective cross-sectional study using samples collected during the conduct of previously reported studies on hospitalized children with SAM at the university teaching hospital (Lusaka, Zambia) and stunted children from a high-density community in Lusaka known to have with a high risk of malnutrition and enteropathy. 
•	A meta-analysis of all publicly available duodenal microbiome 16s data on undernourished children retrieved from SRA Archive.

## Description of Datasets:
### Zambian Datasets:
1. BEECH - Duodenal samples were collected from 53 children with linear growth stunting living in an urban slum in Lusaka, Zambia and used for 16s paired amplicon amplicon sequencing at the WashU lab
2. EEBiomarker - Duodenal samples were collected from 24 children with SAM admitted at UTH in Lusaka, Zambia with no clinical explanation for chronic diarrhoea and used for 16s paired amplicon amplicon sequencing at the WashU lab
### Additional datasets for Meta-Analysis:
3. BEED	Duodenal samples were collected from 36 children with linear growth stunting, living in an urban slum in Dhaka, Bangladesh, who had biopsy-confirmed EED.	      Project ID PRJEB32184	(16s paired amplicon sequenced at the WashU lab)
4. SEEM	Duodenal samples were collected from 43 children with linear growth wasting, living in Pakistan 3-6 months of age.	
    Project ID - PRJEB75181	(16s paired amplicon sequenced at the WashU lab)
5. Afribiota - Duodenal samples were collected from 46 stunted children aged 2-5 years living in Bangui, Central African Republic and in Antananarivo, Madagascar.	     Project ID - PRJEB27868	(single-ended amplicons sequenced by  Microbiome insights)

## SAM vs Stunting Analysis
Almost all children with SAM were stunted but stunted children were not servery malnourished i.e severe wasting (WLZ < -2). The ASV counts and taxonomic identification from fastQ files was done by the WashU team. The code deposited here, details the analysis done to relate this data to metdata (clincal features for these children). A guidance document providing deatiled description of this analysis and data used to implement this is provided: **Stunting_SAM/16S_DuodenalMicrobiome_MalZam_GuidanceDocument.docx**  

