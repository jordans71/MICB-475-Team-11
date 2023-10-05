# MICB-475-Team-11

## Thursday October 5th MICB 475 Meeting Agenda
### Research Question
After reading the PD paper in greater detail, we realized that [Cristea et al. (2020)](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28052)  already did PICRUSt2 analysis [(supplementary figure 4)](https://movementdisorders.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fmds.28052&file=mds28052-sup-0002-FigureS1.pdf). Control microbiota was enriched for dietary carbohydrate degradation, PD microbiota rich in nucleic acid degradation & amino acid degradation 

Below are some of the further research questions we have in mind: 
* Sub-clustering of metabolic pathways using PD microbiota data only. The paper discussed two clusters – cluster A being more abundant, cluster B less abundant. Do these two clusters correlate with the shift towards proteolytic metabolic pathways? Literature suggests the PD microbiota is more “diverse” in that it has a high rarity index, but the presence of certain species does not necessarily translate to having the corresponding functions. Our question is if this functional metabolic shift is due to cluster A & cluster B as identified in microbial composition 
* If we identified sub-clusters based on functional data, are there any correlations between these clusters to other metadata categories e.g. patient age, depression, anxiety, BMI, drug usage? 
* How is lipid metabolism different in PD vs. controls? Is there a difference? PIECRUSt2 figure only displays carbohydrate, nucleic acid, and protein pathways 
* Significance of results – exactly how big of a difference is there?
* What would be our future directions if all the questions above come back with negative answers? I.e. no correlation, insignificant results 

### Serum Metabolomics data format
What can we do with serum metabolomics? Is the format of the data such that functional analysis can be run (would not be possible with 16S region)? 

### Data processing
When reproducing figures for the project proposal, should we use code in the original paper? Do so exactly as the MICB 475 QIIME2 modules have outlined except for DADA2 denoising/clustering parameters? 
We have not started data processing yet 
### Is there any potential to change our meeting time? 
An alternative time all group members are available would be Wednesday 1PM, but we understand if this change cannot be accommodate.

### GitHub and R script
Created R script and stored in GitHub

## October 5 Team Meeting Notes
### Housekeeping
Meeting time changed to Wednesdays @ 1PM

### Research Question and Aims

5 AIMS:
*  Parkinson's vs control looking at lipid pathways
*  PD only, mental health factors (how does microbiome change in PD patients with depression/anxiety?)
*  PD only, drug usage factors (how does microbiome change in PD patients that use drugs?)
*  PD only, BMI factors (how does microbiome change in PD patients with different BMI values?)
*  RNA seq data comparison to metabolic pathways
      * split up aims between group members

If research questions come back with negative answers, not an issue. There will be things to report, especially if we use RNA Seq data.

### Code
Use the researcher's code so that our analysis later on are validated and more concrete
* replicate only to 2 cluster image, also want sufficient amounts of novel investigation/analysis

PiCRUST analysis:
* different pipeline on command line from QIIME2
* different set of commands
* need to outline steps we are going through for data processing

PiCRUST outputs:
* 3 tables given
  * 1 represents metabolic pathways
  * 1 represents enzymes
    * can use this to assign to proteins rather than pathways if we want
  * 1 represents metagenomic data

### Serum metabolomics
We should not work with this data; it is very hard to work with, the dataset is not well annotated, and it requires complicated software

### For Next Week
Start on proposal (main focus for next week)
* split up work between group members to focus on different parts of the proposal
* look at sample on canvas and follow their format
* specify we are doing processing with QIIME and PiCRUST
* Sample Timeline from Chris:
  * Week 1/2 = QIIME and PiCRUST processing
  * Week 3 = PD vs Control lipid processing
  * Week 4/5 = look at diversity metrics
* We can send parts of our proposal to Chris for review before it is due 

Can start data processing with QIIME and PiCRUST (inform Chris which 2 individuals are going to do PiCRUST analysis)
* denoise, cluster, find sampling depth, generate .qza files and change format for transfer to R


