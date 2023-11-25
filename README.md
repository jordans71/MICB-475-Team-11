# MICB-475-Team-11
edit
## Wednesday October 11th MICB 475 Meeting Agenda
### RNA-sequencing: to proceed or not to with RNA-seq data?
 * Access? If no access, will investigate lipid metabolism

### On Canvas, there is a data processing checklist for the proposal. However, most of it concerns QIIME2 processing. Is there an equivalent for what we would need to accomplish in PICRUSt2? 
* Aiming to get PICRUSt figures ready by the next team meeting. We were busy this week and have yet to start PICRUSt2 analysis

### Proposal figures
* When Evelyn mentioned having at least 1 figure in the proposal, did any of these need to be figures on the data? Or just graphical representations of our research direction (as show in the example proposals)
* If we do need figures made from the dataset itself, should they be original (i.e. new investigations) or just reproducing the paper’s figures? 

### Post data-processing
* QIIME2 data processing started this week. Should processed files be put into Github? The documentation video on canvas mentioned keeping a digital notebook documenting what was done, with links to processed data/code, but did not detail how files should be stored besides this 

### Tips for how MICB 475 workload should be split 
Our current plan right now is to split 1) the proposal by section and 2) coding work by QIIME2 processing + R, PICRSt2, and RNA-Seq analysis (if available) 
Wanted to ask how previous MICB 475 teams splitted workload and what you’d recommend for us to be successful 

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
**Meeting time changed to Wednesdays @ 1PM**

### Research Question and Aims

5 AIMS:
*  Parkinson's vs control looking at lipid pathways
*  PD only, mental health factors (how does microbiome change in PD patients with depression/anxiety?)
*  PD only, drug usage factors (how does microbiome change in PD patients that use drugs?)
*  PD only, BMI factors (how does microbiome change in PD patients with different BMI values?)
*  RNA seq data (focus on human gene) comparison to metabolic pathways
      * split up aims between group members

If research questions come back with negative answers, not an issue. There will be things to report, especially if we use RNA Seq data.

### Code
Use the researcher's code so that our analysis later on are validated and more concrete
* replicate only to 2 cluster image, also want sufficient amounts of novel investigation/analysis

**PiCRUST analysis:**
* different pipeline on the command line from QIIME2
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
* look at the sample proposal on Canvas and follow their format
* specify we are doing processing with QIIME and PiCRUST
* Sample Timeline from Chris:
  * Week 1/2 = QIIME and PiCRUST processing
  * Week 3 = PD vs. Control lipid processing
  * Week 4/5 = Look at diversity metrics
* Provide the Gantt Chart to visualize the timeline
* Divide each aim into pieces and specify the corresponding approaches
* We can send parts of our proposal to Chris for review before it is due 

Can start data processing with QIIME and PiCRUST (inform Chris which 2 individuals are going to do PiCRUST analysis)
* denoise, cluster, find sampling depth, generate .qza files, and change format for transfer to R


## October 11 Team Meeting Notes
### General notes
* Don't have RNA seq data that we thought we did
  * If we did use it as a last aim, it wouldn't be so involved - would be taking results of other papers (lit review) and commenting on similar changes.

### For Next Meeting
* Focus on proposal first
* Don't worry about picrust figures

### Proposal
* Make flow chart of workflow
* Dataset overview
   * Table
      * Might need to have 2 figures (1 for QIIME, 1 for PiCRUST
* Don't need to have diversity metric plots
* Checklist done for QIIME
   * For PiCRUST:
      * Count as an additional check in the overview checklist
      * Just need required files (rep seqs, table)

### Preprocessing Steps
* New column = anxiety score binned
   * have categorical groups high or low anxiety
   * 80 or over = high anxiety
* For depression, can look at 'yes/no' column, or categorize numerical ones
   * should do literature review to understand what scores constitue high/medium/low classification
* Take out less than 5 reads from no mitochondria / no chloroplast table for start with PiCRUST analysis
   * or else PiCRUST runs slowly, and they all get filtered out anyways

### Aims Updates
* Looking at links between PD patients' gut microbiome and metabolic function with respect to depression
   * Bin metadata
   * Diversity metrics
   * PiCRUST
* Looking at links between PD patients' gut microbiome and metabolic function with respect to anxiety
   * Bin metadata
   * Diversity metrics
   * PiCRUST
* Looking at links between PD patients' gut microbiome and metabolic function with respect to fatigue
   * Bin metadata
   * Diversity metrics
   * PiCRUST

### Splitting work
* Can split work based on aims (chosen method) or based on analysis
   * Group members will either be doing PiCRUST and R analysis on Depression aim, or Anxiety and Fatigue aim
  
