#Import libraries

2
library(tidyverse)

3
library(phyloseq)

4
library(ape) 

5
library(vegan)

6
library(picante)

7
### Load all of the necessary files to make phyloseq object ###

8
#recreate phlyoseq object with new metadata containing binned depression and anxiety scores

9
metafp <- "parkinsons_metadata_new_edited.csv"

10
meta <- read_delim(metafp, delim=",")

11

12
otufp <- "feature-table.txt"

13
otu <- read_delim(file = otufp, delim="\t", skip=1)

14

15
taxfp <- "taxonomy.tsv"

16
tax <- read_delim(taxfp, delim="\t")

17

18
phylotreefp <- "tree.nwk"

19
phylotree <- read.tree(phylotreefp)

20

21
## Format OTU Table ##

22
otu_mat <- as.matrix(otu[,-1])

23

24
# Make first column (#OTU ID) the rownames of the new matrix

25
rownames(otu_mat) <- otu$`#OTU ID`

26

27
# Use the "otu_table" function to make an OTU table

28
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

29

30
## Format metadata ##

31
# Save everything except sampleid as new data frame

32
samp_df <- as.data.frame(meta[,-2])

33

34
# Make sampleids the rownames

35
rownames(samp_df)<- meta$'X.SampleID'

36

37
# Make phyloseq sample data with sample_data() function

38
SAMP <- sample_data(samp_df)

39

40
## Format taxonomy ##

41
# Convert taxon strings to a table with separate taxa rank columns

42
tax_mat <- tax %>% select(-Confidence)%>%
  
43
separate(col=Taxon, sep=";"
         
         44
         , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  
  45
as.matrix() 

46

47
# Save everything except feature IDs 

48
tax_mat <- tax_mat[,-1]

49

50
# Make sampleids the rownames

51
rownames(tax_mat) <- tax$`Feature ID`

52

53
# Make taxa table

54
TAX <- tax_table(tax_mat)

55

56
## Create phyloseq object ##

57
parkinsons <- phyloseq(OTU, SAMP, TAX, phylotree)

58

59
#check objects

60
otu_table(parkinsons)

61
sample_data(parkinsons)

62
tax_table(parkinsons)

63
phy_tree(parkinsons)

64

65

66
######### ANALYZE ##########

67
# Remove non-bacterial sequences, if any

68
parkinsons_filt <- subset_taxa(parkinsons,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

69

70
# Remove samples with less than 100 reads

71
parkinsons_filt_nolow_samps <- prune_samples(sample_sums(parkinsons_filt)>100, parkinsons_filt)

72

73

74
### RAREFY###

75
rarecurve(t(as.data.frame(otu_table(parkinsons_filt_nolow_samps))), cex=0.1)

76
parkinsons_rare <- rarefy_even_depth(parkinsons_filt_nolow_samps, rngseed = 1, sample.size = 3797)

77

78
#to save rarefied phyloseq object, filtered phyloseq object 

79
save(parkinsons_filt_nolow_samps, file="parkinsons2_filt_nolow_samps.RData")

80
save(parkinsons_rare, file="parkinsons_edited_rare.RData")

81

82
View(sample_data(parkinsons_rare))

83

84
#create specific phyloseq object for each aim

85
# Remove samples where anxiety is na

86
parkinsons_final_anxiety <- subset_samples(parkinsons_rare, !is.na(anxiety_binned))

87
View(sample_data(parkinsons_final_anxiety))

88
save(parkinsons_final_anxiety, file="parkinsons_final_anxiety.RData")

89

90
#remove samples where sleep is na

91
parkinsons_final_sleep <- subset_samples(parkinsons_rare, !is.na(Sleep_problems))

92
View(sample_data(parkinsons_final_sleep))

93
save(parkinsons_final_sleep, file="parkinsons_final_sleep.RData")

94

95
#remove samples where depression is na

96
parkinsons_final_depression <- subset_samples(parkinsons_rare, !is.na(depression_binned))

97
View(sample_data(parkinsons_final_depression))

98
save(parkinsons_final_depression, file="parkinsons_final_depression.RData")