

library(tidyverse)
#Import the metadata file 
sampdatFP  <- "parkinsons_metadata.txt"
sampdat <- read.delim(file = sampdatFP, sep = "\t")

#Binning the anxiety level
metadata_new = sampdat %>%
  mutate(anxiety_binned = ifelse(sampdat$STAI_anxiety_score>=80,"High",              #If higher than 80, make it "High"
                                 ifelse(sampdat$STAI_anxiety_score<80,"Low","NA")))        #Else if its lower than 80, make it "Low", and besides that, make it NA.
metadata_new = metadata_new %>%
  mutate(depression_binned = ifelse(metadata_new$BDI_depression_score>=11, "High", 
                                    ifelse(metadata_new$BDI_depression_score<11, "Low","NA")) )
write.csv(file = "parkinsons_metadata_new_edited.csv", metadata_new)

