# Proteomics processing
# Jessica Ewald
# From the October 12, 2022 script; adapted Feb 13, 2023 for updated dataset

library(WGCNA)
library(dplyr)
library(ggplot2)
library(imp4p)

# RF works well for proteomics data: https://www.nature.com/articles/s41598-021-81279-4#:~:text=Missing%20value%20imputation%20is%20a,and%20disadvantages%20of%20each%20method
# imp4p requires log-normalized data (takes a long time to run)
# imp4p RF just a wrapper for 'MissForest' package: https://academic.oup.com/bioinformatics/article/28/1/112/219101

# get proteomics data
prot <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/abundance.csv")
rownames(prot) <- prot$SampleID
prot <- prot[,-1]

# remove proteins with >50% missing values
prot <- prot[goodSamplesGenes(t(prot))$goodGenes, ]

# normalize
prot <- as.matrix(prot)
rand.inds <- sample(1:dim(prot)[2], 25)
boxplot(prot[,rand.inds], main = "Raw")

prot <- apply(prot, 2, function(x){x/median(x, na.rm=T)});
boxplot(prot[,rand.inds], main = "Sample Median")

prot[!is.na(prot)] <- log10(prot[!is.na(prot)])
boxplot(prot[,rand.inds], main = "Log transformed")

prot <- impute.RF(prot, conditions = factor(rep("1", dim(prot)[2])), verbose = TRUE)
boxplot(prot[,rand.inds], main = "Imputed")

# write out results
write.csv(prot, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/abundance_imputedRF.csv")



