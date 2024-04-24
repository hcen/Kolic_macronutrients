# Process metadata for co-expression analysis
# Jessica Ewald
# From October 12, 2022 script

# note - must be run after the proc_proteomics script: relies on list of donors with proteomics data

library(dplyr)
library(ggplot2)

# read in metadata 
prot <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/abundance_imputedRF.csv")
meta <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/metadata.csv")
pd.ct <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/pd_culture_time.csv")

# merge together
meta <- merge(meta, pd.ct, by = "Donor.ID", all.x = TRUE, all.y = FALSE)

meta <- meta[meta$Donor.ID %in% colnames(prot),]

rm(prot)

# combine KCl
kcl <- meta[ ,c("AUC.30mm.KCl.gluc", "AUC.30mm.KCl.leu", "AUC.30mm.KCl.olp")]
kcl.mean <- apply(kcl, 1, function(x){mean(x, na.rm = TRUE)})
meta$kcl.mean <- kcl.mean

# remove GSIS, peak (re:discussion with Jelena & Jim)
meta <- meta[,!(colnames(meta) %in% c("GSIS", "Peak.15mm.gluc", "Peak.5mm.leu", "Peak.1_5mm.olp", "AUC.30mm.KCl.gluc", "AUC.30mm.KCl.leu", "AUC.30mm.KCl.olp"))]

# add flags for high responders
meta$olp.high.resp <- NA
meta$olp.high.resp[meta$Donor.ID %in% c("R241", "R260", "R353", "R273", "R251", "R279", "R426", "R276", "R396", "R278", "R270")] <- 1
meta$olp.high.resp[!(meta$Donor.ID %in% c("R241", "R260", "R353", "R273", "R251", "R279", "R426", "R276", "R396", "R278", "R270"))] <- 0

# make numerical
rownames(meta) <- meta$Donor.ID
meta <- meta[,-1]
meta$T2D <- meta$T2D %>% as.factor() %>% as.numeric()
meta$Sex <- meta$Sex %>% as.factor() %>% as.numeric()

# re-order columns
col.order <- c("Age", "Sex", "BMI", "T2D", "HbA1c", "pd.culture.time", "Culture.time", "CIT", "Digestion.time", 
               "Purity.percent", "Insulin.content", "IPI", "AUC.3mm.glucose", "AUC.15mm.gluc", "AUC.6mm.gluc", 
               "AUC.5mm.leu", "AUC.6mm.gluc.leu", "AUC.1_5mm.olp", "AUC.6mm.gluc.olp", "kcl.mean", "olp.high.resp")

meta <- meta[,col.order]


# visualize each column looking for roughly normal dist
for(i in c(1:dim(meta)[2])){
  hist(as.numeric(meta[,i]), main = colnames(meta)[i])
}

# transform right-skew variables
meta$HbA1c <- log10(meta$HbA1c)
meta$Insulin.content <- log10(meta$Insulin.content)
meta$IPI <- log10(meta$IPI)
meta$AUC.3mm.glucose <- log10(meta$AUC.3mm.glucose)
meta$AUC.15mm.gluc <- log10(meta$AUC.15mm.gluc)
meta$AUC.6mm.gluc <- log10(meta$AUC.6mm.gluc)
meta$AUC.5mm.leu <- log10(meta$AUC.5mm.leu)
meta$AUC.6mm.gluc.leu <- log10(meta$AUC.6mm.gluc.leu)
meta$AUC.1_5mm.olp <- log10(meta$AUC.1_5mm.olp + abs(min(meta$AUC.1_5mm.olp, na.rm=TRUE)) + 1000)
meta$AUC.6mm.gluc.olp <- log10(meta$AUC.6mm.gluc.olp + abs(min(meta$AUC.6mm.gluc.olp, na.rm=TRUE)) + 1000)
meta$kcl.mean <- log10(meta$kcl.mean + abs(min(meta$kcl.mean, na.rm=TRUE)) + 1000)

# check to make sure I didn't miss any
for(i in c(1:dim(meta)[2])){
  hist(meta[,i], main = colnames(meta)[i])
}

# save as numerical matrix
meta.num <- as.matrix(meta)
saveRDS(meta.num, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/meta_num.rds")







