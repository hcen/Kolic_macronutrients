# Compute co-expression network & characterize modules
# Jessica Ewald
# October 12, 2022

library(WGCNA)
library(dplyr)
library(ggplot2)
library(RSQLite)
library(pheatmap)
library(fgsea)


# read in imputed data
prot <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/abundance_imputedRF.csv",
                 row.names = 1)
meta <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/Data/metadata.csv")

# remove samples with no meta data
prot <- prot[,colnames(prot) %in% meta$Donor.ID]

# plot PCA
pca <- prcomp(t(prot))
pca.dat <- pca$x %>% data.frame()

ggplot(pca.dat, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_classic(base_size = 15) +
  scale_size(guide = 'none') +
  theme(legend.position = 'bottom', 
        legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5))

##### WGCNA #####

dataT <- t(prot)

# set powers to check
powers = c(seq(from = 2, to=40, by=2))

# check each threshold
sft = pickSoftThreshold(dataT, powerVector = powers)

par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sftThr <- 10 # got this from the above plot - it's fine if not >90%, but must be at "elbow"

# compute network
adj = adjacency(dataT, power = sftThr, type = "signed", corFnc = "bicor")
TOM = TOMdist(adj, TOMType = "signed")

# set minimum module size
minModuleSize <- 15

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(TOM), method = "average")

# detect modules
mColor=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                      minClusterSize = minModuleSize, cutHeight = 0.99,
                      deepSplit = ds, distM = TOM)
  mColor=cbind(mColor,labels2colors(tree$labels));
}

plotDendroAndColors(geneTree, mColor, paste("dpSplt =", 0:3), main = "Detected modules", dendroLabels=FALSE);

modules <- mColor[,2]
test <- mergeCloseModules(dataT, modules)
mColor <- cbind(mColor, as.vector(test$colors))
mColor <- mColor[,c(1,2,5,3,4)]
modules <- mColor[,3]
table(modules)

png("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/modules_merge.png", width = 11, height = 8, units = "in", res = 300)
plotDendroAndColors(geneTree, mColor, paste("dpSplt =",c(0,1,1.1,2,3)), main = "Detected modules", dendroLabels=FALSE);
dev.off()

# write out modules & data for reproducibility of downstream steps
mods <- data.frame(protID = colnames(dataT), module = modules)
write.csv(mods, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/modules.csv")
saveRDS(dataT, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/dataT.rds")
saveRDS(mods, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/mods.rds")


#### ORA of modules ####
mods <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/mods.rds")
dataT <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/dataT.rds")

# get ID mapping files
id.map <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics/data/ids.csv")
con <- dbConnect(RSQLite::SQLite(), "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics/hsa_genes.sqlite")
dbListTables(con)
entrez <- dbReadTable(con, "entrez")
uniprot <- dbReadTable(con, "entrez_uniprot")
dbDisconnect(con)

# get pathway files
pw.rds <- c("go_bp", "go_cc", "go_mf", "go_panthbp", "go_panthcc", "go_panthmf", "kegg", "reactome")

# read in all pw files
ListGSC <- list()
gs.names <- data.frame()
for(i in c(1:length(pw.rds))){
  gsl <- readRDS(paste0("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics/functional_libraries/", pw.rds[i], ".rds"))
  ListGSC[[i]] <- gsl[['sets']]
  temp.name <- data.frame(names = gsl$term, IDs = names(gsl[['sets']]))
  gs.names <- rbind(gs.names, temp.name)
}
names(ListGSC) <- pw.rds
gs.names <- distinct(gs.names)
gs.names <- gs.names[!duplicated(gs.names$IDs),]

modules <- unique(mods$module)

# convert lists from uniprot to entrez
mod.entrez <- list()
for(i in c(1:length(modules))){
  ids <- mods$protID[mods$module == modules[i]]
  entrez <- uniprot$gene_id[uniprot$accession %in% ids]
  mod.entrez[[i]] <- entrez
}
names(mod.entrez) <- modules

universe <- unique(uniprot$gene_id[uniprot$accession %in% colnames(dataT)])

# conduct ORA
ora.res <- data.frame()
for(i in c(1:length(modules))){
  
  hits <- mod.entrez[[i]]
  
  # create gsca object
  ora.temp <- data.frame()
  for(j in c(1:length(pw.rds))){
    res <- fora(ListGSC[[pw.rds[j]]], hits, universe)
    res$library <- pw.rds[[j]]
    ora.temp <- rbind(ora.temp, res)
  }
  ora.temp$module <- modules[i]
  ora.res <- rbind(ora.res, ora.temp)
  
}

saveRDS(ora.res, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/ora_all.rds")

# get top results
ora.res <- ora.res[ora.res$overlap > 0,]
ora.res <- ora.res[ora.res$pval < 0.05, ]
ora.res <- merge(gs.names, ora.res, by.x = "IDs", by.y = "pathway")

ora.sig <- ora.res[ora.res$padj < 0.05, ]

# write out results
saveRDS(ora.sig, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/ora_sig.rds")
ora.sig$overlapGenes <- as.character(ora.sig$overlapGenes)
write.csv(ora.sig, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/ora_sig.csv",
          quote = TRUE, row.names = FALSE)
