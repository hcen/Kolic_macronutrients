# Module:Phenotype relationships
# Oct 12, 2022
# Jessica Ewald

library(WGCNA)
library(dplyr)
library(ggplot2)

library(RSQLite)
library(pheatmap)
library(fgsea)

#### Module analysis - read in dataT and mods from previous step. Previous step contains some randomness 
# in color selection (although produces relatively stable functional modules). 
dataT <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/dataT.rds")
mods <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/mods.rds")
meta.num <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/meta_num.rds")
modules <- mods$module

# re-order meta to match
meta.num <- meta.num[match(rownames(dataT),rownames(meta.num)),]

# eigengene analysis
# compute the eigengenes of the fold changes
MEs <- moduleEigengenes(dataT, colors = modules)$eigengenes
MEs = orderMEs(MEs)
moduleCor = WGCNA::cor(MEs, meta.num, use = "p")
modulePvalue = corPvalueStudent(moduleCor, dim(dataT)[1])

# adjust p-values
modulePvalue.adj <- p.adjust(modulePvalue, method = "fdr")
modulePvalue.adj <- matrix(modulePvalue.adj, nrow = dim(modulePvalue)[1], ncol = dim(modulePvalue)[2])
rownames(modulePvalue.adj) <- rownames(modulePvalue)
colnames(modulePvalue.adj) <- colnames(modulePvalue)

# re-order so everything hierarchically clustered
row.col.order <- pheatmap(moduleCor, silent = TRUE)
moduleCor <- moduleCor[row.col.order$tree_row$order, row.col.order$tree_col$order]
modulePvalue.adj <- modulePvalue.adj[row.col.order$tree_row$order, row.col.order$tree_col$order]

# visualize AO results
# Will display correlations and their p-values
textMatrix = paste(signif(moduleCor, 2), "\n(",
                   signif(modulePvalue.adj, 1), ")", sep = "")
dim(textMatrix) = dim(moduleCor)

# Display the correlation values within a heatmap plot
png(filename = "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/results/module_phen_all.png", 
width = 15, height = 10, units = "in", res = 300)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation val
labeledHeatmap(Matrix = moduleCor,
               xLabels = colnames(moduleCor),
               yLabels = rownames(moduleCor),
               ySymbols = rownames(moduleCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("All Modules"))
dev.off()


# look at important module:phenotype relationships
meta <- data.frame(meta.num)

plotEigenPhen <- function(ME, phen, modName, phenName){
  temp <- data.frame(Module = ME, Phenotype = phen)
  if(is.factor(phen)){
    ggplot(temp, aes(x = Module, y = Phenotype)) +
      geom_jitter(width=0.15, color='grey') +
      geom_boxplot(fill=NA) +
      theme_classic(base_size = 15) +
      xlab(modName) +
      ylab(phenName)
  } else {
    ggplot(temp, aes(x = Module, y = Phenotype)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      theme_classic(base_size = 15) +
      xlab(modName) +
      ylab(phenName)
  }
}

plotEigenPhen(MEs$MEgreenyellow, meta$AUC.3mm.glucose, "Light yellow", "AUC 3mm glucose")

### heatmap of eigengenes labelled with meta-data (T2D, olp HR, Purity 'eye')
# clustering
col.ann <- meta.num[,c("T2D", "olp.high.resp")] %>% as.data.frame()
col.ann <- col.ann[order(col.ann$olp.high.resp, decreasing = TRUE),]
ME.hm <- MEs[match(rownames(col.ann), rownames(MEs)), ] %>% t()
rownames(ME.hm) <- gsub("ME", "", rownames(ME.hm))
pheatmap(ME.hm, show_colnames = FALSE, annotation_col = col.ann, cluster_cols = FALSE)

png(filename = "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/results/prot_eigengene_heatmap.png", 
    width = 15, height = 10, units = "in", res = 300)
pheatmap(ME.hm, show_colnames = FALSE, annotation_col = col.ann, cluster_cols = FALSE)
dev.off()

