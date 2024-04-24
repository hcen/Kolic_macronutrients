# Partial correlation analysis: try correst for % purity
# Jessica Ewald
# October 27, 2022

library(dplyr)
library(WGCNA)
library(ppcor)

### Partial correlation analysis
#### Module analysis - read in dataT and mods from previous step. Previous step contains some randomness 
# in color selection (although produces relatively stable functional modules). 
dataT <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/dataT.rds")
mods <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/mods.rds")
meta.num <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/meta_num.rds")
modules <- mods$module

vars.remove <- c("T2D", "pd.culture.time", "olp.high.resp")

meta.num <- meta.num[,!(colnames(meta.num) %in% vars.remove)]

# re-order meta to match
meta.num <- meta.num[match(rownames(dataT),rownames(meta.num)),]

# eigengene analysis
# compute the eigengenes of the fold changes
MEs <- moduleEigengenes(dataT, colors = modules, excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs)

vars.adj <- c("Purity.percent", "Culture.time", "Digestion.time", "CIT")
meta.adj <- meta.num[, vars.adj]
meta.num <- meta.num[,!(colnames(meta.num) %in% vars.adj)]

# function
par.cor <- function(ME, meta){
  temp <- data.frame(ME = ME, meta = meta)
  temp <- cbind(temp, meta.adj)
  temp <- na.omit(temp)
  pcor <- pcor.test(temp$ME,  temp$meta, temp[,-c(1:2)])
  return(pcor)
}

all.r <- data.frame()
all.p <- data.frame()

for(i in c(1:dim(MEs)[2])){
  mod.r <- rep(NA, dim(meta.num)[2])
  mod.p <- rep(NA, dim(meta.num)[2])
  for(k in c(1:dim(meta.num)[2])){
    res <- par.cor(MEs[,i], meta.num[,k])
    mod.r[k] <- res$estimate
    mod.p[k] <- res$p.value
  }
  all.r <- rbind(all.r, mod.r)
  all.p <- rbind(all.p, mod.p)
}

rownames(all.r) <- colnames(MEs)
rownames(all.p) <- colnames(MEs)
colnames(all.r) <- colnames(meta.num)
colnames(all.p) <- colnames(meta.num)

all.r <- as.matrix(all.r)
all.p <- as.matrix(all.p)

# adjust p-values
all.ap <- p.adjust(all.p, method = "fdr")
all.ap <- matrix(all.ap, nrow = dim(all.p)[1], ncol = dim(all.p)[2])
rownames(all.ap) <- rownames(all.p)
colnames(all.ap) <- colnames(all.p)

# re-order so everything hierarchically clustered
row.col.order <- pheatmap(all.r, silent = TRUE)
all.r <- all.r[row.col.order$tree_row$order, ]
all.ap <- all.ap[row.col.order$tree_row$order, ]


# visualize results

# Will display correlations and their p-values
textMatrix = paste(signif(all.r, 2), "\n(",
                   signif(all.ap, 1), ")", sep = "")

# Display the correlation values within a heatmap plot
png(filename = "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/results/module_phen_adj.png", 
width = 6*1.5, height = 5*1.5, units = "in", res = 300)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation val
labeledHeatmap(Matrix = all.r,
               xLabels = colnames(all.r),
               yLabels = rownames(all.r),
               ySymbols = rownames(all.r),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               textAdj = c(0.6, 0.6),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Adjusted for isolation parameters"),
               plotLegend = FALSE)
dev.off()


# make heatmap

# clustering
meta <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/meta_num.rds")
meta <- meta[match(rownames(dataT),rownames(meta)),]

col.ann <- meta[,c("T2D", "olp.high.resp", "Purity.percent")] %>% as.data.frame()
pheatmap(t(MEs), show_colnames = FALSE, annotation_col = col.ann)



