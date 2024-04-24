# Compute Module Membership
# Jessica Ewald
# February 14, 2023

### get hub gene TOM stats
library(WGCNA)
library(dplyr)
library(reshape2)
library(RSQLite)

dataT <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/dataT.rds")
mods <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/mods.rds")
meta.num <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/data/meta_num.rds")
modules <- mods$module

# re-order meta to match
meta.num <- meta.num[match(rownames(dataT),rownames(meta.num)),]

res <- chooseTopHubInEachModule(dataT, mods$module, omitColors = "grey", type = "signed")

MEs <- moduleEigengenes(dataT, colors = modules)$eigengenes
modNames <- substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(dataT, MEs, use = "p"));
names(geneModuleMembership) <- modNames

MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), dim(dataT)[1]));
names(MMPvalue) <- modNames

geneModuleMembership$uniprotID <- rownames(geneModuleMembership)
MMPvalue$uniprotID <- rownames(MMPvalue)

MM <- reshape2::melt(geneModuleMembership, id.vars = "uniprotID")
colnames(MM) <- c("uniprotID", "module", "module.membership")
MMp <- reshape2::melt(MMPvalue, id.vars = "uniprotID")
colnames(MMp) <- c("uniprotID", "module", "p-value")

MM.df <- merge(MM, MMp, by = c("uniprotID", "module"))

# subset to only show proteins in their modules
MM.df <- MM.df[MM.df$module != "grey", ]
uniq.mods <- unique(modules)[-1]
MM.res <- data.frame()
for(i in c(1:length(uniq.mods))){
  mod <- uniq.mods[i]
  mod.genes <- mods$protID[mods$module == mod]
  df.temp <- MM.df[MM.df$module == mod, ]
  df.temp <- df.temp[df.temp$uniprotID %in% mod.genes, ]
  MM.res <- rbind(MM.res, df.temp)
}


# merge with gene description
# get ID mapping files
con <- dbConnect(RSQLite::SQLite(), "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics/hsa_genes.sqlite")
dbListTables(con)
entrez <- dbReadTable(con, "entrez")
uniprot <- dbReadTable(con, "entrez_uniprot")
dbDisconnect(con)

prot.desc <- merge(entrez, uniprot, by = "gene_id")
colnames(prot.desc) <- c("entrezID", "Official.gene.symbol", "Description", "uniprotID")

MM.res <- merge(MM.res, prot.desc, by = "uniprotID", all.x = TRUE, all.y = FALSE)
MM.res <- MM.res[order(MM.res$`p-value`), ]
MM.res <- MM.res[order(MM.res$module), ]
rownames(MM.res) <- NULL

write.csv(MM.res, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Analysis/islet_proteomics_v2/results/module_membership.csv",
          quote = TRUE, row.names = FALSE)
