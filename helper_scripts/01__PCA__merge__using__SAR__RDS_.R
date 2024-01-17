#!/usr/bin/env Rscript

library(DESeq2)
library(SARTools)
library(dplyr)
#``````````````````````````````````````````````````````````````````````````````
# Load in RDS object ====
args <-  commandArgs(trailingOnly = T)

sar_data <- args[1]
top_var <- args[2]


load(sar_data)

#``````````````````````````````````````````````````````````````````````````````
# Get VSD Assay ====

vsd <- varianceStabilizingTransformation(out.DESeq2$dds, blind=T)

# calculate the variance for each gene
rv <- rowVars(assay(vsd))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(top_var, length(rv)))]

#``````````````````````````````````````````````````````````````````````````````
# Get PCA ====

# perform a PCA on the data in assay(vsd) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

#``````````````````````````````````````````````````````````````````````````````
# Get contribution to the total variance for each component ====

percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
pVar.df <- as.data.frame(percentVar)
pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))

pVar.df = pVar.df[ , order(names(pVar.df))]
pVar.df$percentVar = pVar.df$percentVar * 100
pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)

group = target$group
intgroup.df <- as.data.frame(colData(vsd)[, "group", drop=FALSE])

#``````````````````````````````````````````````````````````````````````````````
# Set Up PCA DataFrames ====
df.pca <- data.frame(pca$x, name=rownames(pca$x))
df.pca.annotated <- left_join(target, df.pca, by= c("label"="name"))


pca.explorer = list(PCA.df = df.pca.annotated, Variance.df = pVar.df, prcomp.out = pca)
saveRDS(object = pca.explorer, file = paste0(pin, "_PCA_Matrix_nTop_", top_var, ".RDS"))




