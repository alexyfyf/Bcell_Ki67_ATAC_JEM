#! /usr/bin/env Rscript
library(DiffBind)
library(tidyverse)
library(rtracklayer)
library(plyranges)

#args <- commandArgs(trailingOnly = T)

samples <- read.csv("samplesheet_atac_ZD.csv")

dbObj <- dba(sampleSheet = samples, minOverlap = 2)

olap.rate <- dba.overlap(dbObj, dbObj$masks$All, mode = DBA_OLAP_RATE)
pdf("plot/overlap_rate.pdf", height = 5, width = 5)
plot(olap.rate, type='b', ylab="# peaks",
     xlab="Overlap at least this many peaksets")
dev.off()

pdf("plot/peak_heatmap.pdf", height = 5, width = 5)
dba.plotHeatmap(dbObj)
dev.off()

pdf("plot/peak_pca.pdf", height = 5, width = 5)
dba.plotPCA(dbObj)
dev.off()

register(MulticoreParam())
dbObj.count <- dba.count(dbObj, bUseSummarizeOverlaps = T, minOverlap = 2)

## PCA based on peak score
pdf("plot/counts_heatmap.pdf", height = 5, width = 5)
dba.plotHeatmap(dbObj.count)
dev.off()

pdf("plot/counts_pca.pdf", height = 5, width = 5)
dba.plotPCA(dbObj.count)
dev.off()

saveRDS(dbObj.count, "data/db.rds")
