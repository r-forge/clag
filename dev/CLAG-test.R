# Author: Raphael Champeimont
# UMR 7238 Genomique des Microorganismes

rm(list=ls())

while (sink.number() > 0)
  sink()
#while (length(dev.list()) > 0)
#  dev.off()

library(CLAG)

# Load my functions
MYDOCPATH <- Sys.getenv("MYDOCPATH")[[1]]
R_SCRIPTS_PATH <- paste(MYDOCPATH, "/A/R-scripts", sep="")

setwd("~/Documents/CLAG-tests/")

pdf("plots.pdf")

BREAST <- CLAG.loadExampleData("BREAST")
BREAST.RES <- CLAG.clust(BREAST, delta=0.2, threshold=0, verbose=TRUE)
CLAG.plotClusters(BREAST.RES)

GLOBINE <- CLAG.loadExampleData("GLOBINE")
GLOBINE.RES <- CLAG.clust(GLOBINE, analysisType=3, delta=0.2, threshold=0.5, verbose=TRUE)
CLAG.plotClusters(GLOBINE.RES)

dev.off()


