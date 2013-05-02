# Author: Raphael Champeimont
# UMR 7238 Genomique des Microorganismes

rm(list=ls())

while (sink.number() > 0)
  sink()
#while (length(dev.list()) > 0)
#  dev.off()

if (exists("CLAG.clust")) {
  detach("package:CLAG", unload=TRUE)
}
library(CLAG)

# Load my functions
MYDOCPATH <- Sys.getenv("MYDOCPATH")[[1]]
R_SCRIPTS_PATH <- paste(MYDOCPATH, "/A/R-scripts", sep="")
source(paste(R_SCRIPTS_PATH, "/auxfunctions.R", sep=""))

setwd("~/Documents/CLAG-tests/")

microarrayColors <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)

pdf("plots6.pdf")

DATA <- CLAG.loadExampleData("GLOBINE")
M <- DATA$M

RES <- CLAG.clust(M, delta=0.2, threshold=0.5, analysisType=3, rowIds=DATA$rowIds, colIds=DATA$colIds, verbose=TRUE, keepTempFiles=TRUE)

o <- order(RES$cluster)
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]

M2 <- M[o,o]

heatmap(M2, symm=TRUE, Colv=NA, Rowv=NA, scale="none", col=microarrayColors, ColSideColors=clusterColors, RowSideColors=clusterColors)

dev.off()


