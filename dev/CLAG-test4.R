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

pdf("plots4.pdf")

# First example with real variables
DATA <- CLAG.loadExampleData("DIM128-subset")
RES <- CLAG.clust(DATA, verbose=TRUE, keepTempFiles=TRUE)
# Display points in 2D using a PCA and color them by cluster
# except unclunsted points which are left black.
PCA <- prcomp(DATA)
clusterColors <- c("black", rainbow(RES$ncluster))
plot(PCA$x[,1], PCA$x[,2], col=clusterColors[RES$cluster+1], main=paste(RES$nclusters, "clusters"))

# Second example with boolean variables
# (we replace each variable with TRUE or FALSE
# depending on whether the value is below or above 128)
DATA2 <- DATA > 128
RES <- CLAG.clust(DATA2, analysisType=2, verbose=TRUE, keepTempFiles=TRUE)
PCA <- prcomp(DATA2)
clusterColors <- c("black", rainbow(RES$ncluster))
plot(PCA$x[,1], PCA$x[,2], col=clusterColors[RES$cluster+1], main=paste(RES$nclusters, "clusters"))



dev.off()


