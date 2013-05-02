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


pdf("plots_test_all.pdf")


cat("Test general matrix\n")
# First example with real variables
data(DIM128_subset, package="CLAG")
RES <- CLAG.clust(DIM128_subset)
# Display points in 2D using a PCA and color them by cluster
# except unclunsted points which are left black.
PCA <- prcomp(DIM128_subset)
clusterColors <- c("black", rainbow(RES$ncluster))
plot(PCA$x[,1], PCA$x[,2], col=clusterColors[RES$cluster+1], main=paste(RES$nclusters, "clusters"))


cat("Test binary matrix\n")
# Second example with boolean variables
# (we replace each variable with TRUE or FALSE
# depending on whether the value is below or above 128)
DATA2 <- DIM128_subset > 128
RES <- CLAG.clust(DATA2, analysisType=2)
PCA <- prcomp(DATA2)
clusterColors <- c("black", rainbow(RES$ncluster))
plot(PCA$x[,1], PCA$x[,2], col=clusterColors[RES$cluster+1], main=paste(RES$nclusters, "clusters"))




cat("Test specific matrix\n")
data(GLOBINE, package="CLAG")
M <- GLOBINE$M
RES <- CLAG.clust(M, delta=0.2, threshold=0.5, analysisType=3)
o <- order(RES$cluster)
M2 <- M[o,o]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(M2, symm=TRUE, Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors, RowSideColors=clusterColors)



cat("Test specific matrix\n")
data(GLOBINE, package="CLAG")
M <- GLOBINE$M
M <- M[1:30,]
RES <- CLAG.clust(M, delta=0.2, threshold=0.5, analysisType=3, rowIds=1:30)
o <- order(RES$cluster)
M2 <- M[o,o]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(M2, symm=TRUE, Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors, RowSideColors=clusterColors)

data(GLOBINE, package="CLAG")
M <- GLOBINE$M
M <- M[1:30,1:30]
RES <- CLAG.clust(M, delta=0.2, threshold=0.5, analysisType=3)
o <- order(RES$cluster)
M2 <- M[o,o]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(M2, symm=TRUE, Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors, RowSideColors=clusterColors)


invisible(dev.off())


