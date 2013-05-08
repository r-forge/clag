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
source(paste(R_SCRIPTS_PATH, "/auxali.R", sep=""))

setwd("~/Documents/CLAG-tests/")

microarrayColors <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)

pdf("plots-BREAST.pdf")

data(BREAST, package="CLAG")

# Correct solution
o <- c(1L, 3L, 4L, 6L, 7L, 8L, 10L, 11L, 12L, 2L, 5L, 9L, 13L, 18L, 
       14L, 15L, 16L, 17L, 19L, 20L)

RES <- CLAG.clust(BREAST, delta=0.2, normalization="affine-global")
#o <- order(RES$cluster)
M2 <- RES$A[o,]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(t(M2[,1:50]), Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors)


RES <- CLAG.clust(BREAST, delta=0.2, normalization="rank-column")
#o <- order(RES$cluster)
M2 <- RES$A[o,]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(t(M2[,1:50]), Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors)

RES <- CLAG.clust(BREAST, delta=0.1, normalization="rank-column")
#o <- order(RES$cluster)
M2 <- RES$A[o,]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(t(M2[,1:50]), Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors)

RES <- CLAG.clust(BREAST, delta=0.4, normalization="rank-column")
#o <- order(RES$cluster)
M2 <- RES$A[o,]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(t(M2[,1:50]), Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors)

RES <- CLAG.clust(BREAST, delta=0.5, normalization="rank-column")
#o <- order(RES$cluster)
M2 <- RES$A[o,]
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(t(M2[,1:50]), Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors)



dev.off()


