# Author: Raphael Champeimont
# UMR 7238 Genomique des Microorganismes

rm(list=ls())

while (sink.number() > 0)
  sink()
while (length(dev.list()) > 0)
  dev.off()

library(CLAG)

# Load my functions
MYDOCPATH <- Sys.getenv("MYDOCPATH")[[1]]
R_SCRIPTS_PATH <- paste(MYDOCPATH, "/A/R-scripts", sep="")

source(paste(R_SCRIPTS_PATH, "/auxfunctions.R", sep=""))


pdf.width <- 5
pdf.height <- 5

figureABC <- function(text, right=FALSE) {
  if (right) {
    text(usrFromRelativeX(1), usrFromRelativeY(0.93), labels=text, pos=2, cex=3)
  } else {
    text(usrFromRelativeX(0), usrFromRelativeY(0.93), labels=text, pos=4, cex=3)
  }
}

pdf("~/Documents/workspace/RCLAGarticle/figures/clusters.pdf", width=pdf.width, height=pdf.height)

oldmar <- par("mar")

data(DIM128_subset, package="CLAG")
RES <- CLAG.clust(DIM128_subset)
# Display points in 2D using a PCA and color them by cluster
# except unclunsted points which are left black.
PCA <- prcomp(DIM128_subset)
clusterColors <- c("black", rainbow(RES$ncluster))


par(mar=c(2, 2, 0.5, 0.5))

plot(PCA$x[,1], PCA$x[,2], col=clusterColors[RES$cluster+1], main=NA, xlab=NA, ylab=NA)

#figureABC("A")

dev.off()



pdf("~/Documents/workspace/RCLAGarticle/figures/clusters2.pdf", width=pdf.width, height=pdf.height)



par(mar=c(2, 2, 0.5, 0.5))

data(GLOBINE, package="CLAG")
M <- GLOBINE$M
RES <- CLAG.clust(M, delta=0.2, threshold=0.5, analysisType=3)

o <- order(RES$cluster)
M2 <- M[o,o]

clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]
colorScale <- colorRampPalette(c("blue", "green","yellow","red","darkred"))(1000)
heatmap(M2, symm=TRUE, Colv=NA, Rowv=NA, scale="none", col=colorScale,
        ColSideColors=clusterColors, RowSideColors=clusterColors,
        labRow=NA, labCol=NA, margins=c(0,0), mar=c(2, 4, 0.5, 0.5))

#figureABC("B")



par(mar=oldmar)

dev.off()

