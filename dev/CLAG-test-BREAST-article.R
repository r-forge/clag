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

microarrayColors <- colorRampPalette(c("darkgreen", "green","yellow","red","darkred"))(1000)

pdf("plots-BREAST-article.pdf")

data(BREAST, package="CLAG")

# Correct solution
o <- c(1L, 3L, 4L, 6L, 7L, 8L, 10L, 11L, 12L, 2L, 5L, 9L, 13L, 18L, 
       14L, 15L, 16L, 17L, 19L, 20L)



RES <- CLAG.clust(BREAST, delta=0.2, normalization="affine-global")
clusterColors <- c("black", rainbow(RES$nclusters))[RES$cluster[o]+1]

#RES.NORM <- CLAG.clust(t(BREAST), normalization="affine-column")
#M2 <- RES.NORM$A[,o]
M2 <- t(BREAST[o,1:189])
#a <- quantile(M2, 0.01)
#b <- quantile(M2, 0.99)
#M2[M2 > b] <- b
#M2[M2 < a] <- a
M2 <- matrix(rank(M2), nrow(M2))

heatmap(M2, Colv=NA, Rowv=NA, scale="none", col=microarrayColors,
        ColSideColors=clusterColors)



dev.off()


