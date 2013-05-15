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

setwd("~/Documents/CLAG-tests/")

data(DIM128_subset)
DATA <- DIM128_subset

correctClustering <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 
                       4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 
                       7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 
                       10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 
                       13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 
                       15, 15, 15, 16, 16, 16, 16, 16, 16, 16)


pdf("plots3-250.pdf")

delta <- 0.05

#DATA <- DATA[seq(1, nrow(DATA), by=10),]
PCA <- prcomp(DATA)

for (caseNumber in 1:3) {
  
  M <- DATA
  
  if (caseNumber == 1) {
    caseName <- "original"
  } else if (caseNumber == 2) {
    caseName <- "translated"
  } else {
    caseName <- "other"
  }
  for (j in 1:ncol(M)) {
    if (caseNumber == 1) {
    } else if (caseNumber == 2) {
      b <- j*250
      cat("column", j, " + ", b, "\n")
      M[,j] <- M[,j] + b
    } else if (caseNumber == 3) {
      if (j %% 2 == 0) {
        M[,j] <- exp(M[,j]/10) + j*250
      } else {
        M[,j] <- M[,j] + j*250
      }
    } else {
      stop()
    }
  }
  
  V <- c()
  C <- c()
  for (j in 1:10) {
    V <- c(V, M[,j])
    C <- c(C, rep(j, nrow(M)))
  }
  #rlcCompareDens(V, C, bw=10, main=paste("data=", caseName, sep=""))
  
  for (normalizationMethod in c("affine-global", "affine-column", "rank-column")) {
  
    for (delta in c(0.05)) {
      RES <- CLAG.clust(M, delta=delta, threshold=0.5, normalization=normalizationMethod, verbose=TRUE)
      
      V <- c()
      C <- c()
      for (j in 1:10) {
        V <- c(V, RES$A[,j])
        C <- c(C, rep(j, nrow(M)))
      }
      r <- range(V)
      rlcCompareDens(V, C, bw=(r[2]-r[1])/100, main=paste("data=", caseName, " normalization=", normalizationMethod, sep=""))
      
      clusterRemapped <- rep(-1, length(RES$cluster))
      k <- 1
      for (i in unique(RES$cluster[RES$cluster >= 1])) {
        clusterRemapped[RES$cluster == i] <- k
        k <- k + 1
      }
      
      PCA <- prcomp(M)
      
      mapping <- compareClusterings(RES$cluster, correctClustering, verbose=TRUE)
      plot(-PCA$x[,1], PCA$x[,2], col=mapColors(clusterRemapped), pch=4, cex.main=0.8,
           main=paste("DIM128 ", caseName, " normalization=", normalizationMethod, " delta=", delta, " #clusters=", RES$nclusters, " bad=", mapping$ndiff, sep=""))
      
      points(-PCA$x[mapping$diffclust,1], PCA$x[mapping$diffclust,2])
      
    }
  }
}

dev.off()

