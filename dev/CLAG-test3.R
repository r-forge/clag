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
source(paste(R_SCRIPTS_PATH, "/auxfunctions.R", sep=""))

setwd("~/Documents/CLAG-tests/")

DATA <- CLAG.loadExampleData("DIM128-subset")

pdf("plots3.pdf")

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
      b <- j*1000
      cat("column", j, " + ", b, "\n")
      M[,j] <- M[,j] + b
    } else if (caseNumber == 3) {
      if (j %% 2 == 0) {
        M[,j] <- exp(M[,j]/10) + j*1000
      } else {
        M[,j] <- M[,j] + j*1000
      }
    } else {
      stop()
    }
  }
  
  V <- c()
  C <- c()
  for (j in 1:ncol(M)) {
    V <- c(V, M[,j])
    C <- c(C, rep(j, nrow(M)))
  }
  rlcCompareDens(V, C)
  
  for (normalizationMethod in c("affine-global", "affine-column", "rank-column")) {
  
    for (delta in c(0.05)) {
      RES <- CLAG.clust(M, delta=delta, threshold=0.5, normalization=normalizationMethod, verbose=TRUE)
      
      clusterRemapped <- rep(-1, length(RES$cluster))
      k <- 1
      for (i in unique(RES$cluster[RES$cluster >= 1])) {
        clusterRemapped[RES$cluster == i] <- k
        k <- k + 1
      }
      
      plot(-PCA$x[,1], PCA$x[,2], col=mapColors(clusterRemapped),
           main=paste("DIM128 ", caseName, " normalization=", normalizationMethod, " delta=", delta, " #clusters=", RES$nclusters, sep=""))
    }
  }
}

dev.off()


