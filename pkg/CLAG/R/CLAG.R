# CLAG: an unsupervised non hierarchical clustering algorithm
# Copyright (c) 2013, Linda DIB, Raphael CHAMPEIMONT, Alessandra CARBONE
#                     UMR 7238 Génomique des Microorganismes
#                     CNRS - Université Pierre et Marie Curie
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this list
#   of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author of R package: Raphael Champeimont
# UMR 7238 Genomique des Microorganismes
# This R package allows for easy use of CLAG.
# CLAG: an unsupervised non hierarchical clustering algorithm
# handling biological data, by Linda Dib and Alessandra Carbone

CLAG.path <- system.file("Perl/TOOLCLAG-minimal", package="CLAG")
CLAG.data.path <- system.file("extdata/dataCLAG", package="CLAG")

CLAG.writeInput <- function(A, outfile, rowIds=NULL, colIds=NULL) {
  f <- file(outfile, "w")
  N <- nrow(A)
  M <- ncol(A)
  if (is.null(rowIds)) {
    rowIds <- 1:nrow(A)
  }
  if (is.null(rowIds)) {
    rowIds <- 1:ncol(A)
  }
  for (i in 1:N) {
    for (j in 1:M) {
      cat(rowIds[[i]], " ", colIds[[j]], " ", A[i,j], "\n", sep="", file=f)
    }
  }
  close(f)
}

CLAG.readInput <- function(infile, p=1) {
  rawInput <- read.table(infile, sep=" ", header=FALSE, stringsAsFactors=FALSE)
  rowIds <- sort(unique(rawInput[,1]))
  colIds <- sort(unique(rawInput[,2]))
  destRows <- match(rawInput[,1], rowIds)
  destCols <- match(rawInput[,2], colIds)
  colIndices <- 
  MAT <- matrix(NA, nrow=length(rowIds), ncol=length(colIds))
  for (i in 1:nrow(rawInput)) {
    MAT[destRows[i], destCols[i]] <- rawInput[i,3]
  }
  DATA <- list(M=MAT, rowIds=rowIds, colIds=colIds)
  
  return(DATA)
}

# Run CLAG
# Possibles values for p:
# 1 a general matrix
# 2 a binary matrix
# 3 a matrix where N included in E
CLAG.exec <- function(f, p=1, k=0, d=NULL, verbose=TRUE) {
  exefile <- paste(CLAG.path, "/exe-RCommand.pl", sep="")
  if (!file.exists(exefile)) {
    stop(paste(exefile, "does not exist"))
  }
  cargs <- c(paste("-f=", f, sep=""),
             paste("-p=", p, sep=""))
  if (!is.null(k)) {
    cargs <- c(cargs, paste("-k=", k, sep=""))
  }
  if (!is.null(d)) {
    cargs <- c(cargs, paste("-d=", d, sep=""))
  }
  if (verbose) cat(exefile, paste(cargs, collapse=" "), "\n")
  oldwd <- getwd()
  setwd(CLAG.path)
  if (verbose) {
    stdout <- ""
  } else {
    stdout <- FALSE
  }
  system2(exefile, args=cargs, stdout=stdout)
  setwd(oldwd)
}

CLAG.readClusters <- function(fpath) {
  return(read.table(fpath, sep=":", header=FALSE, fill=TRUE, stringsAsFactors=FALSE, row.names=NULL))
}

CLAG.removeDir <- function(outputDir) {
  unlink(outputDir, recursive=TRUE)
}

# DATA is either a matrix or a list composed of
# - a matrix M
# - a vector colIds with column ids (same length as ncol(M))
# - a vector rowIds with row ids (same length as nrow(M))
# analysisType: type of matrix (-p parameter of CLAG)
#    1 for generic case where lines are elements and columns characters
#      -> only environnment score is computed
#    3 if characters are a subset of elements
#      -> both environnment score and symmetric scores are computed
#      (you can use 1 if you don't want the symmetric score)
# delta: delta parameter (examples: 0.05, 0.1 or 0.2)
# threshold: score threshold in [-1,+1] (examples: 0, 0.5, 1)
CLAG.clust <- function(DATA, delta=0.05, threshold=0, analysisType=1, verbose=FALSE, keepTempFiles=FALSE) {
  p <- analysisType
  k <- threshold
  d <- 100*delta
  if (is.matrix(DATA) || is.data.frame(DATA) || is.null(DATA$M)) {
    M <- DATA
    rowIds <- 1:nrow(DATA)
    colIds <- 1:ncol(DATA)
  } else {
    M <- DATA$M
    rowIds <- DATA$rowIds
    colIds <- DATA$colIds
  }
  stopifnot(nrow(M) == length(rowIds))
  stopifnot(ncol(M) == length(colIds))
  outdir <- tempfile("CLAG_")
  dir.create(outdir)
  if (verbose) cat("Writing to", outdir, "\n")
  fn <- paste(outdir, "/input.txt", sep="")
  CLAG.writeInput(M, fn, rowIds=rowIds, colIds=colIds)
  CLAG.exec(f=outdir, p=p, k=k, d=d, verbose=verbose)
  respath <- paste(outdir, "/aggregation-", d, ".txt", sep="")
  rawClusters <- CLAG.readClusters(respath)
  if (!keepTempFiles) {
    CLAG.removeDir(outdir)
  }
  
  RES <- list()
  RES$cluster <- rep(-1, nrow(M))
  RES$firstScores <- c()
  RES$lastScores <- c()
  if (nrow(rawClusters) > 0) {
    for (i in 1:nrow(rawClusters)) {
      elements <- as.integer(strsplit(rawClusters[i,1], split=" ")[[1]])
      indices <- match(elements, rowIds)
      if (any(RES$cluster[indices] != -1)) {
        stop("element in 2 clusters")
      }
      RES$cluster[indices] <- i
      if (p == 3) {
        RES$firstSymScore <- c(RES$firstScores, as.numeric(rawClusters[i,3]))
        RES$lastSymScore <- c(RES$lastScores, as.numeric(rawClusters[i,4]))
        RES$firstEnvScore <- c(RES$firstScores, as.numeric(rawClusters[i,5]))
        RES$lastEnvScore <- c(RES$lastScores, as.numeric(rawClusters[i,6]))
      } else {
        RES$firstEnvScore <- c(RES$firstScores, as.numeric(rawClusters[i,3]))
        RES$lastEnvScore <- c(RES$lastScores, as.numeric(rawClusters[i,4]))
      }
    }
  }
  #RES$inputNormalized <- CLAG.readInput(paste(outdir, "/input.txt", sep=""), p=p)
  RES$nclusters <- nrow(rawClusters)
  RES$delta <- d/100
  RES$threshold <- k
  RES$analysisType <- p
  RES$M <- M
  RES$rowIds <- rowIds
  RES$colIds <- colIds
  RES$rawClusters <- rawClusters
  return(RES)
}


CLAG.plotClusters <- function(RES, HC=TRUE) {
  library(lattice)
  
  MAT <- RES$M
  colIds <- RES$colIds
  rowIds <- RES$rowIds
  
  if (RES$analysisType == 3) {
    # In this case, columns correspond to elements
    # so order columns according to the clusters.
    rowRanks <- rank(RES$cluster)
    rowRanks <- ifelse(RES$cluster == -1, Inf, rowRanks)
    # Where is each column corresponding row?
    colRows <- match(colIds, rowIds)
    # For each column, rank of the corresponding cluster
    colRanks <- rowRanks[colRows]
    M <- MAT[, order(colRanks)]
  } else {
    if (HC) {
      # cluster with hierarchical clustering the columns
      M <- MAT[, hclust(dist(t(MAT), method="manhattan"))$order]
    } else {
      M <- MAT
    }
  }
  
  potentials <- NULL
  for (cl in 1:RES$nclusters) {
    if (!is.null(potentials)) {
      potentials <- rbind(potentials, rep(-2, ncol(M)))
    }
    potentials <- rbind(potentials, M[which(RES$cluster == cl), ])
  }
  colorFun <- colorRampPalette(c("lightblue","blue","lightgreen","darkgreen","yellow" ,"orange","red", "darkred"))
  x <- seq(0, 1, by=0.04)
  print(levelplot(potentials,
            scales=list(tck=0, x=list(rot=90),cex=0.1),
            col.regions=colorFun,
            main=NULL,
            xlab=list("key aggregates", cex=1),
            ylab=list("Environment", cex=1),
            aspect="iso", at=x, cut=50))
}


# Possible values are:
# BRAIN  BREAST  GLOBINE DIM128
CLAG.loadExampleData <- function(set="BREAST") {
  if (set %in% c("BREAST","BRAIN","GLOBINE")) {
    return(CLAG.readInput(paste(CLAG.data.path, "/", set, "/input.txt", sep="")))
  } else if (set == "DIM128") {
    return(read.table(paste(CLAG.data.path, "/DIM128.txt", sep="")))
  } else {
    stop(paste("unknown CLAG example data set:", set))
  }
}

