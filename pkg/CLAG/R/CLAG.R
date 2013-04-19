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
  if (is.null(colIds)) {
    colIds <- 1:ncol(A)
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
  if (!file.exists(fpath)) {
    stop("No CLAG result file.")
  }
  return(read.table(fpath, sep=":", header=FALSE, fill=TRUE, stringsAsFactors=FALSE, row.names=NULL))
}

CLAG.removeDir <- function(outputDir) {
  unlink(outputDir, recursive=TRUE)
}


CLAG.clust <- function(M,
                       delta=0.05,
                       threshold=0,
                       analysisType=1,
                       normalization="affine-global",
                       rowIds=NULL,
                       colIds=NULL,
                       verbose=FALSE,
                       keepTempFiles=FALSE) {
  
  if (analysisType != 1 && analysisType != 2 && analysisType != 3) {
    stop("analysisType must be 1, 2, or 3")
  } else {
    analysisType <- as.integer(analysisType)
  }
  
  
  if (analysisType != 2) {
    # Variables are real
    if (delta <= 0 || delta > 1) {
      stop("delta must be between 0 and 1")
    }
    d <- 100*delta
  }
  
  if (is.null(rowIds)) {
    rowIds <- 1:nrow(M)
  } else {
    if (nrow(M) != length(rowIds)) {
      stop("rowIds need to contain as many elements as rows in matrix")
    }
    if (anyDuplicated(rowIds)) {
      stop("rowIds need to be unique")
    }
  }
  
  if (is.null(colIds)) {
    colIds <- 1:ncol(M)
  } else {
    if (ncol(M) != length(colIds)) {
      stop("colIds need to contain as many elements as columns in matrix")
    }
    if (anyDuplicated(colIds)) {
      stop("colIds need to be unique")
    }
  }
  
  if (analysisType == 3) {
    if (!all(colIds %in% rowIds)) {
      stop("column ids need to be a subset of row ids when analysisType=3")
    }
  }
  
  if (normalization == "affine-global") {
    # do nothing, it is done automatically by CLAG
  } else if (normalization == "affine-column") {
    for (j in 1:ncol(M)) {
      r <- range(M[,j], na.rm=TRUE)
      M[,j] <- (M[,j] - r[1])/(r[2] - r[1])
    }
  } else if (normalization == "rank-column") {
    for (j in 1:ncol(M)) {
      M[,j] <- rank(M[,j])
    }
  } else {
    stop("invalid normalization method specified")
  }

  outdir <- tempfile("CLAG_")
  dir.create(outdir)
  if (verbose) cat("Writing to", outdir, "\n")
  fn <- paste(outdir, "/input.txt", sep="")
  CLAG.writeInput(M, fn, rowIds=rowIds, colIds=colIds)
  CLAG.exec(f=outdir, p=analysisType, k=threshold, d=d, verbose=verbose)
  if (analysisType != 2) {
    respath <- paste(outdir, "/aggregation-", d, ".txt", sep="")
  } else {
    respath <- paste(outdir, "/aggregation.txt", sep="")
  }
  rawClusters <- CLAG.readClusters(respath)
  if (!keepTempFiles) {
    CLAG.removeDir(outdir)
  }
  
  RES <- list()
  RES$cluster <- rep(0L, nrow(M))
  if (nrow(rawClusters) > 0) {
    for (i in 1:nrow(rawClusters)) {
      elements <- as.integer(strsplit(rawClusters[i,1], split=" ")[[1]])
      indices <- match(elements, rowIds)
      if (any(RES$cluster[indices] != 0L)) {
        stop("element in 2 clusters")
      }
      RES$cluster[indices] <- i
      if (analysisType == 3) {
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
  RES$nclusters <- nrow(rawClusters)
  RES$delta <- delta
  RES$threshold <- threshold
  RES$analysisType <- analysisType
  RES$M <- M
  RES$rowIds <- rowIds
  RES$colIds <- colIds
  RES$rawClusters <- rawClusters
  return(RES)
}



CLAG.loadExampleData <- function(set="BREAST") {
  if (set %in% c("BREAST","BRAIN","GLOBINE")) {
    return(CLAG.readInput(paste(CLAG.data.path, "/", set, "/input.txt", sep="")))
  } else if (set == "DIM128") {
    return(read.table(paste(CLAG.data.path, "/DIM128.txt", sep="")))
  } else if (set == "DIM128-subset") {
    return(read.table(paste(CLAG.data.path, "/DIM128-subset.txt", sep="")))
  } else {
    stop(paste("unknown CLAG example data set:", set))
  }
}

