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

CLAG.normalize <- function(X) {
  r <- range(X, na.rm=TRUE)
  return((X - r[1])/(r[2]-r[1]))
}

# Noramlize each column of a matrix independently.
# This is a simple homothetic transformation to [0,1].
CLAG.normalizeColumns <- function(M) {
  for (j in 1:ncol(M)) {
    M[,j] <- CLAG.normalize(M[,j])
  }
  return(M)
}


# Noramlize each column of a matrix independently.
# This is a simple homothetic transformation to [0,1].
CLAG.rankColumns <- function(M) {
  for (j in 1:ncol(M)) {
    M[,j] <- rank(M[,j])
  }
  return(M)
}



setwd("~/Documents/CLAG-tests/")

pdf("plots2.pdf")

for (f in c(1, 10, 100, 1000)) {
  delta <- 0.05
  M <- iris[,-5]
  M$Sepal.Width <- M$Sepal.Width * f
  M$Petal.Length <- M$Petal.Length / f
  RES <- CLAG.clust(M, delta=delta, verbose=TRUE, keepTempFiles=TRUE)
  plot(cbind(M, iris$Species), col=RES$cluster, main=paste("iris delta=", delta, " f=", f, " #clusters=", RES$nclusters, sep=""))
  
  delta <- 0.1
  M2 <- CLAG.normalizeColumns(M)
  RES <- CLAG.clust(M2, delta=delta, verbose=TRUE)
  plot(cbind(M, iris$Species), col=RES$cluster, main=paste("iris delta=", delta, " f=", f, " normByCol", " #clusters=", RES$nclusters, sep=""))
  
  delta <- 0.1
  M2 <- CLAG.rankColumns(M)
  RES <- CLAG.clust(M2, delta=delta, verbose=TRUE)
  plot(cbind(M, iris$Species), col=RES$cluster, main=paste("iris delta=", delta, " f=", f, " normByCol", " #clusters=", RES$nclusters, sep=""))
  
}

plot(iris, col=iris$Species, main="iris species")

dev.off()


