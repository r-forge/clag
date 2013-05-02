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

microarrayColors <- rev(colorRampPalette(c("darkgreen", "green","yellow","red","darkred"))(1000))
#image(matrix(1:1000), col=microarrayColors)

pdf("plots5.pdf")

# First example with real variables
DATA <- CLAG.loadExampleData("BREAST")
M <- DATA$M

RES <- CLAG.clust(M, delta=0.2, verbose=TRUE)

M2 <- M
o <- order(RES$cluster)
M2 <- M2[o,]
for (j in 1:ncol(M2)) {
  r <- range(M2[,j])
  M2[,j] <- (M2[,j]-r[1])/(r[2]-r[1])
}

heatmap(t(M2), col=microarrayColors, scale="none", ColSideColors=c("red","blue","green")[RES$cluster[o]])


dev.off()


