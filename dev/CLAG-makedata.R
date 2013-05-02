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


setwd("~/Documents/A/R-scripts/CLAG")

dataSetNames <- c("BREAST", "GLOBINE", "DIM128", "DIM128-subset")

BREAST <- CLAG.loadExampleData("BREAST")
GLOBINE <- CLAG.loadExampleData("GLOBINE")
DIM128 <- CLAG.loadExampleData("DIM128")
DIM128_subset <- CLAG.loadExampleData("DIM128-subset")

save(BREAST, file="BREAST.rda", compress=TRUE)
save(GLOBINE, file="GLOBINE.rda", compress=TRUE)
save(DIM128, file="DIM128.rda", compress=TRUE)
save(DIM128_subset, file="DIM128_subset.rda", compress=TRUE)



