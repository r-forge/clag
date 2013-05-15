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

alignment <- readAlignment("CFTR.fasta")
alignment <- alignment[nrow(alignment):1, ]
print(dim(alignment))

pdf("plots8.pdf")

VALS <- rownames(taylorColors)
ZZ <- match(alignment, VALS)
ZZ[ZZ == 1] <- NA
M2 <- matrix(ZZ, ncol=ncol(alignment), nrow=nrow(alignment))
rownames(M2) <- rownames(alignment)

heatmap(M2, col=taylorColors$color, scale="none", Colv=NA, Rowv=NA)

RES <- CLAG.clust(alignment, analysisType=2, threshold=0.95,
                  verbose=TRUE, keepTempFiles=TRUE)

heatmap(M2, col=taylorColors$color, scale="none", Colv=NA, Rowv=NA, RowSideColors=c("black", rainbow(RES$nclusters))[RES$cluster+1])

dev.off()


