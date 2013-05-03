# Author: Raphael Champeimont
# UMR 7238 Genomique des Microorganismes

# A few extra scripts usefull for further CLAG result analysis


# ===== CLUSTER COMPARISON ======

# Count how many elements are differently cluster
# considering the cluster mapping assoc.
# assoc[i] = j means cluser i in cl1 is cluster j in cl2
CLAG.mapClusteringsTest <- function(cl1, cl2, assoc) {
  bad <- 0
  for (i in 1:length(cl1)) {
    if (cl1[[i]] == 0) {
      # i is unclusted in cl1
      if (cl2[[i]] != 0) {
        bad <- bad + 1
      }
    } else {
      # find associated cluster
      if (cl2[[i]] == 0 || assoc[[cl1[[i]]]] != cl2[[i]]) {
        bad <- bad + 1
      }
    }
  }
  return(bad)
}


# Hypothesis: cl1 and cl2 have non zeros (only clustered elements)
CLAG.mapClusteringsAux <- function(cl1, cl2) {
  stopifnot(length(cl1) == length(cl2))
  stopifnot(all(cl1 > 0))
  stopifnot(all(cl2 > 0))
  # number of clusters in cl1 and cl2
  ncl1 <- max(cl1)
  ncl2 <- max(cl2)
  bestForNow <- Inf
  aux <- function(assoc, k, cost) {
    if (k > length(assoc)) {
      cat("assoc=", assoc, " bad=", cost, " bestForNow=", bestForNow, "\n")
      return(list(bad=cost, assoc=assoc))
    } else {
      # Search the partner for cluster k
      best <- Inf
      bestassoc <- NA
      candidates <- unique(cl2[cl1 == k & cl2 != 0])
      candidates <- c(candidates[! (candidates %in% assoc[1:(k-1)])], 0)
      elements <- which(cl1 == k)
      for (j in candidates) {
        assoc[k] <- j
        # How many differently clustered elements does this choice generate?
        # special case: if j = 0 all will be different (hyp: no 0 in cl2)
        morebad <- sum(cl2[elements] != j)
        newcost <- cost + morebad
        # If cost already exceeds (or equals) lowest cost known,
        # it is usless to explore the branch (cut it).
        if (newcost < bestForNow) {
          res <- aux(assoc, k+1, newcost)
          if (res$bad < best) {
            best <- res$bad
            bestassoc <- res$assoc
            if (res$bad < bestForNow) {
              bestForNow <<- res$bad
            }
          }
        }
      }
      return(list(bad=best, assoc=bestassoc))
    }
  }
  return(aux(rep(0L, ncl1), 1, 0))
}


# Given two clusterings, find the best association between clusters
# such that the number of differently clustered element is minimum.
# Assumed:
# 0 = unclustered
# integer >= 1: number of cluster to which the element belong
CLAG.mapClusterings <- function(cl1, cl2) {
  stopifnot(length(cl1) == length(cl2))
  # When an element is unclusted in one clustering only,
  # it will count as "different" whatever the cluster mapping is
  unavoidableCost <- sum(xor(cl1 == 0, cl2 == 0))
  # It is useless to try to map unclunstered elements, se remove then
  keep <- cl1 != 0 & cl2 != 0
  res <- CLAG.mapClusteringsAux(cl1[keep], cl2[keep])
  res$bad <- res$bad + unavoidableCost
  
  # Check the score by an independant method
  badrecount <- CLAG.mapClusteringsTest(cl1, cl2, res$assoc)
  stopifnot(badrecount == res$bad)
  
  return(res)
}


# Assumed:
# 0 = unclustered
# integer >= 1: number of cluster to which the element belong
CLAG.compareClusterings <- function(cl1, cl2) {
  res <- CLAG.mapClusterings(cl1, cl2)
  cat("Cluster mapping:", res$assoc, "\n")
  cl1a <- rep(NA, length(cl1))
  for (i in 1:length(cl1)) {
    if (cl1[[i]] == 0) {
      cl1a[[i]] <- 0
    } else {
      cl1a[[i]] <- res$assoc[[cl1[[i]]]]
    }
  }
  cat("CL1 original:", cl1, "\n")
  cat("CL1 renamed: ", cl1a, "\n")
  cat("CL2 original:", cl2, "\n")
  cat("Differently clustered elements:", res$bad, "\n")
  res$diffclust <- cl1a != cl2
  return(res)
}

cat("\n\n")
CLAG.compareClusterings(c(0,1,1,1,2,2,2,3,3,3,3), c(0,3,3,3,1,1,1,2,2,2,2))

cat("\n\n")
CLAG.compareClusterings(c(0,1,1,1,2,2,2,3,3,3,3,0), c(3,1,3,3,1,1,1,2,0,2,4,0))


