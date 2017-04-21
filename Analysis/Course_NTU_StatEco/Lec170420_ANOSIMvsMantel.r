library(vegan)
library(cluster)
library(fpc)

#################################################
#standardize data
data = envdata[,-1]
data.std = scale(data)

#distance matrix
distEU = dist(data.std,method='euclidean')
matDist = as.matrix(distEU)

#clustering
resTree = agnes(distEU,method='ward')
plot(resTree)
resClus = cutree(resTree, k=5)
names(resClus) <- 1:34



#################################################
# get pairs
n <- 34
idxPair <- t(combn(1:n, 2))

# Distance
distPair < apply(idxPair, 1, function(x){
    idx1 <- x[1]
    idx2 <- x[2]
    res <- matDist[idx1,idx2]
}) # end apply
distRank <- rank(-distPair) # rank from the smallest


# within vs between
isWithin <- apply(idxPair, 1, function(x){
    idx1 <- x[1]
    idx2 <- x[2]
    res <- ifelse(resClus[idx1] == resClus[idx2], 1, 0)
})
isBetween <- ifelse(isWithin, 0, 1)

#combine
idxPair <- cbind(idxPair, distPair, distRank)
idxPair <- cbind(idxPair, isWithin, isBetween)
#convert to data frame
idxPair <- as.data.frame(idxPair)


#################################################
#ANOSIM
rankW <- idxPair$distRank[
    ifelse(idxPair$isWithin, T, F)]
rankB <- idxPair$distRank[
    ifelse(idxPair$isBetween, T, F)]
M <- n * (n-1) / 2

statANOSIM <- (mean(rankW) - mean(rankB)) / (M/2)


#################################################
#Mantel
x <- idxPair$distPair 
y <- idxPair$isBetween
statMantel <- sum(x*y)

