
envdata = read.table('enviANDdensity.txt',header=T)

library(vegan)
library(cluster)
library(fpc)
library(TWIX)


data = envdata[,-1]
data.std = scale(data)

#compute distance matrix
dist.eucl = dist(data.std,method='euclidean')

y.eucl.ward = hclust(dist.eucl,method='ward.D')
y.eucl.ward.cut = cutree(y.eucl.ward, 5)

MRPP <- function(data, group, n = 999, ...){
    col_n = ncol(data)
    group = factor(group)
    k = nlevels(group)
    group_split = split(data, group)
    d = numeric()
    for(i in 1:k){
        d[i] = mean(dist(matrix(group_split[[i]], length(group_split[[i]])/col_n, col_n), ...))
    }
    delta = sum(table(group) / length(group) * d)
    
    delta_perm = numeric()
    for(j in 1:n){
        group_p = sample(group)
        group_p_split = split(data, group_p)
        d_p = numeric()
        for(i in 1:k){
            d_p[i] = mean(dist(matrix(group_p_split[[i]], length(group_p_split[[i]])/col_n, col_n), ...))
        }
        delta_perm[j] = sum(table(group_p) / length(group_p) * d_p)
    }
    A = 1-delta/mean(delta_perm)
    sig = sum(delta_perm <= delta) 
    if(sig == 0){sig = paste("<", format(1/(n+1), scientific = F))}
    return(list(observ_delta = delta, permuted_delta = delta_perm, A = A, significance = sig))
}

test = MRPP(data.std, y.eucl.ward.cut)