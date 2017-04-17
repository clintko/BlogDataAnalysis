pgm_train_innerFun <- function(dat){
    covMat <- cov(dat)
    logDet <- determinant(covMat)
    lst <- list()
    lst$mu1 = apply(dat, 2, mean)
    
    lst$sigma1 = covMat
    lst$prec1  = solve(covMat)
    lst$detsig_log = as.vector(
        logDet$sign * logDet$modulus)
    
    lst$N1 = nrow(dat)
    return(lst)
}

pgm_train <- function(itemName, datTrain){
    lst <- list()
    for(idx in itemName){
        lst[[idx]] <- pgm_train_innerFun(datTrain[[idx]])
    } # end for loop
    return(lst)
}

