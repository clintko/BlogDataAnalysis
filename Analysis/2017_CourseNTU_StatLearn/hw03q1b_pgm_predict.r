myKernal <- function(x, lstModel){
    p <- length(x)
    x <- matrix(as.numeric(x), ncol=1)
    
    mu          <- matrix(lstModel$mu1, ncol=1)
    sigmaInv    <- lstModel$prec1
    sigmaDetlog <- lstModel$detsig_log
    
    valKernal <- -0.5 * t(x - mu) %*% sigmaInv %*% (x - mu)
    res <- as.numeric((-p/2) * sigmaDetlog + valKernal)
    return(res)
}

pgm_predict <- function(lstOfModel, dat){
    # 
    if(ncol(dat) != 3){
        return(NULL)
    } # end if
    
    res <- apply(dat, 1, function(x){
        vec <- lapply(lstOfModel, function(lstModel){
            return(myKernal(x, lstModel))
        }) # end lapply
        
        idx <- order(
            unlist(vec), 
            decreasing = TRUE)[1]
        return(idx)
    }) # end apply
    
    # return
    names(res) <- NULL
    return(res)
}


