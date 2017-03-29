gen_utagmat <- function(utagsvec, y) {
    if (is.null(utagsvec)){
        constant <- rep(1, length(y))
        mat <- as.matrix(constant)
        colnames(mat) <- "constant"
        return(mat)
    } # end if
    
    lst <- strsplit(utagsvec, split=",")
    
    tmp <- unlist(lst)
    tmp <- na.omit(tmp)
    
    threshold <- 5
    count <- table(tmp)
    count <- count[count >= threshold]
    count <- sort(count, decreasing = TRUE)
    
    setAll <- as.numeric(names(count))
    if (length(setAll) == 0){
        constant <- rep(1, length(y))
        mat <- as.matrix(constant)
        colnames(mat) <- "constant"
        return(mat)
    } # end if
    
    mat <- lapply(
        lst, 
        function(x){
            as.numeric(setAll %in% x)
        } # end lapply func
    ) # end lapply
    
    mat <- do.call(rbind, mat)
    colnames(mat) <- paste("user", setAll, sep="_")
    
    val <- apply(
        mat, 2, 
        function(x) {
            xmat = matrix(1, ncol=2, nrow=length(y))
            xmat[,2] = x
            bhead  = solve(t(xmat) %*% xmat, t(xmat) %*% y)
            yhead  = xmat %*% bhead
            e1     = y - yhead
            var1   = sum(e1 * e1) / (length(e1)-2)
            sigma2 = solve(t(xmat) %*% xmat) * var1
            t1     = bhead[2] / sqrt(sigma2[2,2])
            return(t1)
        } # end apply func
    ) # end apply
    
    val <- abs(val)
    val <- sort(val, decreasing = TRUE)
    val <- val[val >= 1]
    
    mat <- mat[,names(val)]
    mat <- cbind(1,mat)
    colnames(mat)[1] <- "constant"
    
    return(mat)
} # end func gen_utagmat

