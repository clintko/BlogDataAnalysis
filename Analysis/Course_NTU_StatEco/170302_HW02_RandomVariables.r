
###########################################################
# create my functions for each random variable generator
###########################################################
# Gaussian Distribution
myGauss <- function(n, mu, sigma2){
    x <- rnorm(n)
    y <- x * (sigma2^0.5) + mu
    return(y)
} # end function myGauss

# Shifting a Gaussian Distribution into desire mu and variance
ShiftDistr <- function(x, mu, sigma2){
    # first transform the distribution into standard norm
    z <- x / sd(x) - mean(x)
    # then shift the distribution into desire one
    y <- z * (sigma2^0.5) + mu
    return(y)
} # end function myGaussShift

# Uniform Distribution [0,1]
myUnif <- function(n, a=0, b=1) {
    x <- runif(n)
    y <- x * (b-a) + a
    return(y)
} # end func myUnif

# Bernoulli Distribution
myBern  <- function(n, prob=0.5){
    res <- ifelse(runif(n) <= prob, 1, 0)
    return(res)
} # end func myBern

# Binomial Distribution
myBinom <- function(n, size, prob=0.5){
    # Generate random numbers from binomial distribution
    # which is the sum of bernoulli random variable
    res <- replicate(n, sum(myBern(size, prob)))
    return(res)
} # end func myBinom

#########################################
# Helper function to plot matrix
#########################################
plotMatrix <- function(dat, isXlab=T, isYlab=T){
    # given a data (matrix or dataframe), plot the
    # value of each element as color
    # ---------------------------------
    # Set color scale
    colPalette.WR <- colorRampPalette(c(
        "#fef0d9", "#fdbb84", "#fc8d59", "#e34a33"))(256)
    
    # Flipping the matrix
    mat <- as.matrix(dat)
    mat <- t(mat)
    mat <- mat[,ncol(mat):1,drop=FALSE]
    
    # Plot the matrix
    image(mat, xaxt= "n", yaxt= "n",
          col = colPalette.WR)
    
    # Show ticks if needed
    if(isXlab){ # if X ticks are added
        axis(1, at=seq(0,1,length.out=ncol(dat)), 
             labels=colnames(dat), las=1)
    } # end if
    if(isYlab){ # if Y ticks are added
        axis(2, at=seq(0,1,length.out=nrow(dat)), 
             labels=rev(rownames(dat)), las=2)
    } # end if
} # end func plotMatrix

###########################################
# 1a Generate 10000 random numbers from 
# Gaussian distribution with mean=20 and 
# variance=10, and plot the distribution
###########################################
# Decide the number of random numbers
N <- 10000

# Generate standard normal distribution
x <- myGauss(N, 0, 1)

# Shift the distribution into mean=20, variance=10
y <- ShiftDistr(x, 20, 10)

# generate the histogram
gp1 <- hist(x, plot=FALSE)
gp2 <- hist(y, plot=FALSE)

# visualize both histogram in one plot
plot(gp1, col=rgb(0,0,1,1/4), 
     xlim=c(-10, 40), 
     ylim=c(  0, N/3))
plot(gp2, col=rgb(1,0,0,1/4), add=T)

###########################################
# 1b. Generate 10000 random numbers from 
# Binomial distribution with p=0.5 and n=40, 
# and plot the distribution.
###########################################
# Decide the number of random numbers
N <- 10000

# Generate random numbers from binomial distribution
x <- myBinom(N, 40)

# plot the results
hist(x, col="grey80", xlim=c(0, 40))
axis(side=1, lwd = 2)
axis(side=2, lwd = 2)

###########################################
# 1c. Compare the distribution of 
#     1a and 1b, what do you find?
###########################################
# increase N
N     <- c(10, 50, 100, 1000)
Nchar <- as.character(N)

# intiate list to store random numbers
listBinom <- list()
listGauss <- list()

# generate binomial and Gaussian random numbers with different N
for (idx in 1:length(N)){
    # Binomial distribution
    listBinom[[Nchar[idx]]] <- 
        myBinom(N[idx], size=40, prob=0.5)
    # Gaussian distribution
    listGauss[[Nchar[idx]]] <- 
        myGauss(N[idx], mean(x), var(x))
} # end for loop

# set plot
par(mfrow=c(2,4))

# plot binomial distribution
for (idx in Nchar){
    hist(listBinom[[idx]], 
         col = "#b3cde3",
         main=paste0("Binomial Distr.\nN=", idx),
         xlab = "Value")
} # end for loop

# plot Gaussian distribution
for (idx in Nchar){
    hist(listGauss[[idx]], 
         col = "#fbb4ae",
         main=paste0("Gauss Distr.\nN=", idx),
         xlab = "Value")
} # end for loop

##################################################################
# 2. Make a program that can select our candidates for 
#    presentation next week. This program should select randomly 
#    but avoid selecting the numbers that had been selected before. 
#    (hint: uniform distribution, loop, and if statements.) 
##################################################################
# create a function to pick two numbers
Get2NumRand <- function(ids){
    # pick a student using uniform distribution
    idx <- myUnif(1, 0, length(ids))
    idx <- ceiling(idx)
    id01 <- ids[idx]
    
    # remove the first chosen student
    ids <- setdiff(ids, id01)
    
    # pick another student from the remain students
    idx <- myUnif(1, 0, length(ids))
    idx <- ceiling(idx)
    id02 <- ids[idx]
    
    # return result
    return(c(id01, id02))
} # end func 

# number of simulation
N   <- 100000 # 100,000

# suppose we have 18 students
ids <- 1:18

# simulation
res <- replicate(N, Get2NumRand(ids))

# initiate a matrix to store the result
mat <- matrix(0, nrow=length(ids), ncol=length(ids))

# convert the result into a matrix
for (idx in 1:N){
    mat[res[1,idx], res[2,idx]] <- 
        mat[res[1,idx], res[2,idx]] + 1
} # end for loop


par(mfrow=c(1,3), 
    #mar=c(5.1,4.1,4.1,2.1))
    mar=c(10,4.1,8,2))

# plot the frequencies of first student
barplot(table(res[1,]) / N, 
        main="First Student", 
        ylim=c(0,0.2))

# plot the frequencies of second student
barplot(table(res[2,]) / N, 
        main="Second Student", 
        ylim=c(0,0.2))

# plot the counting result
plotMatrix(mat, isXlab = F, isYlab = F)
title(
    main="student selection",
    xlab="Second Student",
    ylab="First Student")