#################################################
#####          Set Environemnt              #####
#################################################
library(readxl)

# set work space
workdir <- "C:\\Users\\clint\\Documents\\GitHub\\BlogDataAnalysis"
setwd(file.path(workdir, "Analysis"))

# path for data
filePathData <- file.path(
    workdir,
    "Data\\StatCompEcology\\Copepod")

#################################################
#####    Import and Summary the Data        #####
#################################################
# read in data from excel file
dat <- read_excel(
    file.path(filePathData, "enviANDdensity.xls"),
    sheet=1) 
#head(dat)

# get the density of fish and copepod
fishDens <- dat$`FishDensity (ind./1000m3)`
cpodDens <- dat$`CopepodDensity ind./m3)`

# copepod densities
res <- list()
res$Val    <- cpodDens
res$Mean   <- mean(res$Val)
res$SD     <- sd(res$Val)
res$SE     <- res$SD/length(res$Val)^0.5
res$Median <- median(res$Val)
cpod <- res

# fish densities
res <- list()
res$Val    <- fishDens
res$Mean   <- mean(res$Val)
res$SD     <- sd(res$Val)
res$SE     <- res$SD/length(res$Val)^0.5
res$Median <- median(res$Val)
fish <- res
#-------------------------------
cat("",
    "Mean of Copepod Density", cpod$Mean, "\n",
    "SE   of Copepod Density", cpod$SE,   "\n",
    "Mean of Fish    Density", fish$Mean, "\n",
    "SE   of Fish    Density", fish$SE)


#################################################
######    Setting all the functions        ######
#################################################
# function to perform regression analysis
RegBeta <- function(x, y){
    # Function RegBeta return the coefficient of 
    # linear regression
    # x: independent variable
    # y: dependent variable
    ###################################
    # construct Design matrix
    X <- as.matrix(cbind(1, x)) 
    # solve by projecting to the solution space
    b <- solve(t(X) %*% X) %*% t(X) %*% y
    # return the parameters
    return(b)
} # end func RegBeta
#==================================
# Calculate & Store 
# parameters of linear regression
beta <- list()
beta$Val <- RegBeta(x=cpodDens, y=fishDens)
beta$b0  <- beta$Val[1]
beta$b1  <- beta$Val[2]

# Function for probability distribution
myUnif <- function(n, a=0, b=1) {
    # Continuous Uniform Distribution 
    # Default: R.V. ~ Uniform( [0,1] )
    #====================================
    x <- runif(n)      # R.V. ~~~ Uniform( [0,1] )
    y <- x * (b-a) + a # Shift to Uniform( [a,b] )
    return(y)
} # end func myUnif

# Function for Sampling
mySample <- function(x){
    # The function is based on 
    # return sample of x with replacement
    #===========================================
    # generate random numbers following uniform distribution
    n <- length(x)
    id <- myUnif(n, a=0, b=n)
    # convert continuous value into discrete
    id <- ceiling(id)
    # remove 0 (very rare case, p(X=0) is zero 
    # for any continuous random variable)
    id[id==0] <- 1 
    
    # return the result
    return(x[id])
} # end mySample

# Function for Bootstrapping
myBootstrap <- function(x, f, N){UseMethod("myBootstrap")}

myBootstrap.data.frame <- function(x, f, N){
    print("The input data is a data.frame...")
    # get id, the idx of rows
    id  <- seq_len(nrow(x))
    # bootstrapping
    print("Performing bootstrapping...")
    val <- replicate(
        N, f(x[mySample(id),])
    ) # end replicate
    
    # summary statistics & return the result 
    print("...Finished")
    res <- list()
    res$Val <- val       # bootstrap pedvalue
    res$MU  <- mean(val) # mean of bootstrapped value
    res$SE  <- sd(val)   # standard error of mean
    return(res)
} # end myJackknife.data.frame

myBootstrap.numeric <- function(x, f, N){
    print("The input data is a data.frame...")
    # get id, the idx of rows
    id  <- seq_len(length(x))
    # bootstrapping
    print("Performing bootstrapping...")
    val <- replicate(
        N, f(x[mySample(id)])
    ) # end replicate
    
    # summary statistics & return the result 
    print("...Finished")
    res <- list()
    res$Val <- val       # bootstrapped value
    res$MU  <- mean(val) # mean of bootstrapped value
    res$SE  <- sd(val)   # standard error of mean
    return(res)
} # end myJackknife.numeric


# Function for Jackknife
myJackknife <- function(x, f){UseMethod("myJackknife")}

myJackknife.data.frame <- function(x, f){
    print("Input is data frame")
    n  <- nrow(x)
    id <- 1:n
    matID <- sapply(
        id, function(x){id[-x]}
    ) # end sapply
    colnames(matID) <- paste0("JK", id)
    
    print("Performing Jackknife...")
    val <- apply(
        matID, 2, 
        function(idx){
            f(
                x[idx,]
            )
        }  # end apply func
    ) # end apply
    
    print("...Finish")
    res <- list()
    res$MatID <- matID     
    res$Val   <- val       # Jackknife value
    res$MU    <- mean(val) # mean of Jackknife value
    # standard error of mean
    res$SE    <- ((n-1) / n * sum((val-mean(val))^2))^0.5
    return(res)
} # myJackknife.data.frame

myJackknife.numeric <- function(x, f){
    print("Input is numeric vector")
    n <- length(x)
    id <- 1:n
    matID <- sapply(
        id, function(x){id[-x]}
    ) # end sapply
    colnames(matID) <- paste0("JK", id)
    
    print("Performing Jackknife...")
    val <- apply(
        matID, 2, 
        function(idx){
            f(
                x[idx]
            )
        }  # end apply func
    ) # end apply
    
    print("...Finish")
    res <- list()
    res$MatID <- matID
    res$Val   <- val
    res$MU    <- mean(val)
    res$SE    <- ((n-1) / n * sum((val-mean(val))^2))^0.5
    return(res)
} # myJackknife.numeric


# plot histogram
PlotHist <- function(res, resResampling, ...){
    # res: a number
    # resBStrap: a list
    ##############################
    hist(resResampling$Val,
         col="grey90",
         xlim=c(range(resResampling$Val)[1] * 0.7,
                range(resResampling$Val)[2] * 1.3),
         ...)
    # add verticle lines
    abline(v=resResampling$MU, lwd=5, col="steelblue1")
    abline(v=res,          lwd=2, col="red")
} # end func PlotHist


# Helper Function to plot the statistics
# # function to plot the deviation
PlotArrows <- function(x1, y1, x2, y2, ...){
    for (idx in seq_along(x1)){
        arrows(x1[idx],y1[idx],
               x2[idx],y2[idx],
               ...)
    } # end for loop
} # end func PlotArrows

# function to plot the summary statistics
PlotStat <- function(MU, SE, ...){
    # plot the figure
    plot(MU, xaxt='n', 
         xlab="", ylab="",
         xlim=c(0,length(MU)*1.3),
         ylim=c(min(MU)-max(SE)*2, 
                max(MU)+max(SE)*2),
         ...)
    
    # add standard error
    PlotArrows(
        seq_along(MU), MU-SE, 
        seq_along(MU), MU+SE,
        code=3,
        length=0.1,
        lwd=2,
        angle=90,
        col='red')
    
    # add mean
    points(seq_along(MU), MU, 
           pch=19, cex=1.0, 
           bg="black")
} # end func PlotStat


#################################################
# Question 1
# Compute the mean and standard error of the mean 
# for the fish and copepod density (all data points) 
# respectively using Jackknife. Plot the histogram of 
# Jackknife means.
#################################################

# Perform Jackknife analysis for mean of density
jKnife_FishMean <- myJackknife(fish$Val, f=mean)
jKnife_CpodMean <- myJackknife(cpod$Val, f=mean)

#plot the results
par(mfrow=c(1,2))
hist(jKnife_FishMean$Val, 
     main="Distribution of\nFish Density Mean (Jackknife)",
     xlab="Value")
hist(jKnife_CpodMean$Val,
     main="Distribution of\nCopepod Density Mean (Jackknife)",
     xlab="Value")


#################################################
# Question 2
# Compute the regression coefficients for 
# fish = beta0+beta1 * copepod and Jackknife SE of 
# beta0 and beta1. Plot the histogram of Jackknife 
# beta0 and beta1
#################################################

# first encapsulate the calculation of regression coefficients
# encapsulate the function
myFunc <- function(dat){
    # the function input a data.frame, which contains two column X and Y
    # perform the linear regression and acquire the coefficients
    res <- RegBeta(x=dat$X, y=dat$Y)
    return(as.vector(res))
} # end myFunc

# encapsulate the data
dat <- data.frame(
    X=cpod$Val,
    Y=fish$Val)

# perform Jackknife analysis for regression
jKnife_Beta <- myJackknife(dat, f=myFunc)
n <- nrow(dat)

res <- list()
val <- jKnife_Beta$Val[1,]
res$Val   <- val
res$MU    <- mean(val)
res$SE    <- ((n-1) / n * sum((val-mean(val))^2))^0.5
jKnife_Beta0 <- res

res <- list()
val <- jKnife_Beta$Val[2,]
res$Val   <- val
res$MU    <- mean(val)
res$SE    <- ((n-1) / n * sum((val-mean(val))^2))^0.5
jKnife_Beta1 <- res


#output the result
par(mfrow=c(1,2))
hist(jKnife_Beta0$Val, 
     main="Distribution of Beta0 (Jackknife)",
     xlab="Value")
hist(jKnife_Beta1$Val,
     main="Distribution of Beta1 (Jackknife)",
     xlab="Value")

#################################################
# Question 3 
# Compare the estimates for Q1 and Q2 obtained from 
# normal theory, bootstrap, and jackknife
#################################################

# Perform boostrapping analysis for mean
bStrap_FishMean <- myBootstrap(fish$Val, f=mean, N=1000)
bStrap_CpodMean <- myBootstrap(cpod$Val, f=mean, N=1000)


# Perform boostrapping analysis for regression
# boostrapping
N <- 1000
tmp <- myBootstrap(x=dat, f=myFunc, N=N)

# store the results: beta0
res <- list()
res$Val <- tmp$Val[1,]
res$MU  <- mean(res$Val)
res$SE  <- sd(res$Val)
bStrap_Beta0 <- res

# store the results: beta1
res <- list()
res$Val <- tmp$Val[2,]
res$MU  <- mean(res$Val)
res$SE  <- sd(res$Val)
bStrap_Beta1 <- res


# compare the result of mean
# plot summary statistics of bootstraped mean
par(mfrow=c(1,1))
PlotStat(
    c(cpod$Mean,   
      bStrap_CpodMean$MU, 
      jKnife_CpodMean$MU, 
      fish$Mean,
      bStrap_FishMean$MU,
      jKnife_FishMean$MU),
    c(cpod$SE,
      bStrap_CpodMean$SE, 
      jKnife_CpodMean$SE, 
      fish$SE,
      bStrap_FishMean$SE,
      jKnife_FishMean$SE))

# add labels
title(main="Statistics of Bootstrap Mean")
axis(1, at=c(1,2,3,4,5,6), labels=FALSE)
text(c(1,2,3,4,5,6)-0.2, par("usr")[3] - 330, 
     labels = c(
         "Copepod\n(Sample)", "Copepod\n(Bootstrap)", "Copepod\n(Jackknife)", 
         "Fish\n(Sample)", "Fish\n(Bootstrap)", "Fish\n(Jackknife)"),
     cex = 0.8,
     srt = 60,
     pos = 1, 
     xpd = TRUE)



# compare the result of regression: Beta0
PlotStat(
    c(bStrap_Beta0$MU,
      jKnife_Beta0$MU),
    c(bStrap_Beta0$SE,
      jKnife_Beta0$SE))

# Add sample median
points(x=0.5, y=beta$b0, pch=20, col="#9932cc")
text(  x=0.5, y=beta$b0 * 0.7,
       cex=0.7, col="Purple",
       "Beta0\nSample")

# add labels
title(main="Statistics of Beta0\n(Bootstrap vs Jackknife)")
axis(1, at=c(1,2), labels=c("Beta0\n(Boostrap)", "Beta0\n(Jackknife"))


# compare the result of regression: Beta1
# plot summary statistics of bootstraped median
PlotStat(
    c(bStrap_Beta1$MU,
      jKnife_Beta1$MU),
    c(bStrap_Beta1$SE,
      jKnife_Beta1$SE))

# Add sample median
points(x=0.5, y=beta$b1, pch=20, col="#9932cc")
text(  x=0.5, y=beta$b1 * 0.8,
       cex=0.7, col="Purple",
       "Beta1\nSample")


# add labels
title(main="Statistics of Beta1\n(Bootstrap vs Jackknife)")
axis(1, at=c(1,2), labels=c("Beta1\n(Boostrap)", "Beta1\n(Jackknife"))





