##### Set  Function #####
getDeri <- function(f, tol=0.000001){
    # the function returns a 
    # derivative function of input f
    #######################################
    newFun <- function(x){
        # dy/dx
        delta <- (f(x+tol) - f(x)) / tol
        return(delta)
    } # end inner newFun
    
    return(newFun)
} # end func 

findRoot <- function(f, x, tol=0.000001){
    # the function return the next x 
    # using Newton-Ralphson method
    #########################################
    return(x - f(x) / getDeri(f, tol=tol)(x))
} # end func findRoot

getRoot <- function(f, x0, N=100){
    # f: input function
    # x0: initial value
    # N: number of iterations
    ########################################
    # initialization
    y0 <- f(x0)
    xi <- x0
    yi <- y0
    px <- c()
    py <- c()
    
    # search for the root
    for (dummyNum in 1:N){
        # store each point (x, y)
        px <- c(px, xi)
        py <- c(py, yi)
        
        # update point (x, y)
        xi <- findRoot(myFun, xi)
        yi <- myFun(xi)
    } # end for loop
    
    # store the result
    res <- list()
    res$x <- px
    res$y <- py
    return(res)
} # end func getRoot

getLine <- function(slope, x, y) {
    # from a slope and a point, we
    # are able to define a line
    #######################################
    newFun <- function(px){
        py  <- y + slope * (px-x)
        return(py)
    } # end inner newFun
    
    return(newFun)
} # end fun getLine

plotTangent <- function(f, x, y, tol=0.000001, ...){
    # this function help to plot 
    # the tangent line at the given point
    #########################################
    m <- getDeri(f, tol=tol)(x)
    abline(a=getLine(m, x, y)(0), b=m, ...)
} # end fun plotLine

##### Test Function #####
myFun <- function(x){
    x^3
} # end func myFun

# choose a range
r <- c(-2, 10)

# sketch the function by points
x <- seq(r[1], r[2], by=0.1)
y <- myFun(x)

# choose a initial point
x0 <- 9
y0 <- myFun(x0)

# search for the root for 5 iteration
N <- 5
res <- getRoot(myFun, x0 = x0, N = N)

# plot the function
plot(x, y, type="l", lwd=2, xlim=c(-0.5, 10))
abline(h=0, col="grey50")
abline(v=0, col="grey50")

# initial point
points(x0, y0, pch=20, col="red")

# show the process of searching a root
for(idx in 1:N){
    x <- res$x[idx]
    y <- res$y[idx]
    plotTangent(myFun, x, y, col="#26a5a5", lwd=2)
    segments(x,0, x,y, col="#c0f0f0", lwd=1, lty=2)
} # end for loop

points(res$x, res$y, col="blue")

##### Question 01   #####
# the function in question 01
# f(x) = sin(x) / x - 0.6
myFun <- function(x){
    sin(x) / x - 0.6
} # end func myFun

# choose a range
r <- c(-100, 100)

# sketch the function by points
x <- seq(r[1], r[2], by=0.1)
y <- myFun(x)

# choose a initial point
x0 <- 5
y0 <- myFun(x0)

# search for the root
N <- 200
res <- getRoot(myFun, x0 = x0, N = N)

# plot the function
plot(x, y, type="l", lwd=2)
abline(h=0, col="grey50")
abline(v=0, col="grey50")

# initial point
points(x0, y0, pch=19, col="red")

# process of searching root
points(res$x, res$y, col="blue", pch=20)

# show value of x and y during the process of searching root; 
# as you can see, after around 160 iterations, the value (x,y)
# reach to a point
par(mfrow=c(1,2))
plot(res$x, type="l", xlab="Generation", ylab="Value", main="Point X")
plot(res$y, type="l", xlab="Generation", ylab="Value", main="Point Y")
par(mfrow=c(1,1))

##### Question 02   #####
# read in the data
filePathData <- "C:\\Users\\clint\\Documents\\GitHub\\BlogDataAnalysis\\Data\\StatCompEcology\\VidalTvsDuration.txt"
dat <- read.table(filePathData, head=TRUE)

# setting the Belehradek's equation
# D = a ( T - \alpha)^b
belEq <- function(x, alpha, a, b=-2.05){
    return(a * (x - alpha)^b)
} # end func myFun

# set a function that calculate the residual 
# sum of squares based on a input dataset
myFun <- function(par, x, y){
    # parameters
    alpha <- par[1]
    a     <- par[2]
    # residual sum of squares
    y2  <- sapply(x, function(val){belEq(val, alpha, a)})
    res <- sum((y2 - y)^2)
    
    return(res)
} # end newFun

# fit the models for each dataset
fitC2 <- optim(par=c(0,0), fn=myFun, x=dat$X.tempearture, y=dat$C2)
fitC3 <- optim(par=c(0,0), fn=myFun, x=dat$X.tempearture, y=dat$C3)
fitC4 <- optim(par=c(0,0), fn=myFun, x=dat$X.tempearture, y=dat$C4)
fitC5 <- optim(par=c(0,0), fn=myFun, x=dat$X.tempearture, y=dat$C5)

# plot the result
plot(
    rep(dat$X.tempearture,4), 
    c(dat$C2, dat$C3, dat$C4, dat$C5), 
    col="red", pch=20,
    xlab="Temperature",
    ylab="Stage duration",
    xlim=c(7, 16),
    ylim=c(0, 30))

x <- seq(8, 15.5, by=0.5)
lines(x, belEq(x, fitC2$par[1], fitC2$par[2]))
lines(x, belEq(x, fitC3$par[1], fitC3$par[2]))
lines(x, belEq(x, fitC4$par[1], fitC4$par[2]))
lines(x, belEq(x, fitC5$par[1], fitC5$par[2]))