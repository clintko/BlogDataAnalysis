setwd("C:\\Users\\clint\\Documents\\GitHub\\BlogDataAnalysis\\Data\\StatCompEcology\\moths_data")

##RDA
library(vegan)

moths = read.csv('moths.csv',header=T)
hab.full = read.csv('moths.full.csv',header=T)
hab.plot = read.csv('moths.plot.csv',header=T)
hab.patch = read.csv('moths.patch.csv',header=T)
hab.land = read.csv('moths.land.csv',header=T)
hab.space = read.csv('moths.space.csv',header=T)

y = moths[,-1]       #exclude site id
x.full = hab.full[,-1]
x.plot = hab.plot[,-1]
x.patch = hab.patch[,-1]
x.land = hab.land[,-1]
x.space = hab.space[,-1]

y.log = log(y+1)    #log-transform y
y.chord = data.frame(matrix(0,nrow(y.log),ncol(y.log)))    #row normalization
rsum = apply(y.log,1,sum)
for(i in 1:nrow(y.log)){
	y.chord[i,] = y.log[i,]/rsum[i]
}

#Determine an appropriate species response model 
y.dca = decorana(y.log)
y.dca
plot(y.dca,type='text')


#Check multicollinearity in explanatory variables
round(as.dist(cor(x.plot)),2)    #found that CH.OAK and SO.M are highly correlated
round(as.dist(cor(x.patch)),2)
round(as.dist(cor(x.land)),2)

rda(y.chord~CH.OAK,data=x.plot)   #Keep the variable with the larger constrained eigenvalue
rda(y.chord~SO.M,data=x.plot)

vif.cca(rda(y.chord~.,data=x.plot))
vif.cca(rda(y.chord~.,data=x.patch))
vif.cca(rda(y.chord~.,data=x.land))

#Stepwise variable selection: forward
y.full = rda(y.chord~.,data=x.plot)
y.red = rda(y.chord~1,data=x.plot)
y.rda = step(y.red,scope=list(lower=~1,upper=formula(y.full)))
y.rda$anova

#Stepwise variable selection: backward
y.rda = step(y.full,scope=list(lower=formula(y.red),upper=formula(y.full)))
y.rda

#Finding the subset of environmental variables so that the Euclidean distances of environmental variable have the max correlation with community dissimilarities
summary(bioenv(y.chord,x.plot,index='euclidean'))
summary(bioenv(y.chord,x.patch,index='euclidean'))
summary(bioenv(y.chord,x.land,index='euclidean'))

#Suppose we select these environmental variables
x.new = x.full[,c(19,21,22,10,12,15,1,6,7)]


#Constrained ordination
y.rda = rda(y.chord~.,data=x.new)
summary(y.rda)

#Total inertia = total sum of variance of species columns
sum(apply(y.chord,2,var))

#Compare differences in WA and LC sample scores
plot(y.rda,choices=c(1,2),type='points',display='lc',scaling=1)
points(y.rda,choices=c(1,2),display='wa',pch=19,scaling=1)

#Species and environment correlation for each axis
spenvcor(y.rda)

#IntER- and IntRA-set correlation (structure coefficients)
intersetcor(y.rda)

intrasector = function(object){
	w = weights(object)
	lc = sweep(object$CCA$u,1,sqrt(w),"*")
	cor(qr.X(object$CCA$QR),lc)
}
intrasector(y.rda)

#Triplot
plot(y.rda,choices=c(1,2),display=c('wa','sp','bp'),scaling=2)  #try scaling=1, 2, or 3

plot(y.rda,choices=c(1,2),type='none',scaling=2)   #Use these lines to plot 1 item at a time
points(y.rda,choices=c(1,2),display='wa',pch=19,cex=.5,scaling=2)
text(y.rda,choices=c(1,2),display='sp',col='red',cex=.75,scaling=2)
text(y.rda,choices=c(1,2),display='bp',col='blue')



















