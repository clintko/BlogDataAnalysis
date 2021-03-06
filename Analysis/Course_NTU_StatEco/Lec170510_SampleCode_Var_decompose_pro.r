#Variance partitioning
library(vegan)

moths = read.csv('moths.csv',header=T)
hab.full = read.csv('moths.full.csv',header=T)
hab.plot = read.csv('moths.plot.csv',header=T)
hab.patch = read.csv('moths.patch.csv',header=T)
hab.land = read.csv('moths.land.csv',header=T)
hab.space = read.csv('moths.space.csv',header=T)

y = moths[,-1]       #exclude site id
x.plot = hab.plot[,c(2,9,10)]
x.patch = hab.patch[,c(4,6,7)]
x.land = hab.land[,c(3,8,9)]
x.space = hab.space[,c(3,7,8)]

y.log = log(y+1)    #log-transform y
y.chord = data.frame(matrix(0,nrow(y.log),ncol(y.log)))    #row normalization
rsum = apply(y.log,1,sum)
for(i in 1:nrow(y.log)){
	y.chord[i,] = y.log[i,]/rsum[i]
}

#RDA with plot variables after partialling out patch, landscape, and space effects (modify z matrix for partialling out other effects)
z = cbind(x.patch,x.land,x.space)
z.plot = rda(y.chord,x.plot,z)

summary(z.plot)
plot(z.plot,choices=c(1,2),display=c('wa','sp','bp'),scaling=2) 
