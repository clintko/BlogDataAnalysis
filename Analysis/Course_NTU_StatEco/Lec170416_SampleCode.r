workdir <- "C:\\Users\\clint\\Documents\\GitHub\\BlogDataAnalysis"
setwd(workdir)

filePathData <- file.path(
    workdir,
    "Data/StatCompEcology/Copepod")

envdata = read.table(
    file.path(
        filePathData,
        'enviANDdensity.txt'),header=T)

#########################################################
install.packages("vegan")
install.packages("cluster")
install.packages("fpc")

#install.packages("TWIX")
#Note:170417
#Package ‘TWIX’ was removed from the CRAN repository.

library(vegan)
library(cluster)
library(fpc)
#library(TWIX)

##Conduct a Nonhierarchical Clustering (NHC)

#standardize data
data = envdata[,-1]
data.std = scale(data)

#compute distance matrix
dist.eucl = dist(data.std,method='euclidean')

#Determine number of clusters
##Need to write a program for making scree plot and silhouette width

#compute k-means NHC around 5 mediods (clusters)
y.pam = pam(dist.eucl,k=5) 
y.clara = clara(data.std,k=5) 

summary(y.pam)

#evaluate cluster stability
y.eucl.boot = clusterboot(data.std, B=100, metric='euclidean', bootmethod=c('boot','subset'), clustermethod=claraCBI, usepam=TRUE, k=5, count=FALSE)
print(y.eucl.boot)			


##Conduct a Hierarchical Clustering (HC)

y.eucl.ward = hclust(dist.eucl,method='ward.D')  #two way of doing HC
y.eucl.ward = agnes(dist.eucl,method='ward')
y.eucl.dia = diana(dist.eucl)

#dendrogram
plot(y.eucl.ward,main='Wards-linkage Dendrogram',xlab='Station',labels=envdata[,1])
rect.hclust(y.eucl.ward,k=5) #cut by #clusters
rect.hclust(y.eucl.ward,h=5) #cut by height

#aggolmerative coefficient
#y.eucl.ave = hclust(dist.eucl,method='average') # Edit 170417
y.eucl.ave = agnes(dist.eucl, method='average') # Edit 170417
summary(y.eucl.ave)$ac 

#cophenetic correlation
cor(dist.eucl,cophenetic(y.eucl.ave))
plot(dist.eucl,cophenetic(y.eucl.ave),xlab="Original distances",ylab="Cophenetic distances")

