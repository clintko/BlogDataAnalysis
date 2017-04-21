#Discriminant analysis
herring = read.csv('herring.csv',header=TRUE) 

y = herring[,-c(1:2)] 
grp = herring[,1] 

library(boot) 
library(energy) 
library(MASS) 
library(vegan) 
library(nortest) 

boxplot(SL~grp,data = y)   #check equal variance assumption; decide if need to transform data

par(mfrow=c(1,3))
hist(y$SL[1:134])   #check for normality
hist(y$SL[135:268])
hist(y$SL[269:402])

as.dist(cor(y))   #check for multicollinearity
summary(aov(y$SL~as.factor(grp)))   #go through each variable and drop variables with small F values

pairs(y[3:8])    #check linearity among variables

#Variable selection
y.mat = as.matrix(y) 
summary(lm(y.mat~as.factor(grp)))  #justify variables on statistical and ecological grounds

#Derive discriminant functions
y.lda<-lda(y,grouping=grp) 
y.lda 

#Assessing the importance of the canonical functions
y.lda$svd   #svd^2 gives eigenvalues
y.lda$svd^2/sum(y.lda$svd^2)    #relative magnitude of eigenvalue
y.lda.pred = predict(y.lda)   #classification based on DA
y.lda.pred 
scores = y.lda.pred$x
summary(lm(scores~as.factor(grp)))   #R^2 gives squared canonical correlation coef

x = as.data.frame(cbind(grp,scores))   #visually inspect group difference in canonical scores
boxplot(LD1~grp,data = x,notch=TRUE,varwidth=TRUE)

y.table = table(grp,y.lda.pred$class)   #classification table
y.table
sum(diag(y.table))/sum(y.table)   #correct classification rate


#CART

library(rpart)

z = rpart(grp~.,data=y,method='class',parms=list(split='gini')) 

z = cart(grp~.,data=y,method='class',parms=list(prior=c(.33,.33,.34),split='gini'))  #specify prior

z
summary(z)

plot(z,margin=0.1)

text(z)

1-sum(residuals(z)^2)/sum((grp-mean(grp,na.rm=T))^2,na.rm=T)


#Tree prune

z = rpart(grp~.,data=y,method='class',parms=list(split='gini'),cp=.001) 


cp1 = numeric(50)
for(i in 1:50){
	z = rpart(grp~.,data=y,method='class',parms=list(split='gini'),cp=.001) 
	a = data.frame(printcp(z))
	b = which.min(a[,4])
	cp1[i] = a[b,1]
}
cp1

plotcp(z)

y = numeric(50)+1
t = cbind(cp1,y)
xtabs(y~cp1,t)

zp = prune.rpart(z,0.00373)
plot(zp,margin=0.1)
text(zp)



