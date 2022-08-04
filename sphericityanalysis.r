#functions
library(coda)
library(psych)
library(ggplot2)
library(scales)
Garray.pop<-readRDS("Garray.pop.RDS")
m=4
n=5
MCMCsamp=dim(Garray.pop)[4]


#### Distribution of eigenvalues ####
#Observed values from data
lambda.array<-array(, c(m,n,MCMCsamp))
hpd.lambda<-array(,c(n,2,m))
for (i in 1:m){
  for (j in 1:MCMCsamp){
    lambda.array[i,,j]<-eigen(Garray.pop[,,i,j])$values
  }
}
dimnames(lambda.array)<-list(Gnames, c(paste("e", rep(1:5), sep="")))
                             
                             
for (i in 1:m) {
  for (j in 1:n) {
    hpd.lambda[j,,i]<-HPDinterval(as.mcmc(lambda.array[i,j,]))
    }
  }
mean.lambda<-apply(lambda.array, 1:2, FUN=mean)
hpd.lambdas<-cbind(hpd.lambda[,,1],hpd.lambda[,,2],hpd.lambda[,,3],hpd.lambda[,,4])
#First comparison
#No co-variances and homogenized variances
lambda.r.array<-array(, c(m,n,MCMCsamp))
for (i in 1:m){
  for (j in 1:MCMCsamp){
    null.mat<-diag(tr(Garray.pop[,,i,j])/5,5,5)
    lambda.r.array[i,,j]<-eigen(null.mat)$values
    }
  }
hpd.r.lambda<-array(,c(n,2,m))
for (i in 1:m) {
  for (j in 1:n) {
    hpd.r.lambda[j,,i]<-HPDinterval(as.mcmc(lambda.r.array[i,j,]))
    }
  }
mean.r.lambda<-apply(lambda.r.array, 1:2, FUN=mean)
hpd.r.lambdas<-cbind(hpd.r.lambda[,,1],hpd.r.lambda[,,2],hpd.r.lambda[,,3],hpd.r.lambda[,,4])

#Second comparison
#No co-variances, variances from observed data
lambda.r2.array<-array(, c(m,n,MCMCsamp))
for (i in 1:m){
  for (j in 1:MCMCsamp){
    null.mat<-diag(diag(Garray.pop[,,i,j]))
    lambda.r2.array[i,,j]<-eigen(null.mat)$values
    }
  }
hpd.r2.lambda<-array(,c(n,2,m))
for (i in 1:m) {
  for (j in 1:n) {
    hpd.r2.lambda[j,,i]<-HPDinterval(as.mcmc(lambda.r2.array[i,j,]))
    }
  }
mean.r2.lambda<-apply(lambda.r2.array, 1:2, FUN=mean)
hpd.r2.lambdas<-cbind(hpd.r2.lambda[,,1],hpd.r2.lambda[,,2],hpd.r2.lambda[,,3],hpd.r2.lambda[,,4])

#Third comparison
#Same co-variances as data, homogenized variances
lambda.r3.array<-array(, c(m,n,MCMCsamp))
for (i in 1:m){
  for (j in 1:MCMCsamp){
    null.mat<-Garray.pop[,,i,j]
    for (k in 1:n){
      null.mat[k,k]<-(tr(Garray.pop[,,i,j])/5)
      }
    lambda.r3.array[i,,j]<-eigen(null.mat)$values
    }
  }
hpd.r3.lambda<-array(,c(n,2,m))
for (i in 1:m) {
  for (j in 1:n) {
    hpd.r3.lambda[j,,i]<-HPDinterval(as.mcmc(lambda.r3.array[i,j,]))
    }
  }
mean.r3.lambda<-apply(lambda.r3.array, 1:2, FUN=mean)
hpd.r3.lambdas<-cbind(hpd.r.lambda[,,1],hpd.r.lambda[,,2],hpd.r.lambda[,,3],hpd.r.lambda[,,4])

###plots
#First comparison plot 
plot((1:n)-0.15,   mean.lambda[1,],type="p",xlab="lambda",ylab="value",pch=16,cex=1.2, lwd = 2, col="#311B36",xaxt="n",frame.plot=F, ylim = c(0,1.3), xlim = c(0.8, 5.2), font = 2)
points((1:n)-0.05, mean.lambda[2,],type="p",xlab="",pch=16,cex=1.2, col="#7A6FAE")
points((1:n)+0.05, mean.lambda[3,],type="p",xlab="",pch=16,cex=1.2, col="#C66234")
points((1:n)+0.15, mean.lambda[4,],type="p",xlab="",pch=16,cex=1.2, col="#8F4024")

points((1:n)-0.15, mean.r.lambda[1,],type="p",xlab="",pch=16,cex=1.2, col=alpha("#311B36",0.5))
points((1:n)-0.05, mean.r.lambda[2,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#7A6FAE", 0.5))
points((1:n)+0.05, mean.r.lambda[3,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#C66234", 0.5))
points((1:n)+0.15, mean.r.lambda[4,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#8F4024", 0.5))
axis(1,at=1:n,labels=colnames(mean.lambda), font = 2)
arrows((1:n)-0.15,mean.lambda[1,], (1:n)-0.15, hpd.lambda[,1,1],length=0.05,angle=90,lwd = 2, col="#311B36")
arrows((1:n)-0.05,mean.lambda[2,], (1:n)-0.05, hpd.lambda[,1,2],length=0.05,angle=90,lwd = 2, col="#7A6FAE")
arrows((1:n)+0.05,mean.lambda[3,], (1:n)+0.05, hpd.lambda[,1,3],length=0.05,angle=90,lwd = 2, col="#C66234")
arrows((1:n)+0.15,mean.lambda[4,], (1:n)+0.15, hpd.lambda[,1,4],length=0.05,angle=90,lwd = 2, col="#8F4024")
arrows((1:n)-0.15,mean.lambda[1,], (1:n)-0.15, hpd.lambda[,2,1],length=0.05,angle=90,lwd = 2, col="#311B36")
arrows((1:n)-0.05,mean.lambda[2,], (1:n)-0.05, hpd.lambda[,2,2],length=0.05,angle=90,lwd = 2, col="#7A6FAE")
arrows((1:n)+0.05,mean.lambda[3,], (1:n)+0.05, hpd.lambda[,2,3],length=0.05,angle=90,lwd = 2, col="#C66234")
arrows((1:n)+0.15,mean.lambda[4,], (1:n)+0.15, hpd.lambda[,2,4],length=0.05,angle=90,lwd = 2, col="#8F4024")

arrows((1:n)-0.15,mean.r.lambda[1,], (1:n)-0.15, hpd.r.lambda[,1,1],length=0.05,angle=90,lwd = 2, col=alpha("#311B36",0.5))
arrows((1:n)-0.05,mean.r.lambda[2,], (1:n)-0.05, hpd.r.lambda[,1,2],length=0.05,angle=90,lwd = 2, col=alpha("#7A6FAE",0.5))
arrows((1:n)+0.05,mean.r.lambda[3,], (1:n)+0.05, hpd.r.lambda[,1,3],length=0.05,angle=90,lwd = 2, col=alpha("#C66234",0.5))
arrows((1:n)+0.15,mean.r.lambda[4,], (1:n)+0.15, hpd.r.lambda[,1,4],length=0.05,angle=90,lwd = 2, col=alpha("#8F4024",0.5))
arrows((1:n)-0.15,mean.r.lambda[1,], (1:n)-0.15, hpd.r.lambda[,2,1],length=0.05,angle=90,lwd = 2, col=alpha("#311B36",0.5))
arrows((1:n)-0.05,mean.r.lambda[2,], (1:n)-0.05, hpd.r.lambda[,2,2],length=0.05,angle=90,lwd = 2, col=alpha("#7A6FAE",0.5))
arrows((1:n)+0.05,mean.r.lambda[3,], (1:n)+0.05, hpd.r.lambda[,2,3],length=0.05,angle=90,lwd = 2, col=alpha("#C66234",0.5))
arrows((1:n)+0.15,mean.r.lambda[4,], (1:n)+0.15, hpd.r.lambda[,2,4],length=0.05,angle=90,lwd = 2, col=alpha("#8F4024",0.5))

#second comparison plot
plot((1:n)-0.15,   mean.lambda[1,],type="p",xlab="lambda",ylab="value",pch=16,cex=1.2, lwd = 2, col="#311B36",xaxt="n",frame.plot=F, ylim = c(0,1.3), xlim = c(0.8, 5.2), font = 2)
points((1:n)-0.05, mean.lambda[2,],type="p",xlab="",pch=16,cex=1.2, col="#7A6FAE")
points((1:n)+0.05, mean.lambda[3,],type="p",xlab="",pch=16,cex=1.2, col="#C66234")
points((1:n)+0.15, mean.lambda[4,],type="p",xlab="",pch=16,cex=1.2, col="#8F4024")

points((1:n)-0.15, mean.r2.lambda[1,],type="p",xlab="",pch=16,cex=1.2, col=alpha("#311B36",0.25))
points((1:n)-0.05, mean.r2.lambda[2,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#7A6FAE", 0.25))
points((1:n)+0.05, mean.r2.lambda[3,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#C66234", 0.25))
points((1:n)+0.15, mean.r2.lambda[4,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#8F4024", 0.25))
axis(1,at=1:n,labels=colnames(mean.lambda), font = 2)
arrows((1:n)-0.15,mean.lambda[1,], (1:n)-0.15, hpd.lambda[,1,1],length=0.05,angle=90,lwd = 2, col="#311B36")
arrows((1:n)-0.05,mean.lambda[2,], (1:n)-0.05, hpd.lambda[,1,2],length=0.05,angle=90,lwd = 2, col="#7A6FAE")
arrows((1:n)+0.05,mean.lambda[3,], (1:n)+0.05, hpd.lambda[,1,3],length=0.05,angle=90,lwd = 2, col="#C66234")
arrows((1:n)+0.15,mean.lambda[4,], (1:n)+0.15, hpd.lambda[,1,4],length=0.05,angle=90,lwd = 2, col="#8F4024")
arrows((1:n)-0.15,mean.lambda[1,], (1:n)-0.15, hpd.lambda[,2,1],length=0.05,angle=90,lwd = 2, col="#311B36")
arrows((1:n)-0.05,mean.lambda[2,], (1:n)-0.05, hpd.lambda[,2,2],length=0.05,angle=90,lwd = 2, col="#7A6FAE")
arrows((1:n)+0.05,mean.lambda[3,], (1:n)+0.05, hpd.lambda[,2,3],length=0.05,angle=90,lwd = 2, col="#C66234")
arrows((1:n)+0.15,mean.lambda[4,], (1:n)+0.15, hpd.lambda[,2,4],length=0.05,angle=90,lwd = 2, col="#8F4024")

arrows((1:n)-0.15,mean.r2.lambda[1,], (1:n)-0.15, hpd.r2.lambda[,1,1],length=0.05,angle=90,lwd = 2, col=alpha("#311B36",0.5))
arrows((1:n)-0.05,mean.r2.lambda[2,], (1:n)-0.05, hpd.r2.lambda[,1,2],length=0.05,angle=90,lwd = 2, col=alpha("#7A6FAE",0.5))
arrows((1:n)+0.05,mean.r2.lambda[3,], (1:n)+0.05, hpd.r2.lambda[,1,3],length=0.05,angle=90,lwd = 2, col=alpha("#C66234",0.5))
arrows((1:n)+0.15,mean.r2.lambda[4,], (1:n)+0.15, hpd.r2.lambda[,1,4],length=0.05,angle=90,lwd = 2, col=alpha("#8F4024",0.5))
arrows((1:n)-0.15,mean.r2.lambda[1,], (1:n)-0.15, hpd.r2.lambda[,2,1],length=0.05,angle=90,lwd = 2, col=alpha("#311B36",0.5))
arrows((1:n)-0.05,mean.r2.lambda[2,], (1:n)-0.05, hpd.r2.lambda[,2,2],length=0.05,angle=90,lwd = 2, col=alpha("#7A6FAE",0.5))
arrows((1:n)+0.05,mean.r2.lambda[3,], (1:n)+0.05, hpd.r2.lambda[,2,3],length=0.05,angle=90,lwd = 2, col=alpha("#C66234",0.5))
arrows((1:n)+0.15,mean.r2.lambda[4,], (1:n)+0.15, hpd.r2.lambda[,2,4],length=0.05,angle=90,lwd = 2, col=alpha("#8F4024",0.5))

### Third comparison plot
plot((1:n)-0.15,   mean.lambda[1,],type="p",xlab="lambda",ylab="value",pch=16,cex=1.2, lwd = 2, col="#311B36",xaxt="n",frame.plot=F, ylim = c(0,1.3), xlim = c(0.8, 5.2), font = 2)
points((1:n)-0.05, mean.lambda[2,],type="p",xlab="",pch=16,cex=1.2, col="#7A6FAE")
points((1:n)+0.05, mean.lambda[3,],type="p",xlab="",pch=16,cex=1.2, col="#C66234")
points((1:n)+0.15, mean.lambda[4,],type="p",xlab="",pch=16,cex=1.2, col="#8F4024")

points((1:n)-0.15, mean.r3.lambda[1,],type="p",xlab="",pch=16,cex=1.2, col=alpha("#311B36",0.25))
points((1:n)-0.05, mean.r3.lambda[2,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#7A6FAE", 0.25))
points((1:n)+0.05, mean.r3.lambda[3,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#C66234", 0.25))
points((1:n)+0.15, mean.r3.lambda[4,],type="p",xlab="",pch=16,cex=1.2,alpha = 0.5, col=alpha("#8F4024", 0.25))

axis(1,at=1:n,labels=colnames(mean.lambda), font = 2)

arrows((1:n)-0.15,mean.lambda[1,], (1:n)-0.15, hpd.lambda[,1,1],length=0.05,angle=90,lwd = 2, col="#311B36")
arrows((1:n)-0.05,mean.lambda[2,], (1:n)-0.05, hpd.lambda[,1,2],length=0.05,angle=90,lwd = 2, col="#7A6FAE")
arrows((1:n)+0.05,mean.lambda[3,], (1:n)+0.05, hpd.lambda[,1,3],length=0.05,angle=90,lwd = 2, col="#C66234")
arrows((1:n)+0.15,mean.lambda[4,], (1:n)+0.15, hpd.lambda[,1,4],length=0.05,angle=90,lwd = 2, col="#8F4024")
arrows((1:n)-0.15,mean.lambda[1,], (1:n)-0.15, hpd.lambda[,2,1],length=0.05,angle=90,lwd = 2, col="#311B36")
arrows((1:n)-0.05,mean.lambda[2,], (1:n)-0.05, hpd.lambda[,2,2],length=0.05,angle=90,lwd = 2, col="#7A6FAE")
arrows((1:n)+0.05,mean.lambda[3,], (1:n)+0.05, hpd.lambda[,2,3],length=0.05,angle=90,lwd = 2, col="#C66234")
arrows((1:n)+0.15,mean.lambda[4,], (1:n)+0.15, hpd.lambda[,2,4],length=0.05,angle=90,lwd = 2, col="#8F4024")

arrows((1:n)-0.15,mean.r3.lambda[1,], (1:n)-0.15, hpd.r3.lambda[,1,1],length=0.05,angle=90,lwd = 2, col=alpha("#311B36",0.5))
arrows((1:n)-0.05,mean.r3.lambda[2,], (1:n)-0.05, hpd.r3.lambda[,1,2],length=0.05,angle=90,lwd = 2, col=alpha("#7A6FAE",0.5))
arrows((1:n)+0.05,mean.r3.lambda[3,], (1:n)+0.05, hpd.r3.lambda[,1,3],length=0.05,angle=90,lwd = 2, col=alpha("#C66234",0.5))
arrows((1:n)+0.15,mean.r3.lambda[4,], (1:n)+0.15, hpd.r3.lambda[,1,4],length=0.05,angle=90,lwd = 2, col=alpha("#8F4024",0.5))
arrows((1:n)-0.15,mean.r3.lambda[1,], (1:n)-0.15, hpd.r3.lambda[,2,1],length=0.05,angle=90,lwd = 2, col=alpha("#311B36",0.5))
arrows((1:n)-0.05,mean.r3.lambda[2,], (1:n)-0.05, hpd.r3.lambda[,2,2],length=0.05,angle=90,lwd = 2, col=alpha("#7A6FAE",0.5))
arrows((1:n)+0.05,mean.r3.lambda[3,], (1:n)+0.05, hpd.r3.lambda[,2,3],length=0.05,angle=90,lwd = 2, col=alpha("#C66234",0.5))
arrows((1:n)+0.15,mean.r3.lambda[4,], (1:n)+0.15, hpd.r3.lambda[,2,4],length=0.05,angle=90,lwd = 2, col=alpha("#8F4024",0.5))
                            