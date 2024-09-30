# mixed normal
rm(list=ls())

setwd("~/density_simulation/")

source("cSScGGdensityEstimationFunctions.R")

library(huge);library(gaussDiff);library(MASS);library(ks);library(QUIC);
library(Matrix);library(nortest);library(mvtnorm)

n <- 200; p<-20; px <- 3; maxLambda <- maxTheta <- 0.5; 
iterLambda <- iterTheta <- 5; alpha <- 2;
pb <- 0.2; # determine sparsity of precision matrix
perc <- .5; sdx <- .5 # perc denotes "omega", and sdx denotes "sigma"
gridsize=0.1

nsim <- 3
includeQ <- len <- rep(NA, nsim)
runResult <- 
runResultLookl <- 
cSScGGResultKLX <- 
runResult_nt <- 
runResultLookl_nt <- 
runResultQuic <- 
MleResult <- 
SkdeResult <- 
SkdeResultKLX <- NULL

set.seed(123)
for(index in 1:nsim){
# generate data
Lxy <- huge.generator(n, d=(p+px), graph = "random", prob = pb, verbose = FALSE)
SigmaTotal <- Lxy$sigma
Sigmax <- SigmaTotal[c((p+1):(p+px)),c((p+1):(p+px))]
# X follows Mixed Gaussian distribution
d1  <- 1; d2 <- 1; d3 <- 1
x1 <- rf1(n,sdx=sdx, d=d1, perc=perc)
x2 <- rf2(n,sdx=sdx, d=d2, perc=perc)
x3 <- rf3(n,sdx=sdx, d=d3, perc=perc)
centers <- c(perc*d1, -(1-perc)*d2, (1-perc)*d3-perc*d3)#true centers of x1, x2 and x3
x <- cbind(x1, x2, x3)
x <- scale(x, center=centers, scale=FALSE)

# construct Y
tF <- as.matrix(Lxy$omega)[c(1:p),c((p+1):(p+px))];
tF[abs(tF) < 1e-5] <- 0
Sig <- as.matrix(Lxy$sigma)[c(1:p),c(1:p)]
Sig[abs(Sig) < 1e-5] <- 0

y <- GenerateData(x, n, p, px, Sig, tF)

#true matrices
CTotal <- as.matrix(Lxy$theta)
tC <- CTotal[c(1:p), c(1:p)]; 
tCx <- CTotal[c((p+1):(p+px)),c((p+1):(p+px))]
PTotal <- Lxy$omega
Py <- PTotal[c(1:p), c(1:p)]; Pyx <- PTotal[c(1:p),c((p+1):(p+px))];
Px <- PTotal[c((p+1):(p+px)),c((p+1):(p+px))];


# run each algorithm and get KL distance of f(y|x) and f(x).

#############################
####        cSScGG       ####
#############################

cat("Fitting cSScGG + BIC/CV for y|x... \n")
####        cSScGG_CV       ####
runResult <- c(runResult, list(runCSSCGG(y, x, n, p, px, tF, tC, tCx, Py, Pyx, Px, maxLambda, maxTheta, 
                     iterLambda, iterTheta, method="alasso")))
#runResult # KL distance of y|x

####        cSScGG_LOOKL       ####
runResultLookl <- c(runResultLookl, list(runCSSCGG_Lookl(y, x, n, p, px, tF, tC, tCx, Py, Pyx,maxLambda, maxTheta, 
                           iterLambda, iterTheta, method="alasso", alpha)))
#runResultLookl # KL distance of y|x

cat("Fitting cSScGG: use smoothing Spline for x... \n")
# It's possible (though rare) to see a negative estimated kl distance when grid size is large
KERNELESTIMATIONGRIDSIZE = 0.1

cSScGGResultKLX <- c(cSScGGResultKLX, list(runSs(x,n,centers, Sigmax,gridsize=KERNELESTIMATIONGRIDSIZE,d1, d2, d3, perc, sdx)
))#cSScGGResultKLX # KL distance of x


#############################
####        cSScGG_NT    ####
#############################
cat("Fitting cSScGG + CV + Normality Test for y|x... \n")

# decide non-normal part by normal test
z <- cbind(y,x)
colnames(z)=c(paste0("y",as.character(1:p)), paste0("x",as.character(1:px)))

p1=NULL
for (j in 1:(px+p)) {
  #p1=c(p1, shapiro.test(z[,j])$p.value)
  p1=c(p1, ad.test(z[,j])$p.value)
}
sorted.names <- colnames(z)[order(p1)]
xnames <- sorted.names[sort(p1) < 0.05/(px+p)]
includeQ[index] <- sum(is.element(c("x1", "x2", "x3"), xnames))==3
len[index] <- length(xnames)
print(xnames)


detect.xidx <- (1:(px+p))[p1 < 0.05/(px+p)]#(1:(px+p))[order(p1)][1:3]#
detect.x <- z[,detect.xidx, drop=F]
detect.y <- z[,-detect.xidx]


# Estimating KL without known Gaussian part is time consuming. 
NSAMPLE <- 1e5

runResult_nt <- c(runResult_nt, list(runCSSCGG_nt(detect.y, detect.x, n, p, px,  centers,SigmaTotal, maxLambda, maxTheta, 
                             iterLambda, iterTheta, method="alasso", detect.xidx, NSAMPLE)))
#runResult_nt # KL distance of z

####        cSScGG_LOOKL_NT       ####
runResultLookl_nt <- c(runResultLookl_nt, list(runCSSCGGLookl_nt(detect.y, detect.x, n, p, px,  centers, SigmaTotal, maxLambda, maxTheta, 
                                       iterLambda, iterTheta, method="alasso", alpha, detect.xidx, NSAMPLE)
))#runResultLookl_nt # KL distance of z


#############################
####        QUIC         ####
#############################
cat("Fitting Quic... \n")
# BIC is used to choose hyperparameters
# In practice we did not observe obvious performance difference between BIC and CV results
PARAMETRICESTIMATIONGRIDSIZE = 0.1
runResultQuic <- c(runResultQuic, list(runQuic(y, x, x, n, p, px, tF, tC, tCx, Py, Pyx, Px, centers, Sigmax, maxLambda, 
                                     iterLambda,gridsize=PARAMETRICESTIMATIONGRIDSIZE,d1, d2, d3, perc, sdx)
))#runResultQuic # KL distance of y|x and x

#############################
####        MLE          ####
#############################
cat("Fitting MLE... \n")
MleResult <- c(MleResult, list(runMle(y, x, n, p, px, Py, Pyx, centers, Sigmax)))
#MleResult # KL distance of y|x and x

#############################
####        SDKE         ####
#############################
cat("Fitting skde parametric part: use MLE to estimate y|x...  \n")
SkdeResult <- c(SkdeResult, list(runSkdeParametric(y, x, n, p, px, Py, Pyx)))
#SkdeResult # KL distance of y|x
cat("Fitting skde nonparametric part: use kernel density to estimate x and return KL distance ...  \n")
KERNELESTIMATIONGRIDSIZE = 0.05
SkdeResultKLX <- c(SkdeResultKLX, list(runKde(x,n,centers, Sigmax,gridsize=KERNELESTIMATIONGRIDSIZE,d1, d2, d3, perc, sdx)
))#SkdeResultKLX # KL distance of x
}


#dump(c("runResult", "cSScGGResultKLX", "runResultLookl", "runResult_nt", "runResultLookl_nt", "runResultQuic",
#       "MleResult", "SkdeResult", "SkdeResultKLX"), "caseSDX0.5Perc0.5.R")
#source("caseSDX0.5Perc0.5.R")

###X###
mean(unlist(cSScGGResultKLX))
sqrt(var(unlist(cSScGGResultKLX)))

#SKDE
mean(unlist(SkdeResultKLX))
sqrt(var(unlist(SkdeResultKLX)))

#QUIC
mean(unlist(lapply(runResultQuic, function(x) x[2])))
sqrt(var(unlist(lapply(runResultQuic, function(x) x[2]))))

#MLE
mean(unlist(lapply(MleResult, function(x) x[2])))
sqrt(var(unlist(lapply(MleResult, function(x) x[2]))))


#####Y|X######
mean(unlist(lapply(runResult, function(x) x[2])))
sqrt(var(unlist(lapply(runResult, function(x) x[2]))))

mean(unlist(runResultLookl))
sqrt(var(unlist(runResultLookl)))

#SKDE
mean(unlist(SkdeResult))
sqrt(var(unlist(SkdeResult)))

#QUIC
mean(unlist(lapply(runResultQuic, function(x) x[1])))
sqrt(var(unlist(lapply(runResultQuic, function(x) x[1]))))

#MLE
mean(unlist(lapply(MleResult, function(x) x[1])))
sqrt(var(unlist(lapply(MleResult, function(x) x[1]))))

####TOTAL####
mean(unlist(cSScGGResultKLX))+mean(unlist(lapply(runResult, function(x) x[2])))
sqrt(var(unlist(cSScGGResultKLX)+unlist(lapply(runResult, function(x) x[2]))))


mean(unlist(runResult_nt))
sqrt(var(unlist(runResult_nt)))

mean(unlist(cSScGGResultKLX))+mean(unlist(runResultLookl))
sqrt(var(unlist(cSScGGResultKLX)+unlist(runResultLookl)))

mean(unlist(runResultLookl_nt))
sqrt(var(unlist(runResultLookl_nt)))

#SKDE#
mean(unlist(SkdeResultKLX))+mean(unlist(SkdeResult))
sqrt(var(unlist(SkdeResultKLX)+unlist(SkdeResult)))

#QUIC
mean(unlist(lapply(runResultQuic, function(x) x[2])))+mean(unlist(lapply(runResultQuic, function(x) x[1])))
sqrt(var(unlist(lapply(runResultQuic, function(x) x[2]))+unlist(lapply(runResultQuic, function(x) x[1]))))

#MLE
mean(unlist(lapply(MleResult, function(x) x[2])))+mean(unlist(lapply(MleResult, function(x) x[1])))
sqrt(var(unlist(lapply(MleResult, function(x) x[2]))+unlist(lapply(MleResult, function(x) x[1]))))



