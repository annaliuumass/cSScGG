# npn case --  


rm(list=ls())
library(assist)
library(batch)
library(fuzzySim)
library(invgamma)
library(huge)
library(MVN)
library(nortest)
library(QUIC)

parentDir <- "~/edgeselection"
currentDir <- paste0(parentDir, "/edge_selection_withSelection")
setwd(currentDir)
targetLogDir <- paste0(parentDir, "/model_log")

source("cSScGGEdgeSelectionFunctions.R")
load("quinticfun.RData")


cat("Script starting... \n")
n <- 300; p<-20; px <- 3; maxLambda <- maxTheta <- 0.5; 
iterLambda <- iterTheta <- 5; rep <- 3; 
pb <- 0.2;
option <- "nonNormal"; 

cutoff <- 0.03
method <- "alasso"; 

seed=123
set.seed(seed)

includeQ <- len <- rep(NA, rep)

parseCommandArgs()
cat("n:", n, "p:", p, "px:", px, "pb:", pb, "\n")

outfile_name <- paste(Sys.Date(), "-", format(Sys.time(), "%H"), "_n", n, "p", p, "px", px, "pb", pb,
                      "cutRMVTdf10", cutoff, "Npn_v1.5", "method", method, option, "seed", seed, sep = "")



  targetDir <- paste0(parentDir, "/edge_selection_withSelection/model_result/", outfile_name)
  if(!dir.exists(targetDir)) dir.create(targetDir)   
  # SpGSS_QUIC
  TotalFileName <- paste0(targetDir, "/Alg1QuicResultTotal_BIC_", unclass(Sys.time()))
  LambdaFileName <- paste0(targetDir, "/Alg1QuicResultLambda_BIC_", unclass(Sys.time()))
  ThetaFileName <- paste0(targetDir, "/Alg1QuicResultTheta_BIC_", unclass(Sys.time()))
  XFileNameSemi <- paste0(targetDir, "/Alg1QuicResultXSemi_BIC_", unclass(Sys.time()))
  tmpFileName <- paste0(targetDir, "/Alg1QuicResulttmp_BIC_", unclass(Sys.time()))
  
  # QUIC
  TotalFileQuic <- paste0(targetDir, "/QuicResultTotal_BIC_", unclass(Sys.time()))
  LambdaFileQuic <- paste0(targetDir, "/QuicResultLambda_BIC_", unclass(Sys.time()))
  ThetaFileQuic <- paste0(targetDir, "/QuicResultTheta_BIC_", unclass(Sys.time()))
  XFileQuic <- paste0(targetDir, "/QuicResultX_BIC_", unclass(Sys.time()))
  # Npn
  TotalFileNpn <- paste0(targetDir, "/NpnResultTotal_BIC_", unclass(Sys.time()))
  LambdaFileNpn <- paste0(targetDir, "/NpnResultLambda_BIC_", unclass(Sys.time()))
  ThetaFileNpn <- paste0(targetDir, "/NpnResultTheta_BIC_", unclass(Sys.time()))
  XFileNpn <- paste0(targetDir, "/NpnResultX_BIC_", unclass(Sys.time()))

  index <- 1
  # generate model and data
  while (index <= rep) {
    cat("current index: ", index)
     # nonnormal data generation
   
      Lxy <- huge.generator(n, d=(p+px), graph = "random", prob = .2, 
                              verbose = FALSE,vis = FALSE, v=.6)
        
       # true structure
      tF <- as.matrix(Lxy$omega)[c(1:p),c((p+1):(p+px))];
      tF[abs(tF) < 1e-5] <- 0
      tF[1,] <- 0
      tF[,1] <- 0
        
      tC <- as.matrix(Lxy$omega)[c(1:p),c(1:p)]
      tC[abs(tC) < 1e-5] <- 0
     
      x <- t(generate.x(n, chain.len=2000, c1=0.5, c2=10)$x)
      x.scaled <- scale(x)
      tCx <- matrix(c(0,0,1, 0, 0, 1, 1,1, 0), byrow=T, nrow=3)
      CTotal <- rbind(cbind(tC, tF), cbind(t(tF), tCx))
      y <- GenerateData(x.scaled, n, p, px, solve(tC), tF)
      x <- x.scaled
      y <- scale(y)
  
   
      z <- cbind(y,x)
    
      colnames(z)=c(paste0("y",as.character(1:p)), paste0("x",as.character(1:px)))
    
      # marginal test
      p1=NULL
      for (j in 1:(px+p)) {
        #p1=c(p1, shapiro.test(z[,j])$p.value)
        p1=c(p1, ad.test(z[,j])$p.value)
      }

      sorted.names <- colnames(z)[order(p1)]
      xnames <- sorted.names[sort(p1) < 0.05/(p+px)]
    
    
    
      includeQ[index] <- sum(is.element(c("x1", "x2", "x3"), xnames))==3
      len[index] <- length(xnames)
      print(xnames)
  
      x.new <- z[,is.element(colnames(z), xnames), drop=FALSE]
      y.new <- z[,!is.element(colnames(z), xnames), drop=FALSE]
 
      pos <- sapply(colnames(z), function(xx) (1:(p+px))[is.element(colnames(cbind(y.new,x.new)), xx)])
   
      cat("\n Fitting Quic... \n")
      runResultQuic <- runQuic(y, x, n, p, px, tF, tC, tCx, maxLambda, iterLambda)
      runResultQuic
      tCx
      cat("Fitting SpGSS_QUIC...\n")
      runResult <- try(runAlg2(y.new, x.new, n, p, px, tF, tC, tCx, maxLambda, 
                       maxTheta, iterLambda, iterTheta, method="alasso", cutoff, N, pos=pos))
     
      if(class(runResult) == "try-error") next;
    
      QuicResultLambdaBIC <- unlist(runResultQuic$perfLambda.quic)
      QuicResultThetaBIC <- unlist(runResultQuic$perfTheta.quic)
      QuicResultXBIC <- unlist(runResultQuic$perfX.quic)
      QuicResultTotalBIC <- unlist(runResultQuic$perfTotal.quic)
    
    
      Alg1ResultLambdaBIC <- unlist(runResult$perfLambda)
      Alg1ResultThetaBIC <- unlist(runResult$perfTheta)
      Alg1ResultXsemiBIC <- unlist(runResult$semi.perfX)
      Alg1ResultTotalBIC <- unlist(runResult$perfTotal)

   
      z <- cbind(y,x)
      fit.npn = huge.npn(z)
      lambda_seq <- seq(0,maxLambda, length.out = (iterLambda+1))
      out.npn = huge(fit.npn,lambda=lambda_seq,method = "glasso")
      npn.stars = huge.select(out.npn,criterion="ebic")
      out_npn_bic <- as.matrix(npn.stars$opt.icov)
      out_npn_bic <- out_npn_bic[pos, pos]
      Lambda_npn_bic <- out_npn_bic[1:p,1:p]
      Theta_npn_bic <- out_npn_bic[1:p,(p+1):(p+px)]
      x_npn_bic <- out_npn_bic[(p+1):(p+px),(p+1):(p+px)]
      ROC.total.npn <- ROCval(CTotal,out_npn_bic, 1)
      perfTotal.npn <- Calmcc(ROC.total.npn)
     
      ROCLambda.npn <- ROCval(tC,Lambda_npn_bic, 1)
      perfLambda.npn <- Calmcc(ROCLambda.npn)
      ROCTheta.npn <- ROCval(tF, Theta_npn_bic, 2)
      perfTheta.npn <- Calmcc(ROCTheta.npn)
      ROCX.npn <- ROCval(tCx, x_npn_bic, 1)
      
      print(x_npn_bic)
      perfX.npn <- Calmcc(ROCX.npn)
    # record performancesfor Lambda, Theta and X part
      NpnResultLambdaBIC <- unlist(perfLambda.npn)
      NpnResultThetaBIC <- unlist(perfTheta.npn)
      NpnResultXBIC <- unlist(perfX.npn)
      print(NpnResultXBIC)
      NpnResultTotalBIC <- unlist(perfTotal.npn)
     
    # write results
    # Lambda
    write.table(Alg1ResultLambdaBIC, file = LambdaFileName, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(QuicResultLambdaBIC, file = LambdaFileQuic, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(NpnResultLambdaBIC, file = LambdaFileNpn, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    # Theta
    write.table(Alg1ResultThetaBIC, file = ThetaFileName, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(QuicResultThetaBIC, file = ThetaFileQuic, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(NpnResultThetaBIC, file = ThetaFileNpn, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    # X
    write.table(Alg1ResultXsemiBIC, file = XFileNameSemi, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(QuicResultXBIC, file = XFileQuic, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(NpnResultXBIC, file = XFileNpn, append = TRUE, row.names=FALSE,
                col.names=FALSE)
   
    # Total
    write.table(Alg1ResultTotalBIC, file = TotalFileName, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(QuicResultTotalBIC, file = TotalFileQuic, append = TRUE, row.names=FALSE,
                col.names=FALSE)
    write.table(NpnResultTotalBIC, file = TotalFileNpn, append = TRUE, row.names=FALSE,
                col.names=FALSE)

    index <- index + 1
  }

#summarise results - overall comparison
#need to replace the files names with the generated ones after the loop above
Alg1X <- read.table("model_result/2024-09-28-18_n300p20px3pb0.2cutRMVTdf100.03Npn_v1.5methodalassononNormalseed123/Alg1QuicResultTotal_BIC_1727560947.1393")
Alg1X <- matrix(Alg1X[,1], byrow=T, ncol=4)
apply(na.omit(Alg1X[1:rep,2:4]), 2, mean)
sqrt(apply(na.omit(Alg1X[1:rep,2:4]), 2, var))

QuicX <- read.table("model_result/2024-09-28-18_n300p20px3pb0.2cutRMVTdf100.03Npn_v1.5methodalassononNormalseed123/QuicResultTotal_BIC_1727560947.14005")
QuicX <- matrix(QuicX[,1], byrow=T, ncol=4)
apply(na.omit(QuicX[1:rep,2:4]), 2, mean)
sqrt(apply(na.omit(QuicX[1:rep,2:4]), 2, var))

NpnX <- read.table("model_result/2024-09-28-18_n300p20px3pb0.2cutRMVTdf100.03Npn_v1.5methodalassononNormalseed123/NpnResultTotal_BIC_1727560947.14075")
NpnX <- matrix(NpnX[,1], byrow=T, ncol=4)
apply(na.omit(NpnX[1:rep,2:4]), 2, mean)
sqrt(apply(na.omit(NpnX[1:rep,2:4]), 2, var))

