Alg1_alasso = function(y,x,lam_cov,lam_F,EnoughIter,thr)
{
  eps <- 2e-2; prevfx <- 1e15;n <- dim(y)[1]; p <- dim(y)[2]; px <- dim(x)[2];
  smnum <- 1e-7; Cy <- 1/n*t(y)%*%y; Cyx<-1/n*t(y)%*%x; Cx <- 1/n*t(x)%*%x;
  gamma <- .5
  diff_Theta <- diff_Lambda <- 2*thr + 1
  Ipx <- diag(px);
  S <- Cy - Cyx%*%solve(Cx)%*%t(Cyx);
  # initial values
  if(abs(det(S))<smnum)
  {
    Lambda <- solve(S + smnum*diag(p));
  }else{
    Lambda <- solve(S); #theta
  }
  Lambda_init <- Lambda
  # from MLE
  Theta_init <- -Lambda %*% Cyx %*% solve(Cx)  #theta_xy
  # tC%*% Cyx %*% solve(Cx)  true theta_xy
  iter <- 0; Sigma_old <- S; Theta_old <- Theta_init; Lambda_old <- Lambda_init
  sXmatrix <- Cy%*%Lambda_old + Sigma_old%*%Theta_old%*%Cx%*%t(Theta_old)
  trsX <- sum(diag(sXmatrix))
  l1normLambda <- lam_cov * (sum(Lambda_old != 0)-p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -log(det(Lambda_old)) + trsX + l1normLambda + l1normTheta
  record1 <- record2 <- rep(0, EnoughIter)
  while( (iter<EnoughIter) && (abs((fx - prevfx)/fx) > eps) )
  {
    iter=iter+1;
    # cat("iteration at: ", iter, "\n")
    # cat("change in objective function: ", abs((fx - prevfx)/fx), "\n")
    #print("Updating Lambda...")
    dat <- y + x%*%t(Theta_old)%*%Sigma_old
    S_temp <- 1/n*t(dat)%*%dat
    fit.lasso = glasso::glasso(S_temp, lam_cov,maxit=1,penalize.diagonal=FALSE)
    wi.lasso = solve(fit.lasso$w)
    rhomat = lam_cov/p/2 * matrix(1, p, p)/(pmax(abs(wi.lasso)^gamma, 1e-05))
    fit.alasso = glasso::glasso(S_temp, rhomat,maxit=1,penalize.diagonal=FALSE)
    Lambda_new = solve(fit.alasso$w); Sigma_new <- fit.alasso$w
    obF1 <- cdF_nw(Sigma_new,Cx,Cyx,Theta_old,lam_F,1e-3);
    Theta_new <- obF1$F
    sXmatrix_new <- Cy%*%Lambda_new + Sigma_new%*%Theta_new%*%Cx%*%t(Theta_new)
    trsX <- sum(diag(sXmatrix_new))
    l1normLambda <- lam_cov * (sum(Lambda_new != 0)-p)
    l1normTheta <- lam_F * sum(Theta_new != 0)
    fx1 <- -log(det(Lambda_new)) + trsX + l1normLambda + l1normTheta
    # cat("fx1 value:", fx1, "\n")
    prevfx <- fx
    fx <- fx1
    if (is.na(abs((fx - prevfx)/fx))) break;
    Theta_old <- Theta_new;  Lambda_old <- Lambda_new; Sigma_old <- Sigma_new
  }
  list(Theta=Theta_new, Lambda=Lambda_new)
}


# update Theta with active set
cdF_nw <- function(Sigma, Cx, Cyx, Fini, lam, thr_F) {
  diff_F <- 2*thr_F+1; eps <- 1e-7;
  p <- dim(Fini)[1]; px <- dim(Fini)[2];
  F <- matrix(0, nrow=p, ncol=px);			
  Fini0 <- Fini;
  for(i in 1:p)
  {
    for(j in 1:px)
    {
      ac_set <- Cyx[i,j] + (Cx%*%t(Fini)%*%Sigma)[j,i]
      if (abs(ac_set) <= lam/2 && Fini[i,j] == 0) next;
      gij <- -ac_set + Cx[j,j]*Sigma[i,i]*Fini[i,j];
      if(abs(gij) > lam*0.5){
        preabs <- abs(gij) - lam*0.5;
        coef <- max(eps,Cx[j,j]) * max(eps,Sigma[i,i]);
        absF <- preabs/coef;
        F[i,j] <- sign(gij)*absF;
      }else{
        F[i,j] <- 0;
      }
      Fini[i,j] <- F[i,j];
    }
  }
  list(F=F)
}

CallBIC=function(y,x,Lambda,Theta,opt)
{
  n=dim(y)[1]; p=dim(y)[2]; px=dim(x)[2]
  Sy=1/n*t(y)%*%y;
  Syx=1/n*t(y)%*%x;
  Sx=1/n*t(x)%*%x;
  ebic.gamma <- .5
  if(opt==1)   #for glasso
  {
    z <- cbind(y,x)
    S=cov(z);
  }else if(opt==2){
    dat <- y + x%*%t(Theta)%*%solve(Lambda)
    S <- 1/n*t(dat)%*%dat
  }
  Lambda_offd <- Lambda - diag(diag(Lambda))
  sn <- sum(Lambda_offd != 0)
  if(opt==1) {
    neglik <- -n*log(det(Lambda)) + n*sum(diag(S%*%Lambda)) 
    BIC <- neglik + log(n)*sn/2
    ebic <- BIC + 4*ebic.gamma*log(p)
  }
  else if(opt==2) {
    kn <- sum(Theta != 0)
    prod <- Sy%*%Lambda + 2*Syx%*%t(Theta) + solve(Lambda)%*%Theta%*%Sx%*%t(Theta);
    neglik <- -n*log(det(Lambda)) + n*sum(diag(prod))
    BIC <- neglik + log(n)*(dim(Theta)[1] + sn/2 + kn)
    #BIC <- neglik + log(n)*(dim(Theta)[1] + p*px/total_df*sn/2 + p^2/total_df*kn)
    ebic <- BIC + 4*ebic.gamma*log(p)
  }
  list(BIC = BIC, EBIC = ebic, neglik = neglik)
}

ROCval=function(truePara, estimatedPara, whichPara) {
  p <- dim(truePara)[1]; q <- dim(truePara)[2]
  truePara <- ifelse(truePara!=0, 1, 0)
  estimatedPara <- ifelse(estimatedPara!=0, 1, 0)
  
  # for square matrix Lambda and structure within X
  if (whichPara == 1) {
    diag(truePara) <- diag(estimatedPara) <- 1
    P <- sum(truePara != 0) - p
    TP <- sum(truePara * estimatedPara) - p
    N <- p^2 - p - P;
    allOne <- matrix(1, nrow=p, ncol=p)
    FP <- sum(estimatedPara * (allOne - truePara))
  }
  # for Theta
  else if (whichPara == 2) {
    P <- sum(truePara)
    TP <- sum(truePara * estimatedPara)
    N <- p * q - P
    allOne <- matrix(1, nrow = p, ncol = q)
    FP <- sum(estimatedPara * (allOne - truePara))
  }
  
  # calculate TPR & FPR
  if (P > 0) {
    TPR <- TP / P
  } else {
    TPR <- 0
  }
  if (N > 0) {
    FPR <- FP / N
  } else {
    FPR <- 0
  }
  
  list(P=P, N=N, FPR=FPR, TPR=TPR)
}

Calmcc=function(perf)
{
  FPR <- perf$FPR; TPR <- perf$TPR; P <- perf$P; N <- perf$N;
  FP <- N * FPR
  TP <- P * TPR
  if(P == 0) P <- P + 1e-5
  if(TP == 0) TP <- TP + 1e-5
  FN <- P - TP
  TN <- N - FP
  mcc <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
  Sp <- TN / N
  Sn <- TP / P
  F1 <- 2* TP / (2*TP+FN+FP)
  list(mcc=mcc, Sp=Sp, Sn=Sn, F1=F1)
}

# utility functions to generate X from mixed gaussian distribution
rf1 <- function(n,sdx,d,perc) {
  u <- runif(n); x0 <- rnorm(n)
  ifelse(u>perc,x0*sdx,x0*sdx+d)
}

rf2 <- function(n,sdx,d,perc) {
  u <- runif(n); x0 <- rnorm(n)
  ifelse(u>perc,x0*sdx-d,x0*sdx)
}

rf3 <- function(n,sdx,d,perc) {
  u <- runif(n); x0 <- rnorm(n)
  ifelse(u>perc,x0*sdx+d,x0*sdx-d)
}

Mpower=function(Sig,q)
{
  a=svd(Sig);
  d=a$d^(q);
  m=a$u%*%diag(d)%*%t(a$v);
  m	
}

GenerateData=function(x,n,p,px,Sig,F)
{
  mu <- -x%*%t(F)%*%Sig
  v <- matrix(rnorm(n*p),nrow=n);
  y <- mu + v%*%Mpower(Sig,1/2);
  return(y)
}

selectBICForAlg1 <- function(trainy, trainx, Lam_cov, Lam_F, N, N_j, thr, useEBIC = TRUE, method = "quic")
{
  smallnum=1e-4; 	EnoughIter=20;
  n1=dim(trainy)[1]; 	p=dim(trainy)[2];		px=dim(trainx)[2];
  #trainx_mean <- apply(trainx,2,mean)
  #trainx <- sweep(data.matrix(trainx), 2, trainx_mean)
  if (method=="quic") {
    mod <- Alg1(trainy,trainx,Lam_cov,Lam_F,EnoughIter,thr)
  } else if (method=="alasso") {
    mod <- Alg1_alasso(trainy,trainx,Lam_cov,Lam_F,EnoughIter,thr)
  } else if (method=="scad") {
    mod <- Alg1_scad(trainy,trainx,Lam_cov,Lam_F,EnoughIter,thr)
  } else {
    stop("Inappropriate penalty type")
  }
  if(useEBIC) 
  {
    cvob=CallBIC(trainy,trainx,mod$Lambda,mod$Theta,2)$EBIC
  } else {cvob=CallBIC(trainy,trainx,mod$Lambda,mod$Theta,2)$BIC}
  opt_cv_egec=cvob;
  optlam_egec_cov=Lam_cov;
  optlam_egec_F=Lam_F;
  optLambda_egec=mod$Lambda;
  optTheta_egec=mod$Theta;	
  optSig_egec = solve(optLambda_egec)
  Theta2s <- matrix(as.integer(abs(optLambda_egec)>smallnum),nrow=p);
  optadj_egec=Theta2s - diag(diag(Theta2s));
  
  
  bicmatrix <- matrix(NA,N+1,N_j+1)
  count <- 0
  for(i in 1:(N+1))
  {
    lam_cov=Lam_cov/N*(N+1-i);
    for(j in 1:(N_j+1))
    {
      lam_F=Lam_F/N_j*(N_j+1-j);
      # cat("current i = ", i, ", j = ", j, "\n")
      if (method=="quic") { 
        mod2 <- Alg1(trainy,trainx,lam_cov,lam_F,EnoughIter,thr)
      } else if (method=="alasso") {
        mod2 <- Alg1_alasso(trainy,trainx,lam_cov,lam_F,EnoughIter,thr)
      } else if (method=="scad") {
        mod2 <- Alg1_scad(trainy,trainx,lam_cov,lam_F,EnoughIter,thr)
      } else {
        stop("Inappropriate penalty type")
      }
      if(useEBIC) {cvob2=CallBIC(trainy,trainx,mod2$Lambda,mod2$Theta,2)$EBIC;}
      else {cvob2=CallBIC(trainy,trainx,mod2$Lambda,mod2$Theta,2)$BIC;}
      bicmatrix[i,j] <- cvob2
      Theta2 <- mod2$Lambda;
      Theta2s <- matrix(as.integer(abs(Theta2)>smallnum),nrow=p);
      cv_egec <- cvob2;
      if(cv_egec < opt_cv_egec)
      {
        opt_cv_egec=cv_egec;	
        optlam_egec_cov=lam_cov;
        optlam_egec_F=lam_F;
        optTheta_egec= mod2$Theta
        optadj_egec=Theta2s - diag(diag(Theta2s));
        optLambda_egec <- Theta2
      }
      
    }
  }	
  list( optBIC = opt_cv_egec, optlam_cov = optlam_egec_cov, optlam_F = optlam_egec_F,
        opt_Lam = optLambda_egec, opt_Theta = optTheta_egec,
        opt_adj = optadj_egec,
        bicmatrix = bicmatrix)		
}

selectCVForAlg1_new <- function(trainy, trainx, Lam_cov, Lam_F, N, N_j, k, method="quic") {
  EnoughIter <- 10; smallnum=0;
  n <- dim(trainx)[1]
  Theta.list <- list()
  Lambda.list <- list()
  
  #Get overall solution path
  Lam_cov <- round(seq(0+smallnum, Lam_cov, length = N),1);
  Lam_F <- round(seq(0+smallnum, Lam_F, length = N_j),1);
  thr <- 1e-3
  CV <- matrix(NA, ncol = N_j, nrow = N)
  grps <- cut(1 : n, k, labels = FALSE)[sample(n)]
  loglikFold <- matrix(NA, ncol = k, nrow = 1)
  for(i in 1 : N) {
    for(j in 1 : N_j) {
      # cat("cv for i = ", i, "j= ", j, "\n")
      lambda1 <- Lam_cov[i]
      lambda2 <- Lam_F[j]
      #cat("lambda_1=", lambda1, "lambda_2=", lambda2, "\n")
      l2pe <- matrix(NA, ncol = k, nrow = 1)
      for(kth in 1 : k) {
        # cat("current fold: ", kth)
        omit <- which(grps == kth)
        if (method=="quic") {
          modFold <- Alg1(trainy[-omit, ], trainx[-omit, , drop = FALSE], lambda1, lambda2, EnoughIter,thr)
        } else if (method=="alasso") {
          modFold <- Alg1_alasso(trainy[-omit, ], trainx[-omit, , drop = FALSE], lambda1, lambda2, EnoughIter,thr)
        } else if (method=="scad") {
          modFold <- Alg1_scad(trainy[-omit, ], trainx[-omit, , drop = FALSE], lambda1, lambda2, EnoughIter,thr)
        } else {
          stop("Inappropriate penalty type in function selectCVForAlg1_new")
        }
        Theta.list[[kth]] <- as.matrix(modFold$Theta)
        Lambda.list[[kth]] <- as.matrix(modFold$Lambda)
        tmp <- CallBIC(trainy[omit, ], trainx[omit, , drop = FALSE], Lambda.list[[kth]], Theta.list[[kth]], 2)$neglik
        loglikFold[ ,kth] <- tmp
      }
      ## CV
      CV[i,j] <- (apply(loglikFold, 1, sum)) / k
    }
  }
  minIndex <- which(CV == min(CV), arr.ind = TRUE)
  optCV <- min(CV)
  optlam_cov <- Lam_cov[minIndex[1]]
  optlam_F <- Lam_F[minIndex[2]]
  if (method=="quic") {
    finalMod <- Alg1(trainy, trainx, optlam_cov, optlam_F, EnoughIter,thr) 
  } else if (method=="alasso") {
    finalMod <- Alg1_alasso(trainy, trainx, optlam_cov, optlam_F, EnoughIter,thr)
  } else if (method=="scad") {
    finalMod <- Alg1_scad(trainy, trainx, optlam_cov, optlam_F, EnoughIter,thr)
  }
  optLambda_egec <- finalMod$Lambda
  optTheta_egec <- finalMod$Theta
  list(optCV = optCV, optlam_cov = optlam_cov, optlam_F = optlam_F, cvMatrix = CV, 
       opt_Lam = optLambda_egec, opt_Theta = optTheta_egec)	
}

selectCV_QUIC <- function(trainy, trainx, SigmaTotal, Lam_cov, N, k) {
  EnoughIter <- 10; smallnum=0;
  n <- dim(trainx)[1]
  px <- dim(trainx)[2]
  p <- dim(trainy)[2]
  ztotal <- cbind(trainy,trainx)
  Stotal <- cov(ztotal)
  
  #Get overall solution path
  thr <- 1e-3
  CV <- rep(NA,N+1)
  grps <- cut(1 : n, k, labels = FALSE)[sample(n)]
  loglikFold <- matrix(NA, ncol = k, nrow = 1)
  for(i in 1 : (N+1)) {
    #cat("cv for i = ", i, "\n")
    lam_cov=Lam_cov/N*(N+1-i);
    #cat("lambda = ", lambda1, "lambda2 = ", lambda2, "\n")
    for(kth in 1 : k) {
      omit <- which(grps == kth)
      z <- cbind(trainy[-omit, ],trainx[-omit, , drop = FALSE])
      S <- cov(z)
      modFold <- QUIC(S, rho = lam_cov, path = NULL,  msg = 0)
      tmp <- CallBIC(trainy[omit, ],trainx[omit, , drop = FALSE],modFold$X,0,1)$neglik;
      loglikFold[ ,kth] <- tmp
    }
    ## CV
    CV[i] <- (apply(loglikFold, 1, sum)) / k
  }
  minIndex <- which(CV == min(CV), arr.ind = TRUE)
  optCV <- min(CV)
  optlam_cov <- Lam_cov/N*(N+1-minIndex);
  finalMod <- QUIC(Stotal, rho = optlam_cov, path = NULL,  msg = 0)
  optLambda_egec <- finalMod$X
  optTheta_egec <- finalMod$W
  list(optCV = optCV, optlam_cov = optlam_cov, cvMatrix = CV, 
       opt_Lam = optLambda_egec, opt_Theta = optTheta_egec)	
}

getPerf <- function(Lambda, Theta, tC, tF) {
  ROCLambda <- ROCval(tC, Lambda, 1)
  perfLambda <- Calmcc(ROCLambda)
  ROCTheta <- ROCval(tF, Theta, 2)
  perfTheta <- Calmcc(ROCTheta)
  list(perfLambda=perfLambda, perfTheta=perfTheta)
}

calKL_2 <- function(n,mu1, mu2, sig1, sig2) {
  tmp = 0
  for (i in 1:n) {
    # reverse dist 1 dist 2 location
    tmp = tmp + normdiff(mu1=mu2[,i],sigma1=sig2,mu2=mu1[,i],sigma2=sig1,method="KL")[1]
  }
  return(tmp)
}

calKL_oracle <- function(trainy, trainx, n,xs, Py, Pyx, Lam_cov, Lam_F, N, N_j, thr, method = "glasso")
{
  library(gaussDiff)
  mu1 <- -solve(Py) %*% Pyx %*% t(xs)
  sig1 <- solve(Py)
  smallnum=1e-4; 	EnoughIter=10;
  n1=dim(trainy)[1]; 	p=dim(trainy)[2];		px=dim(trainx)[2];
  #trainx_mean <- apply(trainx,2,mean)
  #trainx <- sweep(data.matrix(trainx), 2, trainx_mean)
  if (method=="glasso") {
    mod <- Alg1(trainy,trainx,Lam_cov,Lam_F,EnoughIter,thr)
  } else if (method=="alasso") {
    mod <- Alg1_alasso(trainy,trainx,Lam_cov,Lam_F,EnoughIter,thr)
  } else if (method=="scad") {
    mod <- Alg1_scad(trainy,trainx,Lam_cov,Lam_F,EnoughIter,thr)
  } else {
    stop("Inappropriate penalty type")
  }
  
  optlam_egec_cov=Lam_cov;
  optlam_egec_F=Lam_F;
  optLambda_egec=mod$Lambda;
  optTheta_egec=mod$Theta;	
  optSig_egec = solve(optLambda_egec)
  
  mu2 <- -optSig_egec %*% optTheta_egec %*% t(xs)
  sig2 <- optSig_egec
  cvob <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  opt_cv_egec=cvob;
  
  count <- 0
  for(i in 1:(N+1))
  {
    lam_cov=Lam_cov/N*(N+1-i);
    for(j in 1:(N_j+1))
    {
      lam_F=Lam_F/N_j*(N_j+1-j);
      if (method=="glasso") { 
        mod2 <- Alg1(trainy,trainx,lam_cov,lam_F,EnoughIter,thr)
      } else if (method=="alasso") {
        mod2 <- Alg1_alasso(trainy,trainx,lam_cov,lam_F,EnoughIter,thr)
      } else if (method=="scad") {
        mod2 <- Alg1_scad(trainy,trainx,lam_cov,lam_F,EnoughIter,thr)
      } else {
        stop("Inappropriate penalty type")
      }
      Lambda2 <- mod2$Lambda;
      Theta2 <- mod2$Theta
      mu2 <- -solve(Lambda2) %*% Theta2 %*% t(xs)
      sig2 <- solve(Lambda2)
      cvob2 <- calKL_2(n,mu1, mu2, sig1, sig2)/n
      cv_egec <- cvob2;
      if(cv_egec < opt_cv_egec)
      {
        opt_cv_egec=cv_egec;	
        optlam_egec_cov=lam_cov;
        optlam_egec_F=lam_F;
        optTheta_egec= Theta2
        optLambda_egec <- Lambda2
      }
      
    }
  }	
  list( KL_oracle = opt_cv_egec, optlam_cov = optlam_egec_cov, optlam_F = optlam_egec_F,
        opt_Lam = optLambda_egec, opt_Theta = optTheta_egec)		
}


runCSSCGG <- function(y, x, n, p, px, tF, tC, tCx, Py, Pyx, Px, maxLambda, maxTheta, 
                      iterLambda, iterTheta, method) {
  library(gaussDiff); library(gss); library(mvtnorm)
  mu1 <- -solve(Py) %*% Pyx %*% t(x)
  sig1 <- solve(Py)
  
  cat("Fitting bic... \n")
  # bic
  out_alg1_bic <- selectBICForAlg1(y, x, maxLambda, maxTheta, iterLambda, iterTheta, 1e-3,
                                   useEBIC = FALSE, method = method)
  Theta_alg1_bic <- out_alg1_bic$opt_Theta
  Lambda_alg1_bic <- out_alg1_bic$opt_Lam
  optlam_cov_alg1_bic <- out_alg1_bic$optlam_cov
  optlam_F_alg1_bic <- out_alg1_bic$optlam_F
  covnorm_alg1_bic <- norm(Lambda_alg1_bic-Py, "F")
  Fnorm_alg1_bic <- norm(Theta_alg1_bic-Pyx, "F")
  covnorm2_alg1_bic <- base::norm(Lambda_alg1_bic-Py, "2")
  Fnorm2_alg1_bic <- base::norm(Theta_alg1_bic-Pyx, "2")
  mu2 <- -solve(Lambda_alg1_bic) %*% Theta_alg1_bic %*% t(x)
  sig2 <- solve(Lambda_alg1_bic)
  kl_alg_bic <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  #kl_alg_bic <- calKL_3(y, mu1, mu2, sig1, sig2)/n
  #kl_alg_bic2 <- calKL(y,x,Lambda_alg1_bic,Py,Theta_alg1_bic,Pyx)$trueKL/n
  
  # cv
  cat("Fitting cv...\n")
  out_alg1_cv <- selectCVForAlg1_new(y, x, maxLambda, maxTheta, 
                                     iterLambda, iterTheta, k = 5, method = method)
  Theta_alg1_cv <- out_alg1_cv$opt_Theta
  Lambda_alg1_cv <- out_alg1_cv$opt_Lam
  optlam_cov_alg1_cv <- out_alg1_cv$optlam_cov
  optlam_F_alg1_cv <- out_alg1_cv$optlam_F
  covnorm_alg1_cv <- base::norm(Lambda_alg1_cv-Py, "F")
  Fnorm_alg1_cv <- base::norm(Theta_alg1_cv-Pyx, "F")
  covnorm2_alg1_cv <- base::norm(Lambda_alg1_cv-Py, "2")
  Fnorm2_alg1_cv <- base::norm(Theta_alg1_cv-Pyx, "2")
  mu2 <- -solve(Lambda_alg1_cv) %*% Theta_alg1_cv %*% t(x)
  sig2 <- solve(Lambda_alg1_cv)
  kl_alg_cv <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  
  # evaluate performance
  cat("Evaluating perfomance...\n")
  # bic
  alg1_bic_perf <- getPerf(Lambda_alg1_bic, Theta_alg1_bic, tC, tF)
  # cv
  alg1_cv_perf <- getPerf(Lambda_alg1_cv, Theta_alg1_cv, tC, tF)
  
  # kl oracle
  if (method=="quic") {
    kl_oracle <- calKL_oracle(y, x, n, x, Py, Pyx, maxLambda, maxTheta, iterLambda, 
                              iterTheta, thr=1e-3, method = "glasso")$KL_oracle
  } else kl_oracle <- 0
  
  list(
    # alg1_bic_perf=alg1_bic_perf, 
    # optlam_cov_alg1_bic=optlam_cov_alg1_bic, 
    # optlam_F_alg1_bic=optlam_F_alg1_bic,
    # covnorm_alg1_bic=covnorm_alg1_bic,
    # Fnorm_alg1_bic=Fnorm_alg1_bic, 
    # covnorm2_alg1_bic=covnorm2_alg1_bic,
    # Fnorm2_alg1_bic=Fnorm2_alg1_bic, 
    kl_alg_bic=kl_alg_bic,
    # alg1_cv_perf=alg1_cv_perf,
    # optlam_cov_alg1_cv=optlam_cov_alg1_cv,
    # optlam_F_alg1_cv=optlam_F_alg1_cv, 
    # covnorm_alg1_cv=covnorm_alg1_cv,
    # Fnorm_alg1_cv=Fnorm_alg1_cv, 
    # covnorm2_alg1_cv=covnorm2_alg1_cv,
    # Fnorm2_alg1_cv=Fnorm2_alg1_cv, 
    kl_alg_cv=kl_alg_cv,
    kl_oracle= kl_oracle)
}

LOOKL <- function(y,x,lam_cov,lam_F,alpha,method) {
  n=dim(y)[1]; p=dim(y)[2]; px=dim(x)[2]; EnoughIter <- 10; thr <- 1e-4
  Sy=1/n*t(y)%*%y; Syx=1/n*t(y)%*%x; Sx=1/n*t(x)%*%x;
  z <- cbind(y,x)
  
  if (method=="quic") {
    mod1 <- Alg1(y,x,lam_cov,lam_F,EnoughIter,thr)
  } else if (method=="alasso") {
    mod1 <- Alg1_alasso(y,x,lam_cov,lam_F,EnoughIter,thr)
  } else if (method=="scad") {
    mod1 <- Alg1_scad(y,x,lam_cov,lam_F,EnoughIter,thr)
  } else {
    stop("Inappropriate penalty type in function selectLOOKL_Alg1")
  }
  Lambda <- mod1$Lambda; Theta <- mod1$Theta
  prod <- Sy%*%Lambda + 2*Syx%*%t(Theta) + solve(Lambda)%*%Theta%*%Sx%*%t(Theta);
  neglik <- -log(det(Lambda)) + sum(diag(prod))
  
  SxInv <- solve(Sx)
  LamInv <- solve(Lambda)
  LamInvF <- LamInv %*% Theta
  tLamInvF <- t(LamInvF)
  FLamInvF <- t(Theta) %*% LamInvF
  A <- -2 * kronecker(LamInv, Sx)
  SxFtLamInv <- Sx %*% t(LamInvF)
  B <- 2 * kronecker(LamInv, SxFtLamInv)
  C_1 <- Theta %*% solve(SxInv + (2*FLamInvF)) %*% t(Theta)
  D <- t(B)
  E <- -2 * kronecker(LamInvF, diag(px))
  FF <- kronecker(LamInvF, LamInvF)
  G <- -diag(p^2)
  H <- -2 * diag(px*p)
  CInv <-  -kronecker(Lambda, Lambda) + 2 * kronecker(Lambda, C_1)
  BCInv <- -2 * kronecker(diag(p), Sx %*% t(Theta)) + 4 * kronecker(diag(p), SxFtLamInv %*% C_1)
  DAInv <- - kronecker(diag(p), LamInvF)
  Inv2 <- solve(-A + BCInv%*%D)
  Inv1 <- -CInv + t(BCInv) %*% Inv2 %*% BCInv
  Coef11 <- FF - DAInv %*% E
  Coef12 <- DAInv %*% H
  Coef21 <- E - BCInv %*% FF
  Coef22 <- BCInv %*% G
  
  BiasLookl <- function(z) {
    yk <- z[1:p]; xk <- z[(p+1):length(z)]; d <- length(xk) 
    Syk <- yk %*% t(yk); Sxk <- xk %*% t(xk); Sxyk <- xk %*% t(yk)
    dSx <- cbind(c(Sxk - Sx))
    dSy <- cbind(c(Syk - Sy))
    dSxy <- cbind(c(Sxyk - t(Syx)))
    # if (abs(det(FLamInvF)) == 0) FLamInvFInv <- solve(FLamInvF+ridge*diag(d))
    d1Lam <- as.vector(LamInv - Syk + LamInvF %*% Sxk %*% tLamInvF)
    d1Theta <- as.vector(-2 * Sxyk - 2 * Sxk %*% tLamInvF)
    bias_1 <- d1Lam %*% Inv1 %*% (Coef11 %*% dSx - Coef12 %*% dSxy + G %*% dSy)
    bias_2 <- d1Theta %*% Inv2 %*% (Coef21 %*% dSx - Coef22 %*% dSy + H %*% dSxy)
    bias <- (bias_1 + bias_2)
    return(bias)
  }
  bias <- 1/(n*(n-1)) * sum(apply(z, MARGIN = 1, BiasLookl))
  lookl <- neglik + alpha * bias
  list(lookl=lookl, neglik=neglik, bias=bias)
  #return(lookl)
}

selectLOOKL_Alg1 <- function(trainy,trainx,Lam_cov,Lam_F,N, N_j, method="quic",alpha) {
  library(glasso)
  thr <- 1e-3; EnoughIter <- 10
  opt_lookl <- LOOKL(trainy,trainx,lam_cov=Lam_cov,lam_F=Lam_F,alpha,method=method)
  optlam_cov <- Lam_cov; optlam_F <- Lam_F
  Lam_cov <- seq(0, Lam_cov, length = N+1);
  Lam_F <- seq(0, Lam_F, length = N_j+1);
  lookl_alpha <- list()
  looklMat <- neglikMat <- biasMat <- matrix(NA,N+1,N_j+1)
  for(i in 1:(N+1))
  {
    # process_bar <- txtProgressBar(min = 0, max = N+1, style = 3)
    # Sys.sleep(0.1)
    # setTxtProgressBar(process_bar, i)
    for(j in 1:(N_j+1)) 
    {
      lam_cov=Lam_cov[i]; lam_F=Lam_F[j];
      modellookl <- LOOKL(trainy, trainx, lam_cov=lam_cov,lam_F=lam_F,alpha,method=method)
      tmp_neglik <- modellookl$neglik
      tmp_bias <- modellookl$bias
      neglikMat[i,j] <- tmp_neglik; biasMat[i,j] <- tmp_bias
      # if(tmp < opt_lookl)
      # {
      #   opt_lookl <- tmp
      #   optlam_cov <- lam_cov
      #   optlam_F <- lam_F
      # }
    }
  }
  for(k in 1:length(alpha)) 
  {
    looklMat <- neglikMat + alpha[k] * biasMat
    minIndex <- which(looklMat == min(looklMat), arr.ind = TRUE)
    optlam_cov <- Lam_cov[minIndex[1]]
    optlam_F <- Lam_F[minIndex[2]]
    if (method=="quic") {
      finalMod <- Alg1(trainy, trainx, optlam_cov, optlam_F, EnoughIter,thr)
    } else if (method=="alasso") {
      finalMod <- Alg1_alasso(trainy, trainx, optlam_cov, optlam_F, EnoughIter,thr)
    } else if (method=="scad") {
      finalMod <- Alg1_scad(trainy, trainx, optlam_cov, optlam_F, EnoughIter,thr)
    } else {
      stop("Inappropriate penalty type in function selectLOOKL_Alg1")
    }
    optLambda_egec <- finalMod$Lambda
    optTheta_egec <- finalMod$Theta
    lookl_alpha[[k]] <- list(looklMat, optLambda_egec, optTheta_egec, optlam_cov, optlam_F)
  }
  return(lookl_alpha)
}

runCSSCGG_Lookl <- function(y, x, n, p, px, tF, tC, tCx, Py, Pyx, maxLambda, maxTheta, 
                            iterLambda, iterTheta, method, alpha) {
  library(gaussDiff)
  mu1 <- -solve(Py) %*% Pyx %*% t(x)
  sig1 <- solve(Py)
  # lookl
  cat("Fitting lookl...\n")
  out_alg1_lookl <- selectLOOKL_Alg1(y, x, maxLambda, maxTheta, iterLambda, 
                                     iterTheta, method = method, alpha)
  perf_alpha_list <- list()
  for (i in 1:length(out_alg1_lookl)) {
    Theta_alg1_lookl <- out_alg1_lookl[[i]][[3]]
    Lambda_alg1_lookl <- out_alg1_lookl[[i]][[2]]
    optlam_cov_alg1_lookl <- out_alg1_lookl[[i]][[4]]
    optlam_F_alg1_lookl <- out_alg1_lookl[[i]][[5]]
    cat("Evaluating lookl...\n")
    alg1_lookl_perf <- getPerf(Lambda_alg1_lookl, Theta_alg1_lookl, tC, tF)
    mu2 <- -solve(Lambda_alg1_lookl) %*% Theta_alg1_lookl %*% t(x)
    sig2 <- solve(Lambda_alg1_lookl)
    mu1 <- -solve(Py) %*% Pyx %*% t(x)
    sig1 <- solve(Py)
    tmp = 0
    for (j in 1:n) {
      tmp = tmp + normdiff(mu1=mu2[,j],sigma1=sig2,mu2=mu1[,j],sigma2=sig1,method="KL")[1]
    }
    kl_alg_lookl <- tmp/n
    covnorm_lookl <- norm(Lambda_alg1_lookl-Py, "F")
    Fnorm_lookl <- norm(Theta_alg1_lookl-Pyx, "F")
    covnorm2_lookl <- base::norm(Lambda_alg1_lookl-Py, "2")
    Fnorm2_lookl <- base::norm(Theta_alg1_lookl-Pyx, "2")
    perf_alpha <- list(
      # alg1_lookl_perf=alg1_lookl_perf, 
      # optlam_cov_alg1_lookl=optlam_cov_alg1_lookl,
      # optlam_F_alg1_lookl=optlam_F_alg1_lookl,
      kl_alg_lookl=kl_alg_lookl
      # covnorm_lookl=covnorm_lookl,
      # Fnorm_lookl=Fnorm_lookl,
      # covnorm2_lookl=covnorm2_lookl,
      # Fnorm2_lookl=Fnorm2_lookl
    )
    perf_alpha_list[[i]] <- perf_alpha
  }
  names(perf_alpha_list) <- round(alpha, digits = 2)
  return(perf_alpha_list)
}

calKL_nt <- function(n, p, px, centers, SigmaTotal, fit.x, detect.xidx, Theta, Lambda) {
  cat("sampling data points from f0...\n")
  
  
  # X follows Mixed Gaussian distribution
  d1  <- 1; d2 <- 1; d3 <- 1
  x1 <- rf1(n,sdx=sdx, d=d1, perc=perc)
  x2 <- rf2(n,sdx=sdx, d=d2, perc=perc)
  x3 <- rf3(n,sdx=sdx, d=d3, perc=perc)
  x <- cbind(x1, x2, x3)
  x <- scale(x, center=centers, scale=FALSE)
  # construct Y
  tF <- as.matrix(Lxy$omega)[c(1:p),c((p+1):(p+px))];
  tF[abs(tF) < 1e-5] <- 0
  Sig <- as.matrix(Lxy$sigma)[c(1:p),c(1:p)]
  Sig[abs(Sig) < 1e-5] <- 0
  y <- GenerateData(x, n, p, px, Sig, tF)
  z <- cbind(y, x)
  colnames(z)=c(paste0("y",as.character(1:p)), paste0("x",as.character(1:px)))
  
  kl_n <- kl_nt <- 0
  cat("calculating kl...\n")
  for (i in 1:n) {
    zpt <- z[i,,drop=F]
    xpt <- zpt[1,(p+1):(p+px)]
    ypt <- zpt[1,1:p]
    den_p1 <- (1-perc)*dnorm(xpt[1]+centers[1], 0, sdx) + perc*dnorm(xpt[1]+centers[1], d1, sdx)
    den_p2 <- (1-perc)*dnorm(xpt[2]+centers[2], -d2, sdx) + perc*dnorm(xpt[2]+centers[2], 0, sdx)
    den_p3 <- (1-perc)*dnorm(xpt[3]+centers[3], d3, sdx) + perc*dnorm(xpt[3]+centers[3], -d3, sdx)
    
    log_den_z <- log(den_p1) + log(den_p2) + log(den_p3)+
      dmvnorm(ypt, mean=-xpt%*%t(tF)%*%Sig, sigma = Sig, log = TRUE)
    # hat(f)(z)
    detect.x <- zpt[1, detect.xidx, drop=F]
    if(all(detect.x <=  fit.x$domain[2,]) && all(detect.x >=  fit.x$domain[1,])) {
      detect.y <- zpt[1, -detect.xidx]
      if(length(detect.xidx)==1) den_xfit <- dssden(fit.x, detect.x)
      else      den_xfit <- dssden(fit.x, as.data.frame(detect.x))
      log_den_xfit <- log(den_xfit)
      sig_fit <- solve(Lambda)
      mu_fit <- -detect.x%*%t(Theta)%*%sig_fit
      log_den_yfit <- dmvnorm(detect.y, mean=mu_fit, sigma = sig_fit, log = TRUE)
      log_den_zfit <- log_den_xfit + log_den_yfit
      kl_nt <-  kl_nt + log_den_z - log_den_zfit
      kl_n <- kl_n + 1
    }
  }
  cat("real sampe size = ", kl_n, "\n")
  return(kl_nt/kl_n)
}

runCSSCGG_nt <- function(y, x, n, p, px, centers, SigmaTotal, maxLambda, maxTheta, 
                         iterLambda, iterTheta, method, detect.xidx, nsample) {
  library(gss);
  out_alg1_cv <- selectCVForAlg1_new(y, x, maxLambda, maxTheta, 
                                     iterLambda, iterTheta, k = 5, method = method)
  Theta_alg1_cv <- out_alg1_cv$opt_Theta
  Lambda_alg1_cv <- out_alg1_cv$opt_Lam
  
  mn <- apply(x,2,min);
  mx <- apply(x,2,max)
  domain <- data.frame(rbind(mn, mx))
  x <- data.frame(x)
  if(ncol(x)==1) fit.x <- ssden(~., domain=domain, data=x)
  else fit.x <- ssden1(~(.)^2, domain=domain, data=x)
  
  cat("calculating KL...")
  kl_nt <- calKL_nt(n=nsample, p, px, centers, SigmaTotal, fit.x, detect.xidx, 
                    Theta=Theta_alg1_cv, Lambda=Lambda_alg1_cv)
  return(kl_nt)
}

runCSSCGGLookl_nt <- function(y, x, n, p, px, centers, SigmaTotal, maxLambda, maxTheta, 
                              iterLambda, iterTheta, method, alpha, detect.xidx, nsample) {
  # lookl
  cat("Fitting lookl...\n")
  out_alg1_lookl <- selectLOOKL_Alg1(y, x, maxLambda, maxTheta, iterLambda, 
                                     iterTheta, method = method, alpha)
  mn <- apply(x,2,min);
  mx <- apply(x,2,max)
  domain <- data.frame(rbind(mn, mx))
  x <- data.frame(x)
  if(ncol(x)==1) fit.x <- ssden(~., domain=domain, data=x)
  else
    fit.x <- ssden1(~(.)^2, domain=domain, data=x)
  perf_alpha_list <- list()
  for (i in 1:length(out_alg1_lookl)) {
    Theta_alg1_lookl <- out_alg1_lookl[[i]][[3]]
    Lambda_alg1_lookl <- out_alg1_lookl[[i]][[2]]
    cat("calculating KL...\n")
    kl_nt <- calKL_nt(n=nsample, p, px, centers, SigmaTotal, fit.x, detect.xidx, 
                      Theta=Theta_alg1_lookl, Lambda=Lambda_alg1_lookl)
    perf_alpha_list[[i]] <- kl_nt
  }
  names(perf_alpha_list) <- round(alpha, digits = 2)
  return(perf_alpha_list)
}

runSs <- function(x, n, centers, Sigmax, gridsize,d1, d2, d3, perc, sdx) {
  cat("Estimating x using smoothing spline...\n")
  mn <- apply(x,2,min);
  mx <- apply(x,2,max)
  domain <- data.frame(rbind(mn, mx))
  x <- data.frame(x)
  fit.x <- ssden1(~(.)^2, domain=domain, data=x)
  
  calKL_ss_mixnormal <- function(x, centers, sig1, fit.x, gridsize, d1, d2, d3, perc) {
    x1.wk <- seq(mn[1],mx[1],by=gridsize)
    x2.wk <- seq(mn[2],mx[2],by=gridsize)
    x3.wk <- seq(mn[3],mx[3],by=gridsize)
    tmp <- 0
    for (i in 1:length(x1.wk)) {
      #if (i%%10 == 0) cat("i = ", i, "\n")
      for (j in 1:length(x2.wk)) {
        for (k in 1:length(x3.wk)) {
          # if (k%%100 == 0) cat("j = ", k, "\n")
          xpt <- data.frame(x1=x1.wk[i], x2=x2.wk[j], x3=x3.wk[k])
          den_ss <- dssden(fit.x, xpt)
          den_p1 <- (1-perc)*dnorm(xpt[1,1]+centers[1], 0, sdx) + perc*dnorm(xpt[1,1]+centers[1], d1, sdx)
          den_p2 <- (1-perc)*dnorm(xpt[1,2]+centers[2], -d2, sdx) + perc*dnorm(xpt[1,2]+centers[2], 0, sdx)
          den_p3 <- (1-perc)*dnorm(xpt[1,3]+centers[3], d3, sdx) + perc*dnorm(xpt[1,3]+centers[3], -d3, sdx)
          den_para <- den_p1 * den_p2 *den_p3
          tmp <- tmp + den_para * (log(den_para) - log(den_ss)) * gridsize^3
        }
      }
    }
    return(tmp)
  }
  
  #ptm <- proc.time()
  kl_ss <- calKL_ss_mixnormal(x, centers=centers, sig1 = Sigmax, fit.x = fit.x, grid=gridsize, d1=d1, d2=d2, d3=d3, perc=perc)
  return(kl_ss)
  #proc.time() - ptm
}


selectBIC_QUIC <- function(trainy,trainx,Lam_cov,N,opt) 
{
  library(QUIC)
  smallnum=1e-4; ridge=1e-2; gamma=gama=0.5
  n1=dim(trainy)[1]; 	p=dim(trainy)[2];		px=dim(trainx)[2];
  z <- cbind(trainy,trainx)
  S <- cov(z)
  pz <- p+px
  F <- 0
  bicmatrix <- rep(NA,N+1)
  obj_quic <- QUIC(S, rho = Lam_cov, 
                   path = NULL,  msg = 0)
  Theta2 = obj_quic$X
  Theta2s=matrix(as.integer(abs(Theta2)>smallnum),nrow=pz);
  cvob2=CallBIC(trainy,trainx,Theta2,0,1)$BIC
  opt_cv_egec=cvob2;
  optlam_egec_cov=Lam_cov;
  optTheta_egec=Theta2
  optadj_egec=Theta2s - diag(diag(Theta2s));
  optSig_egec <- obj_quic$W
  for(i in 1:(N+1))
  {
    lam_cov=Lam_cov/N*(N+1-i);
    if(opt==1){
      if(abs(det(S))<smallnum)
      {
        Wei_S=abs(solve(S+ridge*diag(p)))^(-gama);
      }else{			
        Wei_S=abs(solve(S))^(-gama);
      }
      rho1=lam_cov*Wei_S;
      obj1=QUIC(S, rho = rho1, path = NULL,  msg = 0)
    }else{
      obj1=QUIC(S, rho = lam_cov, path = NULL,  msg = 0)
    }
    cvob2=CallBIC(trainy,trainx,obj1$X,0,1)$BIC;
    bicmatrix[i] = cvob2
    Theta2=obj1$X
    Theta2s=matrix(as.integer(abs(Theta2)>smallnum),nrow=pz);
    cv_egec=cvob2;
    if(cv_egec < opt_cv_egec)
    {
      opt_cv_egec=cv_egec;  
      optlam_egec_cov=lam_cov;
      optTheta_egec=Theta2
      optadj_egec=Theta2s - diag(diag(Theta2s));
      optSig_egec <- obj1$W
    }
  }	
  list( optBIC=opt_cv_egec, optlam_cov=optlam_egec_cov,
        opt_Lam = optTheta_egec,
        opt_adj=optadj_egec, opt_Sigma=optSig_egec, bicmatrix=bicmatrix)
}

calKL_normal_mixnormal <- function(x, centers, sig, gridsize, d1, d2, d3, perc, sdx) {
  n <- dim(x)[1]
  mn <- apply(x,2,min);
  mx <- apply(x,2,max)
  sig <- as.matrix(forceSymmetric(sig,uplo="U"))
  tmp = 0
  x1.wk <- seq(mn[1],mx[1],by=gridsize)
  x2.wk <- seq(mn[2],mx[2],by=gridsize)
  x3.wk <- seq(mn[3],mx[3],by=gridsize)
  tmp <- 0
  for (i in 1:length(x1.wk)) {
    #if (i%%10 == 0) cat("i = ", i, "\n")
    for (j in 1:length(x2.wk)) {
      for (k in 1:length(x3.wk)) {
        # if (k%%100 == 0) cat("j = ", k, "\n")
        xpt <- c(x1.wk[i], x2.wk[j], x3.wk[k])
        den_p1 <- (1-perc)*dnorm(xpt[1]+centers[1], 0, sdx) + perc*dnorm(xpt[1]+centers[1], d1, sdx)
        den_p2 <- (1-perc)*dnorm(xpt[2]+centers[2], -d2, sdx) + perc*dnorm(xpt[2]+centers[2], 0, sdx)
        den_p3 <- (1-perc)*dnorm(xpt[3]+centers[3], d3, sdx) + perc*dnorm(xpt[3]+centers[3], -d3, sdx)
        den_para <- den_p1 * den_p2 *den_p3
        log_den_fit <- dmvnorm(xpt, mean=colMeans(x), sigma = sig, log = TRUE)
        tmp <- tmp + den_para * (log(den_para) - log_den_fit) * gridsize^3
      }
    }
  }
  return(tmp)
}

runQuic <- function(y, x, xs, n, p, px, tF, tC, tCx, Py, Pyx, Px, centers, Sigmax, maxLambda, iterLambda,
                    gridsize, d1, d2, d3, perc, sdx) 
{
  library(gaussDiff); library(gss); library(mvtnorm)
  mu1 <- -solve(Py) %*% Pyx %*% t(xs)
  sig1 <- solve(Py)
  
  # QUIC bic
  cat("Fitting QUIC...\n")
  out_quic_bic <- selectBIC_QUIC(y, x, maxLambda, iterLambda, 2)$opt_Lam
  Lambda_quic_bic <- out_quic_bic[1:p,1:p]
  Theta_quic_bic <- out_quic_bic[1:p,(p+1):(p+px)]
  mu2 <- -solve(Lambda_quic_bic) %*% Theta_quic_bic %*% t(xs)
  sig2 <- solve(Lambda_quic_bic)
  kl_quic_bic <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  
  # x part
  # kl quic
  cat("Calculating KL divergence for quic... \n")
  x_quic_bic <- solve(out_quic_bic)[(p+1):(p+px),(p+1):(p+px)]
  kl_x_quic_bic <- calKL_normal_mixnormal(x, centers, x_quic_bic, gridsize, d1, d2, d3, perc, sdx)
  
  list(
    kl_quic_bic=kl_quic_bic,
    kl_x_quic_bic=kl_x_quic_bic)
}

runMle <- function(y, x, n, p, px, Py, Pyx, centers, Sigmax) 
{
  library(gaussDiff); library(mvtnorm)
  mu1 <- -solve(Py) %*% Pyx %*% t(x)
  sig1 <- solve(Py)
  z <- cbind(y, x)
  sigma_mle <- cov(z)*(n-1)/n
  lambda_mle <- solve(sigma_mle)
  
  # from MLE
  Lambda_mle <- lambda_mle[1:p,1:p]
  Theta_mle <- lambda_mle[1:p,(p+1):(p+px)]
  covnorm_mle <- base::norm(Lambda_mle-Py, "F")
  Fnorm_mle <- base::norm(Theta_mle-Pyx, "F")
  covnorm2_mle <- base::norm(Lambda_mle-Py, "2")
  Fnorm2_mle <- base::norm(Theta_mle-Pyx, "2")
  mu2 <- -solve(Lambda_mle) %*% Theta_mle %*% t(x)
  sig2 <- solve(Lambda_mle)
  kl_mle <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  # x part
  x_mle <- sigma_mle[(p+1):(p+px),(p+1):(p+px)]
  kl_x_mle <- calKL_normal_mixnormal(x, centers, x_mle, gridsize, d1, d2, d3, perc, sdx)
  #normdiff(mu1=rep(0,px),sigma1=Sigmax, mu2=rep(0,px),sigma2=x_mle,method="KL")[1]
  
  list(kl_mle=kl_mle,
       # covnorm_mle=covnorm_mle,
       # Fnorm_mle=Fnorm_mle,
       # covnorm2_mle=covnorm2_mle,
       # Fnorm2_mle=Fnorm2_mle,
       kl_x_mle=kl_x_mle)
}

runSkdeParametric <- function(y, x, n, p, px, Py, Pyx) 
{
  library(gaussDiff); library(mvtnorm)
  mu1 <- -solve(Py) %*% Pyx %*% t(x)
  sig1 <- solve(Py)
  
  
  Theta_skde <- matrix(0, nrow=p, ncol=px)
  Sigma_skde <- 1/n * t(y) %*% y
  mu2 <- - Sigma_skde%*% Theta_skde %*% t(x)
  sig2 <- Sigma_skde
  Lambda_skde <- solve(Sigma_skde)
  kl_skde <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  
  covnorm_skde <- norm(Lambda_skde-Py, "F")
  Fnorm_skde <- norm(Theta_skde-Pyx, "F")
  covnorm2_skde <- base::norm(Lambda_skde-Py, "2")
  Fnorm2_skde <- base::norm(Theta_skde-Pyx, "2")
  
  list(kl_skde=kl_skde
       # covnorm_skde=covnorm_skde,
       # Fnorm_skde=Fnorm_skde,
       # covnorm2_skde=covnorm2_skde,
       # Fnorm2_skde=Fnorm2_skde
  )
}

runKde = function(x, n, centers, Sigmax, gridsize,d1, d2, d3, perc, sdx) {
  library(ks)
  cat("Estimating x using kde and output KL distance...\n")
  mn <- apply(x,2,min);
  mx <- apply(x,2,max)
  
  x1.wk <- seq(mn[1],mx[1],by=gridsize)
  x2.wk <- seq(mn[2],mx[2],by=gridsize)
  x3.wk <- seq(mn[3],mx[3],by=gridsize)
  xrow <- length(x1.wk)*length(x2.wk)*length(x3.wk)
  x.wk <- matrix(NA, nrow= xrow, ncol=px)
  idx <- 0
  for (i in 1:length(x1.wk)) {
    for (j in 1:length(x2.wk)) {
      for (k in 1:length(x3.wk)) {
        #if (idx%%10000 == 0) cat("idx = ", idx, "\n")
        idx <- idx + 1
        x.wk[idx,] <-c(x1.wk[i], x2.wk[j], x3.wk[k])
      }
    }
  }
  ptm <- proc.time()
  # try different bandwidth selection
  # H <- Hscv(x) # much better than Hpi
  H <- Hscv.diag(x) 
  den_dke <- kde(x = x, H = H,eval.points = x.wk)$estimate
  proc.time() - ptm
  # den_dke <- kde(x = x)
  
  ptm <- proc.time()
  tmp <- 0
  for (l in 1:xrow) {
    #if (l%%100000 == 0) cat("l = ", l, "\n")
    den_dke_l <- den_dke[l]
    den_p1 <- (1-perc)*dnorm(x.wk[l,1]+centers[1], 0, sdx) + perc*dnorm(x.wk[l,1]+centers[1], d1, sdx)
    den_p2 <- (1-perc)*dnorm(x.wk[l,2]+centers[2], -d2, sdx) + perc*dnorm(x.wk[l,2]+centers[2], 0, sdx)
    den_p3 <- (1-perc)*dnorm(x.wk[l,3]+centers[3], d3, sdx) + perc*dnorm(x.wk[l,3]+centers[3], -d3, sdx)
    den_para <- den_p1 * den_p2 *den_p3
    tmp <- tmp + den_para * (log(den_para) - log(den_dke_l)) * gridsize^3
  }
  proc.time() - ptm
  return(tmp)
}
