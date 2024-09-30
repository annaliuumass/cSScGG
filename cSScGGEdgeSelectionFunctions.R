SoftThreshold <- function(a, kappa) {
  if (a > 0) {
    return(max(0, a-kappa))
  }
  else {
    return(-max(0, -a-kappa))
  }
}

l1OffDiagNorm <- function(M) {
  res <- 0
  row <- dim(M)[1]
  col <- dim(M)[2]
  for (i in 1:row) {
    for (j in i:col){
      if (i < j) res <- res + 2*abs(M[i,j])
      else next;
    }
  }
  return(res)
}

# update Theta in scggm
cdF_cmu <- function(n,Sigma, X, Cx, Cxy, Theta, lam, Q, R, max_Theta_iters=1) 
{
  p=dim(Theta)[2]; px=dim(Theta)[1];
  V <- Theta%*%Sigma
  n=dim(X)[1]
  iter <- 0
  while (iter < max_Theta_iters) {
    for(i in 1:px){
      for(j in 1:p){
        ac_set <- Cxy[i,j] + (Cx%*%V)[i,j]
        if (abs(ac_set) <= lam/2 && Theta[i,j] == 0) next;
        a <- 2 * Sigma[j,j] * Cx[i,i]
        b <- 2 * Cxy[i,j] + 2 * Cx[,i] %*% V[,j]
        c <- Theta[i,j]
        mu <- -c + SoftThreshold(c - b/a, lam/a)
        Theta[i,j] <- Theta[i,j] + mu
        V[i,] <- V[i,] + as.numeric(mu)*Sigma[j,]
      }
    }
    iter <- iter + 1
  }
  Q <- 1/sqrt(n) * X %*% Theta
  R <- Q %*% Sigma
  list(Theta=Theta, Q=Q, R=R)
}

# update Lambda in scggm
cdL_cmu <- function(X, Cx, Cxy, Cy, Theta_old, Lambda_old, logdetLambda_old, Sigma_old, lam, Q_old, R_old, Psi_old, trRtQ_old,
                    max_Lambda_iters=1,alpha=1,sig=1e-4,beta=0.5,max_ls_iters=10) 
{
  p=dim(Theta_old)[2]; px=dim(Theta_old)[1]; Delta <- matrix(0, ncol = p, nrow = p)
  U <- Delta%*%Sigma_old
  iter <- 0
  while (iter < max_Lambda_iters) {
    for(i in 1:p){
      for(j in 1:p){
        if (i == j) {
          a <- Sigma_old[i,i]^2 + 2 * Sigma_old[i,i] * Psi_old[i,i]
          b <- Cy[i,i] - Sigma_old[i,i] - Psi_old[i,i] + Sigma_old[i,]%*%U[,j] + 2*Psi_old[i,]%*%U[,i]
          mu <- -b/a
          Delta[i,j] <- Delta[i,j] + mu
          U[i,] = U[i,] + as.numeric(mu) * Sigma_old[i,]
        } else {
          ac_set <- Cy[i,j] - Sigma_old[i,j] - Psi_old[i,j]
          if (abs(ac_set) <= lam && Lambda_old[i,j] == 0) next;
          a <- Sigma_old[i,j]^2 + Sigma_old[i,i]*Sigma_old[j,j] + Sigma_old[i,i]*Psi_old[j,j] + 2*Sigma_old[i,j]*Psi_old[i,j] + Sigma_old[j,j]*Psi_old[i,i]
          b <- Cy[i,j] - Sigma_old[i,j] + Sigma_old[i,] %*% U[,j] - Psi_old[i,j] + Psi_old[i,]%*%U[,j] + Psi_old[j,]%*%U[,i]
          c <- Lambda_old[i,j] + Delta[i,j]
          mu <- -c + SoftThreshold(c - b/a, lam/a)
          Delta[i,j] <- Delta[i,j] + mu
          U[i,] <- U[i,] + as.numeric(mu)*Sigma_old[j,]
          U[j,] <- U[j,] + as.numeric(mu)*Sigma_old[i,]
        }
      }
    }
    iter <- iter + 1
  }
  # backtracking line search
  success <- FALSE
  Delta <- as.matrix(forceSymmetric(Delta,uplo="U"))
  trGradDelta <- sum(diag(Delta%*%Cy)) - sum(diag(Delta%*%Sigma_old)) - sum(diag(t(Delta)%*%Psi_old)) 
  LambdaPlusDelta <- Lambda_old + Delta
  RHS <- sig*(trGradDelta + lam*(l1OffDiagNorm(LambdaPlusDelta) - l1OffDiagNorm(Lambda_old)));
  for (lsiter in 1:max_ls_iters) {
    Lambda_alpha = Lambda_old + alpha*Delta
    if (!all(eigen(Lambda_alpha)$values > 0)) {
      #cat("Lambda_alpha not pd")
      alpha <- alpha * beta
      next;
    }
    R_new <- Q_old %*% solve(Lambda_alpha)
    trRtQ_new <- sum(diag(t(R_new)%*%Q_old))
    logdetLambda_new <- log(det(Lambda_alpha))
    LHS <- -logdetLambda_new + logdetLambda_old + alpha*sum(diag(Cy%*%Delta)) + 
      trRtQ_new - trRtQ_old + 
      lam*(l1OffDiagNorm(Lambda_alpha) - l1OffDiagNorm(Lambda_old))
    if (LHS <= alpha*RHS) {
      success <- TRUE
      Lambda_old <- Lambda_old + alpha*Delta
      break;
    }
    alpha <- alpha * beta
  }
  #if (!success)
  #cat("Lambda line search failed. \n")
  list(Lambda=Lambda_old, R=R_new, trRtQ=trRtQ_new, logdetLambda=logdetLambda_new)
}

runAlg2 <- function(y, x, n, p, px, tF, tC, tCx,   maxLambda, maxTheta, 
                    iterLambda, iterTheta, method, cutoff, N, pos) {
  # Parametric Algo to identify Y-Y & X-Y
  paraOut <- selectBICForAlg1(y, x, maxLambda, maxTheta, iterLambda, iterTheta, 1e-3,
                              useEBIC = FALSE, method = method)
  Theta <- paraOut$opt_Theta
  Lambda <- paraOut$opt_Lam
  optlam_cov <- paraOut$optlam_cov
  optlam_F <- paraOut$optlam_F
  Delta <- (-1/2) * t(Theta) %*% Mpower(Lambda, -1) %*% Theta 
  print(Delta)
  TDelta <- (-1/2) * t(tF) %*%  Mpower(tC, -1) %*% tF 
  
  px.new <- ncol(x)
  semi.CxxEst <- diag(px.new)
  
  if(px.new > 1 ) 
  {
    # Nonparametric Algo to indentify X-X
    library(gss)
    #cat("Revealing structure of non-parametric part...\n")
    mn <- apply(x,2,min)
    mx <- apply(x,2,max)
    numterm <- px.new + px.new * (px.new - 1) / 2
    domain <- data.frame(rbind(mn, mx))
    # semi parametric estimate
    ##semi.ssres <- ssX(x=data.frame(x), Delta, semiflag = TRUE, cutoff = cutoff, N=N)
    ##semi.CxxEst <- XandX(semi.ssres,diag(x = .5, px))
    
    ##semi.CxxEst <- semi.CxxEst + t(semi.CxxEst)
    
    colnames(x) <- paste("X", 1:px.new, sep="")
    par <- typelist <- vector("list", px.new)
    r <- min <- rep(NA, px.new)
    names(typelist) <- colnames(x)
    for(i in 1:px.new){
      par[[i]] <- list(nphi=2,mkphi=mkphi.quintic,mkrk=mkrk.quintic,
                       env=c(min(x[,i]),max(x[,i]))+c(-1,1)*(max(x[,i])-min(x[,i]))*.05)
      typelist[[i]] <- list('custom', par[[i]])
      r[i] <- diff(par[[i]]$env)
      min[i] <- par[[i]]$env[1]
    }
    
    
    #order of terms for px=3:
    #S1, S1^2, S2, S2^2, S3, S3^2, S1:S2, S1^2:S2
    #S1:S2^2, S1^2:S2^2, S1:S3, S1^2:S3, S1:S3^2, S1^2:S3^2
    #S2:S3, S2^2:S3, S2:S3^2, S2^2:S3^2
    #where S_i=(X_i-par$env[1])/diff(par$env) - 0.5)
    
    ncoef <- 2*px.new+(2*px.new)*(2*px.new-1)/2-px.new 
    coefs <- rep(0, ncoef)
    
    
    interaction.terms <- NULL 
    cter <- 0
    for(i in 1:px.new){
      coefs[2*i-1] <- 0
      for(j in 1:px.new) coefs[2*i-1] <- coefs[2*i-1] + Delta[i,j]*(r[j]/2+min[j])
      coefs[2*i-1] <- 2*coefs[2*i-1]*r[i]
      coefs[2*i] <- Delta[i,i]*r[i]^2
      index.mat <- triMatInd(Delta, lower=FALSE)
      if(px.new > i) 
        for(j in (i+1):px.new)  {
          coefs[2*px.new+cter*4+1] <- 
            2*Delta[index.mat][cter+1]*r[index.mat[cter+1,1]]*r[index.mat[cter+1,2]]
          interaction.terms <- c(interaction.terms, paste("X", index.mat[cter+1,1], ":X", index.mat[cter+1,2], sep=""))
          cter <- cter +1
        }}
    
    print(cter)
    terms <- c(colnames(x), interaction.terms)
    print(terms)
    print(coefs)
    fit.x <- ssden(~.^2, data=as.data.frame(x),  type=typelist)
    print("density estimation done")
    fit.x$d <- fit.x$d + coefs
    
    
    ratios <- rep(NA, px.new*(px.new-1)/2)
    for(i in 1:length(ratios)){
      ratios[i] <- project(fit.x, include=terms[-(px.new+i)])$ratio
    }
    
    print(ratios)
    candidate.terms <- terms[-c(1:px.new)][order(ratios,  decreasing = TRUE)]
    rlen <- length(candidate.terms)
    r <- rep(NA, rlen)
    r[1] <- project(fit.x, include=terms[1:px.new])$ratio
    iter <- 1
    repeat{
      if(r[iter] < cutoff| rlen < 1 | iter >= rlen) break
      r[iter+1] <- project(fit.x, include=c(terms[1:px.new], candidate.terms[1:iter]))$ratio
      iter <- iter + 1
    }            
    r[rlen] <- min(ratios)
    print(r)
    add.row <- as.numeric(sapply(candidate.terms[r >= cutoff], substring, 2,2))
    add.col <- as.numeric(sapply(candidate.terms[r >= cutoff], substring, 5,5))
    print(cbind(add.row, add.col))
    semi.CxxEst[cbind(add.row, add.col)] <- 1
    semi.CxxEst <- semi.CxxEst+t(semi.CxxEst)-diag(px.new)
  }
  print(semi.CxxEst)
  
  
  # CV
  # k <- 5
  # paraOut.cv <- selectCVForAlg1_new(y, x, Py, Pyx, maxLambda, maxTheta,
  #                                iterLambda, iterTheta, k)
  # Theta.cv <- paraOut.cv$opt_Theta
  # Lambda.cv <- paraOut.cv$opt_Lam
  # optlam_cov.cv <- paraOut.cv$optlam_cov
  # optlam_F.cv <- paraOut.cv$optlam_F
  # Delta.cv <- (-1/2) * t(Theta.cv) %*% Mpower(Lambda.cv, -1) %*% Theta 
  # semi.ssres.cv <- ssX(x=data.frame(x), Delta.cv, semiflag = TRUE, cutoff = cutoff, N=N)
  # semi.CxxEst.cv <- XandX(semi.ssres.cv,diag(x = .5, px))
  # semi.CxxEst.cv <- semi.CxxEst.cv + t(semi.CxxEst.cv)
  
  # evaluate performance
  cat("Evaluating perfomance...\n")
  # BIC
  SemiTotal <- rbind(cbind(Lambda, Theta), cbind(t(Theta), semi.CxxEst))
  SemiTotal <- SemiTotal[pos, pos]
  Lambda <- SemiTotal[1:p, 1:p]
  Theta <- SemiTotal[1:p, -(1:p)]
  semi.CxxEst <- SemiTotal[-(1:p), -(1:p)]
  ROCLambda <- ROCval(tC, Lambda, 1)
  perfLambda <- Calmcc(ROCLambda)
  ROCTheta <- ROCval(tF, Theta, 2)
  perfTheta <- Calmcc(ROCTheta)
  semi.ROCX <- ROCval(tCx, semi.CxxEst, 1)
  print("semi.CxxEst")
  print(semi.CxxEst)
  semi.perfX <- Calmcc(semi.ROCX)
  print(unlist(semi.perfX))
  
  
  
  CTotal <- rbind(cbind(tC, tF), cbind(t(tF), tCx))
  
  semi.total <- ROCval(CTotal, SemiTotal, 1)
  perfTotal <- Calmcc(semi.total)
  
  # CV
  # ROCLambda.cv <- ROCval(tC, Lambda.cv, 1)
  # perfLambda.cv <- Calmcc(ROCLambda.cv)
  # ROCTheta.cv <- ROCval(tF, Theta.cv, 2)
  # perfTheta.cv <- Calmcc(ROCTheta.cv)
  # semi.ROCX.cv <- ROCval(tCx, semi.CxxEst.cv, 1)
  # semi.perfX.cv <- Calmcc(semi.ROCX.cv)
  # SemiTotal.cv <- rbind(cbind(Lambda.cv, Theta.cv), cbind(t(Theta.cv), semi.CxxEst))
  # semi.total.cv <- ROCval(CTotal, SemiTotal.cv, 1)
  # perfTotal.cv <- Calmcc(semi.total.cv)
  
  perf.bic <- list(perfTotal=perfTotal, perfTheta = perfTheta, perfLambda = perfLambda, 
                   optlam_cov = optlam_cov, optlam_F = optlam_F,
                   semi.perfX = semi.perfX, tmp=data.frame(r, candidate.terms))
  return(perf.bic)
}

den <- function(x1,x2,x3, c1=0.5, c2=10){
  if(x1 > 1 | x1 < 0 | x2 > 1 | x2 < 0 | x3 > 1 | x3 < 0) return(0)
  else{
    z1=c1*x1-.2*x3+0.1#x1-0.7*x3+0.1#x1-0.3*x3+0.1
    z2=c1*x2-.4*x3+0.1#x2-0.7*x3+0.1#x2-0.2*x3+0.1
    return((exp(-5*(z1-0.3)^2)+2*exp(-5*(z1-0.7)^2))*
             (exp(-50*(z2-0.3)^2)+2*exp(-50*(z2-0.7)^2))*
             exp(-c2*(x3-0.8)^2))
  }}

generate.x <- function(n, chain.len, c1, c2){
  xs <- matrix(NA, 3, chain.len)
  for(i in 1:chain.len){
    x1=x2=0.5
    repeat{
      fun1 <- function(x3) return(-den(x1,x2,x3, c1, c2))
      max.fun1 <- optimize(fun1, interval=c(0,1))
      U <- runif(1)
      x3.pros <- runif(1)
      ratio1 <- -fun1(x3.pros)/(-max.fun1$objective)
      #print(c(U, ratio1))
      if(U <  ratio1) {x3 <- x3.pros
      break}
    }
    repeat{
      p1 <- rbinom(1, 1, 1/3)
      x1 <-  ifelse(p1==1, rnorm(1, 0.3*x3+0.2, 0.1), rnorm(1, 0.3*x3+0.6, 0.1))
      if(x1 < 1 & x1 > 0) break}
    
    repeat{
      p1 <- rbinom(1, 1, 1/3)
      x2 <- ifelse(p1==1, rnorm(1, 0.2*x3+0.2, 0.1), rnorm(1, 0.2*x3+0.6, 0.1))
      if(x2 < 1 & x2 > 0) break} 
    xs[,i] <- c(x1, x2, x3)
  }
  return(list(x=xs[,(chain.len-n+1):chain.len]))
}

runQuic <- function(y, x, n, p, px, tF, tC, tCx,  maxLambda, iterLambda) 
{
  library(gaussDiff); library(gss); library(mvtnorm)
  CTotal <- rbind(cbind(tC, tF), cbind(t(tF), tCx))
  
  # QUIC bic
  cat("Fitting QUIC...\n")
  out_quic_bic <- selectBIC_QUIC(y, x, maxLambda, iterLambda, 2)$opt_Lam
  Lambda_quic_bic <- out_quic_bic[1:p,1:p]
  Theta_quic_bic <- out_quic_bic[1:p,(p+1):(p+px)]
  x_quic_bic <- out_quic_bic[(p+1):(p+px),(p+1):(p+px)]
  
  # QUIC CV
  # cat("Fitting QUIC...\n")
  # k <- 5
  # out_quic_cv <- selectCV_QUIC(y, x, SigmaTotal, maxLambda, iterLambda, k)$opt_Lam
  # Lambda_quic_cv <- out_quic_cv[1:p,1:p]
  # Theta_quic_cv <- out_quic_cv[1:p,(p+1):(p+px)]
  # x_quic_cv <- out_quic_cv[(p+1):(p+px),(p+1):(p+px)]
  
  # evaluate performance
  cat("Evaluating perfomance...\n")
  ROC.total.quic <- ROCval(CTotal,out_quic_bic, 1)
  perfTotal.quic <- Calmcc(ROC.total.quic)
  ROCLambda.quic <- ROCval(tC,Lambda_quic_bic, 1)
  perfLambda.quic <- Calmcc(ROCLambda.quic)
  ROCTheta.quic <- ROCval(tF, Theta_quic_bic, 2)
  perfTheta.quic <- Calmcc(ROCTheta.quic)
  ROCX.quic <- ROCval(tCx, x_quic_bic, 1)
  perfX.quic <- Calmcc(ROCX.quic)
  
  # ROC.total.quic.cv <- ROCval(CTotal,out_quic_cv, 1)
  # perfTotal.quic.cv <- Calmcc(ROC.total.quic.cv)
  # ROCLambda.quic.cv <- ROCval(tC,Lambda_quic_cv, 1)
  # perfLambda.quic.cv <- Calmcc(ROCLambda.quic.cv)
  # ROCTheta.quic.cv <- ROCval(tF, Theta_quic_cv, 2)
  # perfTheta.quic.cv <- Calmcc(ROCTheta.quic.cv)
  # ROCX.quic.cv <- ROCval(tCx, x_quic_cv, 1)
  # perfX.quic.cv <- Calmcc(ROCX.quic.cv)
  
  perf.bic <- list(perfTotal.quic=perfTotal.quic,perfLambda.quic=perfLambda.quic,
                   perfTheta.quic=perfTheta.quic,perfX.quic=perfX.quic)
  # perf.cv <- list(perfTotal.quic.cv=perfTotal.quic.cv,perfLambda.quic.cv=perfLambda.quic.cv,
  #                 perfTheta.quic.cv=perfTheta.quic.cv,perfX.quic.cv=perfX.quic)
  
  #list(perf.bic=perf.bic, perf.cv=perf.cv)
  return(perf.bic)
}


Alg_cmu = function(y,x,lam_cov,lam_F,EnoughIter)
{
  eps <- 1e-4; prevfx <- 1e15; n <- dim(y)[1]; p <- dim(y)[2]; px <- dim(x)[2];
  Cy <- 1/n*t(y)%*%y; Cxy<-1/n*t(x)%*%y; Cx <- 1/n*t(x)%*%x;
  # initialize 
  Theta_init <- matrix(0, nrow = px, ncol = p); Lambda_init = diag(p)
  iter <- 0; Sigma_old <- diag(p); Theta_old <- Theta_init; Lambda_old <- Lambda_init
  Q_old <- Q_init <- 1/sqrt(n) * x %*% Theta_init
  R_old <- R_init <- Q_init %*% Sigma_old
  Psi_old <- Psi_init <- t(R_init) %*% R_init
  trRtQ_old <- sum(diag(t(R_init)%*%Q_init))
  logdetLambda_old <- log(det(Lambda_old))
  trsX <- sum(diag(Cy%*%Lambda_old)) + trRtQ_old
  l1normLambda <- lam_cov * (sum(Lambda_old != 0)-p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -logdetLambda_old + trsX + l1normLambda + l1normTheta
  while( (iter<EnoughIter) && (abs((fx - prevfx)/fx) > eps))
  {
    iter=iter+1;
    #cat("iteration index: ", iter, "\n")
    # CD for Lambda, also updates Lambda, R 
    Lambda_update <- cdL_cmu(x, Cx, Cxy, Cy, Theta_old, Lambda_old, logdetLambda_old,
                             Sigma_old, lam=lam_cov,
                             Q_old, R_old, Psi_old, trRtQ_old)
    Lambda_new <- Lambda_update$Lambda
    Sigma_new <- solve(Lambda_new)
    R_new <- Lambda_update$R
    trRtQ_new <- Lambda_update$trRtQ
    logdetLambda_new <- Lambda_update$logdetLambda
    # CD for Theta, also updates Q and R
    Theta_update <- cdF_cmu(n,Sigma_new, x, Cx, Cxy, Theta_old, lam = lam_F, Q_old, R_new)
    Theta_new <- Theta_update$Theta
    Q_new <- Theta_update$Q
    R_new2 <- Theta_update$R
    Psi_new <- t(R_new2) %*% R_new2
    # objective function
    trsX <- sum(diag(Cy%*%Lambda_new)) + trRtQ_new
    l1normLambda <- lam_cov * (sum(Lambda_new != 0)-p)
    l1normTheta <- lam_F * sum(Theta_new != 0)
    fx1 <- -logdetLambda_new + trsX + l1normLambda + l1normTheta
    prevfx <- fx
    fx <- fx1
    Theta_old <- Theta_new;  Lambda_old <- Lambda_new; Sigma_old <- Sigma_new
    Q_old <- Q_new; R_old <- R_new2; Psi_old <- Psi_new
    trRtQ_old <- trRtQ_new; logdetLambda_old <- logdetLambda_new
  }
  list(Theta=t(Theta_new), Lambda=Lambda_new)
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

scadrightderv <- function(lamhat, a, lam)
{
  pmax(lam*((lamhat<=lam) + pmax(a*lam-lamhat,0) * (lamhat>lam)/(a-1) ), 1e-10)
}

selectBIC_glasso <- function(trainy,trainx,Lam_cov,N,opt, method="quic") {
  library(glasso)
  smallnum=1e-4; ridge=1e-2; gamma=gama=0.5
  n1=dim(trainy)[1]; 	p=dim(trainy)[2];		px=dim(trainx)[2];
  z <- cbind(trainy,trainx)
  S <- cov(z)
  pz <- p+px
  F <- 0
  bicmatrix <- rep(NA,N+1)
  mod1=glasso(S,rho=Lam_cov,penalize.diagonal=FALSE)
  Theta2 = mod1$wi
  if (method=="alasso") {
    rhomat = Lam_cov/(p+px)/2 * matrix(1, (p+px), (p+px))/(pmax(abs(Theta2)^gamma, 1e-05))
    fit.alasso = glasso(S, rhomat,penalize.diagonal=FALSE)
    Theta2 = fit.alasso$wi
  } else if (method=="scad") {
    wi.scad = Theta2
    epsi = 1
    count = 1
    while (epsi > 1e-04) {
      if (epsi < 0.001 & count > 20) 
        break
      count = count + 1
      wi.scad.old = wi.scad
      rhomat = scadrightderv(abs(wi.scad.old), 3.7, Lam_cov)
      fit.scad = glasso(S, rhomat,penalize.diagonal=FALSE)
      Theta2 = fit.scad$wi
      epsi = mean(abs(wi.scad - wi.scad.old))
      if (count > 50) {
        warning("scad iteration does not converge")
        break
      }
    }
  } else if (method != "quic") {
    stop("Inappropriate penalty type in function selectBIC_glasso")
  }
  Theta2s=matrix(as.integer(abs(Theta2)>smallnum),nrow=pz);
  cvob2=CallBIC(trainy,trainx,mod1$wi,0,1)$BIC
  opt_cv_egec=cvob2;
  optlam_egec_cov=Lam_cov;
  optTheta_egec=Theta2
  optadj_egec=Theta2s - diag(diag(Theta2s));
  optSig_egec <- mod1$wi
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
      obj1=glasso(S,rho=rho1,penalize.diagonal=FALSE);
    }else{
      
      obj1=glasso(S,rho=lam_cov,penalize.diagonal=FALSE);
    }
    Theta2=obj1$wi
    if (method=="alasso") {
      rhomat = Lam_cov/(p+px)/2 * matrix(1, (p+px), (p+px))/(pmax(abs(Theta2)^gamma, 1e-05))
      obj1 = glasso(S, rhomat, penalize.diagonal=FALSE)
      Theta2 = obj1$wi
    } else if (method=="scad") {
      wi.scad = Theta2
      epsi = 1
      count = 1
      while (epsi > 1e-04) {
        if (epsi < 0.001 & count > 20) 
          break
        count = count + 1
        wi.scad.old = wi.scad
        rhomat = scadrightderv(abs(wi.scad.old), 3.7, Lam_cov)
        obj1 = glasso(S, rhomat,penalize.diagonal=FALSE)
        Theta2 = obj1$wi
        epsi = mean(abs(wi.scad - wi.scad.old))
        if (count > 50) {
          warning("scad iteration does not converge")
          break
        }
      }
    }
    
    cvob2=CallBIC(trainy,trainx,obj1$wi,0,1)$BIC;
    bicmatrix[i] = cvob2
    Theta2s=matrix(as.integer(abs(Theta2)>smallnum),nrow=pz);
    cv_egec=cvob2;
    if(cv_egec < opt_cv_egec)
    {
      opt_cv_egec=cv_egec;  
      optlam_egec_cov=lam_cov;
      optTheta_egec=Theta2
      optadj_egec=Theta2s - diag(diag(Theta2s));
      optSig_egec <- obj1$w
    }
  }	
  list( optBIC=opt_cv_egec, optlam_cov=optlam_egec_cov,
        opt_Lam = optTheta_egec,
        opt_adj=optadj_egec, opt_Sigma=optSig_egec, bicmatrix=bicmatrix)
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

# maxit = 1 in QUIC
# use difference in objective function as convergence criterion 
Alg1 = function(y,x,lam_cov,lam_F,EnoughIter,thr)
{
  library(QUIC)
  eps <- 1e-4; prevfx <- 1e15; n <- dim(y)[1]; p <- dim(y)[2]; px <- dim(x)[2];
  smnum <- 1e-7; Cy <- 1/n*t(y)%*%y; Cyx<-1/n*t(y)%*%x; Cx <- 1/n*t(x)%*%x;
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
  # code start 
  # Lambda_init <- diag(p)
  # Theta_init <- matrix(0,p,px)
  # tC%*% Cyx %*% solve(Cx)  true theta_xy
  iter <- 0; Sigma_old <- S; Theta_old <- Theta_init; Lambda_old <- Lambda_init
  sXmatrix <- Cy%*%Lambda_old + Sigma_old%*%Theta_old%*%Cx%*%t(Theta_old)
  trsX <- sum(diag(sXmatrix))
  l1normLambda <- lam_cov * (sum(Lambda_old != 0)-p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -log(det(Lambda_old)) + trsX + l1normLambda + l1normTheta
  record1 <- record2 <- rep(0, EnoughIter)
  # while( (iter<EnoughIter) && ((diff_Theta > thr) ||
  #                              (diff_Lambda > thr))
  #        && (abs((fx - prevfx)/fx) > eps))
  while( (iter<EnoughIter) && (abs((fx - prevfx)/fx) > eps))
    # while( (iter<EnoughIter) && ((diff_Theta > thr) ||
    #                                (diff_Lambda > thr)))
  {
    iter=iter+1;
    #print("Updating Lambda...")
    dat <- y + x%*%t(Theta_old)%*%Sigma_old
    S_temp <- 1/n*t(dat)%*%dat
    #S_temp <- cov2cor(S_temp)
    obj_quic <- QUIC(S_temp, rho = lam_cov, 
                     path = NULL,  msg = 0,
                     X.init = Lambda_old, W.init = Sigma_old,
                     maxIter = 1)
    Sigma_new <- obj_quic$W; Lambda_new <- obj_quic$X
    obF1 <- cdF_nw(Sigma_new,Cx,Cyx,Theta_old,lam_F,1e-3);
    Theta_new <- obF1$F
    diff_Theta <- mean(abs(Theta_new - Theta_old))
    diff_Lambda <- mean(abs(Lambda_new - Lambda_old))
    record1[iter] <- mean(abs(Theta_new - Theta_old))
    record2[iter] <- mean(abs(Lambda_new - Lambda_old))
    sXmatrix_new <- Cy%*%Lambda_new + Sigma_new%*%Theta_new%*%Cx%*%t(Theta_new)
    trsX <- sum(diag(sXmatrix_new))
    l1normLambda <- lam_cov * (sum(Lambda_new != 0)-p)
    l1normTheta <- lam_F * sum(Theta_new != 0)
    fx1 <- -log(det(Lambda_new)) + trsX + l1normLambda + l1normTheta
    prevfx <- fx
    fx <- fx1
    Theta_old <- Theta_new;  Lambda_old <- Lambda_new; Sigma_old <- Sigma_new
    #cat("Stopped at iteration: ", iter, "\n")
  }
  list(Theta=Theta_new, Lambda=Lambda_new)
}

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
    wi.lasso = fit.lasso$wi
    rhomat = lam_cov/p/2 * matrix(1, p, p)/(pmax(abs(wi.lasso)^gamma, 1e-05))
    fit.alasso = glasso::glasso(S_temp, rhomat,maxit=1,penalize.diagonal=FALSE)
    Lambda_new = fit.alasso$wi; Sigma_new <- fit.alasso$w
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

Alg1_scad = function(y,x,lam_cov,lam_F,EnoughIter,thr)
{
  eps <- 1e-4; prevfx <- 1e15; n <- dim(y)[1]; p <- dim(y)[2]; px <- dim(x)[2];
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
  record1 <- record2 <- rep(0, EnoughIter)
  sXmatrix <- Cy%*%Lambda_old + Sigma_old%*%Theta_old%*%Cx%*%t(Theta_old)
  trsX <- sum(diag(sXmatrix))
  l1normLambda <- lam_cov * (sum(Lambda_old != 0)-p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -log(det(Lambda_old)) + trsX + l1normLambda + l1normTheta
  while((iter<EnoughIter) && (abs((fx - prevfx)/fx) > eps)  )
  {
    iter=iter+1;
    #print("Updating Lambda...")
    dat <- y + x%*%t(Theta_old)%*%Sigma_old
    S_temp <- 1/n*t(dat)%*%dat
    
    fit.lasso = glasso::glasso(S_temp, lam_cov,maxit = 1,penalize.diagonal=FALSE)
    wi.lasso = fit.lasso$wi
    wi.scad = wi.lasso
    epsi = 1
    count = 1
    while (epsi > 1e-04) {
      if (epsi < 0.001 & count > 20) 
        break
      count = count + 1
      wi.scad.old = wi.scad
      rhomat = scadrightderv(abs(wi.scad.old), 3.7, lam_cov)
      fit.scad = glasso::glasso(S_temp, rhomat, maxit = 1,penalize.diagonal=FALSE)
      wi.scad = fit.scad$wi
      epsi = mean(abs(wi.scad - wi.scad.old))
      if (count > 50) {
        warning("scad iteration does not converge")
        break
      }
    }
    Lambda_new = fit.scad$wi; Sigma_new <- fit.scad$w
    
    obF1 <- cdF_nw(Sigma_new,Cx,Cyx,Theta_old,lam_F,1e-3);
    Theta_new <- obF1$F
    diff_Theta <- mean(abs(Theta_new - Theta_old))
    diff_Lambda <- mean(abs(Lambda_new - Lambda_old))
    record1[iter] <- mean(abs(Theta_new - Theta_old))
    record2[iter] <- mean(abs(Lambda_new - Lambda_old))
    Theta_old <- Theta_new;  Lambda_old <- Lambda_new; Sigma_old <- Sigma_new
  }
  list(Theta=Theta_new, Lambda=Lambda_new)
}


selectBICForAlg1 <- function(trainy, trainx, Lam_cov, Lam_F, N, N_j, thr, useEBIC = TRUE, method = "quic")
{
  smallnum=1e-4; 	EnoughIter=20;
  n1=dim(trainy)[1]; 	p=dim(trainy)[2];		px=dim(trainx)[2];
  #trainx_mean <- apply(trainx,2,mean)
  #trainx <- sweep(data.matrix(trainx), 2, trainx_mean)
  if (method=="quic") {
    mod <- Alg1(trainy,trainx,smallnum,smallnum,EnoughIter,thr)
  } else if (method=="alasso") {
    mod <- Alg1_alasso(trainy,trainx,smallnum,smallnum,EnoughIter,thr)
  } else if (method=="scad") {
    mod <- Alg1_scad(trainy,trainx,smallnum,smallnum,EnoughIter,thr)
  } else {
    stop("Inappropriate penalty type")
  }
  if(useEBIC) 
  {
    cvob=CallBIC(trainy,trainx,mod$Lambda,mod$Theta,2)$EBIC
  } else {cvob=CallBIC(trainy,trainx,mod$Lambda,mod$Theta,2)$BIC}
  opt_cv_egec=cvob;
  optlam_egec_cov=smallnum;
  optlam_egec_F=smallnum;
  optLambda_egec=mod$Lambda;
  optTheta_egec=mod$Theta;	
  optSig_egec = solve(optLambda_egec)
  bicmatrix <- matrix(NA,N+1,N_j+1)
  count <- 0
  for(i in (N+1):1)
  {
    lam_cov=Lam_cov/N*(N+1-i);
    for(j in (N_j+1):1)
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

# use likelihood as loss function
selectCVForAlg1_new <- function(trainy, trainx, Py, Pyx, Lam_cov, Lam_F, N, N_j, k, method="quic") {
  EnoughIter <- 10; smallnum=0;
  n <- dim(trainx)[1]
  Theta.list <- list()
  Lambda.list <- list()
  
  #Get overall solution path
  Lam_cov <- seq(0+smallnum, Lam_cov, length = N+1);
  Lam_F <- seq(0+smallnum, Lam_F, length = N_j+1);
  thr <- 1e-3
  CV <- matrix(NA, ncol = N_j+1, nrow = N+1)
  grps <- cut(1 : n, k, labels = FALSE)[sample(n)]
  loglikFold <- matrix(NA, ncol = k, nrow = 1)
  for(i in 1 : (N+1)) {
    #cat("cv for i = ", i, "\n")
    for(j in 1 : (N_j+1)) {
      lambda1 <- Lam_cov[i]
      lambda2 <- Lam_F[j]
      #cat("lambda = ", lambda1, "lambda2 = ", lambda2, "\n")
      l2pe <- matrix(NA, ncol = k, nrow = 1)
      for(kth in 1 : k) {
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
  } else {
    stop("Inappropriate penalty type in function selectCVForAlg1_new")
  }
  optLambda_egec <- finalMod$Lambda
  optTheta_egec <- finalMod$Theta
  list(optCV = optCV, optlam_cov = optlam_cov, optlam_F = optlam_F, cvMatrix = CV, 
       opt_Lam = optLambda_egec, opt_Theta = optTheta_egec)	
}

selectCVForCmu <- function(trainy, trainx, Lam_cov, Lam_F, N, N_j, k) {
  EnoughIter <- 50; smallnum=0;
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
    #cat("cv for i = ", i, "\n")
    for(j in 1 : N_j) {
      lambda1 <- Lam_cov[i]
      lambda2 <- Lam_F[j]
      l2pe <- matrix(NA, ncol = k, nrow = 1)
      for(kth in 1 : k) {
        omit <- which(grps == kth)
        modFold <- Alg_cmu(trainy[-omit, ], trainx[-omit, , drop = FALSE], lambda1, lambda2, EnoughIter)
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
  finalMod <- Alg_cmu(trainy, trainx, optlam_cov, optlam_F, EnoughIter)
  optLambda_egec <- finalMod$Lambda
  optTheta_egec <- finalMod$Theta
  list(optCV = optCV, optlam_cov = optlam_cov, optlam_F = optlam_F, cvMatrix = CV, 
       opt_Lam = optLambda_egec, opt_Theta = optTheta_egec)	
}

LOOKL <- function(y,x,lam_cov,lam_F,alpha) {
  n=dim(y)[1]; p=dim(y)[2]; px=dim(x)[2]; EnoughIter <- 10; thr <- 1e-4
  Sy=1/n*t(y)%*%y; Syx=1/n*t(y)%*%x; Sx=1/n*t(x)%*%x;
  z <- cbind(y,x)
  mod1 <- Alg1(y,x,lam_cov,lam_F,EnoughIter,thr)
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
  Inv1 <- -CInv - t(BCInv) %*% Inv2 %*% BCInv
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
  opt_lookl <- LOOKL(trainy,trainx,lam_cov=Lam_cov,lam_F=Lam_F,alpha)
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
      tmp_neglik <- LOOKL(trainy, trainx, lam_cov=lam_cov,lam_F=lam_F,alpha)$neglik
      tmp_bias <- LOOKL(trainy, trainx, lam_cov=lam_cov,lam_F=lam_F,alpha)$bias
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
    if (method=="glasso") {
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

getPerf <- function(Lambda, Theta, tC, tF) {
  ROCLambda <- ROCval(tC, Lambda, 1)
  perfLambda <- Calmcc(ROCLambda)
  ROCTheta <- ROCval(tF, Theta, 2)
  perfTheta <- Calmcc(ROCTheta)
  list(perfLambda=perfLambda, perfTheta=perfTheta)
}

newFileName <- function(targetDir, alg1Flag="alg1", paraOpt="quic", crtr, alpha=2) {
  if (alg1Flag == "alg1") {
    if (crtr != "LOOKL") {
      LambdaFileName <- paste0(targetDir, "/Alg1",paraOpt, "ResultLambda_", crtr, "_", unclass(Sys.time()))
      ThetaFileName <- paste0(targetDir, "/Alg1",paraOpt, "ResultTheta_", crtr, "_", unclass(Sys.time()))
      LambdaNormFileName <- paste0(targetDir, "/Alg1",paraOpt, "LambdaNorm_", crtr, "_", unclass(Sys.time()))
      ThetaNormFileName <- paste0(targetDir, "/Alg1",paraOpt, "ThetaNorm_", crtr, "_", unclass(Sys.time()))
      KlFileName <- paste0(targetDir, "/Alg1",paraOpt, "_KL_", crtr, "_", unclass(Sys.time()))
    } else { 
      LambdaFileName <- ThetaFileName <- KlFileName <- rep(NA, length(alpha))
      LambdaNormFileName <- ThetaNormFileName <- rep(NA, length(alpha))
      for (i in 1:length(alpha)) {
        LambdaFileName[i] <- paste0(targetDir, "/Alg1",paraOpt, "ResultLambda_", crtr, "_", round(alpha[i],2),
                                    "_", unclass(Sys.time()))
        LambdaNormFileName[i] <- paste0(targetDir, "/Alg1",paraOpt, "LambdaNorm_", crtr, "_", round(alpha[i],2),
                                        "_", unclass(Sys.time()))
        ThetaNormFileName[i] <- paste0(targetDir, "/Alg1",paraOpt, "ThetaNorm_", crtr, "_", round(alpha[i],2),
                                       "_", unclass(Sys.time()))
        
        ThetaFileName[i] <- paste0(targetDir, "/Alg1",paraOpt, "ResultTheta_", crtr, "_", round(alpha[i],2),
                                   "_", unclass(Sys.time()))
        KlFileName[i] <- paste0(targetDir, "/Alg1",paraOpt, "_KL_", crtr, "_", round(alpha[i],2),
                                "_", unclass(Sys.time()))
      }
    }
  } else if (alg1Flag == "glasso") {
    LambdaFileName <- paste0(targetDir, "/GlassoResultLambda_", crtr, "_", unclass(Sys.time()))
    ThetaFileName <- paste0(targetDir, "/GlassoResultTheta_", crtr, "_", unclass(Sys.time()))
    LambdaNormFileName <- paste0(targetDir, "/GlassoLambdaNorm_", crtr, "_", unclass(Sys.time()))
    ThetaNormFileName <- paste0(targetDir, "/GlassoThetaNorm_", crtr, "_", unclass(Sys.time()))
    KlFileName <- paste0(targetDir, "/Glasso_KL_", crtr, "_", unclass(Sys.time()))
  } else if (alg1Flag == "quic") {
    LambdaFileName <- paste0(targetDir, "/QuicResultLambda_", crtr, "_", unclass(Sys.time()))
    ThetaFileName <- paste0(targetDir, "/QuicResultTheta_", crtr, "_", unclass(Sys.time()))
    LambdaNormFileName <- paste0(targetDir, "/QuicLambdaNorm_", crtr, "_", unclass(Sys.time()))
    ThetaNormFileName <- paste0(targetDir, "/QuicThetaNorm_", crtr, "_", unclass(Sys.time()))
    KlFileName <- paste0(targetDir, "/Quic_KL_", crtr, "_", unclass(Sys.time()))
  } else if (alg1Flag == "cmu") {
    LambdaFileName <- paste0(targetDir, "/CmuResultLambda_", crtr, "_", unclass(Sys.time()))
    ThetaFileName <- paste0(targetDir, "/CmuResultTheta_", crtr, "_", unclass(Sys.time()))
    LambdaNormFileName <- paste0(targetDir, "/CmuLambdaNorm_", crtr, "_", unclass(Sys.time()))
    ThetaNormFileName <- paste0(targetDir, "/CmuThetaNorm_", crtr, "_", unclass(Sys.time()))
    KlFileName <- paste0(targetDir, "/Cmu_KL_", crtr, "_", unclass(Sys.time()))
  } else break
  list(LambdaFileName=LambdaFileName, ThetaFileName=ThetaFileName, 
       LambdaNormFileName=LambdaNormFileName, ThetaNormFileName=ThetaNormFileName,
       KlFileName=KlFileName)
}

newFileNameX <- function(targetDir, method="ss") {
  if (method == "ss") {
    KlFileName <- paste0(targetDir, "/Ss_x_KL_", unclass(Sys.time()))
  } else if (method == "glasso") {
    KlFileName <- paste0(targetDir, "/Glasso_x_KL_", unclass(Sys.time()))
  } else if (method == "quic") {
    KlFileName <- paste0(targetDir, "/Quic_x_KL_", unclass(Sys.time()))
  } else if (method == "kde") {
    KlFileName <- paste0(targetDir, "/Kde_x_KL_", unclass(Sys.time()))
  } else break
  return(KlFileName=KlFileName)
}

calKL <- function(y, x, Lambda, Py, Theta, Pyx) {
  n=dim(y)[1]; p=dim(y)[2]; px=dim(x)[2]; EnoughIter <- 10; thr <- 1e-4
  Py[abs(Py) < 1e-5] <- 0;  Pyx[abs(Pyx) < 1e-5] <- 0; 
  Sy=1/n*t(y)%*%y; Syx=1/n*t(y)%*%x; Sx=1/n*t(x)%*%x;
  prod <- Sy%*%Lambda + 2*Syx%*%t(Theta) + solve(Lambda)%*%Theta%*%Sx%*%t(Theta);
  neglik <- -log(det(Lambda)) + sum(diag(prod))
  bias1 <- log(det(Py)) - sum(diag(Lambda%*%solve(Py)-prod)) 
  mat1 <- t(Pyx)%*%solve(Py)-t(Theta)%*%solve(Lambda)
  mat2 <- mat1 %*% Lambda %*% t(mat1)
  calbias2 <- function(x) {
    res <- t(x) %*% mat2 %*% x
    return(res)
  }
  bias2 <- 1/n * sum(apply(x, MARGIN = 1, calbias2))
  bias <- bias1 + bias2
  trueKL <- neglik + bias 
  list(trueKL=trueKL, neglik=neglik, bias=bias)
}

# reverse location of n0 and n1 
calKL_2 <- function(n,mu1, mu2, sig1, sig2) {
  tmp = 0
  for (i in 1:n) {
    # reverse dist 1 dist 2 location
    tmp = tmp + normdiff(mu1=mu2[,i],sigma1=sig2,mu2=mu1[,i],sigma2=sig1,method="KL")[1]
  }
  return(tmp)
}

calKL_3 <- function(x, mu1, mu2, sig1, sig2) {
  tmp = 0
  for (i in (1:dim(x)[1])) {
    log_den_para1 <- dmvnorm(x[i,], mean = mu1[,i], sigma = sig1, log = TRUE)
    log_den_para2 <- dmvnorm(x[i,], mean = mu2[,i], sigma = sig2, log = TRUE)
    tmp = tmp + log_den_para1 - log_den_para2
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

calKL_kde_normal <- function(x, mu1, sig1) {
  library(ks)
  kl <- kde(x = x, eval.points = x)$estimate
  record <- rep(NA, length(kl))
  tmp <- 0
  for (i in 1:length(kl)) {
    den_para <- dmvnorm(x[i,], mean = mu1, sigma = sig1, log = FALSE)
    tmp = tmp +  log(den_para) - log(kl[i])
    record[i] = kl[i]
  }
  return(tmp)
}

# allow input x be vector
d.ssden1.v2 <- ## Evaluate density estimate
  function (object,x) {
    if (!("ssden1"%in%class(object))) stop("gss error in d.ssden1: not a ssden1 object")
    if (is.vector(x)) {
      x <- data.frame(t(x))
      colnames(x) <- colnames(object$mf)
    }
    ## rho
    rho <- 1
    for (xlab in names(object$mf)) {
      xx <- x[[xlab]]
      rho.wk <- object$rho[[xlab]]
      if (is.factor(xx)) rho <- rho*rho.wk[xx]
      if (is.vector(xx)&!is.factor(xx)) rho <- rho*dssden(rho.wk,xx)
      if (is.matrix(xx)) rho <- rho*dssden(rho.wk,xx)
    }
    ## exp(eta)
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
      xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
      x.new <- x[,object$terms[[label]]$vlist]
      nphi <- object$terms[[label]]$nphi
      nrk <-  object$terms[[label]]$nrk
      if (nphi) {
        phi <-  object$terms[[label]]$phi
        for (i in 1:nphi) {
          s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
        }
      }
      if (nrk) {
        rk <- object$terms[[label]]$rk
        for (i in 1:nrk) {
          nq <- nq + 1
          r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
        }
      }
    }
    as.vector(rho*exp(cbind(s,r)%*%c(object$d,object$c))*object$scal)
  }

calKL_ss_normal <- function(xs, mu1, sig1, fit.x) {
  tmp = 0
  n=dim(xs)[1]
  for (i in 1:n) {
    # cat("i= ", i, "\n")
    xpt <- xs[i,]
    den_ss <- d.ssden1.v2(fit.x, xpt)
    den_para <- dmvnorm(xs[i,], mean = mu1, sigma = sig1, log = FALSE)
    tmp = tmp + (log(den_para) - log(den_ss))
  }
  return(tmp)
}



I.3d.kde <- function(t,x,Sigmax) {
  library(ks)
  den_kde <- kde(x = x, eval.points = t)$estimate
  den_para <- dmvnorm(t, mean = rep(0,px), sigma = Sigmax, log = FALSE)
  tmp <- den_para *(log(den_para) - log(den_kde))
  return(tmp)
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

runLookl <- function(y, x, xs, n, p, px, tF, tC, tCx, Py, Pyx, maxLambda, maxTheta, 
                     iterLambda, iterTheta, method, alpha) {
  library(gaussDiff)
  mu1 <- -solve(Py) %*% Pyx %*% t(xs)
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
    mu2 <- -solve(Lambda_alg1_lookl) %*% Theta_alg1_lookl %*% t(xs)
    sig2 <- solve(Lambda_alg1_lookl)
    tmp = 0
    for (j in 1:n) {
      tmp = tmp + normdiff(mu1=mu1[,j],sigma1=sig1,mu2=mu2[,j],sigma2=sig2,method="KL")[1]
    }
    kl_alg_lookl <- tmp/n
    perf_alpha <- list(alg1_lookl_perf=alg1_lookl_perf, optlam_cov_alg1_lookl=optlam_cov_alg1_lookl,
                       optlam_F_alg1_lookl=optlam_F_alg1_lookl,kl_alg_lookl=kl_alg_lookl)
    perf_alpha_list[[i]] <- perf_alpha
  }
  names(perf_alpha_list) <- round(alpha, digits = 2)
  return(perf_alpha_list)
}

runMle <- function(y, x, n, p, px, Py, Pyx) 
{
  library(gaussDiff); library(mvtnorm)
  mu1 <- -solve(Py) %*% Pyx %*% t(x)
  sig1 <- solve(Py)
  
  
  Cy <- 1/n*t(y)%*%y; Cyx<-1/n*t(y)%*%x; Cx <- 1/n*t(x)%*%x;
  Ipx <- diag(px);
  Sigma_mle <- Cy - Cyx%*%solve(Cx)%*%t(Cyx);
  Lambda_mle <- solve(Sigma_mle)
  # from MLE
  Theta_mle <- -Lambda_mle %*% Cyx %*% solve(Cx)  #theta_xy
  mu2 <- - Sigma_mle%*% Theta_mle %*% t(x)
  sig2 <- Sigma_mle
  kl_mle <- calKL_2(n,mu1, mu2, sig1, sig2)/n
  
  covnorm_mle <- norm(Lambda_mle-Py, "F")
  Fnorm_mle <- norm(Theta_mle-Pyx, "F")
  covnorm2_mle <- base::norm(Lambda_mle-Py, "2")
  Fnorm2_mle <- base::norm(Theta_mle-Pyx, "2")
  
  list(kl_mle=kl_mle,covnorm_mle=covnorm_mle,Fnorm_mle=Fnorm_mle,covnorm2_mle=covnorm2_mle,Fnorm2_mle=Fnorm2_mle)
}

runSdke <- function(y, x, n, p, px, Py, Pyx) 
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
  
  list(kl_skde=kl_skde,covnorm_skde=covnorm_skde,Fnorm_skde=Fnorm_skde,covnorm2_skde=covnorm2_skde,Fnorm2_skde=Fnorm2_skde)
}


## Calculate integrals of phi and rk for ssden1
mkint2 <- function(mf,type,id.basis,quad,term)
{
  ## Obtain model terms
  mt <- attr(mf,"terms")
  xvars <- as.character(attr(mt,"variables"))[-1]
  xfacs <- attr(mt,"factors")
  term.labels <- labels(mt)
  vlist <- xvars[as.logical(apply(xfacs,1,sum))]
  ## Create phi and rk
  nbasis <- length(id.basis)
  phi.term <- rk.term <- list(NULL)
  nvar <- length(names(mf))
  ns <- nq <- 0
  for (label in term.labels) {
    ns <- ns+term[[label]]$nphi
    nq <- nq+term[[label]]$nrk
    phi.term[[label]] <- rk.term[[label]] <- list(NULL)
    vlist <- xvars[as.logical(xfacs[,label])]
    x <- mf[,vlist]
    dm <- length(vlist)
    phi <- rk <- NULL
    if (dm==1) {
      type.wk <- type[[vlist]][[1]]
      xx <- mf[id.basis,vlist]
      xmesh <- quad[[vlist]]$pt
      if (type.wk%in%c("nominal","ordinal")) {
        ## factor variable
        if (type.wk=="nominal") fun <- gss::mkrk.nominal(levels(x))
        else fun <- gss::mkrk.ordinal(levels(x))
        if (nlevels(x)>2) {
          ## rk
          rk <- fun$fun(xmesh,xx,fun$env,TRUE)
        }
        else {
          ## phi
          wk <- as.factor(names(fun$env$code)[1])
          phi <- fun$fun(xmesh,wk,fun$env)
        }
      }
      if (type.wk=="cubic") {
        ## cubic splines
        range <- type[[vlist]][[2]]
        ## phi
        phi.fun <- gss::mkphi.cubic(range)
        phi <- phi.fun$fun(xmesh,1,phi.fun$env)
        ## rk
        rk.fun <-  gss::mkrk.cubic(range)
        rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
      }
      if (type.wk%in%c("cubic.per","linear","linear.per","sphere")) {
        ## cubic periodic, linear, and linear periodic splines
        range <- type[[vlist]][[2]]
        ## rk
        if (type.wk=="cubic.per") rk.fun <-  gss::mkrk.cubic.per(range)
        if (type.wk=="linear") rk.fun <-  gss::mkrk.linear(range)
        if (type.wk=="linear.per") rk.fun <-  gss::mkrk.linear.per(range)
        if (type.wk=="sphere") rk.fun <-  gss::mkrk.sphere(range)
        rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
      }
      if (type.wk=="tp") {
        ## thin-plate splines
        par <- type[[vlist]][[2]]
        order <- par$order
        mesh <- par$mesh
        weight <- par$weight
        if (is.vector(x)) xdim <- 1
        else xdim <- dim(x)[2]
        ## phi
        phi.fun <-  gss::mkphi.tp(xdim,order,mesh,weight)
        nphi <- choose(xdim+order-1,xdim)-1
        if (nphi>0) {
          for (nu in 1:nphi) {
            phi <- cbind(phi,phi.fun$fun(xmesh,nu,phi.fun$env))
          }
        }
        ## rk
        rk.fun <-  gss::mkrk.tp(xdim,order,mesh,weight)
        rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
      }
      if (type.wk=="custom") {
        ## user-defined
        par <- type[[vlist]][[2]]
        nphi <- par$nphi
        if (nphi>0) {
          phi.fun <- par$mkphi(par$env)
          for (nu in 1:nphi) {
            phi <- cbind(phi,phi.fun$fun(xmesh,nu,phi.fun$env))
          }
        }
        rk.fun <- par$mkrk(par$env)
        rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
      }
      phi.term[[label]][[vlist]] <- phi
      if (is.null(rk)) rk.term[[label]][[vlist]] <- rk
      else {
        nmesh <- length(quad[[vlist]]$wt)
        rk.term[[label]][[vlist]] <- array(rk,c(nmesh,nbasis,1))
      }
    }
    else {
      bin.fac <- n.phi <- phi.list <- rk.list <- NULL
      for (i in 1:dm) {
        type.wk <- type[[vlist[i]]][[1]]
        if (type.wk%in%c("nominal","ordinal")) {
          ## factor variable
          if (type.wk=="nominal")
            rk.wk <-  gss::mkrk.nominal(levels(x[[i]]))
          else rk.wk <-  gss::mkrk.ordinal(levels(x[[i]]))
          phi.wk <- rk.wk
          n.phi <- c(n.phi,0)
          bin.fac <- c(bin.fac,!(nlevels(x[[i]])>2))
        }
        if (type.wk=="cubic") {
          ## cubic or linear splines
          range <- type[[vlist[i]]][[2]]
          ## phi
          phi.wk <- gss::mkphi.cubic(range)
          n.phi <- c(n.phi,1)
          ## rk
          rk.wk <-  gss::mkrk.cubic(range)
          bin.fac <- c(bin.fac,0)
        }
        if (type.wk%in%c("cubic.per","linear","linear.per","sphere")) {
          ## cubic periodic, linear, or linear periodic splines
          range <- type[[vlist[i]]][[2]]
          n.phi <- c(n.phi,0)
          phi.wk <- NULL
          if (type.wk=="cubic.per") rk.wk <- gss::mkrk.cubic.per(range)
          if (type.wk=="linear") rk.wk <- gss::mkrk.linear(range)
          if (type.wk=="linear.per") rk.wk <- gss::mkrk.linear.per(range)
          if (type.wk=="sphere") rk.wk <- gss::mkrk.sphere(range)
          bin.fac <- c(bin.fac,0)
        }
        if (type.wk=="tp") {
          ## thin-plate splines
          par <- type[[vlist[i]]][[2]]
          order <- par$order
          mesh <- par$mesh
          weight <- par$weight
          if (is.vector(x[[i]])) xdim <- 1
          else xdim <- dim(x[[i]])[2]
          phi.wk <- gss::mkphi.tp(xdim,order,mesh,weight)
          n.phi <- c(n.phi,choose(xdim+order-1,xdim)-1)
          rk.wk <- gss::mkrk.tp(xdim,order,mesh,weight)
          bin.fac <- c(bin.fac,0)
        }
        if (type.wk=="custom") {
          ## user-defined
          par <- type[[vlist[i]]][[2]]
          n.phi <- c(n.phi,par$nphi)
          if (par$nphi>0) phi.wk <- par$mkphi(par$env)
          else phi.wk <- NULL
          rk.wk <- par$mkrk(par$env)
          bin.fac <- c(bin.fac,0)
        }
        phi.list <- c(phi.list,list(phi.wk))
        rk.list <- c(rk.list,list(rk.wk))
      }
      ## phi
      id0 <- names(mf)%in%vlist
      nphi <- term[[label]]$nphi
      iphi <- term[[label]]$iphi
      if (nphi>0) {
        for (nu in 1:nphi) {
          ind <- nu - 1
          for (i in 1:dm) {
            phi.wk <- phi.list[[i]]
            xmesh <- quad[[vlist[i]]]$pt
            if (bin.fac[i]) {
              wk <- as.factor(names(phi.wk$env$code)[1])
              phi <- phi.wk$fun(xmesh,wk,phi.wk$env)
            }
            else {
              code <- ind%%n.phi[i] + 1
              ind <- ind%/%n.phi[i]
              phi <- phi.wk$fun(xmesh,code,phi.wk$env)
            }
            phi.term[[label]][[vlist[i]]] <-
              cbind(phi.term[[label]][[vlist[i]]],phi)
          }
        }
      }
      ## rk
      n.rk <- ifelse(n.phi,2,1)
      nrk <- prod(n.rk) - as.logical(nphi)
      if (nrk>0) {
        for (nu in 1:nrk) {
          ind <- nu - !nphi
          for (i in 1:dm) {
            code <- ind%%n.rk[i] + 1
            ind <- ind%/%n.rk[i]
            xx <- mf[id.basis,vlist[[i]]]
            xmesh <- quad[[vlist[i]]]$pt
            if (code==n.rk[i]) {
              rk.wk <- rk.list[[i]]
              rk <- rk.wk$fun(xmesh,xx,rk.wk$env,TRUE)
            }
            else {
              rk <- 0
              phi.wk <- phi.list[[i]]
              for (j in 1:n.phi[i]) {
                phix <- phi.wk$fun(xmesh,j,phi.wk$env)
                phiy <- phi.wk$fun(xx,j,phi.wk$env)
                rk <- rk + outer(phix,phiy)
              }
            }
            nmesh <- length(quad[[vlist[i]]]$wt)
            rk.term[[label]][[vlist[i]]] <-
              array(c(rk.term[[label]][[vlist[i]]],rk),
                    c(nmesh,nbasis,nu))
          }
        }
      }
    }
  }
  ## create arrays
  ss <- matrix(1,ns,ns)
  sr <- array(1,c(ns,nbasis,nq))
  rr <- array(1,c(nbasis,nbasis,nq,nq))
  for (label1 in term.labels) {
    if (!term[[label1]]$nphi) id.s1 <- NULL
    else id.s1 <- term[[label1]]$iphi+(1:term[[label1]]$nphi)-2
    if (!term[[label1]]$nrk) id.r1 <- NULL
    else id.r1 <- term[[label1]]$irk+(1:term[[label1]]$nrk)-1
    irk1 <- term[[label1]]$irk
    for (label2 in term.labels) {
      if (!term[[label2]]$nphi) id.s2 <- NULL
      else id.s2 <- term[[label2]]$iphi+(1:term[[label2]]$nphi)-2
      if (!term[[label2]]$nrk) id.r2 <- NULL
      else id.r2 <- term[[label2]]$irk+(1:term[[label2]]$nrk)-1
      irk2 <- term[[label2]]$irk
      for (xlab in names(mf)) {
        wmesh <- quad[[xlab]]$wt
        phi1 <- phi.term[[label1]][[xlab]]
        phi2 <- phi.term[[label2]][[xlab]]
        rk1 <- rk.term[[label1]][[xlab]]
        rk2 <- rk.term[[label2]][[xlab]]
        ## ss
        if (!is.null(id.s1)&!is.null(id.s2)) {
          if ((!is.null(phi1))&(!is.null(phi2))) {
            ss[id.s1,id.s2] <- ss[id.s1,id.s2]*(t(wmesh*phi1)%*%phi2)
          }
          else {
            if (!is.null(phi1)) {
              ss[id.s1,id.s2] <- ss[id.s1,id.s2]*apply(wmesh*matrix(phi1),2,sum)
            }
            else {
              if (!is.null(phi2)) {
                ss[id.s1,id.s2] <- t(t(ss[id.s1,id.s2])*
                                       apply(wmesh*matrix(phi2),2,sum))
              }
            }
          }
        }
        ## sr
        if (!is.null(id.s1)&!is.null(id.r2)) {
          if ((!is.null(phi1))&(!is.null(rk2))) {
            for (i in id.r2) {
              sr[id.s1,,i] <- sr[id.s1,,i]*(t(wmesh*phi1)%*%rk2[,,i-irk2+1])
            }
          }
          else {
            if (!is.null(phi1)) {
              sr[id.s1,,id.r2] <- sr[id.s1,,id.r2]*apply(wmesh*matrix(phi1),2,sum)
            }
            else {
              if (!is.null(rk2)) {
                for (i in id.r2) {
                  sr[id.s1,,i] <- t(t(sr[id.s1,,i])*
                                      apply(wmesh*rk2[,,i-irk2+1],2,sum))
                }
              }
            }
          }
        }
        ## rr
        if (!is.null(id.r1)&!is.null(id.r2)) {
          if ((!is.null(rk1))&(!is.null(rk2))) {
            for (i in id.r1) {
              for (j in id.r2) {
                rr[,,i,j] <- rr[,,i,j]*(t(wmesh*rk1[,,i-irk1+1])%*%rk2[,,j-irk2+1])
              }
            }
          }
          else {
            if (!is.null(rk1)) {
              for (i in id.r1) {
                rr[,,i,id.r2] <- rr[,,i,id.r2]*apply(wmesh*rk1[,,i-irk1+1],2,sum)
              }
            }
            else {
              if (!is.null(rk2)) {
                for (j in id.r2) {
                  rr[,,id.r1,j] <-
                    aperm(aperm(rr[,,id.r1,j,drop=FALSE],c(2,1,3,4))*
                            apply(wmesh*rk2[,,j-irk2+1],2,sum),c(2,1,3,4))
                }
              }
            }
          }
        }
      }
    }
  }
  list(ss=ss,sr=sr,rr=rr)
}


d.ssden <- ## Evaluate density estimate
  function (object,x) {
    if (class(object)!="ssden") stop("gss error in d.ssden: not a ssden object")
    if (dim(object$mf)[2]==1&is.vector(x)) {
      x <- data.frame(x)
      colnames(x) <- colnames(object$mf)
    }
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
      xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
      x.new <- x[,object$terms[[label]]$vlist]
      nphi <- object$terms[[label]]$nphi
      nrk <-  object$terms[[label]]$nrk
      if (nphi) {
        phi <-  object$terms[[label]]$phi
        for (i in 1:nphi) {
          s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
        }
      }
      if (nrk) {
        rk <- object$terms[[label]]$rk
        for (i in 1:nrk) {
          nq <- nq + 1
          r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
        }
      }
    }
    as.vector(exp(cbind(s,r)%*%c(object$d,object$c))/object$int)
  }

my.d.ssden1 <- ## Evaluate density estimate
  function (object,x) {
    if (!("ssden1"%in%class(object))) stop("gss error in d.ssden1: not a ssden1 object")
    ## rho
    rho <- 1
    for (xlab in names(object$mf)) {
      xx <- x[[xlab]]
      rho.wk <- object$rho[[xlab]]
      if (is.factor(xx)) rho <- rho*rho.wk[xx]
      if (is.vector(xx)&!is.factor(xx)) rho <- rho*dssden(rho.wk,xx)
      if (is.matrix(xx)) rho <- rho*dssden(rho.wk,xx)
    }
    ## exp(eta)
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
      xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
      x.new <- x[,object$terms[[label]]$vlist]
      nphi <- object$terms[[label]]$nphi
      nrk <-  object$terms[[label]]$nrk
      if (nphi) {
        phi <-  object$terms[[label]]$phi
        for (i in 1:nphi) {
          s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
        }
      }
      if (nrk) {
        rk <- object$terms[[label]]$rk
        for (i in 1:nrk) {
          nq <- nq + 1
          r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
        }
      }
    }
    # unscaled version
    pdf.unscaled <- as.vector(cbind(s,r)%*%c(object$d,object$c))
    list(s=s, r=r, pdf.unscaled=pdf.unscaled)
  }

my.dssden <- ## Evaluate density estimate
  function (object,x) {
    ## check input
    if (!("ssden"%in%class(object))) stop("gss error in dssden: not a ssden object")
    if ("ssden1"%in%class(object)) return(my.d.ssden1(object,x))
    else return(d.ssden(object,x))
  }


# generate data samples from rho(x)
gen.data.sp <- function(fit, N=3e3) {
  fit.dm <- fit$domain
  sample.list <- NULL
  for (label in names(fit$mf)) {
    fit.pho <- fit$rho[[label]]
    xx <- seq(fit.dm[1,label], fit.dm[2,label], by=0.01)
    tmp <- my.dssden(fit.pho, xx)
    xx.wt <- tmp / sum(tmp)
    xx.sp <- sample(xx, prob = xx.wt, replace = TRUE, size = N)
    # xx <- sample(x[, label], N, replace = TRUE)
    sample.list[[label]] <- xx.sp
  }
  return(sample.list)
}

## semi-parametric adjustment

cal.semi.mse <- function(object, delta, N=3e3, semiflag=TRUE) { 
  cat("calculating semi-mse. \n")
  fit <- object
  rho1 <- sum(object$rho.int)
  rho2 <- rho1^2-sum(object$rho.int^2)+sum(object$rho.int2)
  s <- object$int$s
  r <- object$int$r
  s.rho <- object$int$s.rho - s*rho1
  r.rho <- object$int$r.rho - r*rho1
  int2 <- mkint2(object$mf,object$int$var.type,
                 object$id.basis,object$quad,object$terms)
  ss <- int2$ss
  sr <- int2$sr
  rr <- int2$rr
  d <- object$d
  c <- object$c
  theta <- object$theta
  nq <- length(theta)
  s.eta <- ss%*%d
  r.eta <- tmp <- NULL
  r.wk <- r.rho.wk <- sr.wk <- rr.wk <- 0
  for (i in 1:nq) {
    tmp <- c(tmp,10^(2*theta[i])*sum(diag(rr[,,i,i])))
    s.eta <- s.eta + 10^theta[i]*sr[,,i]%*%c
    if (length(d)==1) r.eta.wk <- sr[,,i]*d
    else r.eta.wk <- t(sr[,,i])%*%d
    r.wk <- r.wk + 10^theta[i]*r[,i]
    r.rho.wk <- r.rho.wk + 10^theta[i]*r.rho[,i]
    sr.wk <- sr.wk + 10^theta[i]*sr[,,i]
    for (j in 1:nq) {
      r.eta.wk <- r.eta.wk + 10^theta[j]*rr[,,i,j]%*%c
      rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
    }
    r.eta <- cbind(r.eta,r.eta.wk)
  }
  s.eta <- s.eta - s*(sum(s*d)+sum(r.wk*c))
  r.eta <- r.eta - r*(sum(s*d)+sum(r.wk*c))
  ss <- ss - outer(s,s,"*")
  sr.wk <- sr.wk - outer(s,r.wk,"*")
  rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
  rho.eta <- sum(s.rho*d) + sum(r.rho.wk*c)
  eta2 <- sum(c*(rr.wk%*%c)) + sum(d*(ss%*%d)) + 2*sum(d*(sr.wk%*%c))
  old.mse <- eta2 + rho2-rho1^2 + 2*rho.eta
  #cat("old mse is", old.mse, "\n")
  
  
  if (semiflag){
    # generate data frame for monte-carlo simulation
    sample.list <- gen.data.sp(fit, N=N)
    # generate components needed to calculate mse and se
    delta.x <- delta.x2 <- B_21 <- B_22 <- A_31 <- A_32 <- eta.hat <- eta.hat2 <- rep(NaN, N)
    eta.mu.record <- rep(NaN, N)
    s.x.full <- r.x.full <- NULL
    s.eta.mc.full <- r.eta.mc.full <- NULL
    
    for(i in 1:N) {
      data.pt <- NULL
      eta.mu.tmp <- eta.mu <- NULL
      for (label in names(sample.list)) {
        fit.rho <- fit$rho[[label]]
        data.pt <- c(data.pt, sample.list[[label]][i])
        eta.mu.tmp <- c(eta.mu.tmp, my.dssden(fit.rho, sample.list[[label]][i]))
        eta.mu <- sum(log(eta.mu.tmp))
        eta.mu.record[[i]] <- eta.mu
      }
      delta.x[i] <- t(data.pt) %*% delta %*% data.pt
      delta.x2[i] <- delta.x[i]^2
      data.df <- as.data.frame(t(data.pt))
      colnames(data.df) <- names(sample.list)
      multi_density.info <- my.dssden(fit, data.df)
      eta.hat[i] <- multi_density.info$pdf.unscaled
      eta.hat2[i] <- eta.hat[i]^2
      # s.x.full[[i]] <- multi_density.info$s
      # r.x.full[[i]] <- multi_density.info$r
      # s.eta.mc.full[[i]] <- s.x.full[[i]] * eta.hat[[i]]
      # r.eta.mc.full[[i]] <- r.x.full[[i]] * eta.hat[[i]]
      A_31[i] <- (eta.hat[i] + delta.x[i])
      A_32[i] <- (eta.hat[i] + delta.x[i])^2
      B_21[i] <- (eta.hat[i] + eta.mu) * delta.x[i]
      B_22[i] <- (eta.hat[i] + eta.mu)
    }
    # A-B verifies the correctness of this mc approximation. it's close to old.mse
    A <- 1/N * (sum(eta.hat2) + sum(eta.mu.record^2) + 2*sum(eta.mu.record * eta.hat))
    B <- (1/N * sum(eta.hat + eta.mu.record))^2
    A-B
    A_2 <- 1/N * sum(delta.x2) - (1/N * sum(delta.x))^2
    B_2 <- 1/N * sum(B_21) - (1/N * sum(B_22)) * (1/N * sum(delta.x))
    semi.mse <- old.mse + A_2 + 2 * B_2
    zeta2.int <- 1/N * sum(A_32) - (1/N * sum(A_31))^2
    list(old.mse=old.mse,semi.mse=semi.mse,zeta2.int=zeta2.int,sample.list=sample.list,
         eta.hat=eta.hat, delta.x=delta.x)
  }
}

my.project.mc.v2 <- function(object, delta, semiinfo, include, drop1=FALSE, N=3e3, semiflag=TRUE) {
  old.mse <- semiinfo$old.mse; semi.mse <- semiinfo$semi.mse; zeta2.int <- semiinfo$zeta2.int; 
  sample.list <- semiinfo$sample.list; eta.hat <- semiinfo$eta.hat; delta.x <- semiinfo$delta.x; 
  s.delta <- r.delta <- NULL
  fit <- object
  rho1 <- sum(object$rho.int)
  rho2 <- rho1^2-sum(object$rho.int^2)+sum(object$rho.int2)
  s <- object$int$s
  r <- object$int$r
  s.rho <- object$int$s.rho - s*rho1
  r.rho <- object$int$r.rho - r*rho1
  int2 <- mkint2(object$mf,object$int$var.type,
                 object$id.basis,object$quad,object$terms)
  ss <- int2$ss
  sr <- int2$sr
  rr <- int2$rr
  d <- object$d
  c <- object$c
  theta <- object$theta
  nq <- length(theta)
  s.eta <- ss%*%d
  r.eta <- tmp <- NULL
  r.wk <- r.rho.wk <- sr.wk <- rr.wk <- 0
  for (i in 1:nq) {
    tmp <- c(tmp,10^(2*theta[i])*sum(diag(rr[,,i,i])))
    s.eta <- s.eta + 10^theta[i]*sr[,,i]%*%c
    if (length(d)==1) r.eta.wk <- sr[,,i]*d
    else r.eta.wk <- t(sr[,,i])%*%d
    r.wk <- r.wk + 10^theta[i]*r[,i]
    r.rho.wk <- r.rho.wk + 10^theta[i]*r.rho[,i]
    sr.wk <- sr.wk + 10^theta[i]*sr[,,i]
    for (j in 1:nq) {
      r.eta.wk <- r.eta.wk + 10^theta[j]*rr[,,i,j]%*%c
      rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
    }
    r.eta <- cbind(r.eta,r.eta.wk)
  }
  s.eta <- s.eta - s*(sum(s*d)+sum(r.wk*c))
  r.eta <- r.eta - r*(sum(s*d)+sum(r.wk*c))
  ss <- ss - outer(s,s,"*")
  sr.wk <- sr.wk - outer(s,r.wk,"*")
  rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
  rho.eta <- sum(s.rho*d) + sum(r.rho.wk*c)
  eta2 <- sum(c*(rr.wk%*%c)) + sum(d*(ss%*%d)) + 2*sum(d*(sr.wk%*%c))
  
  # se
  ## calculate projection
  rkl <- function(include) {
    cat("include: ", include, "\n")
    cat("calculating zeta \n")
    # verify correctness of mu using mc approximation
    inc.wk <- union(names(object$mf),include)
    id.s <- id.q <- NULL
    for (label in inc.wk) {
      if (!any(label==object$terms$labels)) next
      term <- object$terms[[label]]
      if (term$nphi>0) id.s <- c(id.s,term$iphi+(1:term$nphi)-2)
      if (term$nrk>0) id.q <- c(id.q,term$irk+(1:term$nrk)-1)
    }
    ss.wk <- ss[id.s,id.s]
    r.eta.wk <- r.wk <- sr.wk <- rr.wk <- 0
    for (i in id.q) {
      r.eta.wk <- r.eta.wk + 10^theta[i]*r.eta[,i]
      r.wk <- r.wk + 10^theta[i]*r[,i]
      sr.wk <- sr.wk + 10^theta[i]*sr[id.s,,i]
      for (j in id.q) {
        rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
      }
    }
    sr.wk <- sr.wk - outer(s[id.s],r.wk,"*")
    rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
    v <- cbind(rbind(ss.wk,t(sr.wk)),rbind(sr.wk,rr.wk))
    if (semiflag) {
      form <- as.formula(paste("~",paste(inc.wk, collapse=" + ")))
      fit.reduced <- ssden1(form, domain=object$domain, 
                            id.basis=object$id.basis, data=object$mf)
      fit.reduced$theta <- object$theta[id.q]
      s.eta.mc <- r.eta.mc <- s.x <- r.x <- NULL
      for(i in 1:N) {
        data.pt <- NULL
        for (label in names(sample.list)) {
          data.pt <- c(data.pt, sample.list[[label]][i])
        }
        data.df <- as.data.frame(t(data.pt))
        colnames(data.df) <- names(sample.list)
        reduced.info <- my.dssden(fit.reduced, data.df)
        s.x[[i]] <- reduced.info$s
        r.x[[i]] <- reduced.info$r
        s.eta.mc[[i]] <- s.x[[i]] * eta.hat[[i]]
        r.eta.mc[[i]] <- r.x[[i]] * eta.hat[[i]]
        s.delta[[i]] <- s.x[[i]] * delta.x[i]
        r.delta[[i]] <- r.x[[i]] * delta.x[i]
      }
      mu_1 <- 1/N * Reduce("+",s.eta.mc) - (1/N * sum(eta.hat)) * (1/N * Reduce("+",s.x))
      mu_2 <- 1/N * Reduce("+",r.eta.mc) - (1/N * sum(eta.hat)) * (1/N * Reduce("+",r.x))
      s.delta.int <- 1/N * Reduce("+",s.delta) - (1/N * sum(delta.x)) * (1/N * Reduce("+",s.x))
      r.delta.int <- 1/N * Reduce("+",r.delta) - (1/N * sum(delta.x)) * (1/N * Reduce("+",r.x))
      mu_mc <- c(mu_1 + s.delta.int, mu_2 + r.delta.int)
      mu <- mu_mc
      nn <- length(mu)
      # z <- chol(v,pivot=TRUE)
      # v.ori <-v
      # v <- z
      # rkv <- attr(z,"rank")
      # m.eps <- .Machine$double.eps
      # while (v[rkv,rkv]<2*sqrt(m.eps)*v[1,1]) rkv <- rkv - 1
      # if (rkv<nn) v[(1:nn)>rkv,(1:nn)>rkv] <- diag(v[1,1],nn-rkv)
      # mu <- backsolve(v,mu[attr(z,"pivot")],transpose=TRUE)
      # se <- eta2 - sum(mu[1:rkv]^2)
      dc.tilda <- solve(v) %*% mu
      d.tilda <- dc.tilda[1:length(fit.reduced$d)]
      c.tilda <- dc.tilda[(length(fit.reduced$d)+1) :length(dc.tilda)]
      zeta.tilda <- sum(c.tilda*(rr.wk%*%c.tilda)) + sum(d.tilda*(ss.wk%*%d.tilda)) + 2*sum(d.tilda*(sr.wk%*%c.tilda))
      se <- zeta2.int - zeta.tilda
    }
    # A_3 <- 1/N * sum(A_31) - (1/N * sum(eta.hat + delta.x))^2
    # A_3 <- 1/N * sum(A_31) - (1/N * sum(eta.hat + delta.x))^2
    else {
      mu <- c(s.eta[id.s], r.eta.wk)
      nn <- length(mu)
      z <- chol(v,pivot=TRUE)
      v.ori <-v
      v <- z
      rkv <- attr(z,"rank")
      m.eps <- .Machine$double.eps
      #cat("rkv value is: ", rkv, "\n")
      while (v[rkv,rkv]<2*sqrt(m.eps)*v[1,1]) rkv <- rkv - 1
      #cat("now rkv value is: ", rkv, "\n")
      if (rkv<nn) v[(1:nn)>rkv,(1:nn)>rkv] <- diag(v[1,1],nn-rkv)
      mu <- backsolve(v,mu[attr(z,"pivot")],transpose=TRUE)
      #print("print here")
      se <- eta2 - sum(mu[1:rkv]^2)
    }
  }
  
  ## projection
  # cat("returning ratios")
  if (drop1) {
    se <- NULL
    for (i in 1:length(include)) se <- c(se, rkl(include[-i]))
    if (semiflag) ratio <- se / semi.mse
    else  {
      ratio <- se / old.mse
      cat("ratio: ", ratio)
    }
    names(se) <- names(ratio) <- include
  }
  else se <- rkl(include)
  if (semiflag) ratio <- se / semi.mse
  else ratio <- ratio <- se / old.mse
  list(ratio=ratio,se=se)
}

# Non-parametric
ssX <- function(x, delta, semiflag=TRUE, cutoff=0.05, N) {
  library(gss)
  px <- dim(x)[2]
  n <- dim(x)[1]
  mn <- apply(x,2,min);
  mx <- apply(x,2,max)
  numterm <- px + px * (px - 1) / 2
  domain <- data.frame(rbind(mn, mx))
  x <- data.frame(x)
  fit.x <- ssden1(~(.)^2, domain=domain, data=x)
  #print("Full ssden1 model fitted")
  label <- fit.x$terms$labels[(px+1):numterm]
  
  
  ## use semi-parametric projection
  semi.info <- cal.semi.mse(fit.x, delta, N=N, semiflag = semiflag)
  ratio <- my.project.mc.v2(fit.x,delta,semi.info, include=label,drop1=TRUE,N=N,semiflag=semiflag)$ratio
  #ratio <- project(fit.x,include=label,drop1=TRUE)$ratio
  cat("Projection compared to only main effects model calculated \n")
  ord <- rev(order(ratio))
  po <- length(ord)
  if (semiflag) {
    temp <- my.project.mc.v2(fit.x,delta, semi.info, include=0, drop1=FALSE,N=N,semiflag=semiflag)$ratio
    if (temp < cutoff) {
      cat("ratio with only main effects: ", temp, "\n")
      return("No interaction")
    }
    for(i in 1: po) {
      cat("Doing projection at i = ", i, "\n")
      lis <- ord[1:i]
      ratio <- my.project.mc.v2(fit.x,delta, semi.info, include=label[lis],drop1=FALSE,N=N,semiflag=semiflag)$ratio
      cat("current term: ", label[lis], ", current ratio: ", ratio, "\n")
      if(ratio < cutoff) {
        return(label[lis])
      }
    }
  }
  else {
    if(project(fit.x,include=0)$ratio < cutoff) {
      return("No interaction")
    }
    for(i in 1: po) {
      print(paste("Doing projection at i = ", i, "\n"))
      lis <- ord[1:i]
      ratio <- project(fit.x,include=label[lis])$ratio
      if(ratio < cutoff) {
        return(label[lis])
      }
    }
  }
}

# translate ssden1 result into precision matrix of x
XandX <- function(Xr,CXX) {
  px=dim(CXX)[2];
  CXXr <- CXX
  l <- length(Xr)
  if(Xr[1] != "No interaction") {
    for(k in 1:l) {
      for(i in 1:px) {
        for(j in 1:px) {
          quo <- paste("X",i,":X",j,sep="")
          t <- as.character(Xr[k])
          quo <- as.character(quo)
          if(quo == t) {
            CXXr[i,j] <- 1
          }
          quo2 <- paste("X",j,":X",i,sep="")
          quo2 <- as.character(quo2)
          if(quo2 == t) {
            CXXr[j,i] <- 1
          }
        }
      }
    }
  }
  CXXr
}
