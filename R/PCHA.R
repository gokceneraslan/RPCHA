PCHA <- function(X, noc, I, U, delta, ...) {

  opts <- list(...)
  stopifnot(is.matrix(X))

  if(missing(I)) I <- seq_len(ncol(X))
  if(missing(U)) U <- seq_len(ncol(X))
  if(missing(delta)) delta <- 0

  SST <- sum(X[,U]^2)

  if(exists('conv_crit', where=opts)){
    conv_crit <- opts$conv_crit
  } else {
    conv_crit <- 1e-6
  }

  if(exists('maxiter', where=opts)){
    maxiter <- opts$maxiter
  } else {
    maxiter <- 500
  }

  if(exists('C', where=opts)){
    C <- opts$C
  } else {
   i <- FurthestSum(X[,I], noc, ceiling(length(I) * runif(1)))
   #i <- FurthestSum(X[,I], noc, ceiling(length(I) * 0.7))
   ########################## FIXMEEEEEEEEEEEEEEEEEEEEEEEE ########################
   C <- sparseMatrix(i, seq_len(noc), x=1, dims=c(length(I), noc))
  }

  XC <- X[,I] %*% C

  muS <- 1
  muC <- 1
  mualpha <- 1

  #initialize S
  if(exists('S', where=opts)){
    S <- opts$S
    CtXtXC <- t(XC)%*%XC
    XSt <- X[,U]%*%t(S)
    SSt <- S %*% t(S)
    SSE <- SST - 2*sum(XC*XSt) + sum(CtXtXC*SSt)
  } else {
    XCtX <- t(XC) %*% X[,U]
    CtXtXC <- t(XC) %*% XC
    #f <- as.matrix(read.csv('f.csv', header = F))
    #S <- -log(f)
    ########################## FIXMEEEEEEEEEEEEEEEEEEEEEEEE ########################
    S <- -log(matrix(runif(noc*length(U)), nrow=noc, ncol=length(U)))
    S <- S / (rep(1,noc) %o% colSums(S))
    SSt <- S %*% t(S)
    SSE <- SST - 2*sum(XCtX*S) + sum(CtXtXC*SSt)
    Sup <- Supdate(S, XCtX, CtXtXC, muS, SST, SSE, 25)
    S<-Sup$S ; SSE<-Sup$SSE ; muS <- Sup$muS ; SSt<-Sup$SSt
  }

  iter <- 0
  dSSE <- Inf
  varexpl <- (SST-SSE)/SST

  cat(sprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s ','Iteration','Expl. var.','Cost func.','Delta SSEf.','muC','mualpha','muS'))
  cat(sprintf('-------------+--------------+--------------+--------------+--------------+--------------+--------------+'))

  while(abs(dSSE)>=conv_crit*abs(SSE) && iter < maxiter && varexpl < 0.9999) {
    if (iter %% 100) print(paste0('Iter: ', iter))

    iter <- iter+1
    SSE_old <- SSE

    # C (and alpha) update
    XSt <- X[,U] %*% t(S)
    #[C,SSE,muC,mualpha,CtXtXC,XC]=Cupdate(X(:,I),XSt,XC,SSt,C,delta,muC,mualpha,SST,SSE,10);
    Cup <- Cupdate(X[,I], XSt, XC, SSt, C, delta, muC, mualpha, SST, SSE, 10)
    C <- Cup$C; SSE <- Cup$SSE; muC <- Cup$muC; mualpha <- Cup$mualpha; CtXtXC <- Cup$CtXtXC; XC <- Cup$XC

    # S update
    XCtX <- t(XC) %*% X[,U]
    Sup <- Supdate(S, XCtX, CtXtXC, muS, SST, SSE, 10)
    S<-Sup$S ; SSE<-Sup$SSE ; muS <- Sup$muS ; SSt<-Sup$SSt

    # Evaluate and display iteration
    dSSE <- SSE_old - SSE
    varexpl <- (SST-SSE)/SST
    cat(sprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4e\n',iter,varexpl,SSE,dSSE/abs(SSE),muC,mualpha,muS))
  }

  varexpl <- (SST-SSE) / SST
  cat(sprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4e \n',iter,varexpl,SSE,dSSE/abs(SSE),muC,mualpha,muS))

  # sort components according to importance
  ind <- order(rowSums(S), decreasing = T)
  S <- S[ind, ]
  C <- C[, ind]
  XC <- XC[, ind]
  return(list(XC=XC,S=S,C=C,SSE=SSE,varexpl=varexpl))
}

Supdate <- function(S,XCtX,CtXtXC,muS,SST,SSE,niter) {
  noc <- nrow(S)
  J <- ncol(S)
  e <- rep(1, noc)
  for (k in seq_len(niter)) {
   SSE_old <- SSE
   g <- (CtXtXC %*%S - XCtX) / (SST/J)
   g <- g - e %o% colSums(g * S)
   stp <- 0
   Sold <- S
   while(!stp) {
    S <- Sold - g*muS
    S[S<0] <- 0
    S <- S / (e %o% colSums(S))
    SSt <- S %*% t(S)
    SSE <- SST - 2*sum(XCtX*S) + sum(CtXtXC*SSt)
    if (SSE <= SSE_old*(1+1e-9)) {
      muS <- muS * 1.2
      stp <- 1
    }
    else {
      muS <- muS /2
    }
   }
  }

  return(list(S=S,SSE=SSE,muS=muS,SSt=SSt))
}


Cupdate <- function(X,XSt,XC,SSt,C,delta,muC,mualpha,SST,SSE,niter) {

  J <- nrow(C)
  noc <- ncol(C)
  if (missing(niter)) niter <- 1
  if (delta != 0) {
    alphaC <- colSums(C)
    C <- C %*% diag(1/alphaC)
  }

  e <- rep(1, J) #col vector
  XtXSt <- t(X) %*% XSt
  for (k in seq_len(niter)) {
    SSE_old <- SSE
    g <- (t(X) %*% (XC%*%SSt) - XtXSt) / SST

    if (delta!=0) {
      g <- g %*% diag(alphaC)
    }

    g <- g - e %o% colSums(g*C)
    stp <- 0
    Cold <- C
    while(!stp) {
      C <- Cold - muC*g
      C[C<0] <- 0
      nC <- colSums(C) + .Machine$double.eps
      C <- C %*% sparseMatrix(seq_len(noc), seq_len(noc), x=1/nC)
      if (delta != 0) {
        Ct <- Matrix(C %*% diag(alphaC), sparse = T)
      } else {
        Ct <- Matrix(C, sparse = T)
      }
      XC <- X %*% Ct
      CtXtXC <- t(XC) %*% XC
      SSE <- SST - 2*sum(XC * XSt) + sum(CtXtXC * SSt)
      if (SSE <= SSE_old*(1+1e-9)) {
          muC <- muC*1.2
          stp <- 1
      } else {
          muC <- muC / 2
      }
    }

    # Update alphaC
    SSE_old <- SSE
    if (delta!=0) {
        g <- (diag(CtXtXC %*% SSt) / alphaC - colSums(C*XtXSt)) / (SST*J)
        stp <- 0
        alphaCold <- alphaC
        while (!stp) {
            alphaC <- alphaCold - mualpha*g
            alphaC[alphaC < 1-delta] <- 1-delta
            alphaC[alphaC > 1+delta] <- 1+delta
            XCt <- XC %*% diag(alphaC/alphaCold)
            CtXtXC <- t(XCt) %*% XCt
            SSE <- SST - 2*sum(XCt*XSt) + sum(CtXtXC*SSt)
            if (SSE <= SSE_old*(1+1e-9)) {
                mualpha <- mualpha*1.2
                stp <- 1
                XC <- XCt
            } else {
                mualpha <- mualpha / 2
            }
        }
    }
  }

  if (delta!=0) {
    C <- C %*% diag(alphaC)
  }
  return(list(C=C,SSE=SSE,muC=muC,mualpha=mualpha,CtXtXC=CtXtXC,XC=XC))
}


test.PCHA <- function() {
  X = as.matrix(read.csv('X.csv', header=F, stringsAsFactors = F))
  noc=3 # Number of archetypes

  N=50
  I=FurthestSum(X, N, ceiling(0.2*ncol(X)))
  res <- PCHA(X, noc, I, delta=0.1, maxiter = 1000)
  res
}
