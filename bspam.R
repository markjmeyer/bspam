############# Bayesian Multivariate Matched Proportions #############
# R functions to run Bayesian Multivariate Matched Proportions		
#	from Meyer, Li, and Knutson 2021									
#																	
# By: Mark J Meyer													
#																	
#####################################################################

#### require libraries ####
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(MCMCpack))
suppressPackageStartupMessages(require(mvtnorm))
suppressPackageStartupMessages(require(numbers))
suppressPackageStartupMessages(require(msm))
suppressPackageStartupMessages(require(geepack))
suppressPackageStartupMessages(require(statmod))

#### BSpAM ####

bspam	<- function(X1, X2, B, burnin = NULL, penalty = c('ri', 't', 'la'), ptype = c('l2', 'pe', 'tn'), 
                  pvar = c('fl', 'hc', 'ig'), pvargroup = c('kspec', 'global'), prior = NULL, mu = NULL, up = 100, dots = up/10,
                  verbose = TRUE){
  
  if(ncol(X1) != ncol(X2)){
    stop('column dimensions of X1 and X2 must agree')
  }
  
  if(nrow(X1) != nrow(X2)){
    stop('row dimensions of X1 and X2 must agree')
  }
  
  if(B < 1){
    stop('B must be at least 1')
  }
  
  if(is.null(burnin)){
    burnin    <- B
  }
  
  if(length(penalty) > 1){
    penalty <- penalty[1]
  }
  
  if(penalty != 'ri' & penalty != 't' & penalty != 'la'){
    stop('penalty must be either ri (ridge), t, or la (laplace)')
  }
  
  if(length(ptype) > 1){
    ptype <- ptype[1]
  }
  
  if(ptype != 'l2' & ptype != 'pe' & ptype != 'tn'){
    stop('penalty type must be either l2 (classic), pe (parameter expanded), or tn (truncated normal)')
  }
  
  if(length(pvar) > 1){
    pvar <- pvar[1]
  }
  
  if(pvar != 'fl' & pvar != 'hc' & pvar != 'ig'){
    stop('prior variance must be either flat (fl), half-cauchy (hc), or inverse gamma (ig')
  }
  
  if(is.null(prior)){
    if(penalty == 'ri'){
      prior	<- list()
      if(ptype == 'tn'){
        sigb  <- prior$sigb <- 1/3
      } else{
        at		<- prior$at	<- 1
        bt		<- prior$bt <- 1
        if(ptype == 'pe'){
          as		<- prior$as <- 1
          bs		<- prior$bs <- 1
        }
      }
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          A   <- prior$A <- rep(0.5, ncol(X1))
        } else {
          A   <- prior$A <- 0.5
        }
      }
      if(pvar == 'ig'){
        ae  <- prior$ae <- 1
        be  <- prior$be <- 1
      }
    } else if(penalty == 't'){
      prior	<- list()
      nu		<- prior$nu <- 5
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          A   <- prior$A <- rep(0.5, ncol(X1))
        } else {
          A   <- prior$A <- 0.5
        }
      }
      if(pvar == 'ig'){
        ae  <- prior$ae <- 1
        be  <- prior$be <- 1
      }
    } else {
      prior	<- list()
      at    <- prior$at <- 1
      bt    <- prior$bt <- 1
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          A   <- prior$A <- rep(0.5, ncol(X1))
        } else {
          A   <- prior$A <- 0.5
        }
      }
      if(pvar == 'ig'){
        ae  <- prior$ae <- 1
        be  <- prior$be <- 1
      }
    }
  } else {
    if(penalty == 'ri'){
      if(ptype == 'tn'){
        sigb  <- prior$sigb
      } else{
        at		<- prior$at
        bt		<- prior$bt
        if(ptype == 'pe'){
          as		<- prior$as
          bs		<- prior$bs
        }
      }
      if(pvargroup == 'kspec'){
        A   <- prior$A
        if(length(A) != ncol(X1)){stop('kspec pvargroup requires K elements for A')}
      } else {
        A   <- prior$A
      }
      if(pvar == 'ig'){
        ae  <- prior$ae
        be  <- prior$be
      }
    } else if(penalty == 't'){
      nu		<- prior$nu
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          A   <- prior$A
          if(length(A) != ncol(X1)){stop('kspec pvargroup requires K elements for A')}
        } else {
          A   <- prior$A
        }
      }
      if(pvar == 'ig'){
        ae  <- prior$ae
        be  <- prior$be
      }
    }	else {
      at    <- prior$at
      bt    <- prior$bt
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          A   <- prior$A
          if(length(A) != ncol(X1)){stop('kspec pvargroup requires K elements for A')}
        } else {
          A   <- prior$A
        }
      }
      if(pvar == 'ig'){
        ae  <- prior$ae
        be  <- prior$be
      }
    }
  }
  
  if(pvargroup == 'kspec'){
    if(is.null(prior$ae)){
      ae <- prior$ae <- 1
    } else {
      ae <- prior$ae
    }
    if(is.null(prior$be)){
      be <- prior$be <- 1
    } else {
      be <- prior$be
    }
  }
  
  if(penalty == 't'){
    prior$type	<- 't'
  } else if(penalty == 'la'){
    prior$type	<- 'la'
  } else {
    prior$type	<- ptype
  }
  
  if(is.null(mu)){
    mu   <- 1/2
  }
  
  prior$pvar        <- pvar
  prior$pvargroup   <- pvargroup
  
  if(pvar == 'hc'){
    v <- 1
  } else{
    v <- 0
  }
  
  if(!is.logical(verbose)){
    verbose <- TRUE
  }
  
  n		  <- nrow(X1)
  K		  <- ncol(X1)
  Ym		<- as.matrix(cbind(X1, X2))
  ids		<- rep(1:n, each = 2*K)
  outs	<- rep(rep(1:K, 2), n)
  gDat  <- data.frame(y	= c(t(Ym)), s = rep(1:(2*K), n), id = ids, outcome = outs, 
                      idk = paste(ids, outs, sep = ''))
  L		  <- cbind(diag(K), -diag(K))
  
  Y		  <- gDat$y
  N		  <- length(Y)
  X		  <- model.matrix(y ~ factor(s) - 1, data = gDat)
  
  sumX	  <- c(apply(X1, 2, sum), apply(X2, 2, sum))
  Xs1Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), 
                                  nrow = 4, byrow = TRUE), 2, sum) > 0)
  Xs2Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), 
                                  nrow = 4, byrow = TRUE), 2, sum) > 0) + K
  XsCols	<- c(Xs1Cols, Xs2Cols)
  kids    <- matrix(0, nrow = 2*n, ncol = K)
  for(k in 1:K){
    kids[,k] <- sort(c(k + 2*K*(0:(n-1)), k + K + 2*K*(0:(n-1))))
  }
  
  if(penalty == 'ri'){
    xtx			    <- t(X)%*%X
    xty         <- t(X)%*%Y
    beb         <- rep(mu, 2*K)
    bebf			  <- (sumX + 1)/(n + 2)
    prior$beb   <- beb
    
    beta		    <- matrix(0, nrow = B + burnin, ncol = ncol(X))
    if(pvargroup == 'kspec'){
      lambda      <- matrix(1, nrow = B + burnin, ncol = K)
      zeta        <- matrix(1, nrow = B + burnin, ncol = K)
    } else {
      lambda      <- vector('numeric', length = B + burnin)
      zeta        <- rep(1, length = B + burnin)
    }
    sigk		    <- matrix(1, nrow = B + burnin, ncol = ncol(X))
    
    if(pvargroup == 'kspec'){
      for(k in 1:K){
        lambda[1,k] <- (sum((Y[kids[,k]]-X[kids[,k], c(k, k + K)]%*%bebf[c(k, k + K)])^2))/N
      }
    } else {
      lambda[1]	  <- (sum((Y-X%*%bebf)^2))/N
    }
    if(ptype == 'tn'){
      tau			    <- rep(sigb, length = B + burnin)
    } else {
      tau			    <- vector('numeric', length = B + burnin)
      tau[1]		  <- (2*K)/(t(bebf)%*%(bebf))
    }
    
    
    for(b in 2:(B + burnin)){
      ## sample betas ##
      sigt        <- tau[b-1]*diag(sigk[b-1,])
      if(pvargroup == 'kspec'){
        lambk       <- c(lambda[b-1,], lambda[b-1,])
        precb		    <- (1/lambk)*xtx + sigt
      } else {
        precb		    <- (1/lambda[b-1])*xtx + sigt
      }
      sigb		    <- 1/diag(precb)
      if(pvargroup == 'kspec'){
        mub			    <- sigb*((1/lambk)*xty + sigt%*%beb)
      } else {
        mub			    <- sigb*((1/lambda[b-1])*xty + sigt%*%beb)
      }
      if(ptype == 'tn'){
        beta[b,]	  <- rtnorm(ncol(X), mean = mub, sd = sqrt(sigb), lower = 0, upper = 1)
      } else{
        beta[b,]	  <- rnorm(ncol(X), mean = mub, sd = sqrt(sigb))
      }
      
      ## sample zeta (if HC only) ##
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          for(k in 1:K){
            if(sum(k == Xs1Cols) > 0){
              zeta[b,k]   <- 1
            } else {
              b_alp       <- 1/(A[k]^2) + v/lambda[b-1, k]
              zeta[b,k]   <- rinvgamma(1, v/2 + 1/2, b_alp)
            }
          }
        } else {
          b_alp     <- 1/(A^2) + v/lambda[b-1]
          zeta[b]   <- rinvgamma(1, v/2 + 1/2, b_alp)
        }
      }
      
      ## sample lambda ##
      if(pvargroup == 'kspec'){
        for(k in 1:K){
          if(sum(k == Xs1Cols) > 0){
            Yk            <- Y[kids[,k]]
            Xk            <- X[kids[,k], c(k, k + K)]
            betak         <- beta[b-1,c(k, k + K)]
            b_lam		      <- be + (1/2)*(t(Yk - Xk%*%betak)%*%(Yk - Xk%*%betak))
            lambda[b,k]   <- rinvgamma(1, n + ae, b_lam)
          } else {
            al            <- switch(pvar, fl = 0, hc = v/2, ig = ae)
            bl            <- switch(pvar, fl = 0, hc = v/zeta[b-1, k], ig = be)
            Yk            <- Y[kids[,k]]
            Xk            <- X[kids[,k], c(k, k + K)]
            betak         <- beta[b-1,c(k, k + K)]
            b_lam		      <- bl + (1/2)*(t(Yk - Xk%*%betak)%*%(Yk - Xk%*%betak))
            lambda[b,k]   <- rinvgamma(1, n + al, b_lam)
          }
        }
      } else {
        al          <- switch(pvar, fl = 0, hc = v/2, ig = ae)
        bl          <- switch(pvar, fl = 0, hc = v/zeta[b-1], ig = be)
        b_lam		    <- bl + (1/2)*(t(Y - X%*%beta[b-1,])%*%(Y - X%*%beta[b-1,]))
        lambda[b]   <- rinvgamma(1, N/2 + al, b_lam)
      }
      
      ## sample tau ##
      if(ptype != 'tn'){
        b_tau		    <- bt + (1/2)*sum(sigk[b-1,]*((beta[b-1,] - beb)^2))
        tau[b]		  <- rgamma(1, K + at, b_tau)
      }
      
      ## sample sigmas (if PE only) ##
      if(ptype == 'pe'){
        b_sigk		<- bs + (tau[b]/2)*(beta[b,] - beb)^2
        sigk[b,]	<- rgamma(length(b_sigk), (1/2) + as, b_sigk)
      }
      
      if(verbose){
        if(mod(b, dots) == 0){
          cat('.')
        }
        
        if(mod(b, up) == 0){
          cat(paste("\n",b,"samples completed\n"))
        }
      }
    }
    
    rho			    <- t(L%*%t(beta))
    var.comp	  <- list(lambda = lambda, tau = tau, sigk = sigk, zeta = zeta)
    mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior)
    data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
    
    out			    <- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
    
  } else if(penalty == 't'){
    xtx			    <- t(X)%*%X
    xty         <- t(X)%*%Y
    beb         <- rep(mu, 2*K)
    bebf			  <- (sumX + 1)/(n + 2)
    prior$beb   <- beb
    
    beta		    <- matrix(0, nrow = B + burnin, ncol = ncol(X))
    if(pvargroup == 'kspec'){
      lambda      <- matrix(1, nrow = B + burnin, ncol = K)
      zeta        <- matrix(1, nrow = B + burnin, ncol = K)
    } else {
      lambda      <- rep(1, length = B + burnin)
      zeta        <- rep(1, length = B + burnin)
    }
    tau			    <- rep(1, length = B + burnin)
    sigk		    <- matrix(1, nrow = B + burnin, ncol = ncol(X))
    
    beta[1,]	  <- bebf
    
    for (b in 2:(B + burnin)) { 
      ## sample betas ##
      sigt        <- diag(1/sigk[b-1,])
      if(pvargroup == 'kspec'){
        lambk       <- c(lambda[b-1,], lambda[b-1,])
        precb		    <- (1/lambk)*xtx + sigt
      } else {
        precb		    <- (1/lambda[b-1])*xtx + sigt
      }
      sigb		    <- 1/diag(precb)
      if(pvargroup == 'kspec'){
        mub			    <- sigb*((1/lambk)*xty + sigt%*%beb)
      } else {
        mub			    <- sigb*((1/lambda[b-1])*xty + sigt%*%beb)
      }
      beta[b,]	  <- rnorm(ncol(X), mean = mub, sd = sqrt(sigb))
      
      ## sample zeta (if HC only) ##
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          for(k in 1:K){
            if(sum(k == Xs1Cols) > 0){
              zeta[b,k]   <- 1
            } else {
              b_alp       <- 1/(A[k]^2) + v/lambda[b-1, k]
              zeta[b,k]   <- rinvgamma(1, v/2 + 1/2, b_alp)
            }
          }
        } else {
          b_alp     <- 1/(A^2) + v/lambda[b-1]
          zeta[b]   <- rinvgamma(1, v/2 + 1/2, b_alp)
        }
      }
      
      ## sample lambda ##
      if(pvargroup == 'kspec'){
        for(k in 1:K){
          if(sum(k == Xs1Cols) > 0){
            Yk            <- Y[kids[,k]]
            Xk            <- X[kids[,k], c(k, k + K)]
            betak         <- beta[b-1,c(k, k + K)]
            b_lam		      <- be + (1/2)*(t(Yk - Xk%*%betak)%*%(Yk - Xk%*%betak))
            lambda[b,k]   <- rinvgamma(1, n + ae, b_lam)
          } else {
            al            <- switch(pvar, fl = 0, hc = v/2, ig = ae)
            bl            <- switch(pvar, fl = 0, hc = v/zeta[b-1, k], ig = be)
            Yk            <- Y[kids[,k]]
            Xk            <- X[kids[,k], c(k, k + K)]
            betak         <- beta[b-1,c(k, k + K)]
            b_lam		      <- bl + (1/2)*(t(Yk - Xk%*%betak)%*%(Yk - Xk%*%betak))
            lambda[b,k]   <- rinvgamma(1, n + al, b_lam)
          }
        }
      } else {
        al          <- switch(pvar, fl = 0, hc = v/2, ig = ae)
        bl          <- switch(pvar, fl = 0, hc = v/zeta[b-1], ig = be)
        b_lam		    <- bl + (1/2)*(t(Y - X%*%beta[b-1,])%*%(Y - X%*%beta[b-1,]))
        lambda[b]   <- rinvgamma(1, N/2 + al, b_lam)
      }
      
      ## sample sigmas ##
      b_sigk		  <- (nu*tau[b-1]/2) + (1/2)*(beta[b-1,] - beb)^2
      sigk[b,]	  <- rgamma(length(b_sigk), (nu+1)/2, b_sigk)
      
      ## sample tau ##
      b_tau		    <- (nu/2)*sum(1/sigk[b,])
      tau[b]		  <- rgamma(1, (2*K*nu)/2, b_tau)
      
      if(verbose){
        if(mod(b, dots) == 0){
          cat('.')
        }
        
        if(mod(b, up) == 0){
          cat(paste("\n",b,"samples completed\n"))
        }
      }
    } 
    
    rho			    <- t(L%*%t(beta))
    var.comp	  <- list(lambda = lambda, sigk = sigk, tau = tau, zeta = zeta)
    mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior)
    data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
    
    out			    <- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
    
  } else {
    xtx			    <- t(X)%*%X
    xty         <- t(X)%*%Y
    beb         <- rep(mu, 2*K)
    bebf			  <- (sumX + 1)/(n + 2)
    prior$beb   <- beb
    
    beta		    <- matrix(0, nrow = B + burnin, ncol = ncol(X))
    if(pvargroup == 'kspec'){
      lambda      <- matrix(1, nrow = B + burnin, ncol = K)
      zeta        <- matrix(1, nrow = B + burnin, ncol = K)
    } else {
      lambda      <- rep(1, length = B + burnin)
      zeta        <- rep(1, length = B + burnin)
    }
    tau			    <- rep(1, length = B + burnin)
    rhokj		    <- matrix(1, nrow = B + burnin, ncol = ncol(X))
    
    beta[1,]	  <- bebf
    
    for (b in 2:(B + burnin)) { 
      ## sample betas ##
      sigt        <- (tau[b-1]/4)*diag(rhokj[b-1,])
      if(pvargroup == 'kspec'){
        lambk       <- c(lambda[b-1,], lambda[b-1,])
        precb		    <- (1/lambk)*xtx + sigt
      } else {
        precb		    <- (1/lambda[b-1])*xtx + sigt
      }
      sigb		    <- 1/diag(precb)
      if(pvargroup == 'kspec'){
        mub			    <- sigb*((1/lambk)*xty + sigt%*%beb)
      } else {
        mub			    <- sigb*((1/lambda[b-1])*xty + sigt%*%beb)
      }
      beta[b,]	  <- rnorm(ncol(X), mean = mub, sd = sqrt(sigb))
      
      ## sample zeta (if HC only) ##
      if(pvar == 'hc'){
        if(pvargroup == 'kspec'){
          for(k in 1:K){
            if(sum(k == Xs1Cols) > 0){
              zeta[b,k]   <- 1
            } else {
              b_alp       <- 1/(A[k]^2) + v/lambda[b-1, k]
              zeta[b,k]   <- rinvgamma(1, v/2 + 1/2, b_alp)
            }
          }
        } else {
          b_alp     <- 1/(A^2) + v/lambda[b-1]
          zeta[b]   <- rinvgamma(1, v/2 + 1/2, b_alp)
        }
      }
      
      ## sample lambda ##
      if(pvargroup == 'kspec'){
        for(k in 1:K){
          if(sum(k == Xs1Cols) > 0){
            Yk            <- Y[kids[,k]]
            Xk            <- X[kids[,k], c(k, k + K)]
            betak         <- beta[b-1,c(k, k + K)]
            b_lam		      <- be + (1/2)*(t(Yk - Xk%*%betak)%*%(Yk - Xk%*%betak))
            lambda[b,k]   <- rinvgamma(1, n + ae, b_lam)
          } else {
            al            <- switch(pvar, fl = 0, hc = v/2, ig = ae)
            bl            <- switch(pvar, fl = 0, hc = v/zeta[b-1, k], ig = be)
            Yk            <- Y[kids[,k]]
            Xk            <- X[kids[,k], c(k, k + K)]
            betak         <- beta[b-1,c(k, k + K)]
            b_lam		      <- bl + (1/2)*(t(Yk - Xk%*%betak)%*%(Yk - Xk%*%betak))
            lambda[b,k]   <- rinvgamma(1, n + al, b_lam)
          }
        }
      } else {
        al          <- switch(pvar, fl = 0, hc = v/2, ig = ae)
        bl          <- switch(pvar, fl = 0, hc = v/zeta[b-1], ig = be)
        b_lam		    <- bl + (1/2)*(t(Y - X%*%beta[b-1,])%*%(Y - X%*%beta[b-1,]))
        lambda[b]   <- rinvgamma(1, N/2 + al, b_lam)
      }
      
      ## sample rhos ##
      mu_rhokj    <- 2/(sqrt(tau[b-1])*abs(beta[b-1,]-beb))
      rhokj[b,]    <- rinvgauss(ncol(X), mean = mu_rhokj, shape = 1)
      
      ## sample tau ##
      b_tau		    <- bt + (1/8)*sum(rhokj[b,]*(beta[b,]-beb)^2)
      tau[b]		  <- rgamma(1, ncol(X)/2 + at, b_tau)
      
      if(verbose){
        if(mod(b, dots) == 0){
          cat('.')
        }
        
        if(mod(b, up) == 0){
          cat(paste("\n",b,"samples completed\n"))
        }
      }
    }
    
    rho			    <- t(L%*%t(beta))
    var.comp	  <- list(lambda = lambda, rhokj = rhokj, tau = tau, zeta = zeta)
    mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior)
    data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
    
    out			    <- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
  }
  class(out)	<- 'bspam'
  
  return(out)
  
}

print.bspam		<- function(mod){
  K		<- ncol(mod$data$X)/2
  n		<- nrow(mod$data$X)/(2*K)
  B		<- mod$mcmc.specs$B
  burnin	<- mod$mcmc.specs$burnin
  pvarm   <- switch(mod$mcmc.specs$prior$pvar, hc = 'Half Cauchy', fl = 'Flat', ig = 'Inverse Gamma')
  typem   <- switch(mod$mcmc.specs$prior$type, l2 = 'Ridge', pe = 'PE Ridge', t = 't', la = 'Laplace', tn = 'TN Ridge')
  cat("Bayesian sparsity adjusted matched-proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  cat(paste("Prior variance:", pvarm,"\n"))
  cat(paste("Penalty:", typem,"\n"))
  
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(apply(mod$rho[-(1:burnin),], 2, median), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nMedian difference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.bspam	<- function(mod, ci.alpha = 0.05, chains = 4){
  K		<- ncol(mod$data$X)/2
  n		<- nrow(mod$data$X)/(2*K)
  B		<- mod$mcmc.specs$B
  burnin	<- mod$mcmc.specs$burnin
  pvarm   <- switch(mod$mcmc.specs$prior$pvar, hc = 'Half Cauchy', fl = 'Flat', ig = 'Inverse Gamma')
  typem   <- switch(mod$mcmc.specs$prior$type, l2 = 'Ridge', pe = 'PE Ridge', t = 't', la = 'Laplace', tn = 'TN Ridge')
  cat("Bayesian sparsity adjusted matched-proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  cat(paste("Prior variance:", pvarm,"\n"))
  cat(paste("Penalty:", typem,"\n"))
  
  ## difference in marginal probability ##
  rho		<- mod$rho
  low		<- ci.alpha/2
  high	<- 1-(ci.alpha/2)
  tab		<- t(apply(rho[-c(1:burnin),], 2, quantile, probs = c(0.5, low, high)))
  prob	<- apply(rho[-c(1:burnin),] > 0, 2, mean)
  GR		<- getGelman(rho[-c(1:burnin),], chains)
  tab1	<- cbind(tab, prob, GR)
  colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)", 'GR Est.', 'Upper GR')
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  cat('Notes: P(r > 0) is posterior probability difference > 0.\nGR denotes Gelman-Rubin, < 1.1 suggests convergence.\n')
  
  cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
  cat(paste(burnin, "(discarded). "))
  
  # sparseTables	<- apply(mod$MCMCsp$prior$sparsityCheck, 2, sum)
  if(length(mod$data$SparseCols) > 0){
    cat(paste('Sparse tables detected: ', paste(mod$data$SparseCols[1:(length(mod$data$SparseCols)/2)], collapse = ", "), '.', sep = '')) 
  }
}

coef.bspam		<- function(mod, CI = FALSE, ci.alpha = 0.05){
  burnin		<- mod$mcmc.specs$burnin
  
  dmpTable	<- apply(mod$rho[-(1:burnin),], 2, quantile, probs = c(0.5, ci.alpha/2, (1 - ci.alpha/2)))
  
  if(CI){
    dmpt	<- dmpTable
  } else {
    dmpt	<- dmpTable[1, ]
  }
  
  return(dmpt)
}

pprob	<- function(mod, ...){
  UseMethod("pprob")
}


pprob.bspam	<- function(mod, rho0 = 0){
  burnin	<- mod$mcmc.specs$burnin
  probs	<- apply(mod$rho[-c(1:burnin),] > rho0, 2, mean)
  
  return(probs)
}

#### Klingenberg and Agresti (2006) GEE ####
gbmp	<- function(X1, X2, family = gaussian(link = 'identity'), corstr = 'independence'){
  n		<- nrow(X1)
  K		<- ncol(X1)
  Y		<- as.matrix(cbind(X1, X2))
  gDat	<- data.frame(y	= c(t(Y)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))
  
  model	<- geeglm(y ~ factor(s) - 1, family = family, data = gDat, id = id, corstr = corstr)
  
  beta	<- coef(model)
  covb	<- summary(model)$cov.scaled
  # covb	<- (t(Y - beta)%*%(Y - beta))/(n^2)
  
  
  # tests of differences #
  L	      <- cbind(diag(K), -diag(K))
  rho		  <- L%*%beta
  seRho   <- sqrt(diag(L%*%covb%*%t(L)))
  lower	  <- rho - qnorm(0.975)*seRho
  upper	  <- rho + qnorm(0.975)*seRho
  wald	  <- (rho/seRho)^2
  pval	  <- 1 - pchisq(wald, 1)
  
  ## Simultaneous Marginal Homogeneity ##
  # KA Generalized Score Test #
  zse     <- which(diag(covb) == 0)
  # 	if(length(zse) == 0){
  #   	Tgs		<- (n^2)*t(rho)%*%solve(L%*%(t(Y)%*%Y)%*%t(L))%*%rho
  #   	pTgs	<- 1 - pchisq(Tgs, nrow(L))
  # 	} else {
  Tgs		<- NaN
  pTgs	<- NaN
  # }
  
  ## if estimating 
  if(corstr == 'unstructured'){
    cory	<- matrix(0, 2*K, 2*K)
    cory[lower.tri(cory, diag = FALSE)]	<- summary(model)$corr[,1]
    cory[upper.tri(cory, diag = FALSE)]	<- t(cory)[upper.tri(t(cory), diag = FALSE)]		
    diag(cory)	<- rep(1, 2*K)
  } else {
    cory	<- NULL
  }
  
  l		<- list(rho = rho, beta = beta, sdv = seRho, cov = list(covb = covb, cory = cory, zse = zse), lower = lower, upper = upper, wald = wald, pval = pval, model = model, smh = list(Tgs = Tgs, pTgs = pTgs), data = list(X1 = X1, X2 = X2))
  
  class(l)	<- 'gbmp'
  
  # 	if(zse > 0){
  #   	warning('Sparse columns detected, variance estimates may be singular')
  # 	}
  
  return(l)
}

print.gbmp	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("GEE-based multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(c(mod$rho), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nDifference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.gbmp	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("GEE-based multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## odds ratio ##
  rho		<- mod$rho
  sdv		<- mod$sdv
  wald	<- mod$wald
  pval	<- mod$pval
  tab1	<- cbind(rho, sdv, wald, pval)
  colnames(tab1)	<- c('Estimate', 'Std.err', 'Wald', 'Pr(>|W|)')
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  cat(paste('Notes: Estimates based on GEE with family: ', mod$model$family$family, ', link: ', mod$model$family$link, '.\nStd.err estimated using robust estimator.\n', sep = ''))
  
  cat("\nTest of simultaneous marginal homogeneity\n")
  cat(paste("Score statistic:", round(mod$smh$Tgs, 3), "on", K, "df, p-value =", round(mod$smh$pTgs, 3)))
  
}

coef.gbmp		<- function(mod){
  dmpt	<- c(mod$rho)
  
  return(dmpt)
}

#### Westfall, Troendle, and Pennello (2010) ####

bootmmp	<- function(X1, X2, B, int.hyp = FALSE){
  X		<- cbind(X1, X2)
  K		<- ncol(X1)
  n		<- nrow(X1)
  
  L		<- cbind(diag(K), -diag(K))
  
  D		<- t(L%*%t(X))
  Dbar	<- apply(D, 2, mean)
  S		<- diag(apply((t(t(D) - Dbar))^2, 2, mean))
  if(length(which(diag(S) == 0)) > 0){
    Z <- NaN
  } else{
    Z		<- sqrt(n)*solve(sqrt(S))%*%Dbar
  }
  
  Dboot	<- matrix(0, nrow = B, ncol = K)
  if(int.hyp){
    Zboot	<- matrix(0, nrow = B, ncol = K)
    Zmax	<- vector('numeric', length = B)
  }
  
  for(b in 1:B){
    idb			<- sample(1:n, n, replace = TRUE)
    Ds			<- D[idb,]
    Dbb			<- Dboot[b,] <- apply(Ds, 2, mean)
    if(int.hyp){
      Sb			<- sqrt(apply((t(t(Ds) - Dbb))^2, 2, mean))
      Zboot[b,]	<- c(sqrt(n)*ifelse(Sb > 0, Dbb/Sb, 0))
      Zmax[b]		<- max(abs(Zboot[b,]))
    }	
  }
  
  resMat	<- cbind(t(apply(Dboot, 2, quantile, probs = c(0.5, 0.025, 0.975))), apply(Dboot > 0, 2, mean))
  
  if(int.hyp){
    l	<- list(resMat = resMat, Dboot = Dboot, Zboot = Zboot, Zmax = Zmax, Z = Z, S = S, Dbar = Dbar, data = list(X1 = X1, X2 = X2))	
  } else {
    l	<- list(resMat = resMat, Dboot = Dboot, data = list(X1 = X1, X2 = X2))
  }
  
  class(l)	<- 'bootmmp'
  
  return(l)
}

print.bootmmp	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("Bootstraped multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(c(mod$resMat[,1]), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nDifference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.bootmmp	<- function(mod, ci.alpha = 0.05){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("Bootstraped multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  low				<- ci.alpha/2
  high			<- 1-(ci.alpha/2)
  tab1			<- mod$resMat
  colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)")
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  
}

coef.bootmmp	<- function(mod){
  dmpt	<- c(mod$resMat[,1])
  
  return(dmpt)
}

#### Misc. Functions ####
getCounts	<- function(X1, X2, K){
  n21		<- vector('numeric', length = K)
  n12		<- vector('numeric', length = K)
  n11		<- vector('numeric', length = K)
  n22		<- vector('numeric', length = K)
  
  for(j in 1:K){
    rawtab	<- table(X2[,j], X1[,j])
    temp	<- matrix(0, nrow = 2, ncol = 2)
    temp[1:dim(rawtab)[1], 1:dim(rawtab)[2]]	<- rawtab
    
    patab	<- matrix(c(temp[2,2], temp[1,2], temp[2,1], temp[1,1]), nrow = 2)
    
    n21[j]	<- patab[2,1]
    n12[j]	<- patab[1,2]
    
    n11[j]	<- patab[1,1]
    n22[j]	<- patab[2,2]
  }
  
  l	<- list(n21 = n21, n12 = n12, n11 = n11, n22 = n22)
  return(l)
}

getGelman	<- function(mcmcOut, chains = 4){
  if(chains != 4 & chains != 2){
    stop('Split chains in half or quarters only')
  }
  K		<- ncol(mcmcOut)
  B		<- nrow(mcmcOut)
  GRR		<- matrix(0, ncol = 2, nrow = K)
  cuts		<- split(1:B, rep(1:chains, each = B/chains))
  
  for(k in 1:K){
    if(chains == 4){
      mcmc1	<- mcmc(mcmcOut[cuts[[1]],k])
      mcmc2	<- mcmc(mcmcOut[cuts[[2]],k])
      mcmc3	<- mcmc(mcmcOut[cuts[[3]],k])
      mcmc4	<- mcmc(mcmcOut[cuts[[4]],k])
      gdb		<- gelman.diag(mcmc.list(mcmc1, mcmc2, mcmc3, mcmc4))
      
      GRR[k,]	<- gdb$psrf
    } else {
      mcmc1	<- mcmc(mcmcOut[cuts[[1]],k])
      mcmc2	<- mcmc(mcmcOut[cuts[[2]],k])
      gdb		<- gelman.diag(mcmc.list(mcmc1, mcmc2))
      
      GRR[k,]	<- gdb$psrf		
    }
  }
  colnames(GRR)	<- c('Point est.', 'Upper C.I.')
  rownames(GRR)	<- 1:K
  
  return(GRR)
  
}

#### Sample data ####

soc <- rbind(matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1), times = 17), nrow = 17, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), times = 16), nrow = 16, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0), times = 5), nrow = 5, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1), times = 5), nrow = 5, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 1), times = 3), nrow = 3, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), times = 3), nrow = 3, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 1), times = 3), nrow = 3, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 1, 0, 0, 0, 0, 0, 0), times = 2), nrow = 2, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0), times = 2), nrow = 2, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1), times = 2), nrow = 2, byrow = TRUE),
             matrix(rep(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 1, 0, 0, 0, 1, 1, 1, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 1, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 1, 0, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 1, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(1, 0, 0, 0, 1, 0, 0, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 1, 1, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 1, 0, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 1, 0, 0, 1, 1, 1, 1), times = 1), nrow = 1, byrow = TRUE))


