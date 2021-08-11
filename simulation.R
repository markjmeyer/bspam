#### Source Code ####
source('bspam.R')

#### Manuscript Simulation ####

### Sparse Settings ###
# simualtion settings #
n       <- 75
K       <- 2
M       <- 1000
mu      <- 0.5 #1/(n+2)
r       <- 0 # 0, 1/3, 2/3
sCov    <- (0.01^2)*matrix(c(1, r, 0, 0,
                             r, 1, 0, 0,
                             0, 0, 1, r,
                             0, 0, r, 1),
                           nrow = 2*K, byrow = TRUE)
theta12  <- c(0.03, 0.05, 0.07, 0.09)

# BSpAM specs #
prt     <- list(nu = 5, A = rep(0.5, K), ae = 2.8, be = 0.45)
B       <- 10000
burnin  <- B

# storage #
br1i75s  <- array(0, dim = c(M, 4, length(theta12)))
br2i75s  <- array(0, dim = c(M, 4, length(theta12)))

pw1i75s  <- array(0, dim = c(M, 4, length(theta12)))
pw2i75s  <- array(0, dim = c(M, 4, length(theta12)))

cr1i75s  <- array(0, dim = c(M, 4, length(theta12)))
cr2i75s  <- array(0, dim = c(M, 4, length(theta12)))

spci75s  <- array(0, dim = c(M, 4, length(theta12)))

for(d in 1:length(theta12)){
  tb      <- c(0.05, theta12[d],
               0.25, 0.001)
  for(m in 1:M){
    set.seed(m)
    Ymc   <- rmvnorm(n, mean = tb, sigma = sCov)
    Ymcl  <- c(t(Ymc))
    Yb    <- rbinom(length(Ymcl), 1, prob = Ymcl*(Ymcl > 0)) # abs or set to zero?
    Ys    <- matrix(Yb, nrow = n, byrow = TRUE)
    X1s   <- Ys[,1:K]
    X2s   <- Ys[,(K+1):(2*K)]
    X 	  <- matrix(unlist(getCounts(X1s, X2s, K)), ncol = K, byrow = TRUE)
    
    if(apply(X1s, 2, sum)[2] == 0 & apply(X2s, 2, sum)[2] == 0){
      next
    }
    Xs1Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0)
    Xs2Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0) + K
    XsCols	<- c(Xs1Cols, Xs2Cols)
    
    ## run models ##
    bspamt  <- bspam(X1s, X2s, B, burnin = burnin, penalty = 't', pvar = 'hc', pvargroup = 'kspec', 
                     prior = prt, mu = mu, up = 2000)
    geek    <- gbmp(X1s, X2s)
    boot    <- bootmmp(X1s, X2s, B = 10000)
    lcr121	<- ifelse(X[1,1] == 0, (X[1,1] + 0.5)/(sum(X[,1] + 0.5)), X[1,1]/sum(X[,1]))
    lcr112	<- ifelse(X[2,1] == 0, (X[2,1] + 0.5)/(sum(X[,1] + 0.5)), X[2,1]/sum(X[,1]))
    lcr221	<- ifelse(X[1,2] == 0, (X[1,2] + 0.5)/(sum(X[,2] + 0.5)), X[1,2]/sum(X[,2]))
    lcr212	<- ifelse(X[2,2] == 0, (X[2,2] + 0.5)/(sum(X[,2] + 0.5)), X[2,2]/sum(X[,2]))
    lcRho1 	<- lcr121 - lcr112
    lcRho2 	<- lcr221 - lcr212
    
    ## extract results ###
    ## bias ##
    # rho #
    br1i75s[m,,d]   <- (tb[1] - tb[3]) - c(median(bspamt$rho[-c(1:burnin),1]), lcRho1, geek$rho[1], boot$resMat[1,1])
    br2i75s[m,,d]   <- (tb[2] - tb[4]) - c(median(bspamt$rho[-c(1:burnin),2]), lcRho2, geek$rho[2], boot$resMat[2,1])
    
    ## intervals ##
    # estimates #
    bspamti1    <- quantile(bspamt$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    geeki1      <- c(geek$lower[1], geek$upper[1])
    booti1      <- boot$resMat[1,2:3]
    
    bspamti2    <- quantile(bspamt$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    geeki2      <- c(geek$lower[2], geek$upper[2])
    booti2      <- boot$resMat[2,2:3]
    
    oldw 			<- getOption("warn")
    options(warn = -1)
    lcri1 		<- prop.test(n*c(lcr121, lcr112), c(n,n), correct = FALSE)
    
    lcri2 		<- prop.test(n*c(lcr221, lcr212), c(n,n), correct = FALSE)
    options(warn = oldw)
    
    # power #
    pw1i75s[m,,d]   <- 1*c((prod(bspamti1) > 0), (lcri1$p.val < 0.05), (geek$pval[1] < 0.05), (max(booti1) < 0))
    pw2i75s[m,,d]   <- 1*c((prod(bspamti2) > 0), (lcri2$p.val < 0.05), (geek$pval[2] < 0.05), (max(booti2) > 0))
    
    # coverage #
    cr1i75s[m,,d]   <- 1*c((bspamti1[1] < tb[1] - tb[3] & bspamti1[2] > tb[1] - tb[3]),
                           ((lcri1$conf[1] < (tb[1] - tb[3])) & (lcri1$conf[2] > (tb[1] - tb[3]))),
                           (geeki1[1] < tb[1] - tb[3] & geeki1[2] > tb[1] - tb[3]),
                           (booti1[1] < tb[1] - tb[3] & booti1[2] > tb[1] - tb[3]))
    cr2i75s[m,,d]   <- 1*c((bspamti2[1] < tb[2] - tb[4] & bspamti2[2] > tb[2] - tb[4]),
                           ((lcri2$conf[1] < (tb[2] - tb[4])) & (lcri2$conf[2] > (tb[2] - tb[4]))),
                           (geeki2[1] < tb[2] - tb[4] & geeki2[2] > tb[2] - tb[4]),
                           (booti2[1] < tb[2] - tb[4] & booti2[2] > tb[2] - tb[4]))
    
    ## sparsity checks ##
    spci75s[m,1,d]  <- length(XsCols) > 0
    spci75s[m,2,d]  <- sum(apply(Ys, 2, sum) == 0) > 0
  }
}

bias1   <- matrix(0, nrow = 4, ncol = length(theta12))
bias2   <- matrix(0, nrow = 4, ncol = length(theta12))

power1  <- matrix(0, nrow = 4, ncol = length(theta12))
power2  <- matrix(0, nrow = 4, ncol = length(theta12))

cover1  <- matrix(0, nrow = 4, ncol = length(theta12))
cover2  <- matrix(0, nrow = 4, ncol = length(theta12))

for(i in 1:length(theta12)){
  b1i75s   <- abs(apply(br1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  b2i75s   <- abs(apply(br2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  p1i75s   <- abs(apply(pw1i75s[spci75s[,2,i] == 1,,i], 2, mean, na.rm = TRUE))
  p2i75s   <- abs(apply(pw2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  c1i75s   <- abs(apply(cr1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  c2i75s   <- abs(apply(cr2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  bias1[,i]   <- b1i75s
  bias2[,i]   <- b2i75s
  
  power1[,i]  <- p1i75s
  power2[,i]  <- p2i75s
  
  cover1[,i]  <- c1i75s
  cover2[,i]  <- c2i75s
}

bias2

cover2

power2

### Non-sparse ###
# simualtion settings #
n       <- 75
K       <- 2
M       <- 1000
mu      <- 0.5 #1/(n+2)
r       <- 0 # 0, 1/3, 2/3
sCov    <- (0.01^2)*matrix(c(1, r, 0, 0,
                             r, 1, 0, 0,
                             0, 0, 1, r,
                             0, 0, r, 1),
                           nrow = 2*K, byrow = TRUE)
theta21  <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

# t prior #
prt     <- list(nu = 5, A = rep(0.5, K), ae = 2.8, be = 0.45)

# all models #
br1i75s  <- array(0, dim = c(M, 4, length(theta21)))
br2i75s  <- array(0, dim = c(M, 4, length(theta21)))

pw1i75s  <- array(0, dim = c(M, 4, length(theta21)))
pw2i75s  <- array(0, dim = c(M, 4, length(theta21)))

cr1i75s  <- array(0, dim = c(M, 4, length(theta21)))
cr2i75s  <- array(0, dim = c(M, 4, length(theta21)))

spci75s  <- array(0, dim = c(M, 4, length(theta21)))

for(d in 1:length(theta21)){
  tb      <- c(0.05, 0.05,
               theta21[d], 0.001)
  for(m in 1:M){
    set.seed(m)
    Ymc   <- rmvnorm(n, mean = tb, sigma = sCov)
    Ymcl  <- c(t(Ymc))
    Yb    <- rbinom(length(Ymcl), 1, prob = Ymcl*(Ymcl > 0)) # abs or set to zero?
    Ys    <- matrix(Yb, nrow = n, byrow = TRUE)
    X1s   <- Ys[,1:K]
    X2s   <- Ys[,(K+1):(2*K)]
    X 	  <- matrix(unlist(getCounts(X1s, X2s, K)), ncol = K, byrow = TRUE)
    
    if(apply(X1s, 2, sum)[2] == 0 & apply(X2s, 2, sum)[2] == 0){
      next
    }
    Xs1Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0)
    Xs2Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0) + K
    XsCols	<- c(Xs1Cols, Xs2Cols)
    
    ## run models ##
    bspamt  <- bspam(X1s, X2s, B, burnin = burnin, penalty = 't', pvar = 'hc', pvargroup = 'kspec', 
                     prior = prt, mu = mu, up = 2000)
    geek    <- gbmp(X1s, X2s)
    boot    <- bootmmp(X1s, X2s, B = 10000)
    lcr121	<- ifelse(X[1,1] == 0, (X[1,1] + 0.5)/(sum(X[,1] + 0.5)), X[1,1]/sum(X[,1]))
    lcr112	<- ifelse(X[2,1] == 0, (X[2,1] + 0.5)/(sum(X[,1] + 0.5)), X[2,1]/sum(X[,1]))
    lcr221	<- ifelse(X[1,2] == 0, (X[1,2] + 0.5)/(sum(X[,2] + 0.5)), X[1,2]/sum(X[,2]))
    lcr212	<- ifelse(X[2,2] == 0, (X[2,2] + 0.5)/(sum(X[,2] + 0.5)), X[2,2]/sum(X[,2]))
    lcRho1 	<- lcr121 - lcr112
    lcRho2 	<- lcr221 - lcr212
    
    ## extract results ###
    ## bias ##
    # rho #
    br1i75s[m,,d]   <- (tb[1] - tb[3]) - c(median(bspamt$rho[-c(1:burnin),1]), lcRho1, geek$rho[1], boot$resMat[1,1])
    br2i75s[m,,d]   <- (tb[2] - tb[4]) - c(median(bspamt$rho[-c(1:burnin),2]), lcRho2, geek$rho[2], boot$resMat[2,1])
    
    ## intervals ##
    # estimates #
    bspamti1    <- quantile(bspamt$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    geeki1      <- c(geek$lower[1], geek$upper[1])
    booti1      <- boot$resMat[1,2:3]
    
    bspamti2    <- quantile(bspamt$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    geeki2      <- c(geek$lower[2], geek$upper[2])
    booti2      <- boot$resMat[2,2:3]
    
    oldw 			<- getOption("warn")
    options(warn = -1)
    lcri1 		<- prop.test(n*c(lcr121, lcr112), c(n,n), correct = FALSE)
    
    lcri2 		<- prop.test(n*c(lcr221, lcr212), c(n,n), correct = FALSE)
    options(warn = oldw)
    
    # power #
    pw1i75s[m,,d]   <- 1*c((prod(bspamti1) > 0), (lcri1$p.val < 0.05), (geek$pval[1] < 0.05), (max(booti1) < 0))
    pw2i75s[m,,d]   <- 1*c((prod(bspamti2) > 0), (lcri2$p.val < 0.05), (geek$pval[2] < 0.05), (max(booti2) > 0))
    
    # coverage #
    cr1i75s[m,,d]   <- 1*c((bspamti1[1] < tb[1] - tb[3] & bspamti1[2] > tb[1] - tb[3]),
                           ((lcri1$conf[1] < (tb[1] - tb[3])) & (lcri1$conf[2] > (tb[1] - tb[3]))),
                           (geeki1[1] < tb[1] - tb[3] & geeki1[2] > tb[1] - tb[3]),
                           (booti1[1] < tb[1] - tb[3] & booti1[2] > tb[1] - tb[3]))
    cr2i75s[m,,d]   <- 1*c((bspamti2[1] < tb[2] - tb[4] & bspamti2[2] > tb[2] - tb[4]),
                           ((lcri2$conf[1] < (tb[2] - tb[4])) & (lcri2$conf[2] > (tb[2] - tb[4]))),
                           (geeki2[1] < tb[2] - tb[4] & geeki2[2] > tb[2] - tb[4]),
                           (booti2[1] < tb[2] - tb[4] & booti2[2] > tb[2] - tb[4]))
    
    ## sparsity checks ##
    spci75s[m,1,d]  <- length(XsCols) > 0
    spci75s[m,2,d]  <- sum(apply(Ys, 2, sum) == 0) > 0
  }
}

bias1   <- matrix(0, nrow = 4, ncol = length(theta21))
bias2   <- matrix(0, nrow = 4, ncol = length(theta21))

power1  <- matrix(0, nrow = 4, ncol = length(theta21))
power2  <- matrix(0, nrow = 4, ncol = length(theta21))

cover1  <- matrix(0, nrow = 4, ncol = length(theta21))
cover2  <- matrix(0, nrow = 4, ncol = length(theta21))

for(i in 1:length(theta21)){
  b1i75s   <- abs(apply(br1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  b2i75s   <- abs(apply(br2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  p1i75s   <- abs(apply(pw1i75s[spci75s[,2,i] == 1,,i], 2, mean, na.rm = TRUE))
  p2i75s   <- abs(apply(pw2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  c1i75s   <- abs(apply(cr1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  c2i75s   <- abs(apply(cr2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  bias1[,i]   <- b1i75s
  bias2[,i]   <- b2i75s
  
  power1[,i]  <- p1i75s
  power2[,i]  <- p2i75s
  
  cover1[,i]  <- c1i75s
  cover2[,i]  <- c2i75s
}

bias1

cover1

power1


#### Supplement Simulation ####
### Sparse ###
# simualtion settings #
n       <- 75
K       <- 2
M       <- 1000
mu      <- 0.5 #1/(n+2)
r       <- 0
sCov    <- (0.01^2)*matrix(c(1, r, 0, 0,
                             r, 1, 0, 0,
                             0, 0, 1, r,
                             0, 0, r, 1),
                           nrow = 2*K, byrow = TRUE)
theta12  <- c(0.03, 0.05, 0.07, 0.09)

# BSpAM priors #
prpe    <- list(at = 1, bt = 1, as = 1, bs = 1, A = rep(0.5, K), ae = aer, be = ber)
prri    <- list(at = 1, bt = 1, A = rep(0.5, K), ae = aer, be = ber)
prt     <- list(nu = 5, A = rep(0.5, K), ae = aer, be = ber)
prla    <- list(at = 1, bt = 1, A = rep(0.5, K), ae = aer, be = ber)

# storage #
br1i75s  <- array(0, dim = c(M, 4, length(theta12)))
br2i75s  <- array(0, dim = c(M, 4, length(theta12)))

pw1i75s  <- array(0, dim = c(M, 4, length(theta12)))
pw2i75s  <- array(0, dim = c(M, 4, length(theta12)))

cr1i75s  <- array(0, dim = c(M, 4, length(theta12)))
cr2i75s  <- array(0, dim = c(M, 4, length(theta12)))

spci75s  <- array(0, dim = c(M, 4, length(theta12)))

for(d in 1:length(theta12)){
  tb      <- c(0.05, theta12[d],
               0.25, 0.001)
  for(m in 1:M){
    set.seed(m)
    Ymc   <- rmvnorm(n, mean = tb, sigma = sCov)
    Ymcl  <- c(t(Ymc))
    Yb    <- rbinom(length(Ymcl), 1, prob = Ymcl*(Ymcl > 0)) # abs or set to zero?
    Ys    <- matrix(Yb, nrow = n, byrow = TRUE)
    X1s   <- Ys[,1:K]
    X2s   <- Ys[,(K+1):(2*K)]
    X 	  <- matrix(unlist(getCounts(X1s, X2s, K)), ncol = K, byrow = TRUE)
    
    if(apply(X1s, 2, sum)[2] == 0 & apply(X2s, 2, sum)[2] == 0){
      next
    }
    Xs1Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0)
    Xs2Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0) + K
    XsCols	<- c(Xs1Cols, Xs2Cols)
    
    ## run models ##
    bspamp    <- bspam(X1s, X2s, B, burnin = burnin, penalty = 'ri', ptype = 'pe', pvar = 'hc', pvargroup = 'kspec', 
                     prior = prpe, mu = mu, up = 2000)
    bspamr    <- bspam(X1s, X2s, B, burnin = burnin, penalty = 'ri', ptype = 'l2', pvar = 'hc', pvargroup = 'kspec', 
                     prior = prri, mu = mu, up = 2000)
    bspamt   <- bspam(X1s, X2s, B, burnin = burnin, penalty = 't', pvar = 'hc', pvargroup = 'kspec', 
                     prior = prt, mu = mu, up = 2000)
    bspaml    <- bspam(X1s, X2s, B, burnin = burnin, penalty = 'la', pvar = 'hc', pvargroup = 'kspec', 
                     prior = prla, mu = mu, up = 2000)
    
    ## extract results ###
    ## bias ##
    # rho #
    br1i75s[m,]  <- (tb[1] - tb[3]) - c(median(bspamp$rho[-c(1:burnin),1]), median(bspamr$rho[-c(1:burnin),1]),
                                        median(bspamt$rho[-c(1:burnin),1]), median(bspaml$rho[-c(1:burnin),1]))
    br2i75s[m,]  <- (tb[2] - tb[4]) - c(median(bspamp$rho[-c(1:burnin),2]), median(bspamr$rho[-c(1:burnin),2]),
                                        median(bspamt$rho[-c(1:burnin),2]), median(bspaml$rho[-c(1:burnin),2]))
    
    ## intervals ##
    # estimates #
    bspampi1      <- quantile(bspamp$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    bspamri1      <- quantile(bspamr$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    bspamti1      <- quantile(bspamt$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    bspamli1      <- quantile(bspaml$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    
    bspampi2      <- quantile(bspamp$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    bspamri2      <- quantile(bspamr$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    bspamti2      <- quantile(bspamt$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    bspamli2      <- quantile(bspaml$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    
    # int width #
    wr1i75s[m,]  <- c(abs(diff(bspampi1)), abs(diff(bspamri1)), abs(diff(bspamti1)), abs(diff(bspamli1)))
    wr2i75s[m,]  <- c(abs(diff(bspampi2)), abs(diff(bspamri2)), abs(diff(bspamti2)), abs(diff(bspamli2)))
    
    # power #
    pw1i75s[m,]  <- 1*c((prod(bspampi1) > 0), (prod(bspamri1) > 0), (prod(bspamti1) > 0), (prod(bspamli1) > 0))
    pw2i75s[m,]  <- 1*c((prod(bspampi2) > 0), (prod(bspamri2) > 0), (prod(bspamti2) > 0), (prod(bspamli2) > 0))
    
    # coverage #
    cr1i75s[m,]  <- 1*c((bspampi1[1] < tb[1] - tb[3] & bspampi1[2] > tb[1] - tb[3]),
                        (bspamri1[1] < tb[1] - tb[3] & bspamri1[2] > tb[1] - tb[3]),
                        (bspamti1[1] < tb[1] - tb[3] & bspamti1[2] > tb[1] - tb[3]),
                        (bspamli1[1] < tb[1] - tb[3] & bspamli1[2] > tb[1] - tb[3]))
    cr2i75s[m,]  <- 1*c((bspampi2[1] < tb[2] - tb[4] & bspampi2[2] > tb[2] - tb[4]),
                        (bspamri2[1] < tb[2] - tb[4] & bspamri2[2] > tb[2] - tb[4]),
                        (bspamti2[1] < tb[2] - tb[4] & bspamti2[2] > tb[2] - tb[4]),
                        (bspamli2[1] < tb[2] - tb[4] & bspamli2[2] > tb[2] - tb[4]))
    
    ## sparsity checks ##
    spci75s[m,1,d]  <- length(XsCols) > 0
    spci75s[m,2,d]  <- sum(apply(Ys, 2, sum) == 0) > 0
  }
}

bias1   <- matrix(0, nrow = 4, ncol = length(theta12))
bias2   <- matrix(0, nrow = 4, ncol = length(theta12))

power1  <- matrix(0, nrow = 4, ncol = length(theta12))
power2  <- matrix(0, nrow = 4, ncol = length(theta12))

cover1  <- matrix(0, nrow = 4, ncol = length(theta12))
cover2  <- matrix(0, nrow = 4, ncol = length(theta12))

for(i in 1:length(theta12)){
  b1i75s   <- abs(apply(br1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  b2i75s   <- abs(apply(br2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  p1i75s   <- abs(apply(pw1i75s[spci75s[,2,i] == 1,,i], 2, mean, na.rm = TRUE))
  p2i75s   <- abs(apply(pw2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  c1i75s   <- abs(apply(cr1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  c2i75s   <- abs(apply(cr2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  bias1[,i]   <- b1i75s
  bias2[,i]   <- b2i75s
  
  power1[,i]  <- p1i75s
  power2[,i]  <- p2i75s
  
  cover1[,i]  <- c1i75s
  cover2[,i]  <- c2i75s
}

bias2

cover2

power2


### Non-sparse ###
# simualtion settings #
n       <- 75
K       <- 2
M       <- 1000
mu      <- 0.5 #1/(n+2)
r       <- 0
sCov    <- (0.01^2)*matrix(c(1, r, 0, 0,
                             r, 1, 0, 0,
                             0, 0, 1, r,
                             0, 0, r, 1),
                           nrow = 2*K, byrow = TRUE)
theta21  <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

# BSpAM priors #
prpe    <- list(at = 1, bt = 1, as = 1, bs = 1, A = rep(0.5, K), ae = aer, be = ber)
prri    <- list(at = 1, bt = 1, A = rep(0.5, K), ae = aer, be = ber)
prt     <- list(nu = 5, A = rep(0.5, K), ae = aer, be = ber)
prla    <- list(at = 1, bt = 1, A = rep(0.5, K), ae = aer, be = ber)

# all models #
br1i75s  <- array(0, dim = c(M, 4, length(theta21)))
br2i75s  <- array(0, dim = c(M, 4, length(theta21)))

pw1i75s  <- array(0, dim = c(M, 4, length(theta21)))
pw2i75s  <- array(0, dim = c(M, 4, length(theta21)))

cr1i75s  <- array(0, dim = c(M, 4, length(theta21)))
cr2i75s  <- array(0, dim = c(M, 4, length(theta21)))

spci75s  <- array(0, dim = c(M, 4, length(theta21)))

for(d in 1:length(theta21)){
  tb      <- c(0.05, 0.05,
               theta21[d], 0.001)
  for(m in 1:M){
    set.seed(m)
    Ymc   <- rmvnorm(n, mean = tb, sigma = sCov)
    Ymcl  <- c(t(Ymc))
    Yb    <- rbinom(length(Ymcl), 1, prob = Ymcl*(Ymcl > 0)) # abs or set to zero?
    Ys    <- matrix(Yb, nrow = n, byrow = TRUE)
    X1s   <- Ys[,1:K]
    X2s   <- Ys[,(K+1):(2*K)]
    X 	  <- matrix(unlist(getCounts(X1s, X2s, K)), ncol = K, byrow = TRUE)
    
    if(apply(X1s, 2, sum)[2] == 0 & apply(X2s, 2, sum)[2] == 0){
      next
    }
    Xs1Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0)
    Xs2Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0) + K
    XsCols	<- c(Xs1Cols, Xs2Cols)
    
    ## run models ##
    bspamp    <- bspam(X1s, X2s, B, burnin = burnin, penalty = 'ri', ptype = 'pe', pvar = 'hc', pvargroup = 'kspec', 
                       prior = prpe, mu = mu, up = 2000)
    bspamr    <- bspam(X1s, X2s, B, burnin = burnin, penalty = 'ri', ptype = 'l2', pvar = 'hc', pvargroup = 'kspec', 
                       prior = prri, mu = mu, up = 2000)
    bspamt   <- bspam(X1s, X2s, B, burnin = burnin, penalty = 't', pvar = 'hc', pvargroup = 'kspec', 
                      prior = prt, mu = mu, up = 2000)
    bspaml    <- bspam(X1s, X2s, B, burnin = burnin, penalty = 'la', pvar = 'hc', pvargroup = 'kspec', 
                       prior = prla, mu = mu, up = 2000)
    
    ## extract results ###
    ## bias ##
    # rho #
    br1i75s[m,]  <- (tb[1] - tb[3]) - c(median(bspamp$rho[-c(1:burnin),1]), median(bspamr$rho[-c(1:burnin),1]),
                                        median(bspamt$rho[-c(1:burnin),1]), median(bspaml$rho[-c(1:burnin),1]))
    br2i75s[m,]  <- (tb[2] - tb[4]) - c(median(bspamp$rho[-c(1:burnin),2]), median(bspamr$rho[-c(1:burnin),2]),
                                        median(bspamt$rho[-c(1:burnin),2]), median(bspaml$rho[-c(1:burnin),2]))
    
    ## intervals ##
    # estimates #
    bspampi1      <- quantile(bspamp$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    bspamri1      <- quantile(bspamr$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    bspamti1      <- quantile(bspamt$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    bspamli1      <- quantile(bspaml$rho[-c(1:burnin),1], probs = c(0.025, 0.975))
    
    bspampi2      <- quantile(bspamp$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    bspamri2      <- quantile(bspamr$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    bspamti2      <- quantile(bspamt$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    bspamli2      <- quantile(bspaml$rho[-c(1:burnin),2], probs = c(0.025, 0.975))
    
    # int width #
    wr1i75s[m,]  <- c(abs(diff(bspampi1)), abs(diff(bspamri1)), abs(diff(bspamti1)), abs(diff(bspamli1)))
    wr2i75s[m,]  <- c(abs(diff(bspampi2)), abs(diff(bspamri2)), abs(diff(bspamti2)), abs(diff(bspamli2)))
    
    # power #
    pw1i75s[m,]  <- 1*c((prod(bspampi1) > 0), (prod(bspamri1) > 0), (prod(bspamti1) > 0), (prod(bspamli1) > 0))
    pw2i75s[m,]  <- 1*c((prod(bspampi2) > 0), (prod(bspamri2) > 0), (prod(bspamti2) > 0), (prod(bspamli2) > 0))
    
    # coverage #
    cr1i75s[m,]  <- 1*c((bspampi1[1] < tb[1] - tb[3] & bspampi1[2] > tb[1] - tb[3]),
                        (bspamri1[1] < tb[1] - tb[3] & bspamri1[2] > tb[1] - tb[3]),
                        (bspamti1[1] < tb[1] - tb[3] & bspamti1[2] > tb[1] - tb[3]),
                        (bspamli1[1] < tb[1] - tb[3] & bspamli1[2] > tb[1] - tb[3]))
    cr2i75s[m,]  <- 1*c((bspampi2[1] < tb[2] - tb[4] & bspampi2[2] > tb[2] - tb[4]),
                        (bspamri2[1] < tb[2] - tb[4] & bspamri2[2] > tb[2] - tb[4]),
                        (bspamti2[1] < tb[2] - tb[4] & bspamti2[2] > tb[2] - tb[4]),
                        (bspamli2[1] < tb[2] - tb[4] & bspamli2[2] > tb[2] - tb[4]))
    
    ## sparsity checks ##
    spci75s[m,1,d]  <- length(XsCols) > 0
    spci75s[m,2,d]  <- sum(apply(Ys, 2, sum) == 0) > 0
  }
}

bias1   <- matrix(0, nrow = 4, ncol = length(theta21))
bias2   <- matrix(0, nrow = 4, ncol = length(theta21))

power1  <- matrix(0, nrow = 4, ncol = length(theta21))
power2  <- matrix(0, nrow = 4, ncol = length(theta21))

cover1  <- matrix(0, nrow = 4, ncol = length(theta21))
cover2  <- matrix(0, nrow = 4, ncol = length(theta21))

for(i in 1:length(theta21)){
  b1i75s   <- abs(apply(br1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  b2i75s   <- abs(apply(br2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  p1i75s   <- abs(apply(pw1i75s[spci75s[,2,i] == 1,,i], 2, mean, na.rm = TRUE))
  p2i75s   <- abs(apply(pw2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  c1i75s   <- abs(apply(cr1i75s[spci75s[,2,i] == 1,,i], 2, mean))
  c2i75s   <- abs(apply(cr2i75s[spci75s[,2,i] == 1,,i], 2, mean))
  
  bias1[,i]   <- b1i75s
  bias2[,i]   <- b2i75s
  
  power1[,i]  <- p1i75s
  power2[,i]  <- p2i75s
  
  cover1[,i]  <- c1i75s
  cover2[,i]  <- c2i75s
}

bias1

cover1

power1
