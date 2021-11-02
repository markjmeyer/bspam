#### Source Code ####
source('bspam.R')

#### Data and Specs ####
B     <- 10000
X1		<- soc[,1:5]
X2		<- soc[,6:10]
K     <- ncol(X1)

#### BSpAM ####
### t_5 ###
prt     <- list(nu = 5, A = rep(0.5, K), ae = 2.8, be = 0.45)
model_t <- bspam(X1, X2, B, burnin = NULL, penalty = 't', pvar = 'hc', pvargroup = 'kspec', prior = prt, mu = 0.5,
                 up = 2000, verbose = TRUE)
summary(model_t)

### ridge ###
prri    <- list(at = 1, bt = 1, A = rep(0.5, K), ae = 2.8, be = 0.45)
model_r <- bspam(X1, X2, B, burnin = NULL, penalty = 'ri', pvar = 'hc', pvargroup = 'kspec', prior = prri, mu = 0.5,
                 up = 2000, verbose = TRUE)
summary(model_r)

### Laplace ###
prla    <- list(at = 1, bt = 1, A = rep(0.5, K), ae = 2.8, be = 0.45)
model_l <- bspam(X1, X2, B, burnin = NULL, penalty = 'la', pvar = 'hc', pvargroup = 'kspec', prior = prla, mu = 0.5,
                 up = 2000, verbose = TRUE)
summary(model_l)
