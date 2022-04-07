# bspam
This repository contains code for running Bayesian Sparse-response Adjusted Marginal (BSpAM) models for Multivariate Matched Proportions data as described by Meyer, Li, and Knutson (2022), currently under review, but available at https://arxiv.org/abs/2108.04096.

# File details
The file bspam.R is the source file for all functions related to running the models. Code used to run our simulations is in simulation.R while the a sample script is contained in sample_script.R.

# Brief Function and Object Description
All functions two matrices, X1 and X2, as arguments. X1 is an n x K matrix of the first set of binary multivariate outcomes. X2 is an n x K matrix of the second or paired set of binary multivariate outcomes. Each function has additional arguments for prior specification and output controls. Priors default to those described in the manuscript.

BSpAM model:

> bspam(X1, X2, B, burnin = NULL, penalty = c('ri', 't', 'la'), ptype = c('l2', 'pe', 'tn'), pvar = c('fl', 'hc', 'ig'), 
      pvargroup = c('kspec', 'global'), prior = NULL, mu = NULL, up = 100, dots = up/10, verbose = TRUE)

We also prepare code for two non-Bayesian methods compared in the manuscript.

Klingenberg & Agresti (2006) GEE-based model:

> gbmp(X1, X2, family = gaussian(link = 'identity'), corstr = 'independence')

Westfall, Troendle, and Pennello (2010) Bootstrap-based model:

> bootmmp(X1, X2, B, int.hyp = FALSE)

The data, in the form of Table 1 from the manuscript, is in the object soc:

> head(soc)

Notes:
- gbmp produces the test for simultaneious homogeneity described by Klingenberg & Agresti (2006), by default
- int.hyp = TRUE implements the intersection hypothesis described by Westfall, Troendle, and Pennello (2010)
