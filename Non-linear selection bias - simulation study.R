
##########   NON-LINEAR INTERACTIONS AND COLLIDER BIAS   ##########

## This file contains R code in support of the paper "XXXX".
## Adapted from https://github.com/agkatzionis/Interactions-and-collider-bias

## Author: Ellie Curnow
## Date: 13/12/2024
## Updated:

## Set up R.
setwd("...")
options(digits = 4)
library(xtable)
library(mfp2)

##################################################

##########  R FUNCTIONS   ##########

## Standard expit and logit functions.
expit <- function (x) exp(x) / (1 + exp(x))
logit <- function (x) log(x / (1 - x))

## Auxiliary function that calculates the intercept delta0
## to yield a specified proportion of individuals selected.
## This can be done analytically for the log-additive
## collider model but only numerically for the other models.
compute.d0 <- function(x.distr, y.model, y.family, beta0, beta1, beta2, s.model, 
                       s.model.nl, sel.prob, delta1, delta2, delta3, d4, #latent.sigma,#
                       n.d0) {
  
  ## Simulate exposure data
  ## Simulate exposure data.
  
  if (x.distr == "binary") {
    X <- rbinom(n.d0, 1, 0.3) 
  } else if (x.distr == "unif") {
    X <- runif(n.d0, 0.1, 0.9)
  } else if (x.distr == "normal") {
    X <- rnorm(n.d0, 0.5, 1)
  }
  
  ## Simulate outcome data
  if (y.model == "linear") {
    linpred <- beta0 + beta1 * X 
  } else if (y.model == "quadratic") {
    linpred <- beta0 + beta1 * X + beta2 * X^2 
  }
  
  if (y.family == "gaussian") {
    Y <- linpred + rnorm(n.d0, 0, 0.25)
  } else if (y.family == "binomial") {
    Y.pr <- expit(linpred)
    Y <- rbinom(n.d0, 1, prob = Y.pr)
  } else if (y.family == "poisson") {
    Y.lambda <- exp(linpred)
    Y <- rpois(n.d0, lambda = Y.lambda)
  }
  ## Compute delta0.
  if (s.model == "logadd") {
    if (s.model.nl == "None"){
      d0 <- log(sel.prob) - log( mean( exp(delta1 * X + delta2 * Y + delta3 * X * Y) ) )
    } else if (s.model.nl == "X^2"){
      d0 <- log(sel.prob) - log( mean( exp(delta1 * X + delta2 * Y + delta3 * X * Y + d4 * X^2) ) )
    } else if (s.model.nl == "X^2Y"){
      d0 <- log(sel.prob) - log( mean( exp(delta1 * X + delta2 * Y + delta3 * X * Y + d4 * X^2 * Y) ) )
    } else if (s.model.nl == "Y^2"){
      d0 <- log(sel.prob) - log( mean( exp(delta1 * X + delta2 * Y + delta3 * X * Y + d4 * Y^2) ) )
    }
  }
  ## Out of scope?
  #else if (s.model == "logistic") {
  #damnit <- function (d0) sum((expit(d0 + delta1 * X + delta2 * Y + delta3 * X * Y) - sel.prob)^2)
  #d0 <- optimize(damnit, interval = c(-10, 10))$minimum
  # find.d0 <- function (d0) mean(expit(d0 + delta1 * X + delta2 * Y + delta3 * X * Y)) - sel.prob
  #  d0 <- uniroot(find.d0, interval = c(-10, 10))$root
  #} else if (s.model == "probit") {
  # S.pr <- delta1 * X + delta2 * Y + delta3 * X * Y + rnorm(n.d0, 0, latent.sigma)
  #d0 <- - sort(S.pr)[floor((1 - sel.prob) * n.d0)]
  #} else if (s.model == "double") {
  # d0 <- 0
  #S.pr <- delta1 * X + delta2 * Y + delta3 * X * Y + rnorm(n.d0, 0, latent.sigma)
  #r1 <- sort(S.pr)[floor(sel.prob / 2 * n.d0)]
  #r2 <- sort(S.pr)[floor((1 - sel.prob / 2) * n.d0)]
  #}
  
  ## Return.
  #if (s.model == "double") {
  # return (c(d0, r1, r2))
  #} else {
  return (d0)
  #}
}

## ---------- SIMULATION GENERATOR ---------- ##

## Here is the main simulation function ##
interaction.sims <- function (iter = 1, n = 1e6, x.distr = "unif", y.family = "gaussian",
                              y.model = "linear", xy.fit = "linear",
                              beta0 = 0.1, beta1 = 0.5, beta2 = 0.25,
                              s.model = "logadd", s.model.nl = "X^2Y",
                              sel.prob = 0.5, 
                              delta1 = 0.25, delta2 = 0.25, delta3 = -0.5, 
                              inter.range = seq(-0.5, 0.5, length = 11), #latent.sigma = 1.6, 
                              n.d0 = 1e6) {
  
  ## Store main results for the X-Y associations here.
  k <- length(inter.range)

  xy.lincoef <- matrix(0, iter, k); colnames(xy.lincoef) <- paste("d4 = ", inter.range, sep = "")
  xy.lincoef.se <- matrix(0, iter, k); colnames(xy.lincoef.se) <- paste("d4 = ", inter.range, sep = "")
  xy.int <- matrix(0, iter, k); colnames(xy.int) <- paste("d4 = ", inter.range, sep = "")
  xy.int.se <- matrix(0, iter, k); colnames(xy.int.se) <- paste("d4 = ", inter.range, sep = "")
  xy.quadcoef <- matrix(0, iter, k); colnames(xy.quadcoef) <- paste("d4 = ", inter.range, sep = "")
  xy.quadcoef.se <- matrix(0, iter, k); colnames(xy.quadcoef.se) <- paste("d4 = ", inter.range, sep = "") 
  xy.cubcoef <- matrix(0, iter, k); colnames(xy.cubcoef) <- paste("d4 = ", inter.range, sep = "")
  xy.cubcoef.se <- matrix(0, iter, k); colnames(xy.cubcoef.se) <- paste("d4 = ", inter.range, sep = "") 
  
  fp.powers.sel.X1 <- matrix(0, iter, k); colnames(fp.powers.sel.X1) <- paste("d4 = ", inter.range, sep = "") 
  fp.powers.sel.X2 <- matrix(0, iter, k); colnames(fp.powers.sel.X2) <- paste("d4 = ", inter.range, sep = "") 
  fp.powers.sel.X3 <- matrix(0, iter, k); colnames(fp.powers.sel.X3) <- paste("d4 = ", inter.range, sep = "") 
  fp.powers.full.X1 <- matrix(0, iter, k); colnames(fp.powers.full.X1) <- paste("d4 = ", inter.range, sep = "") 
  fp.powers.full.X2 <- matrix(0, iter, k); colnames(fp.powers.full.X2) <- paste("d4 = ", inter.range, sep = "") 
  fp.powers.full.X3 <- matrix(0, iter, k); colnames(fp.powers.full.X3) <- paste("d4 = ", inter.range, sep = "") 
  
  # Not considered in this study?
  ## Store results from fitting selection models here.
  #m1.delta <- matrix(0, iter, K); colnames(m1.delta) <- paste("I = ", inter.range, sep = "")
  #m1.se <- matrix(0, iter, K); colnames(m1.se) <- paste("I = ", inter.range, sep = "")
  #m2.delta <- matrix(0, iter, K); colnames(m2.delta) <- paste("I = ", inter.range, sep = "")
  #m2.se <- matrix(0, iter, K); colnames(m2.se) <- paste("I = ", inter.range, sep = "")
  #m3.delta <- matrix(0, iter, K); colnames(m3.delta) <- paste("I = ", inter.range, sep = "")
  #m3.se <- matrix(0, iter, K); colnames(m3.se) <- paste("I = ", inter.range, sep = "")
  
  ## Store additional diagnostics here
  min.s.pr <- matrix(0, iter, k); colnames(min.s.pr) <- paste("d4 = ", inter.range, sep = "")
  max.s.pr <- matrix(0, iter, k); colnames(max.s.pr) <- paste("d4 = ", inter.range, sep = "")
  prop.selected <- matrix(0, iter, k); colnames(prop.selected) <- paste("d4 = ", inter.range, sep = "")
  mean.y <- matrix(0, iter, k); colnames(mean.y) <- paste("d4 = ", inter.range, sep = "")
  cutoff <- matrix(0, iter, k); colnames(cutoff) <- paste("d4 = ", inter.range, sep = "")
  delta0 <- rep(0, k);  names(delta0) <- paste("d4 = ", inter.range, sep = "")
  delta4 <- rep(0, k);  names(delta4) <- paste("d4 = ", inter.range, sep = "")
  #if (s.model == "double") {
  #  rho1 <- rep(0, K);  names(rho1) <- paste("I = ", inter.range, sep = "")
  #  rho2 <- rep(0, K);  names(rho2) <- paste("I = ", inter.range, sep = "")
  #  prop.r1 <- matrix(0, iter, K); colnames(prop.r1) <- paste("I = ", inter.range, sep = "")
  #  prop.r2 <- matrix(0, iter, K); colnames(prop.r2) <- paste("I = ", inter.range, sep = "")
  #}
  
  ## Start the loop across delta4 values.
  for (j in 1:k) {
    
    ## Define the seed - EC changed to set seed prior to starting loop (best practice)
    #set.seed(seed + j * 10000)
    
    ## Specify the value of the intercept in the collider model.
    ## This is done by calling a separate function
    d0 <- compute.d0(x.distr = x.distr, y.model = y.model, y.family=y.family,
                     beta0 = beta0, beta1 = beta1, beta2 = beta2, s.model = s.model, 
                     s.model.nl = s.model.nl, sel.prob = sel.prob,
                     delta1 = delta1, delta2 = delta2, delta3 = delta3,
                     d4 = inter.range[j], #latent.sigma = latent.sigma, 
                     n.d0 = n.d0)
    
    
    ## Store results.
    #if (s.model == "double") {
    #  delta0[j] <- d0[1]
    #  rho1[j] <- d0[2]
    #  rho2[j] <- d0[3]
    #} else {
      delta0[j] <- d0
      delta4[j] <- inter.range[j]
    #}
    
    ## Run the simulation.
    for (i in 1:iter) {
      
      ## ----- Simulate data. ------ ##
      
      ## Simulate exposure data.
      if (x.distr == "binary") {
        X <- rbinom(n, 1, 0.3) 
      } else if (x.distr == "unif") {
        X <- runif(n, 0.1, 0.9)
      } else if (x.distr == "normal") {
        X <- rnorm(n, 0.5, 1)
      }
     
      ## Simulate outcome data
      if (y.model == "linear") {
        linpred <- beta0 + beta1 * X 
      } else if (y.model == "quadratic") {
        linpred <- beta0 + beta1 * X + beta2 * X^2 
      }
      if (y.family == "gaussian") {
          Y <- linpred + rnorm(n, 0, 0.25)
        } else if (y.family == "binomial") {
          Y.pr <- expit(linpred)
          Y <- rbinom(n, 1, prob = Y.pr)
      } else if (y.family == "poisson") {
        Y.lambda <- exp(linpred)
        Y <- rpois(n, lambda = Y.lambda)
      }
      
      ## Simulate data for the collider.
      if (s.model == "logadd") {
        
        ## For a log-additive collider.
        if (s.model.nl == "None"){
          S.pr <- exp(d0 + delta1 * X + delta2 * Y + delta3 * X * Y)
        } else if (s.model.nl == "X^2"){
          S.pr <- exp(d0 + delta1 * X + delta2 * Y + delta3 * X * Y + inter.range[j] * X^2)
        } else if (s.model.nl == "X^2Y"){
          S.pr <- exp(d0 + delta1 * X + delta2 * Y + delta3 * X * Y + inter.range[j] * X^2 * Y)
        } else if (s.model.nl == "Y^2"){
          S.pr <- exp(d0 + delta1 * X + delta2 * Y + delta3 * X * Y + inter.range[j] * Y^2)
        }
        
        # If any obs have a calculated selection probability > 1, skip to the next iteration
        # Note that including these obs, or truncating their selection probabilities to 1, 
        # will give misleading (biased) results
        cutoff[i, j] <- sum(S.pr > 1)
        if (cutoff[i,j]>0) next
        
        #Generate selection indicator based on the calculated selection probability
        S <- rbinom(n, 1, prob = S.pr)
        
      } 
      
      # Not considered in this study
      #else if (s.model == "logistic") {
        
        ## Second, for a logistic collider.
        #S.pr <- expit(d0 + delta1 * X + delta2 * Y + inter.range[j] * X * Y)
        #S <- rbinom(n, 1, prob = S.pr)
        
      #} else if (s.model == "probit") {
        
        ## Third, for a probit collider.
       # S.pr <- d0 + delta1 * X + delta2 * Y + inter.range[j] * X * Y + rnorm(n, 0, latent.sigma)
        #S <- as.numeric(S.pr > 0)
        
      #} else if (s.model == "double") {
        
        ## And fourth, for a "double-threshold" collider.
       # S.pr <- d0[1] + delta1 * X + delta2 * Y + inter.range[j] * X * Y + rnorm(n, 0, latent.sigma)
      #  S <- rep(0, n)
       # S[S.pr < d0[2]] <- 1
      #  S[S.pr > d0[3]] <- 1
        
      #}
      
      ## ----- Fit the models. ----- ##
      
      ## Store full data and selected data in data frames
      data.full <- data.frame(Y=Y,X=X)
      data.sel <- data.frame(Y=Y[S == 1],X=X[S == 1])
      
      ## Fit the X-Y model for inference on selected data
      if (xy.fit == "linear"){
        if (y.family == "gaussian") {
          fit1 <- summary(lm(Y ~ X, data=data.sel))$coefficients
        } else if (y.family == "binomial") {
          fit1 <- summary(glm(Y ~ X, data=data.sel, family = binomial(link = "logit")))$coefficients
        } else if (y.family == "poisson") {
          fit1 <- summary(glm(Y ~ X, data=data.sel, family = poisson(link = "log")))$coefficients
        }
      } else if (xy.fit == "quadratic") {
        if (y.family == "gaussian"){
          fit1 <- summary(lm(Y ~ X + I(X^2), data=data.sel))$coefficients
        } else if (y.family == "binomial") {
          fit1 <- summary(glm(Y ~ X + I(X^2), data=data.sel, family = binomial(link = "logit")))$coefficients
        } else if (y.family == "poisson") {
          fit1 <- summary(glm(Y ~ X + I(X^2), data=data.sel, family = poisson(link = "log")))$coefficients
        }
      } else if (xy.fit == "cubic") {
        if (y.family == "gaussian"){
          fit1 <- summary(lm(Y ~ X + I(X^2) + I(X^3), data=data.sel))$coefficients
        } else if (y.family == "binomial") {
          fit1 <- summary(glm(Y ~ X + I(X^2) + I(X^3), data=data.sel, family = binomial(link = "logit")))$coefficients
        } else if (y.family == "poisson") {
          fit1 <- summary(glm(Y ~ X + I(X^2) + I(X^3), data=data.sel, family = poisson(link = "log")))$coefficients
        }
      }
      xy.lincoef[i, j] <- fit1[2, 1]
      xy.lincoef.se[i, j] <- fit1[2, 2]
      xy.int[i, j] <- fit1[1, 1]
      xy.int.se[i, j] <- fit1[1, 2]
      
      if (xy.fit == "quadratic" | xy.fit == "cubic"){
        xy.quadcoef[i, j] <- fit1[3, 1]
        xy.quadcoef.se[i, j] <- fit1[3, 2]  
      }
      
      if (xy.fit == "cubic"){
        xy.cubcoef[i, j] <- fit1[4, 1]
        xy.cubcoef.se[i, j] <- fit1[4, 2]  
      }
      
      # Not performed #
      ## Fit a log-probability model for selection.
      #fit2 <- summary(glm(S ~ X + Y + X * Y, family = poisson(link = "log")))$coefficients
      #m1.delta[I, j] <- fit2[4, 1]
      #m1.se[I, j] <- fit2[4, 2]
      
      ## Fit a logistic model for selection.
      #fit2 <- summary(glm(S ~ X + Y + X * Y, family = binomial(link = "logit")))$coefficients
      #m2.delta[I, j] <- fit2[4, 1]
      #m2.se[I, j] <- fit2[4, 2]
      
      ## Fit a probit model for selection.
      #fit2 <- summary(glm(S ~ X + Y + X * Y, family = binomial(link = "probit")))$coefficients
      #m3.delta[I, j] <- fit2[4, 1]
      #m3.se[I, j] <- fit2[4, 2]
      
      ## Fit model using fractional polynomial selection, allowing up to 3 terms (df=6)
      ## Fit to both selected and full data, keeping at least one X term
      fit.fp.sel <- mfp2(Y ~ fp(X, df=6), data=data.sel, family = y.family,
                           verbose=FALSE, keep="X")
      fit.fp.full <- mfp2(Y ~ fp(X, df=6), data=data.full, family = y.family,
                         verbose=FALSE, keep="X")
      ## Extract selected fractional powers of X
      fp.powers.sel.X1[i, j] <- fit.fp.sel$fp_powers$X[1]
      if (length(fit.fp.sel$fp_powers$X)>=2) {fp.powers.sel.X2[i, j] <- fit.fp.sel$fp_powers$X[2]}
      if (length(fit.fp.sel$fp_powers$X)==3) {fp.powers.sel.X3[i, j] <- fit.fp.sel$fp_powers$X[3]}
      
      fp.powers.full.X1[i, j] <- fit.fp.full$fp_powers$X[1]
      if (length(fit.fp.full$fp_powers$X)>=2) {fp.powers.full.X2[i, j] <- fit.fp.full$fp_powers$X[2]}
      if (length(fit.fp.full$fp_powers$X)==3) {fp.powers.full.X3[i, j] <- fit.fp.full$fp_powers$X[3]}
      
      ## ----- Additional diagnostics. ----- ##
        
      ## Store additional diagnostics.
      min.s.pr[i, j] <- min(S.pr)
      max.s.pr[i, j] <- max(S.pr)
      prop.selected[i, j] <- mean(S)
      mean.y[i, j] <- mean(Y)
      
      ## Store further diagnostics for "double-threshold" models.
      #if (s.model == "double") {
      #  prop.r1[I, j] <- sum(S.pr < d0[2])
      #  prop.r2[I, j] <- sum(S.pr > d0[3])
      #}
      
    }
    
    ## Progress report.
    print(paste("Interaction value ", j, " done.", sep = ""))
    
  }
  
  ## Put everything together in a list.
  res.list <- list("lincoef.est" = xy.lincoef, "lincoef.se" = xy.lincoef.se, 
                   "quadcoef.est" = NA, "quadcoef.se" = NA,
                   "cubcoef.est" = NA, "cubcoef.se" = NA,
                   "int.est" = xy.int, "int.se" = xy.int.se, 
                   "fp.powers.sel.X1" = fp.powers.sel.X1,
                   "fp.powers.sel.X2" = fp.powers.sel.X2,
                   "fp.powers.sel.X3" = fp.powers.sel.X3,
                   "fp.powers.full.X1" = fp.powers.full.X1,
                   "fp.powers.full.X2" = fp.powers.full.X2,
                   "fp.powers.full.X3" = fp.powers.full.X3,
                   #"logadd.est" = m1.delta, "logadd.se" = m1.se, "logit.est" = m2.delta, "logit.se" = m2.se, 
                   #"probit.est" = m3.delta, "probit.se" = m3.se, 
                   "cutoff" = cutoff, "selection.prob" = prop.selected, 
                   "mean.y" = mean.y, "min.sel.probs" = min.s.pr, "max.sel.probs" = max.s.pr, 
                   "delta0" = delta0, "delta4" = delta4)
  #if (s.model == "logadd") res.list$cutoff <- cutoff
  #if (s.model == "double") res.list$rho1 <- rho1
  #if (s.model == "double") res.list$rho2 <- rho2
  #if (s.model == "double") res.list$prop.rho1 <- prop.r1
  #if (s.model == "double") res.list$prop.rho2 <- prop.r2
  if (xy.fit == "quadratic" | xy.fit == "cubic"){
    res.list$quadcoef.est <- xy.quadcoef
    res.list$quadcoef.se <- xy.quadcoef.se
  }
  if (xy.fit == "cubic"){
    res.list$cubcoef.est <- xy.cubcoef
    res.list$cubcoef.se <- xy.cubcoef.se
  }
  #if (length(fit.fp.sel$fp_powers$X)>=2) {res.list$fp.powers.sel.X2 <- fp.powers.sel.X2}
  #if (length(fit.fp.sel$fp_powers$X)==3) {res.list$fp.powers.sel.X3 <- fp.powers.sel.X3}
  #if (length(fit.fp.full$fp_powers$X)>=2) {res.list$fp.powers.full.X2 <- fp.powers.full.X2}
  #if (length(fit.fp.full$fp_powers$X)==3) {res.list$fp.powers.full.X3 <- fp.powers.full.X3}
  
  ## Goodbye.
  return(res.list)
  
}





##################################################

## ---------- DEFAULT (50% SELECTED) ---------- ##

#Start with some basic checks before running the full sim study
# 1. Check that range of Y is generally between -1 and 1 
# Or ideally between 0 and 1 - this will be make it easier 
# to construct a selection model such that all selection probabilities are < 1
# If necessary, change beta0 and beta1 to ensure this
n <- 1e7
beta0 <- 0
beta1 <- 0.5
X <- runif(n, 0, 0.5)
Y <- beta0 + beta1 * X + rnorm(n, 0, 0.15)
summary(Y)

# 2. Run 1 sim with a very large sample size and check selection probabilities 
# are all < 1 and results are as expected i.e. no non-linear r'ship in settings
# where this is not expected, plus bias of X coef of expected magnitude 

## This can be checked by printing 'cutoff', 'lincoef' and 'quadcoef' for each run.
## May need to alter delta1-3 to ensure selection probs < 1

# Using default settings of function unless specified:
#(iter = 1, n = 1e7, x.distr = "unif", y.model = "linear", 
#xy.fit = "linear",
#beta0 = 0, beta1 = 0.5, s.model = "logadd", s.model.nl = "X^2",
#sel.prob = 0.5, 
#delta1 = 0.25, delta2 = 0.25, delta3 = -0.5, 
#inter.range = seq(-0.5, 0.5, length = 11), n.d0 = 1e6)

## We use a sample size of 1e7 to approximate an infinite
## sample. We run each asymptotic experiment once.

## Each run takes about 1 min on EC's laptop (Processor:11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz, 16GB RAM).

## Set seed
set.seed(684632)

ASY0.1 <- interaction.sims(xy.fit="linear", s.model.nl = "None", inter.range = 0,
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
ASY0.2 <- interaction.sims(xy.fit="quadratic", s.model.nl = "None", inter.range = 0,
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
#As expected, no evidence of quad relationship
ASY1.1 <- interaction.sims(xy.fit="linear", s.model.nl = "X^2",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
ASY1.2 <- interaction.sims(xy.fit="quadratic", s.model.nl = "X^2",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
#As expected, no evidence of quad relationship
ASY2.1 <- interaction.sims(xy.fit="linear", s.model.nl = "X^2Y",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
ASY2.2 <- interaction.sims(xy.fit="quadratic", s.model.nl = "X^2Y",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
#Evidence of quad relationship, which increases in (absolute) magnitude 
#with magnitude of delta4
ASY3.1 <- interaction.sims(xy.fit="linear", s.model.nl = "Y^2",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)
#Sims with pos values of delta4 did not run - because these give some sel probs > 1
ASY3.2 <- interaction.sims(xy.fit="quadratic", s.model.nl = "Y^2",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5)

# Some evidence of a quad relationship for some values of delta4
# Try for bigger range of delta4
ASY3.2.1 <- interaction.sims(xy.fit="quadratic", s.model.nl = "Y^2",
                           delta1 = 0.25, delta2 = 0.25, delta3 = -0.5,
                           inter.range = seq(-1, 1, length = 11))
# Largest delta4 values don't run
# Perhaps light can be shed on the reason for this by calculating bias algebraically?

#Try reducing prop selected to try to allow pos values of delta4
ASY3.2.2 <- interaction.sims(xy.fit="quadratic", s.model.nl = "Y^2",
                             delta1 = 0.25, delta2 = 0.25, delta3 = -0.5,
                             sel.prob = 0.2,
                             inter.range = seq(-1, 1, length = 11))
summary(ASY3.2.2$cutoff)
summary(ASY3.2.2$lincoef.est)
summary(ASY3.2.2$lincoef.se)
summary(ASY3.2.2$quadcoef.est)
summary(ASY3.2.2$quadcoef.se)
## ---------- SAVE RESULTS ---------- ##

## Don't save asymptotic runs - purely exploratory
#save(ASY0.1, ASY0.2, ASY1.1, ASY1.2,ASY2.1, ASY2.2, ASY3.1, ASY3.2,
#     file = "Asymptotic_results.RData")

## --------- SAMPLE SIZE CALCULATION -----------##
# Calculate 'worst-case' MCSE for coverage and Type 1 error
# Note this applies for both sample sizes
sqrt(0.5*0.5/1000)
# mcse=0.01581

# Alternatively, if aim for MCSE of 0.02, need the following nsim:
0.5*0.5/(0.02)^2
#625

0.5*0.5/(0.01)^2
#2500 for MCSE of 0.01

# Check for one setting for sample size calculation:
## Set seed
set.seed(16346)

# For n=1000
ssize_calc <- interaction.sims(iter = 1000, n = 1000, inter.range=-0.5) # all other settings as default
summary(ssize_calc$lincoef.est)
summary(ssize_calc$lincoef.se)

# Calculate MCSE of bias of each coefficient estimate based on 1000 sims, using i=-0.5 (largest range)
sqrt(var(ssize_calc$lincoef.est[,1])/1000)
# mcse=0.001531
#sqrt(var(ssize_calc$int.est)/1000)

# Calculate MCSE of average model-based SE of each coefficient estimate based on 1000 sims
sqrt(var((ssize_calc$lincoef.se[,1])^2)/(4*1000*mean((ssize_calc$lincoef.se[,1])^2)))
# mcse=6.464e-05
#sqrt(var((ssize_calc$int.se)^2)/(4*1000*mean((ssize_calc$int.se)^2)))


# For n=1000000 (10 iterations)

## Set seed
set.seed(53)

ssize_calc3 <- interaction.sims(iter = 10, inter.range=-0.5) # all other settings as default
summary(ssize_calc3$lincoef.est)
summary(ssize_calc3$lincoef.se)

# Calculate MCSE of bias of each coefficient estimate based on 10 sims, using i=-0.5 (largest range)
sqrt(var(ssize_calc3$lincoef.est[,1])/10)
# mcse=0.000661
#sqrt(var(ssize_calc$int.est)/1000)

# Calculate MCSE of average model-based SE of each coefficient estimate based on 10 sims
sqrt(var((ssize_calc3$lincoef.se[,1])^2)/(4*10*mean((ssize_calc3$lincoef.se[,1])^2)))
# mcse=4.809e-07
#sqrt(var((ssize_calc$int.se)^2)/(4*1000*mean((ssize_calc$int.se)^2)))

