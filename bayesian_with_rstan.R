############################## Bayesian modeling with STAN #################### 
##
## This is an example script for Bayesian modeling using STAN to estimate parameters.
##
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 10/08/2018
## AUTHORS: Benoit Parmeniter
## Version: 1
## PROJECT: Training for Bayesian modeling
## ISSUE: 
## TO DO:
##
## COMMIT: clean up 
##

## Very good reference:
#http://rpsychologist.com/r-guide-longitudinal-lme-lmer

###################################################
#

#https://ourcodingclub.github.io/2018/04/17/stan-intro.html

in_dir <- "/nfs/bparmentier-data/Data/projects/soilsesfeedback-data/data/CC-Stan-intro-master"
lm1 <- lm(extent_north ~ year, data = seaice)
summary(lm1)
seaice <- read.csv(file.path(in_dir,"seaice.csv"), stringsAsFactors = F)
head(seaice)
plot(extent_north ~ year, pch = 20, data = seaice)

lm1 <- lm(extent_north ~ year, data = seaice)
summary(lm1)
abline(lm1, col = 2, lty = 2, lw = 3)
lm1 <- lm(y ~ x)
summary(lm1)

x <- I(seaice$year - 1978)
y <- seaice$extent_north
N <- length(seaice$year)

lm_alpha <- summary(lm1)$coeff[1]  # the intercept
lm_beta <- summary(lm1)$coeff[2]  # the slope
lm_sigma <- sigma(lm1)  # the residual error


library(rstan)
library(gdata)
library(bayesplot)

#Now let’s turn that into a dataframe for inputting into a Stan model. Data passed to Stan needs to be a list of named objects. The names given here need to match the variable names used in the models (see the model code below).

stan_data <- list(N = N, x = x, y = y)


write("// Stan model for simple linear regression

data {
 int < lower = 1 > N; // Sample size
 vector[N] x; // Predictor
 vector[N] y; // Outcome
}

parameters {
 real alpha; // Intercept
 real beta; // Slope (regression coefficients)
 real < lower = 0 > sigma; // Error SD
}

model {
 y ~ normal(alpha + x * beta , sigma);
}

generated quantities {
} // The posterior predictive distribution",

"stan_model1.stan")


#First, we should check our Stan model to make sure we wrote a file.

stanc("stan_model1.stan")
#Now let’s save that file path.

stan_model1 <- "stan_model1.stan"

fit <- stan(file = stan_model1, data = stan_data, warmup = 500, iter = 1000, chains = 4, cores = 2, thin = 1)

### now use rstanarm for the linear model


################# Part 2: Ordinal logistic model example with rstanarm


#https://www.analyticsvidhya.com/blog/2016/02/multinomial-ordinal-logistic-regression/
#https://stats.idre.ucla.edu/r/dae/ordinal-logistic-regression/
  
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)

dat <- read.dta("http://www.ats.ucla.edu/stat/data/ologit.dta")
dat <- read.dta("https://stats.idre.ucla.edu/stat/data/ologit.dta")
head(dat)

## one at a time, table apply, pared, and public
lapply(dat[, c("apply", "pared", "public")], table)

model_formula <- "apply ~ pared + public + gpa"

mod_polr <- polr(model_formula, data = dat, Hess=TRUE)
summary(mod_polr)

#data,prior = NULL 
#prior_counts = dirichlet(1)
#shape = NULL
#chains = 4
#num_cores = NUL
#seed_val = 100
#iter_val = 200

num_cores <-4
seed_val <- 100

mod_stan_polr <- try(stan_polr(formula=model_formula,
                     data = dat, 
                     prior = NULL, 
                     prior_counts = dirichlet(1),
                     shape = NULL,
                     chains = 4, 
                     cores = num_cores, 
                     seed = 100, 
                     iter = 200))


######################### End of script ##################################