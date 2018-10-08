#### Bayesian stat


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