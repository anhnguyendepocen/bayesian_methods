############### SESYNC Research Support: bayesian_Drosophila_PERA_analysis ########## 
## Exploration of solutions for glmer with bayesian approach.
## 
## DATE CREATED: 04/24/2018
## DATE MODIFIED: 04/24/2018
## AUTHORS: Benoit Parmentier 
## PROJECT: bayesian_Drosophila_PERA_analysis/
## ISSUE: 
## TO DO:
##
## COMMIT: testing simple regression with bayesion framework
##
## Links to investigate:
#http://julianfaraway.github.io/brinla/

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(raster)
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)
library(lubridate)
library(dplyr)
library(rowr)                                # Contains cbind.fill
library(car)
library(sf)
library(brinla)
library(INLA)
library(MCMCglmm)

###### Functions used in this script and sourced from other files

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

### Other functions ####

#function_processing_data <- ".R" #PARAM 1
#script_path <- "/nfs/bparmentier-data/Data/projects/bayesian_Drosophila_PERA_analysis/scripts" #path to script #PARAM 
#source(file.path(script_path,function_processing_data)) #source all functions used in this script 1.


############################################################################
#####  Parameters and argument set up ###########

out_suffix <- "bayesian_exploration_04242018" #output suffix for the files and ouptut folder #param 12

in_dir <- "/nfs/bparmentier-data/Data/projects/bayesian_Drosophila_PERA_analysis/data"
out_dir <- "/nfs/bparmentier-data/Data/projects/bayesian_Drosophila_PERA_analysis/data/outputs"

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
create_out_dir_param=TRUE #PARAM9

############## START SCRIPT ############################

######### PART 0: Set up the output dir ################

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

###################### PART 1: ###########

# Quick Start

data(hubble, package = "brinla")
plot(y ~ x, xlab = "Distance(Mpc)", ylab = "Velocity(km/s)",
     data = hubble)
lmod <- lm(y ~ x - 1, data = hubble) #no intercept hence we use -1
coef(lmod)
summary(lmod)

### check confidence interval
bci <- confint(lmod)
print(bci)

hubtoage <- function(x) 3.09e+19/(x * 60^2 * 24 * 365.25 * 1e+09)
hubtoage(coef(lmod))
(bci <- confint(lmod))

hubtoage(bci)

#library(INLA)

## Let's use an uninformative prior
# set the precision to a very small number. The precision is inverse of variance.
imod <- inla(y ~ x - 1, 
             family = "gaussian",
             control.fixed = list(prec = 1e-09), # set control on prior
             data = hubble)

summary(imod)
(ibci <- imod$summary.fixed)

plot(imod$marginals.fixed$x, type = "l", xlab = "beta",
     ylab = "density", xlim = c(60, 100))

abline(v = ibci[c(3, 5)], lty = 2)

hubtoage(ibci[c(1,3,4,5,6)])
ageden <- inla.tmarginal(hubtoage, imod$marginals.fixed$x)
plot(ageden, type = "l", xlab = "Age in billions of years",
     ylab = "density")
abline(v = hubtoage(ibci[c(3, 5)]), lty = 2)

###### Use a different prior: with prior information

#Before Hubble Space Telescpe age was estimated between 10-20 billion years.
#Transform tis guess using the hubtoage
hubtoage(c(10, 15, 20))

#this suggest a mean prior or 65

var_val <- hubtoage(c(10, 15, 20))[1] - hubtoage(c(10, 15, 20))[3]

variance_prior <- round(var_val/4) # 95% interval goes from -1.96 to 1.96 or about 4 
print(variance_prior)

mean_prior <- 65
sd_prior <- 12
var_prior <- 12^2
precision_prior <- 1/(sd_prior^2)

imod2 <- inla(y ~ x - 1, 
             family = "gaussian",
             control.fixed = list(mean = 65, prec = 1/(12^2)), 
             data = hubble)

(ibci <- imod2$summary.fixed)
plot(imod$summary.fixed)

plot(imod2$marginals.fixed$x,
     main=paste0("Distribution of beta with informed prior(",c(mean_prior,var_prior),")"),
     type = "l", xlab = "beta",
     ylab = "density", xlim = c(60, 100))

########################
##### Use MCMC to sample the Posterior:

trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

# create independent x-values 
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

plot(x,y, main="Test Data")

#https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
  
#a linear relationship y = a*x + b 

likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
    
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

# Example: plot the likelihood profile of the slopmdAnn2016
<- function(x){return(likelihood(c(x, trueB, trueSd)))}
slopelikelihoods <- lapply(seq(3, 7, by=.05), slopevalues )
plot (seq(3, 7, by=.05), slopelikelihoods , type="l", 
      xlab = "values of slope parameter a", ylab = "Log likelihood")

likelihood_fun = dnorm(y, mean = pred, sd = sd, log = T)

#https://stats.idre.ucla.edu/r/faq/normal/
  
https://www4.stat.ncsu.edu/~reich/ST740/code/
  
# fplot <- function(x, params, pdf) {
#   
#   params <- do.call(expand.grid, params)
#     
#   k <- nrow(params)
#     
#   nc <- ceiling(sqrt(k))
#   nr <- ceiling(k / nc)
#   stopifnot(k <= nc * nr)
#     
#   labeller <- function(params) {
#     style <- paste(paste0(names(params), "=%2.2f"), collapse = ", ")
#     do.call(sprintf, c(fmt = style, params))
#   }
#     
#   par(mfrow = c(nr, nc))
#   lapply(1:k, function(i) {
#     y <- do.call(pdf, c(list(x), params[i, ]))
#     plot(x = x, y = y, ylab = "Density", xlab = labeller(params[i, ]), type = "l")
#   })
#   return(invisible(NULL))
# }
# 
# ### Normal distribution
# 
# fplot(
#   x = seq(from = -50, to = 50, by = .1),
#   params = list(mean = c(-10, -5, 5, 10), sd = c(.1, 1, 10, 100)),
#   pdf = dnorm)
# 
# 
# 
# 
# require(foreign)
# 
# require(MCMCglmm)
# 
# dat <- read.spss("http://statistics.ats.ucla.edu/stat/data/hsbdemo.sav", 
#                  to.data.frame=TRUE)
# 
# m1 <- MCMCglmm(write ~ read, data = dat,
#                prior = list(
#                  B = list(mu = c(0, 0), V = diag(2) * 1e5),
#                  R = list(V = 1, nu = .002)),
#                nitt = 11000, thin = 1, burnin = 1000, verbose=FALSE)
# 
# autocorr(m1$VCV)

############################## End of script #############################