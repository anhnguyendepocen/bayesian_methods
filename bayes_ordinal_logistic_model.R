############### Climate perception project ########## 
##
##
#Data used in this project comes from three sources: (1) A survey of NRCS and FSA employees conducted by the USDA Climate Hubs in 2016/2017, (2) Crop idemnity payments made by FSA to farmers for weather-related 
#crop loss between 2013-2016, and (3) mean and standard deviations of drought data over 1, 3, 5, 10, and 15 year periods.
##
## This is an example script using mixed models to test differences in time series
## for different locations/subjects. Autocorrelation structures is taken into account
## in the methods.
##
## DATE CREATED: 09/27/2018
## DATE MODIFIED: 10/08/2018
## AUTHORS: Rachel Schattman, Benoit Parmentier  
## Version: 1
## PROJECT: Climate Percecption
## ISSUE: 
## TO DO:
##
## COMMIT: modifying code
##

## Very good reference:
#http://rpsychologist.com/r-guide-longitudinal-lme-lmer

###################################################
#

###### Library used

## ------------------------------------------------------------------------
library(MASS)
library(lme4)
library(rstanarm)
library("bayesplot")
library("ggplot2")
library("loo")
library("parallel")

##########

###### Functions used in this script and sourced from other files

run_model_assessment <- function(mod,model_type="bayes_stan",k_threshold = 0.7){
  
  if (model_type!="bayes_stan") {
    stop("Model type not implemented")
  }
  
  if(model_type=="bayes_stan"){
    debug(loo)
    mod$formula
    loo_obj <- try(loo(mod, k_threshold = k_threshold))
  }
  
  return(loo_obj)
}

run_model_ordinal_logistic <- function(model_formula,model_type="bayes_stan",data,prior = NULL, prior_counts = dirichlet(1),
                                       shape = NULL,chains = 4, num_cores = NUL, seed_val = 1234, iter_val = 200){
  
  if (model_type!="bayes_stan") {
    stop("Model type not implemented")
  }
  
  if(model_type=="bayes_stan"){
    mod <- try(stan_polr(formula=model_formula,
                         data = data_subset, 
                         prior = prior, 
                         prior_counts = dirichlet(1),
                         shape = shape,
                         chains = chains, 
                         cores = num_cores, 
                         seed = seed_val, 
                         iter = iter_val))
  }
  
  return(mod)
}

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

################### Start script ###################

#Benoit setup
script_path <- "/nfs/bparmentier-data/Data/projects/soilsesfeedback-data/scripts"
#mosaicing_functions <- "weighted_mosaicing_functions_07252018.R"
#source(file.path(script_path,mosaicing_functions))

#########cd ###################################################################
#####  Parameters and argument set up ########### 

#ARGS 1
in_dir <- "/nfs/bparmentier-data/Data/projects/soilsesfeedback-data/data"
#ARGS 2
out_dir <- "/nfs/bparmentier-data/Data/projects/soilsesfeedback-data/outputs"
#ARGS 3:
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 7
out_suffix <-"_10052018" #output suffix for the files and ouptut folder
#ARGS 8
num_cores <- 2 # number of cores

in_filename <- "NRCS_FSAMergeDataset_w_PDSI2_7_28_18.csv"
model_type <- "bayes_stan"
y_var_name <- "Concern_DryDrought"


################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

options(scipen=999)

#set up the working directory
#Create output directory

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

#######################################
### PART 1: Read in DATA #######

#data_df <- read.table(file.path(in_dir,in_filename),
#                      sep=",",
#                      header=T)

#rm(list=ls())            # clear

dataDR <- read.csv(file.path(in_dir,in_filename), #"/nfs/bparmentier-data/Data/projects/soilsesfeedback-data/data/NRCS_FSAMergeDataset_w_PDSI2_7_28_18.csv", 
                   header = TRUE)


dataDR$y_var <- dataDR[[y_var_name]]

dataDR$stdiv <- factor(dataDR$stdiv)
dataDR$Agency <- factor(dataDR$Agency)
dataDR$y_var <- factor(dataDR$y_var)

#### Setting up models var inputs
x_var_clean <- c("PercLossDrought", 
                 "stdiv", 
                 "Agency", 
                 "PDSI_MEAN_2016",
                 "PDSI_MEAN_2014",
                 "PDSI_MEAN_2012",
                 "PDSI_MEAN_2007",
                 "PDSI_MEAN_2002",
                 "PDSI_STD_2016",
                 "PDSI_STD_2014",
                 "PDSI_STD_2012",
                 "PDSI_STD_2007",
                 "PDSI_STD_2002")

y_var_clean <- y_var_name

variables_used <- c(y_var_clean, x_var_clean)

##subset dataset for variables

data_subset <- dataDR [,variables_used]
data_subset$y_var <- factor(data_subset[[y_var_name]])

data_subset <- na.omit(data_subset)

data_subset$y_var <- factor(data_subset[[y_var_name]])

## ------------------------------------------------------------------------

mod_noPDSI <- "y_var ~ PercLossDrought + stdiv"
mod_mean2016 <- "y_var ~ PercLossDrought + PDSI_MEAN_2016"
mod_mean2014 <- "y_var ~ PercLossDrought + PDSI_MEAN_2014"
mod_mean2012 <- "y_var ~ PercLossDrought + PDSI_MEAN_2012"
mod_mean2007 <- "y_var ~ PercLossDrought + PDSI_MEAN_2007"
mod_mean2002 <- "y_var ~ PercLossDrought + PDSI_MEAN_2002"
mod_STD2016 <- "y_var ~ PercLossDrought + PDSI_STD_2016"
mod_STD2014 <- "y_var ~ PercLossDrought + PDSI_STD_2014"
mod_STD2012 <- "y_var ~ PercLossDrought + PDSI_STD_2012"
mod_STD2007 <- "y_var ~ PercLossDrought + PDSI_STD_2007"
mod_STD2002 <-  "y_var ~ PercLossDrought + PDSI_STD_2002"
               
list_model_formulas <- list(mod_noPDSI,mod_mean2016,mod_mean2014,mod_mean2012,mod_mean2007,mod_mean2002,
                    mod_STD2016,mod_STD2014,mod_STD2012,mod_STD2007,mod_STD2002)

############ PART 2: Run model with option for bayesian ordinal logistic


mod2 <- run_model_ordinal_logistic(list_model_formulas[[2]],
                           data = data_subset, 
                           prior = NULL,
                           prior_counts = dirichlet(1),
                           shape = NULL,
                           chains = 4, 
                           num_cores = 4, 
                           seed_val = 1234, 
                           iter_val = 200)
## does not work
#list_mod <- mclapply(list_model_formulas,
#                     FUN=run_model_ordinal_logistic,
#                     data = data_subset, 
#                     prior = NULL,
#                     prior_counts = dirichlet(1),
#                     shape = NULL,
#                     chains = 4, 
#                     num_cores = 4, 
#                     seed_val = 1234, 
#                     iter_val = 200,
#                    mc.preschedule = F,
#                    mc.cores = 1)
         
list_mod <- lapply(list_model_formulas[1:3],
                     FUN=run_model_ordinal_logistic,
                     data = data_subset, 
                     prior = NULL,
                     prior_counts = dirichlet(1),
                     shape = NULL,
                     chains = 4, 
                     num_cores = 4, 
                     seed_val = 1234, 
                     iter_val = 200)

list_mod[[3]]
#save(mod,file= paste("C:\\Users\\rschattman\\Documents\\Research\\climate-drivers\\model",i,"output.rdata", sep ="")) # This save would be useful if you wanted to save each of the 11 models as their own file
  
mod_outfilename <- paste0("list_mod_",out_suffix,".RData")
save(list_mod, 
     file = mod_outfilename)

############# PART 23: Model assessment 

#debug(run_model_assessment)

loo_mod2 <- run_model_assessment(mod2)
loo_mod <- mclapply(list_mod,
                    FUN=run_model_assessment,
                    k_threshold=0.7,
                    mc.preschedule = FALSE,
                    mc.cores=1)

compare_models(list_mod)


#### Collect information in table

## ------------------------------------------------------------------------
str(list_mod)
str(summary(list_mod))

## ------------------------------------------------------------------------

print(list_mod, digets = 3)
round(apply(rstan:: extract(list_mod$stanfit, pars = "drought") [[1]], 2, median), digits = 3)

loo(x, list_mod, cores = 1)
  #save_psis = FALSE, K_threshold = NULL)

table1 <- data.frame (list_mod = 1, 
                    intercept = 1, 
                    intercept.se = 1, 
                    slope = 1, 
                    slope.se = 1, 
                    r.squared = 1, 
                    p.value = 1)

for(i in 1:11){
  x<-rnorm(11)
  y<-rnorm(11)
  list_mod <- list_mod
  summary(list_mod)
}

for (i in 1:11){
  print(i)
}

table1[i,] <- c(i,
                summary(list_mod)[['coefficients']]['(Intercept)','Estimate'],
                     summary(list_mod)[['coefficients']]['(Intercept)','Std. Error'],
                summary(list_mod)[['coefficients']]['x','Estimate'],
                     summary(list_mod)[['coefficients']]['x','Std. Error'],
                summary(list_mod)[['r.squared']],
                     summary(list_mod)[['coefficients']]['x','Pr(>|t|)'])

write.csv(table1, file = 'table1.csv')


################################# End of script ######################################
