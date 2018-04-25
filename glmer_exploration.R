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
#script_path <- "/nfs/bparmentier-data/Data/projects/bayesian_Drosophila_PERA_analysis/R_scripts" #path to script #PARAM 
#source(file.path(script_path,function_processing_data)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

out_suffix <- "connectivity_example_04242018" #output suffix for the files and ouptut folder #param 12

in_dir <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden"
out_dir <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/outputs"
in_dir_grass <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden" 
#in_dir_grass <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/grass_data_urban_garden"

# background reading:
# https://grass.osgeo.org/grass72/manuals/grass_database.html

gisBase <- '/usr/lib/grass72'
#gisDbase <- '/nfs/urbangi-data/grassdata'
gisDbase <- in_dir_grass #should be the same as in_dir

#location <- 'DEM_LiDAR_1ft_2010_Improved_NYC_int'
location <- 'NYC_example'
location <- 'connectivy_example'

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"ny_example_04232018" #output suffix for the files and ouptu folder #PARAM 8
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







#https://stats.idre.ucla.edu/r/faq/normal/
  
  
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
