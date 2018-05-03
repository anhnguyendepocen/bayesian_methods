####################################    Spatio-Temporal Analysis Malaria  #######################################
######################################   Malaria Indian State Analysis #######################################
#This script produces provides an analysis of Malaria case in India.       
#
#AUTHORS: Benoit Parmentier, Neeti Neeti                                             
#DATE CREATED: 04/25/2018 
#DATE MODIFIED: 05/02/2018
#Version: 1
#PROJECT: India research from Neeti            

#COMMENTS: - adding spatial data.frame

#TO DO:
# - run analyses with all states at once
#
#COMMIT: clean up and adding spatial boundaries
#
#################################################################################################

###Loading R library and packages        
### poisson regression
## 9 Chattisgarh, 15 Jharkhand, 30 Uttarakhand


###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(zoo)  #original times series function functionalies and objects/classes
library(xts)  #extension of time series objects and functionalities
library(forecast) #arima and other time series methods
library(lubridate) #date and time handling tools
library(colorRamps) #contains matlab.like color palette
library(rgeos) #spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(sphet) #spatial analyis, regression eg.contains spreg for gmm estimation
library(reshape2)
library(lme4)

###### Functions used in this script

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

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_04072017b.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir <- "/Users/neeti/Documents/neeti/TERIwork/aditi/"
#in_dir <- "/home/bparmentier/Google Drive/Data/India_Research/malaria_study_by_state/data"
#out_dir <- "/home/bparmentier/Google Drive/Data/India_Research/malaria_study_by_state/outputs"

in_dir <- "/home/benoit/Data/India_Research/malaria_study_by_state/data"
out_dir <- "/home/benoit/Data/India_Research/malaria_study_by_state/outputs"

#proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
#proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"malaria_india_04302018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#data_fname <- file.path(in_dir,"dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt")

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######

#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#data_tb <-read.table(data_fname,sep=",",header=T)

#data_fname <- file.path(in_dir,"malaria_incidence_yrs.csv")
data_fname <- file.path(in_dir,"malaria_incidence_yrs_94dmi_Mjo.csv")

#mal_inc <- read.table('malaria_incidence_yrs.csv', header=TRUE, sep=',')
#mal_inc <- read.table('malaria_incidence_yrs_94dmi_Mjo.csv', header=TRUE, sep=',')

mal_inc <- read.table(data_fname, header=TRUE, sep=',')
#rain_index <- read.table('rainfall_index.csv', header=TRUE, sep=',')
rain_fall <- read.table(file.path(in_dir,'rainfall_data.csv'), header=TRUE, sep=',')

View(mal_inc)
View(rain_fall)
dim(rain_fall)
View(rain_fall)
dim(mal_inc)

###### PART I: Reformat data and link to shapefile #################
### Need to combine both together.

names(mal_inc)[1] <- "year"
View(mal_inc)
names(mal_inc)[1:5]

n <- ncol(mal_inc)
n

test <- as.data.frame(t(mal_inc[6:n]))
names(mal_inc)[2:5]

names(test) <- 1994:2017
test$state <- rownames(test)
mdata <- melt(test, 
              variable.name =c("state"),
              value.names = c("state","year")
              #id=c("state","year")
              )
names(mdata)[2] <- "year"
View(mdata)
names(mdata)[3] <- "mal_inc" 
View(mdata)
#### join climate indices by dates
data_df <- merge(mdata,mal_inc[,1:5],by="year")

View(data_df)

#mod_glm_poisson <- glm(mal_inc ~ ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS + rain_fall[,j], 
#                       data=data_df,
#                       family=poisson(link=log))

data_df$year <- as.numeric(as.character(data_df$year))
class(data_df$state)
class(data_df$ONI_DJF)


mod_glm_poisson <- glm(mal_inc ~ state + year + ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS, 
                       data=data_df,
                       family=poisson(link=log))
mod_glm_poisson
summary(mod_glm_poisson)

plot(mal_inc ~ year,data=data_df)
xyplot(mal_inc ~ year | state,data=data_df)

plot(mal_inc ~ year, subset(data_df,state=="GOA"))
plot(mal_inc ~ year, subset(data_df,state=="CHHATTISGARH")) #ok problem in the data, there are zero
                                                            #when it should be NA

unique(data_df$state)
#?over dispersion?
#### Do spatial poisson

#Maybe mixed effect is needed?
mod_glm_poisson <- glm(mal_inc ~ state + year + ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS, 
                       data=data_df,
                       family=poisson(link=log))
summary(mod_glm_poisson)
### Mixed effect model
#https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

#https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
#https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/
#https://rpubs.com/wsundstrom/t_panel

mod_glmer_poisson <- glmer(mal_inc ~ year + ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS + (1|state), 
      data = data_df, family = poisson(link=log))

#mod_glmer_poisson <- glmer(mal_inc ~ year + ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS |state, 
#                           data = data_df, family = poisson(link=log))

## GLM mixed effect model
### here we model the slope and intercept per state for the year variable
#as random effect
### we also have the overall slope for year (fixed effect) and slopes for other var
mod_glmer_poisson <- glmer(mal_inc ~ year + (1+ year| state) + 
                             ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS , 
                           data = data_df, family = poisson(link=log))

summary(mod_glmer_poisson) #did not converge...

summary(mod_glmer_poisson)
# should we normalize by area?

#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)



############# Spatial model



##################################  END OF SCRIPT #####################################


#for (i in 6:35){
#  j <- i -4
#res <- glm(mal_inc1[,i] ~mal_inc1$ONI+mal_inc1$DMI+rain_index[,j], family=poisson(link=log))
#res <- glm(mal_inc1[,i] ~mal_inc1$ONI+mal_inc1$DMI, family=poisson(link=log))
#res <- glm(mal_inc1[,i] ~mal_inc1$ONI_DJF+mal_inc1$DMI_ASO+mal_inc1$MJO_DJFM+mal_inc1$MJO_JJAS + rain_fall[,j], family=poisson(link=log))
#res <- glm(mal_inc1[,i] ~mal_inc1$ONI_DJF+mal_inc1$DMI_ASO+ rain_fall[,j], family=poisson(link=log))

#chaatisg_inc <- mal_inc1[,9]
#chaatisg_inc <- chaatisg_inc[-c(1:6)]

#res <- glm(chaatisg_inc ~mal_inc1[(7:24),2]+mal_inc1[(7:24),3]+rain_fall[(7:24),6], family=poisson(link=log))

#res_sum <- summary(res)
# coef_res <- res_sum$coefficients[,1]
# p_res <- res_sum$coefficients[,4]
# stderror_coef <- res_sum$coefficients[,2]
# aic_model <- res_sum$aic
# aic_total[(i-5),] <- aic_model
# lmp_arr[(i-5),] <- p_res
# lmcoef_arr[(i-5),] <- exp(coef_res)
# lmstd_arr[(i-5),] <- stderror_coef
# lmcilow_arr[(i-5),] <- exp(coef_res - 1.96*stderror_coef)
# lmciabove_arr[(i-5),] <- exp(coef_res + 1.96*stderror_coef)
# #out <- capture.output(summary(res))
# 
# #cat("poisson regression analysis", out, file="summary_poissonregression.txt", sep="n", append=TRUE)
# }
#write.table(lmp_arr, 'pvalue_malpoisson_oni_iod94_rf_mjo.txt')
#write.table(lmcoef_arr, 'coefvalue_malpoisson_oni_iod94_rf_mjo.txt')
#write.table(lmstd_arr, 'stderror_malpoisson_oni_iod94_rf_mjo.txt')
#write.table(lmcilow_arr, 'CIlow_malpoisson_oni_iod94_rf_mjo.txt')
#write.table(lmciabove_arr, 'CIabove_malpoisson_oni_iod94_rf_mjo.txt')

#write.table(aic_total, 'aic_malpoisson_oni_iod_rf_mjo.txt')
