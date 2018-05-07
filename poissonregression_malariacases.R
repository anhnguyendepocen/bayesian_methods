####################################    Spatio-Temporal Analysis Malaria  #######################################
######################################   Malaria Indian State Analysis #######################################
#This script produces provides an analysis of Malaria case in India.       
#
#AUTHORS: Benoit Parmentier, Neeti Neeti                                             
#DATE CREATED: 04/25/2018 
#DATE MODIFIED: 05/04/2018
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
library(sf)
library(INLA)

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
#View(mal_inc)
#View(rain_fall)
dim(rain_fall)
dim(mal_inc)

lf_admin_units_fname <- list.files(path=in_dir,
                                   pattern="IND_adm.*.shp",
                                   full.names = T)

region_sf <- st_read(lf_admin_units_fname[2])
dim(region_sf)
head(region_sf)
View(region_sf)
plot(region_sf)

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


mod_glm_poisson <- glm(mal_inc ~ state + year + 
                         ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS, 
                       data=data_df,
                       family=poisson(link=log))
mod_glm_poisson
summary(mod_glm_poisson)

plot(mal_inc ~ year,data=data_df)
xyplot(mal_inc ~ year | state,data=data_df)

state_selected <- "GOA"
plot(mal_inc ~ year, subset(data_df,state==state_selected),main=state_selected)
state_selected <- "CHHATTISGARH"
plot(mal_inc ~ year, subset(data_df,state==state_selected),main=state_selected)
 #ok problem in the data, there are zero
                                                            #when it should be NA
names(data_df$state)

unique(data_df$state)
#?over dispersion?
#### Do spatial poisson with Mixed effect model to take into account malaria incidence by year and state


### Mixed effect model
#https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

#https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
#https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/
#https://rpubs.com/wsundstrom/t_panel

### Intercept model!!! not what we want
mod_glmer_poisson <- glmer(mal_inc ~ year + 
                             ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS + (1|state), 
      data = data_df, family = poisson(link=log))
summary(mod_glmer_poisson)

ranef(mod_glmer_poisson)$year
ranef(mod_glmer_poisson)

#mod_glmer_poisson <- glmer(mal_inc ~ year + ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS |state, 
#                           data = data_df, family = poisson(link=log))

## GLM mixed effect model
### here we model the slope and intercept per state for the year variable
#as random effect
### we also have the overall slope for year (fixed effect) and slopes for other var
mod_glmer_poisson_year <- glmer(mal_inc ~ year + (1+ year| state) + 
                             ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS , 
                           data = data_df, family = poisson(link=log))
#refit(mod_glmer_poisson_year)

summary(mod_glmer_poisson_year) #did not converge...

# should we normalize by area?

############# INLA model

mod_inla_poisson <- inla(mal_inc ~ year + f(state,model="iid") + 
                             ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS , 
                           data = data_df, family = "poisson")

mod_inla_poisson

### spcefication is wrong
#formula <- mal_inc ~ year + f()
data_df$state <- as.factctor(data_df$state)
mod_inla_poisson <- inla(mal_inc ~ year + f(year,state,model="iid") + 
                           ONI_DJF + DMI_ASO + MJO_DJFM + MJO_JJAS , 
                         data = data_df, family = "poisson")

mod_inla_poisson
summary(mod_inla_poisson)

############# Spatial model

cat(as.character(region_sf$NAME_1))
cat(as.character(data_df$state)[1:3])

test <- agrep(as.character(region_sf$NAME_1)[3],
      as.character(data_df$state))             
#agrep("lasy", "1 lazy 2")
             
             
### join data with names of states
region_sf$state <- unique(as.character(data_df$state))
View(region_sf)

as.character(region_sf$NAME_1)
as.character(region_sf$state)

state_df <- read.table(file.path(in_dir,"state_labeling_matching.csv"),sep=",")
names(state_df) <- c("Names","state")
region_sf$state <- state_df$state
View(region_sf)

data_sf <- merge(region_sf, data_df,by="state")
dim(data_sf)
class(data_sf)
#missing states
plot(data_sf$geometry)

## keep all
data_sf <- merge(region_sf,data_df,by="state",all=T)
data_sf <- merge(data_df,region_sf,by="state",all=T)

plot(data_sf$geometry)
plot(data_df$mal_inc)
dim(data_sf)
dim(data_df)
#problem here, we loose some rows

###################
### Spatial model specification

### First get a neighbour or adjacency matrix object compatible with INLA

region_sp <- as(region_sf,"Spatial")
temp <- poly2nb(region_sp)
nb2INLA("LDN.graph", temp)
#LDN.graph is stored in in the current dir and can be read by INLA
class(temp)
LDN.adj <- paste(getwd(),"/LDN.graph",sep="") # full path to the INLA adjacency/graph file 

H <- inla.read.graph(filename="LDN.graph")

## Neighbours for the states: adjacency matrix
image(inla.graph2matrix(H),xlab="",ylab="")

region_sf$ID_1

names(data_df)
dim(data_sf)
class(data_sf)
View(state_df)

state_df$ID <- 1:36

#Use iCAR specification
####keep all
data_df_test <- merge(data_df,state_df,by="state",all=T)
dim(data_df_test)
dim(data_df)
#lossing 8 rows
names(data_df_test)
names(data_df)
formula.par <- y ~ 1 + f(ID,model="bym",graph=LDN.adj)

## use informative priors:
#1) log of the unstructure effect precision log(tauv) ~ logGamma(1,0.0005)

data_df$y <- data_df$mal_inc
data_df_test$y <- data_df_test$mal_inc

formula.par <- y ~ 1 + f(ID,model="bym",graph=LDN.adj,
                         scale.model=TRUE,
                         hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),
                                    prec.spatial=list(
                                    prior="loggamma",param=c(1,0.001))))

## In bayesian model, priors impact the results so should do a sensitivity analysis

names(data_df)
mod_spatial_inla <- inla(formula.par,
                         family = "poisson",
                         data=data_df_test,
                         E=E, #not found because not defined.
                         control.compute=list(dic=T))

mod_spatial_inla <- inla(formula.par,
                         family = "poisson",
                         data=data_df_test,
#                         E=E,
                         control.compute=list(dic=T))

summary(mod_spatial_inla)

    
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
