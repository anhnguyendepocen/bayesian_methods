
###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)

### poisson regression
## 9 Chattisgarh, 15 Jharkhand, 30 Uttarakhand

in_dir <- "/home/bparmentier/Google Drive/Data/India_Research/malaria_study_by_state"
out_dir <- "/home/bparmentier/Google Drive/Data/India_Research/malaria_study_by_state"
setwd(out_dir)

#setwd("/Users/neeti/Documents/neeti/TERIwork/aditi/")
#mal_inc <- read.table('malaria_incidence_yrs.csv', header=TRUE, sep=',')
mal_inc <- read.table('malaria_incidence_yrs_94dmi_Mjo.csv', header=TRUE, sep=',')

#rain_index <- read.table('rainfall_index.csv', header=TRUE, sep=',')
rain_fall <- read.table('rainfall_data.csv', header=TRUE, sep=',')

View(mal_inc)
dim(mal_inc)
names(mal_inc)

lmp_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 6))
lmcoef_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 6))
lmstd_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 6))
lmcilow_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 6))
lmciabove_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 6))
aic_total <- array(data=NA, dim=c((ncol(rain_index)-1), 1))
#lmp_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 3))
#lmcoef_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 3))
#lmstd_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 3))
#lmcilow_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 3))
#lmciabove_arr <- array(data=NA, dim=c((ncol(rain_index)-1), 3))

#colnames(lmp_arr) <- c('intercept', 'oni', 'iod', 'ri')
#colnames(lmcoef_arr) <- c('intercept', 'oni', 'iod', 'ri')
#colnames(lmstd_arr) <- c('intercept', 'oni', 'iod', 'ri')
#colnames(lmcilow_arr) <- c('intercept', 'oni', 'iod', 'ri')
#colnames(lmciabove_arr) <- c('intercept', 'oni', 'iod', 'ri')

colnames(lmp_arr) <- c('intercept', 'oni', 'iod','mjo_djfm', 'mjo_jjas', 'ri')
colnames(lmcoef_arr) <- c('intercept', 'oni', 'iod', 'mjo_djfm', 'mjo_jjas','ri')
colnames(lmstd_arr) <- c('intercept', 'oni', 'iod','mjo_djfm', 'mjo_jjas', 'ri')
colnames(lmcilow_arr) <- c('intercept', 'oni', 'iod','mjo_djfm', 'mjo_jjas', 'ri')
colnames(lmciabove_arr) <- c('intercept', 'oni', 'iod','mjo_djfm', 'mjo_jjas', 'ri')

#colnames(lmp_arr) <- c('intercept', 'oni', 'iod',)
#colnames(lmcoef_arr) <- c('intercept', 'oni', 'iod')
#colnames(lmstd_arr) <- c('intercept', 'oni', 'iod')
#colnames(lmcilow_arr) <- c('intercept', 'oni', 'iod')
#colnames(lmciabove_arr) <- c('intercept', 'oni', 'iod')
state_names <- colnames(rain_index)
state_names <- state_names[-1]
row.names(lmp_arr) <-state_names
row.names(lmcoef_arr) <-state_names
row.names(lmstd_arr) <-state_names
row.names(lmcilow_arr) <-state_names
row.names(lmciabove_arr) <-state_names
row.names(aic_total) <-state_names
mal_inc1 <- mal_inc[-c(30, 35,36,37,38,40)]
for (i in 6:35){
  j <- i -4
#res <- glm(mal_inc1[,i] ~mal_inc1$ONI+mal_inc1$DMI+rain_index[,j], family=poisson(link=log))
  #res <- glm(mal_inc1[,i] ~mal_inc1$ONI+mal_inc1$DMI, family=poisson(link=log))
res <- glm(mal_inc1[,i] ~mal_inc1$ONI_DJF+mal_inc1$DMI_ASO+mal_inc1$MJO_DJFM+mal_inc1$MJO_JJAS + rain_fall[,j], family=poisson(link=log))
#res <- glm(mal_inc1[,i] ~mal_inc1$ONI_DJF+mal_inc1$DMI_ASO+ rain_fall[,j], family=poisson(link=log))

  #chaatisg_inc <- mal_inc1[,9]
  #chaatisg_inc <- chaatisg_inc[-c(1:6)]
  
#res <- glm(chaatisg_inc ~mal_inc1[(7:24),2]+mal_inc1[(7:24),3]+rain_fall[(7:24),6], family=poisson(link=log))

res_sum <- summary(res)
coef_res <- res_sum$coefficients[,1]
p_res <- res_sum$coefficients[,4]
stderror_coef <- res_sum$coefficients[,2]
aic_model <- res_sum$aic
aic_total[(i-5),] <- aic_model
lmp_arr[(i-5),] <- p_res
lmcoef_arr[(i-5),] <- exp(coef_res)
lmstd_arr[(i-5),] <- stderror_coef
lmcilow_arr[(i-5),] <- exp(coef_res - 1.96*stderror_coef)
lmciabove_arr[(i-5),] <- exp(coef_res + 1.96*stderror_coef)
#out <- capture.output(summary(res))

#cat("poisson regression analysis", out, file="summary_poissonregression.txt", sep="n", append=TRUE)
}
write.table(lmp_arr, 'pvalue_malpoisson_oni_iod94_rf_mjo.txt')
write.table(lmcoef_arr, 'coefvalue_malpoisson_oni_iod94_rf_mjo.txt')
write.table(lmstd_arr, 'stderror_malpoisson_oni_iod94_rf_mjo.txt')
write.table(lmcilow_arr, 'CIlow_malpoisson_oni_iod94_rf_mjo.txt')
write.table(lmciabove_arr, 'CIabove_malpoisson_oni_iod94_rf_mjo.txt')

write.table(aic_total, 'aic_malpoisson_oni_iod_rf_mjo.txt')