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