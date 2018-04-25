#######







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
