rm(list = ls())

library(MCMCpack)

#############################################################################
samps <- 1000     #<-- generate a large number of densities for a smooth curve
xseq <- seq(1e-16, 5, length = samps) #<-- range of x-axis (e.g., variance) to
                                      ## visualize over
nus <- c(0.002, 0.02, 0.2, 1, 10)  #<-- values of `MCMCglmm "nu" prior parameter
  clrs <- c("red", "black", "grey40", "grey60", "blue")  #<-- should = length(nus)
V <- 1  #<-- for a univariate/inverse Gamma prior V=1

# Now generate prior probability density for each value along x-axis
## do this for each separate value in `nus` so will end up with a matrix
dy <- sapply(nus, FUN = function(nu){
  dinvgamma(xseq, shape = nu / 2, scale = (nu * V) / 2)})

# plot these
## setup plot with first prior (`nus[1]`)
plot(dy[, 1] ~ xseq,
        type = "n",   #<-- don't plot data (do that below), just setup region
	main = "Inverse Gamma\nV = 1",
	xlab = "Variance", ylab = "Pobability density",
	xlim = c(0, max(xseq)), ylim = c(0, max(dy)))
  # now plot lines for each prior	
  ## use `rev()` so nu=0.002 plotted on top of all lines (not buried underneath)
  invisible(sapply(rev(seq(length(nus))),
	FUN = function(nu){
	  lines(dy[, nu] ~ xseq, lwd = 2, col = clrs[nu])}))
	        
  legend("topright", lwd = 2, col = clrs,
	title = "nu", legend = as.character(nus), inset = 0.01)


