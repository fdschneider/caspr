# This code loads a vector of patch size, fits several functions (MLE fitting) to the patch size distribution, compares the fits and plots the data + fitted functions.
# Most of the code is based on Cosma Shalizi's code: tuvalu.santafe.edu/~aaronc/powerlaws/
# Dec 2015

# Load the loadfunction.Rdata file (just open the file)
# Load package 'poweRlaw'
library(poweRlaw)

setwd("C:/Users/Sonia/Documents/GitHub/caspr/draft/mledistfit")

# Load data
data<- read.csv("exdata.csv", header=FALSE)
data <- as.vector(t(data))

# PL fit (analytical)
# Note: this is for discrete data (change 'dis' to 'con' for continuous data)
pow = displ$new(data)        # stores data in prep for parameter estimation
minx = estimate_xmin(pow)    #use if you need xmin to be estimated
pow$setXmin(minx)
est1 = estimate_pars(pow)    # estimates parameters of the PL using MLE
pow$setPars(est1$pars)       # stores estimated pars in the appropriate way

# Exp fit (analytical)
# Note: this is for discrete data (change 'dis' to 'con' for continuous data)
exp = disexp$new(data)
#exp$setXmin(minx)
# For comparing distributions, keep same xmin as estimated for the powerlaw fit
# In that case, use this command below in place of the one above: 
exp$setXmin(pow$getXmin())
est2 = estimate_pars(exp)
exp$setPars(est2$pars)

# Lognorm fit (analytical)
# Note: this is for discrete data (change 'dis' to 'con' for continuous data)
lnorm = dislnorm$new(data)
#lnorm$setXmin(minx)
lnorm$setXmin(pow$getXmin())
est3 = estimate_pars(lnorm)
lnorm$setPars(est3$pars)

# Compare fits
com = compare_distributions(pow,exp) # Compare all required combinations (by pair)
com$test_statistic           # vuong test statistic (pvals calculated from this)
com$p_two_sided 
# p = 1 implies that both distributions are equally far from the data
# p = 0 implies that one of the distributions fits the data significantly better
com$p_one_sided
# p = 1 implies that the first distribution is the better fit
# p = 0 implies that the first distribution can be rejected in favour of the second

# p-2 sided = 0; p-1 sided = 1 => 1st distribution wins
# p-2 sided = 0; p-1 sided = 0 => 2nd distribution wins

# Compare these distributions with power law with exp cutoff
# All the distrib have to be numerically fit again for the comparison
exp1 = exp.fit(data,pow$getXmin())        # exponential
pow1 = pareto.fit(data,pow$getXmin())     # power law
lnorm1 = lnorm.fit(data,pow$getXmin())    # log-normal
powexp = powerexp.fit(data,pow$getXmin()) # power law with exponential cut-off

power.powerexp.lrt(pow1,powexp)           # compare pareto and powerexp
# gives only a single p-value. If p=1, pow is better,  If p = 0, powerexp is better
exp.powerexp.lrt(exp1,powexp)
# gives only a single p-value. If p=1, exp is better, If p = 0, powerexp is better
# gives only a single p-value. There is an issue here. 
# If it is an underlying exponential distribution (extreme case), powerexp will be chosen for
# large sample sizes. In such a case, the rate of decline will be similar
vuong(lnorm.powerexp.llr(data, lnorm1, powexp, pow$getXmin())) # compare lnorm and powerexp
# 2 p values similar to the previous analysis
#If p.two.sided = 0, one of the distributions fits better, if p.one.sided = 1, the second distribution is a better fit.


# To evaluate how good the fit of the best distribution is: perform bootstrapping
# (see the documentation of the R 'poweRlaw package)


# Plot disctrete fits
  dat = plot(pow)
  po = lines(pow)
  ex = lines(exp)
  lno = lines(lnorm)
  
  fin = rbind(po,ex,lno)
  fin$type = c(rep("Power law",length(po$x)),rep("Exponential",length(ex$x)),rep("Log-normal",length(lno$x)))
  
  
  ## some problem with the ggplot command below. Will fix and send
  ggp = ggplot(dat, aes(x=x, y=y))  +
    geom_point(size = 1, col = "black") +
    #stat_smooth(data = fin, aes(x = x, y = y, col = type), se = F) +
    xlab("Inter Patch Distance") +
    ylab("Log Inverse CDF") +
    scale_colour_hue(name="Distribution", 
                     breaks=c("Exponential","Log-normal","Power law"),
                     labels=c("Exponential", "Log-normal", "Power law"),
                     l=40) +
    theme_bw() 
  ggp + 
    theme(axis.title.x = element_text(vjust = 0.3, size = 16, face = "bold"), axis.text.x = element_text(size = 12), axis.title.y = element_text(face = "bold", vjust = 0.3, angle = 90, size = 16), axis.text.y = element_text(size = 12)) +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 12))+
    theme(legend.justification=c(1,1), legend.position=c(0.3,0.5)) +
    scale_x_log10()+
    scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1), limits = c(0.0001,1.01))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      ##panel.border = theme_blank(),
      ##panel.background = theme_blank()
    )

