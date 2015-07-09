# 
# A file for testing the use of the forest gap model
# 

initial <- init_landscape(states = c("f","s","0"),
                          cover  = c(.4, .1, .5),
                          width=200, height=200)

parms <- predprey$parms # Import defaults
parms$betaf <- 1/5
parms$betas <- 1/5
parms$delta <- 1/2

output <- ca(initial, parms, model=predprey, subs = 1000)

plot(output, plotstates=seq_along(predprey$states))

library(ggplot2)
output_fmt <- fortify(output)

ggplot(output_fmt) + 
  geom_raster(aes(x,y,fill=state)) + 
  facet_wrap( ~ time ) + 
  scale_fill_manual(values=c('#393939','#E04B29', '#EDEDED'))


# Test over range of deltas
parms$delta <- seq(0, .05, l=23*30)
parjob()
outarr <- carray(predprey, c(.4,.1,.5), parms = parms)

library(ggplot2)
ggplot(outarr) + 
  geom_point(aes(delta, mean_cover.f), alpha=.8)

