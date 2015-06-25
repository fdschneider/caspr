# 
# A file for testing the use of the forest gap model
# 

initial <- init_landscape(states = c("f","s","0"),
                          cover  = c(.3, .3, .3))

parms <- predprey$parms # Import defaults
output <- ca(initial, parms, model=predprey, t_max=1000)

library(ggplot2)
output_fmt <- fortify(output)

ggplot(output_fmt) + 
  geom_raster(aes(x,y,fill=state)) + 
  facet_wrap( ~ time ) + 
  scale_fill_manual(values=c('#111111','#FFFFFF'))
  
