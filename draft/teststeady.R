l_50 <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = 50) # create initial landscape
r_50 <- runca(l_50, model = grazing, t_max = 5000) 

t_eval = 400
steadyval <- sapply(t_eval:5001, function(i) {
  t_1 <- (i-as.integer(t_eval)):(i-0.5*as.integer(t_eval)) 
  # get vector of the last t_eval timesteps 
  t_r <- sample((i-as.integer(t_eval)):i, t_eval/2)
  t_2 <- (i-0.5*as.integer(t_eval):i) 
  
  (mean(r_50$cover[[1]][t_2]) - mean(r_50$cover[[1]][c(t_1,t_2)]))
  #abs(sd(r_50$cover[[1]][t_1])- sd(r_50$cover[[1]][t_2]))
})


plot(steadyval ~ r_50$time[t_eval:5001],  type = "l")


i = 700
interval <- rep(as.factor(1:10), each = t_eval/10)
anova(lm(r_50$cover[[1]][(i-as.integer(t_eval-1)):i]~interval))
boxplot(r_50$cover[[1]][(i-as.integer(t_eval-1)):i]~interval)


l_100 <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = 100) # create initial landscape
r_100 <- runca(l_100, model = grazing, t_max = 500) 
abs(mean(r_100$cover[[1]][401:500]) - mean(r_100$cover[[1]][301:400]))

l_200 <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = 200) # create initial landscape
r_200 <- runca(l_200, model = grazing, t_max = 500) 

par(mfrow = c(1,1))
plot(r_50$cover[[1]] ~ r_50$time, ylim = c(0.5,0.7), type = "l")
lines(r_100$cover[[1]] ~ r_100$time)
lines(r_200$cover[[1]] ~ r_200$time)

c(
abs(mean(r_50$cover[[1]][401:500]) - mean(r_50$cover[[1]][301:400])) ,
abs(mean(r_100$cover[[1]][401:500]) - mean(r_100$cover[[1]][301:400])),
abs(mean(r_200$cover[[1]][401:500]) - mean(r_200$cover[[1]][301:400]))
)

0.01 / c(50*50, 100*100, 200*200)
