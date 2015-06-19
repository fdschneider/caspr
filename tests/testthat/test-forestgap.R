# 
# Automated tests for model implementations
# 

context('Test forest gap model')

test_that('Forestgap is well-specified', { 
  expect_true(all(c("name","ref","states",
                    "cols","parms","update") %in% names(forestgap)))
})

# Helper function
meancover <- function(x, of) summary(x)$mean_cover[of]


landscapes <- list(empty_landscape = init_landscape(c('0','+'), cover=c(1,0)),
                   full_landscape  = init_landscape(c('0','+'), cover=c(0,1)),
                   only_one_cell   = init_landscape(c('0','+'), cover=c(1,0)))
landscapes$only_one_cell$cells[50] <- "+"

test_that("Test the forest gap model behavior", { 
  
  for (landscape in landscapes) { 
    # No alpha (recovery) should produce an empty landscape @ equilibrium
    parms <- forestgap$parms
    parms$alpha <- 0
    expect_true(meancover(ca(landscape, forestgap, parms), of='0') == 1)    
    
    # No spatial process should produce a random pattern in cells
    # <!TODO:Alex!>
}
    
})


# Plot bifurcation diagram and see if it corresponds to (Kubo et al. 1996) 
# This is disabled by default as it requires manual intervention
if (FALSE && require(plyr) && require(ggplot2)) { 
  test_dat <- ddply(data.frame(delta=c(runif(50,0,.3), seq(0.18,0.22, l=100))), 
                ~ delta, 
                function(df) { 
                  # Modify parms and do test run
                  forestgap$parms$delta <- df$delta
                  landscape <- init_landscape(c('0','+'), 
                                              list(c(.98,.02), 
                                                   c(.25,.75))[[sample(c(1,2),1)]],
                                               width=100, height=100) # as in paper
                  result <- ca(landscape, forestgap)
                  as.data.frame(c(forestgap$parms, 
                                  summary(result)$mean_cover))
                }, .progress='time', .parallel=FALSE)
  
  # Compare with Fig (3)/b in (Kubo et al. 1996) 
  ggplot(test_dat) + 
    geom_line(aes(delta, X0, group=X0>.75), color='red') + 
    geom_point(aes(delta, X0), size=2, alpha=.7) + 
    ylab("rho_0")
}
