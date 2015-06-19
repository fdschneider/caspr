
context('Test the initialization of landscapes')

test_that('Landscape init works', { 
  
  for (width in c(1, 100, 1000)) { # [Alex] Maybe add another size ? 
    for (height in c(1, 100, 1000)) { 
      for (states in list(c('-','0','+'), 
                          c('q','p'), 
                          letters, # 26 states
                          c('pred','prey','empty'))) { 
                            
        covers <- runif(length(states), 0, 1)
        covers <- covers/sum(covers)
        names(covers) <- states
        test_landscape <- init_landscape(states, covers, width, height)
        
        # Correct number of cells ?
        expect_equal(length(test_landscape$cells), height*width)
        # Does the returned object has the right attributes ? 
        expect_equal(c(width = width, height = height), test_landscape$dim)
        
        # Is the distribution of cells coherent with the given probabilities ?
        # (works only with large lattices)
        if (width*height >= 10000) {
          expect_equal(round(summary(test_landscape)$cover, digits=2), round(covers, digits = 2) )
        }
        
      }
    } 
  }
})
