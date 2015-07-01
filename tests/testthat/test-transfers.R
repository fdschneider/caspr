

context('Test the transfer functions as.matrix and as.landscape')

test_that('Landscape init works', { 
  
  states <- c("+","0","-")
  cover  <- c(0.5,0.3,0.2) 
  width  <- 50
  height <- 50
  l      <- init_landscape(states, cover, width, height)
  
  # test for matching dimensions
  expect_equal(l$dim, as.landscape(as.matrix(l), states)$dim)
  # test for matching content of vector cells
  expect_equal(l$cells, as.landscape(as.matrix(l), states)$cells)
  # test for matching sequence of levels
  expect_equal(levels(l$cells), 
               levels(as.landscape(as.matrix(l), states)$cells))
  
})
