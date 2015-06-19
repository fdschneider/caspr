

context('Test the transfer functions as.matrix and as.landscape')

test_that('Landscape init works', { 
  
  states = c("+","0","-")
  cover = c(0.5,0.3,0.2) 
  width = 50
  height = 50
  
  l <- init_landscape(states, cover, width, height)
  
  # test for matching dimensions
  expect_equal(l$dim, as.landscape(as.matrix(l))$dim)
  # test for matching content of vector cells
  expect_equal(l$cells, as.landscape(as.matrix(l))$cells)
  # test for matching sequence of levels
  expect_equal(levels(l$cells), levels(as.landscape(as.matrix(l))$cells))

  m <- matrix(sample(states, width*height, replace = TRUE, prob = cover ), ncol = w, byrow = TRUE  )
  
  # test for matching dimensions
  expect_equal(dim(m)[1], as.matrix(as.landscape(m))$dim$width)
  expect_equal(dim(m)[2], as.matrix(as.landscape(m))$dim$height)
  # test for matching content of vector cells
  expect_equal(m, as.matrix(as.landscape(m)))
  
})
