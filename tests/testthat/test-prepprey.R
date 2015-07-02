context('Test predator/prey model')

test_that('Pred/prey is well-specified', { 
  expect_true(all(c("name","ref","states",
                    "cols","parms","update") %in% names(forestgap)))
})


test_that('Conversion of pred/prey states works correctly', { 
  
  x <- init_landscape(states = predprey$states,
                      cover  = rep(1/3,3))
  x_mat <- as.matrix(x, as = 'integer')
  
  # Test whether the back and forth conversion to matrix works
  expect_equal(x$cells, 
               .to_vector(x_mat, predprey$states))
  
  # Check if predprey_core preserves states correctly
  x_new <- with(predprey$parms, 
                predprey_core(x_mat, 1000, 4, betaf, betas, delta))
  expect_equal(levels(x$cells),
               levels(.to_vector(x_new, predprey$states)))
  
})

#   # Create a shark-only landscape -> they all die eventually
#   x_full <- x
#   x_full$cells[] <- rep('s', length(x$cells)) # the [] preserve levels
#   x_full_result <- ca(x_full, model = predprey)
#   ggplot(fortify(x_full_result)) + 
#     geom_raster(aes(x, y, fill=state)) +
#     facet_wrap( ~ time )