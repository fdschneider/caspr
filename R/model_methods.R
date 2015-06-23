# 
# 
# Defines a summary function for class "ca_model" so we can have a quick 
# reminder of the parameters it has
# 

summary.ca_model <- function(x, with_ref=FALSE) { 
  
  
  cat(x$name, "\n")
  cat("--------------------\n")
  
  # Print default parameters
  cat("Possible cell states: ", 
      paste(x$states, collapse=", "), "\n")
  
  
  cat("Default parameters: \n")
  for (i in seq_along(x$parms)) { 
    cat(paste0(names(x$parms)[i], " = ", x$parms[[i]], "\n"))
  }
  
  # Print ref
  if (with_ref) { 
    cat("Original reference:\n")
    cat(x$ref,"\n")
  }
}
