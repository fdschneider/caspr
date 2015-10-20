
#' print function for model objects
#'
#' @param x 
#'
#' @export
#'

print.ca_model <- function(x) { 
  
  cat("\n", x$name, "\n")
  cat(" ", rep("-", nchar(x$name)), "\n\n", sep = "" ) 

  # Print default parameters
  cat(" Possible cell states:")
  cat(" ", paste(x$states, collapse=", "), "\n\n")
  
  cat(" Required parameters and default values: \n")
  for (i in seq_along(x$parms)) { 
    cat("  ", paste0(names(x$parms)[i], " = ", round(x$parms[[i]], 5), "\n"))
  }
  
  # Print ref
    cat("\n Original reference: ")
    cat(x$ref,"\n")
}
