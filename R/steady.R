
#' Test for steady state dynamics of timeseries.
#' 
#' @param i
#' @param result
#' @param steadyparms
#' @param returnvalue
#'   
#' @return a binary value. TRUE, if steady state is reached. FALSE, if still in
#'   transient dynamics.
#' @export
#' 

steady <- function(i, result, 
                   steadyparms = list(t_eval = 400, accept = 0.001),
                   returnvalue = FALSE) {
  
  if( i-as.integer(steadyparms$t_eval) < 0 ) {
    return(FALSE)
  } else {
    
    t_1 <- (i-as.integer(steadyparms$t_eval)):(i-0.5*as.integer(steadyparms$t_eval)) 
    # get vector of the last t_eval timesteps 
    t_2 <- (i-0.5*as.integer(steadyparms$t_eval)+1):(i) 
    
    if(result$cover[[1]][i] > 0) {
      # calculate steadiness, i.e. difference between the mean cover
      #    in the two evaluation periods t_1 and t_2
      
      steadiness <- abs(mean(result$cover[[1]][t_2]) - (mean(result$cover[[1]][c(t_1,t_2)])) )
      #steadiness <- t.test(result$cover[[1]][t_1],result$cover[[1]][t_2])[]$p.value
      
    } else {
      # set steadiness to 0 if cover is 0, immediate stop of simulation
      steadiness <- 0 
    }
    
    if(returnvalue) {steadiness} else {
      return(steadiness <= steadyparms$accept) 
    }
  }
  
  
}

