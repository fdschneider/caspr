#' Count patch sizes
#'
#' @param x A landscape object.
#' @param state A character value. The state of interest to be counted.
#'
#' @return Returns a vector of patch sizes (number of cells) occuring in the landscape.
#' 
#' @details Might be replaced by a faster version implemented in C++. 
#' 

patches <- function(x, state = levels(x$cells)[1]) {
  pattern <- x$cells
  pattern <- pattern %in% state
  map <- rep(NA, times = prod(x$dim))
  old <- rep(99, times = prod(x$dim)) 
  
  while(!identical(old[pattern], map[pattern])) {
    old <- map
    count = as.integer(1)
    for(i in which(pattern)) {
      neighbors <- map[x_with_border][x_to_evaluate[i]+interact]
      if(all(is.na(neighbors)) ) { 
        map[i] <- count
      } else {
        map[i] <- min(neighbors, na.rm = TRUE)
      }
      count <- count +1
    }
    
  }
  
  map <- as.factor(map)
  patchvec <- as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
  
  out <- vector()
  if(length(patchvec) > 0) out <- sort(patchvec) else out <- NA
  #out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
  return(out)
  
} 



###################################
## fit power laws on patch count ##
###################################



fitPL <- function(psd, p_spanning, n = NULL) {
  
  # code of fitted classes
  
  n_plants <- sum(psd$size * psd$n)/n
  
  out <- list()
  out$best <- NA
  out$AIC <- vector("numeric", length = 3)
  out$dAIC <- vector("numeric", length = 3)
  
  # criteria for vegetated state & desert state
  
  ##### linear power law model for parameter estimation
  PLlm <- lm(I(log(p)) ~  1 - I(log(size)) , data = psd) 
  
  ###########
  
  try( {out$TPLdown <- nls(I(log(p)) ~ I( alpha*log(size)-size*Sx ),
                           data = psd,
                           start = list(alpha =  PLlm$coefficients, Sx = 1/200),
                           #algorithm = "port",
                           trace = FALSE
  )}, silent = TRUE
  )    
  
  if(!is.null(out$TPLdown) & !coefficients(out$TPLdown)["Sx"] <= 0 ) {
    out$AIC[1] <- AIC(out$TPLdown) 
  } else {
    out$TPLdown <- list(NA)
    out$AIC[1] <- NA
  }
  
  #####
  
  try({out$PL <- nls(I(log(p)) ~ alpha * log(size), 
                     data = psd,
                     start = list( alpha =  PLlm$coefficients ),
                     trace = FALSE,
                     nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  if(!is.null(out$PL)) {
    out$AIC[2] <- AIC(out$PL)
  } else {
    out$PL  <- list(NA)
    out$AIC[2] <- NA
  }
  
  ###########
  
  try({out$TPLup <- nls(I(log(p)) ~  log(b) + log(1+(size^(alpha))/b ) , 
                        data = psd,
                        start = list( alpha =  PLlm$coefficients, b = p_spanning ) , 
                        nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  
  if(!is.null(out$TPLup)) {
    out$AIC[3] <- AIC(out$TPLup) 
  } else { 
    #result$fit$summary$TPLup  <- list(NA)
    out$TPLup  <- list(NA)
    out$AIC[3] <- NA
  }
  
  ###########
  
  out$dAIC <-   out$AIC -min(out$AIC, na.rm = TRUE)
  
  out$best <- which.min(out$AIC)+1
  
  return(out)
} 

