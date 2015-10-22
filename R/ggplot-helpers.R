#' Fortify method for `ca_result`
#'
#' Transform the snapshots in a ca_result object into a data.frame so it can 
#' be plotted using ggplot.
#' 
#' @param x a ca_result object
#' 
#' @examples
#'
#' # Using the forest gap model
#' initial <- init_landscape(states=c("+","0"), cover=c(.5,.5))
#' parms <- forestgap$parms # Import defaults
#' parms$delta <- .15
#' output <- ca(initial, parms, model=forestgap, t_max=1000)
#' 
#' if (require(ggplot2)) { 
#'   ggplot(fortify(output)) + 
#'     geom_raster(aes(x,y,fill=state)) + 
#'     facet_wrap( ~ time ) + 
#'     scale_fill_manual(values=c('#111111','#FFFFFF'))
#' }
#' @export

fortify.ca_result <- function(x) { 
  output <- lapply(as.list(seq.int(length(x$landscapes))), # for all snapshots
              function(n) { 
                data.frame(time  = x$snaps[n],
                           fortify.landscape(x$landscapes[[n]]))
            })
                      
  do.call(rbind, output)
}

#' #' Fortify method for `landscape`
#'
#' Transform a `landscape` object into a data.frame so it can 
#' be plotted using ggplot.
#' 
#' @param x A `landscape` object
#' 
#' @examples
#'
#' # Using the forest gap model
#' initial <- init_landscape(states=c("+","0"), cover=c(.5,.5))
#' 
#' if (require(ggplot2)) { 
#'   ggplot(fortify(initial)) + 
#'     geom_raster(aes(x,y,fill=state)) + 
#'     scale_fill_manual(values=c('#111111','#FFFFFF'))
#' }
#' @export

fortify.landscape <- function(x) { 
  data.frame(expand.grid(x     = seq.int(x$dim[1]), 
                         y     = seq.int(x$dim[2]), 
                         KEEP.OUT.ATTRS = FALSE),
             state = x$cells)
}

