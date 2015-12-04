#' @title Forest gap model
#'
#' @description Dynamic model for forest pattern after recurring wind
#'   disturbance with two cell states: empty (\code{"0"}) and vegetated
#'   (\code{"+"}).
#'
#' @usage ca(l, forestgap, parms = p)
#'
#' @param alpha A numerical value. Reproduction rate per year.
#' @param delta A numerical value. Death rate depending on local gap density.
#' @param d A numerical value. intrinsic death rate.
#'
#' @author Kubo, T., Y. Iwasa and N. Furumoto.
#' @references
#'
#' \strong{Kubo, T. et al. (1996) Forest spatial dynamics with gap expansion:
#' total gap area and gap size distribution. J. Theor. Biol. 180, 229–246}
#'
#' Kizaki, S. and Katori, M. (1999) Analysis of canopy-gap structures of forests
#' by Ising–Gibbs states – Equilibrium and scaling property of real forests. J.
#' Phys. Soc. Jpn 68, 2553–2560
#'
#' Katori, M. (1998) Forest dynamics with canopy gap expansion and stochastic
#' Ising model. Fractals 6, 81–86
#'
#' Solé, R.V. and Manrubia, S.C. (1995) Self-similarity in rain forests:
#' evidence for a critical state. Phys. Rev. E. 51, 6250–6253
#'
#' Solé, R.V. and Manrubia, S.C. (1995) Are rainforests self-organized in a
#' critical state? J. Theor. Biol. 173, 31–40
#' @family models
#' @details This model can use either an explicit height for trees, in which
#' case states can be anywhere in a range [Smin...Smax] (Solé et al., 1995), or
#' use only two states, vegetated (non-gap) and empty (gap) (Kubo et al., 1996).
#' Here we focus on the version that uses only two states: gap (0) and non-gap
#' (+). Without spatial spreading of disturbance (all cells are independent), a
#' cell transitions from empty to vegetated with a birth probability b and from
#' vegetated to empty with death probability d.
#'
#' However, gap expansion occurs in nature as trees having empty (non-vegetated)
#' surroundings are more likely to fall due to disturbance (e.g. wind blows).
#' Let p(0) be the proportion of neighbouring sites that are gaps. We can
#' implement this expansion effect by modifying the death rate into \code{d +
#' delta p(0)}. Since 0 <= p(0) <= 1, delta represents the maximal added death
#' rate due to gap expansion (i.e. the spatial component intensity).
#'
#' In their simulations, the authors (Kubo et al. 1996) use a 100x100 torus-type
#' lattice (with random initial covers?).
#'
#' The authors consider two cases: one in which the recovery of trees is
#' proportional to the global density of vegetated sites, and one where the
#' recovery is proportional to the local density of vegetation. We use only the
#' first case as it the only one producing bistability.
#'
#' The birth rate b is replaced with alpha rho+ where rho+ represents the global
#' density of non-gap sites and alpha is a positive constant. This can produce
#' alternative stable states over a range of delta values within 0.15-0.2(alpha
#' is fixed to 0.20 and d to 0.01).
#'
#' The state transition probabilities thus become:
#'
#' b = alpha rho+
#'
#' and
#'
#' d = d_0 + delta p_0
#'
#' @examples
#'
#' l <- init_landscape(c("+","0"), c(0.6,0.4), width = 100) # create initial landscape
#' p <- list(alpha = 0.2, delta = 0.17, d = 0.01)   # set parameters
#' r <- ca(l, model = forestgap, parms = p, t_max = 100)    # run simulation
#'
#' r
#' plot(r)
#'
#' @export

"forestgap"

forestgap <- list()
class(forestgap) <- "ca_model"
forestgap$name   <- "Forest Gap Model"
forestgap$ref    <- paste0("Kubo, T., Y. Iwasa and N. Furumoto. 1996. Forest ",
                           "spatial dynamics with gap expansion: Total gap ",
                           "area and gap size distribution. Journal of ",
                           "Theoretical Biology 180:229–246.")
forestgap$states <- c("+", "0")
forestgap$cols   <- grayscale(2)
forestgap$parms  <- list(
  # Default values are taken from the original publication and produce
  # bistability (untested) for \delta within 0.15-0.2.
  alpha = 0.20, # recovery strength
  delta = 0.17, # strength of gap expansion
  d = 0.01      # intrinsic death rate
)

forestgap$update <- function(x,               # landscape object
                             parms,           # set of parms
                             subs = 10) {     # number of iterations/time step ?

  for (s in 1:subs) {

    # Get the local proportion of empty neighbors
    localempty <- neighbors(x, "0") / length(interact)
    rho_plus <- sum(x$cells == "+") / prod(x$dim)

    # 2 - drawing random numbers
    # one random number between 0 and 1 for each cell
    rnum <- runif(prod(x$dim), 0, 1)

    # 3 - set transition probabilities
    p_to_vegetated <- with(parms, (alpha * rho_plus) / subs
    p_to_empty     <- with(parms, (d + delta * localempty) / subs )

    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(p_to_empty, p_to_vegetated) > 1 )) {
      warning("a set probability is exceeding 1, please balance parameters!")
    }

    # 4 - apply transition probabilities
    x_new <- x
    x_new$cells[x$cells == "+" & rnum <= p_to_empty]     <- "0"
    x_new$cells[x$cells == "0" & rnum <= p_to_vegetated] <- "+"

    # 5- store x_new as next x
    x <- x_new

  }

  ## end of single update call
  return(x)
}
