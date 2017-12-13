#' generate_founders function
#'
#' This function generates the first generation of models.
#' @param x(data set)
#' @keywords founders
#' @export
#' @examples
#' generate_founders(mtcars)

generate_founders <- function(X) {
  
  # number of predictors ---------------
  C <- dim(X)[2]
  
  # number of founders ----------------
  P <- 2 * C
  if (P > 200) {P <- 200
  } else if (P < 10) {
    P <- 10}
  
  #randomly generate founders ----------------
  geneSample <- sample(c(0, 1), replace = TRUE,
                       size = ceiling(1.2 * C * P))
  
  #update geneSample to make sure that each gene/variable
  #will exist in at least one chrome, but not all.
  #geneSample <- c(rep(0, C - 1), rep(1, C), 0, geneSample)
  
  #create a first generation
  indices <- seq_along(geneSample)
  first_gen_raw <- split(geneSample, ceiling(indices / C))
  
  generation_t0 <- matrix(unlist(unique(first_gen_raw)[1:P]),
                          ncol = C, byrow = TRUE)
  
  generation_t0 <- generation_t0[apply(generation_t0[, -1], 1,
                                       function(x) !all(x == 0)), ]
  
  
  return(generation_t0)
}