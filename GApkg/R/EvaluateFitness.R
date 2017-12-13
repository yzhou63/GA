#' evaluate_fitness function
#'
#' This function evaluates the goodness of fit with user specified criteria
#' @param genereation_t0, Y, X
#' @keywords fitness
#' @export
#' @import doParallel
#' @examples
#' function(generation_t0, Y, X)

evaluate_fitness <- function(generation_t0, Y, X) {
  
  # get objects from their parent frames----------------
  family <- get("family", mode = "any", envir = parent.frame())
  objFun <- get("objFun", mode = "any", envir= parent.frame())
  parallel <- get("parallel", mode = "any", envir= parent.frame())
  minimize <- get("minimize", mode = "any", envir= parent.frame())
  
  #test parms
  #family <- "gaussian"
  ##objFun <- "AIC"
  #minimize <- TRUE
  #parallel <- FALSE
  
  # number of parent chromosomes or models----------------
  P <- dim(generation_t0)[1]
  
  # calculate objective function ----------------
  calc_objective_function <- function(mod) {
    
    objFun <- get("objFun", mode = "any", envir = parent.frame())
    
    if (objFun == "AIC") {
      extractAIC(mod)[2]
    } else if (objFun == "BIC") {
      BIC(mod)
    } else if (objFun == "logLik") {
      logLik(mod)
    } else if (objFun == "user") {
      user_func(mod)
    }
  }
  
  # rank obj function output ----------------
  rank_objective_function <- function(objFunOutput){
    
    objFun <- get("objFun", mode = "any", envir = parent.frame())
    minimize <- get("minimize", mode = "any", envir = parent.frame())
    
    if (objFun %in% c("AIC", "BIC") | minimize == TRUE) {
      r <- rank(-objFunOutput, na.last = TRUE, ties.method = "first")
    } else if (objFun == "logLike" | minimize == FALSE) {
      r <- rank(objFunOutput, na.last = TRUE, ties.method = "first")
    }
    return(cbind(chr = 1:P, rank = r, objFunOutput))
  }
  
  
  ######
  #evaluate and rank each chromosome with selected objective function
  ######
  
  # serial ----------------
  if (parallel == FALSE) {
    
    # lm ----------------
    if (family == "gaussian") {
      objFunOutput <- sapply(1:P, function(i) {
        mod <- lm(Y ~ X[, generation_t0[i, ] == 1])
        calc_objective_function(mod)
      })
      # glm ----------------
    } else if(family != "gaussian") {
      objFunOutput <- sapply(1:P, function(i) {
        mod <- glm(Y ~ X[, generation_t0[i, ] == 1], family = family)
        calc_objective_function(mod)
      })
    }
    
    # parallel ----------------
  } else if (parallel == TRUE) {
    
    # mclapply options ----------------
    nCores <- detectCores() - 1
    if(dim(X)[1] < 1000) {
      preschedule <- FALSE
    } else {
      preschedule <- TRUE
    }
    
    # lm ----------------
    if (family == "gaussian") {
      objFunOutput <- unlist(mclapply(1:P, function(i) {
        mod <- lm(Y ~ X[, generation_t0[i, ] == 1])
        calc_objective_function(mod)
      }, mc.preschedule = preschedule, mc.cores = nCores))
      # glm ----------------
    } else if(family != "gaussian") {
      objFunOutput <- unlist(mclapply(1:P, function(i) {
        mod <- glm(Y ~ X[, generation_t0[i, ] == 1], family = family)
        calc_objective_function(mod)
      }, mc.preschedule = preschedule, mc.cores = nCores))
    }
  }
  
  # rank ----------------
  rank <- rank_objective_function(objFunOutput)
  
  #return rankings ----------------
  return(rank)
}