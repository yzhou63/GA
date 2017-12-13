#' select function
#'
#' This function is the primary function that do variable selection.
#' @param Y, X, iter, objFun = c("AIC", "BIC", "logLik", "user"),
#'  family = "gaussian",
#'  crossMeth = c("method1", "method2", "method3"),
#'  numChromosomes = NULL,
#'  pCrossover = 1,
#'  user_objFun = NULL,
#'  converge = TRUE, minimize = TRUE, parallel = FALSE
#' @keywords cats
#' @export
#' @import doParallel abind scales
#' @examples
#' select()

select <- function(Y, X, iter, objFun = c("AIC", "BIC", "logLik", "user"),
                   family = "gaussian",
                   crossMeth = c("method1", "method2", "method3"),
                   numChromosomes = NULL,
                   pCrossover = 1,
                   user_objFun = NULL,
                   converge = TRUE, minimize = TRUE, parallel = FALSE) {
  
  #input objects
  if(missing(iter)) {
    iter <- 100
  }
  if(missing(objFun)) {
    objFun <- "AIC"
  }
  if(missing(family)) {
    family <- "gaussian"
  }
  if(missing(converge)) {
    converge <- FALSE
  }
  if(missing(minimize)) {
    minimize <- TRUE
  }
  if(missing(parallel)) {
    parallel <- FALSE
  }
  if(missing(user_objFun)) {
    user_objFun <- NULL
  }
  require(doParallel)
  require(abind)
  
  #error checking
  
  #function objects
  
  # print settings
  
  ##########
  # Perform genetic algorithm
  #########
  
  #step 1: Generate founders ----------------
  generation_t0 <- generate_founders(X)
  P <- nrow(generation_t0) #num chromosomes
  cat("1. Generate founders: ", P, "chromosomes")
  t1 <- Sys.time()
  
  #Step 2. Evaluate founder fitness Fitness of inital pop  ----------------
  objFunOutput_t0 <- evaluate_fitness(generation_t0, Y, X)
  cat("\n2. Evaluate founder fitness")
  t1 <- c(t1, Sys.time())
  
  #create array to store fitness data for each iteration
  convergeData <- array(dim = c(P, 2, 1)) #P x 2 x iter
  
  #save founder fitness evaluation data
  convergeData[, , 1] <- objFunOutput_t0[
    order(objFunOutput_t0[, 2], decreasing = T),
    c(1, 3)]
  
  
  
  #Step 3. loop through successive generations until either:
  #1. finish default or user specified number of iterations
  #2. convergence to specific optimum
  tol = sqrt(.Machine$double.eps) #tolerance for converage check
  converged <- 0
  
  cat("\n3. Begin breeding \n Generations: ")
  
  for (i in 1:iter) {
    
    #breeding: create children
    if (i == 1) {
      generation_t1 <- generation_t0
      objFunOutput_t1 <- objFunOutput_t0
    }
    
    generation_t1 <- create_next_generation(generation_t1,
                                            objFunOutput_t1, iter)
    
    # eval fitness ----------------
    objFunOutput_t1 <- evaluate_fitness(generation_t1, Y, X)
    
    # store fitness data ----------------
    convergeData <- abind(convergeData,
                          objFunOutput_t1[order(objFunOutput_t1[, 2],
                                                decreasing = T), c(1, 3)])
    
    # cat generation and save timing ----------------
    cat(i,"-", sep = "")
    t1 <- c(t1, Sys.time())
    
    # check convergence ----------------
    if (i > 10 & converge == TRUE) {
      if(isTRUE(all.equal(mean(convergeData[1:(P * 0.5), 2, i]),
                          convergeData[1, 2, i],
                          check.names = F,
                          tolerance = tol))) {
        converged <- converged + 1
        if (converged >= 1) {
          cat("\n#### Converged! ####")
          break
        }
      }
    }
  }
  
  
  #process output ----------------
  bestModel <- generation_t1[convergeData[1, 1, i], ]
  value <- convergeData[1, 2, dim(convergeData)[3]]
  if(dim(convergeData)[3] < iter) {converged <- "Yes"
  } else {converged <- "No"}
  
  output <- list("BestModel" = bestModel,
                 objFun = c(objFun, value),
                 iter = dim(convergeData)[3],
                 converged = converged,
                 convergeData = convergeData,
                 timing = t1)
  return(output)
}