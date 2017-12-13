#' create_next_generation function
#'
#' This function creates new generations of model
#' @param generation_t0, objFunOutput_t0, iter
#' @keywords next generation
#' @export
#' @import scales
#' @examples
#' create_next_generation(generation_t0, objFunOutput_t0, iter)

create_next_generation <- function(generation_t0, objFunOutput_t0, iter) {
  #will create children, selecting mate pairs based upon
  #with higher objectFun scores
  #nextGen will be of same size of prev generation
  #will then perform mutation on resulting generation
  
  #phi: prob of selection, softmax using objFun output
  #P: # of chromosomes/models
  #C: # of genes/vars
  #r: object function
  #nextGeneration: stores new generation
  #P1, P2: selected parents
  
  #inherit variables
  crossMeth <- get("crossMeth", mode = "any", envir = parent.frame())
  
  #set variables
  P <- dim(generation_t0)[1]
  C <- dim(generation_t0)[2]
  rank <- objFunOutput_t0[, 2]
  score <- objFunOutput_t0[, 3]
  
  #set probability of selection
  phi <- (2 * rank) / (P * (P + 1))
  
  #set up plot to view sampling in each generation
  #require(scales)
  #plot(phi, score)
  #legend("bottomleft", c("parents", "sampled parents"), pch = c(21, 19))
  
  #Create matrix for next generation
  generation_t1 <- matrix(NA, dim(generation_t0)[1], dim(generation_t0)[2])
  
  #keep high rank parents and put in new generation
  #goodGeneration_t0 <- generation_t0[rank > quantile(rank, probs = 0.95), ]
  #generation_t1[1:dim(goodGeneration_t0)[1], ] <- goodGeneration_t0
  
  #########
  #Selection, Crossover, and Mutation
  #########
  
  i <- 1 #dim(goodGeneration_t0)[1] + 1 #initialize while loop
  while(i <= dim(generation_t1)[1]) {
    
    #SELECTION: select Parents
    #method 1: both parents selected by rank
    #parentInd <- sample(1:P, 2, prob=phi, replace = F) #this is better than method 2
    
    #method 2: first parent by rank, second random
    parentInd <- c(sample(1:P, 1, prob=phi, replace = F),
                   sample(1:P, 1, replace = F))
    
    parents <- generation_t0[parentInd, ]
    parentRank <- rank[parentInd]
    parentScore <- score[parentInd]
    
    #update plot with sampled parents
    #points(phi[parentInd], score[parentInd], pch = 19,
    #       col = alpha("black", 0.5))
    
    #CROSSOVER and MUTATION parent to  create child
    
    if (crossMeth == "method1") {
      
      #METHOD 1:
      #crossover: two crossover points, non-overlapping
      #creates two children
      cross1 <- sample(2:(C/2-1), 1)
      cross2 <- sample((C/2+1):(C-1), 1)
      child1 <- c(parents[1, 1:cross1],
                  parents[2, (cross1+1):cross2],
                  parents[1, (cross2+1):C])
      
      child2 <- c(parents[1, 1:cross1],
                  parents[2, (cross1+1):cross2],
                  parents[1, (cross2+1):C])
      
      #mutation:
      child1 <- abs(round(child1, 0) -
                      rbinom(C, 1, prob = 1 / (P * sqrt(C))))
      child2 <- abs(round(child2, 0) -
                      rbinom(C, 1, prob = 1 / (P * sqrt(C))))
      
      
    } else if (crossMeth == "method2") {
      
      #METHOD 2:
      #crossover: method upweights parent with higher rank high
      #create one child
      childProb <- parents[1, ] * parentRank[1] /
        (parentRank[1] + parentRank[2]) +
        parents[2, ] * parentRank[2]  /
        (parentRank[1] + parentRank[2])
      #mutation: this method has A LOT of mutation
      child1 <- rbinom(C, 1, prob = childProb)
      child2 <- rbinom(C, 1, prob = childProb)
      
    } else if (crossMeth == "method3") {
      
      #METHOD 3
      #randomly samples non-concordant variables
      #between parents, slightly upweights parent selected
      #by prob. proportional to rank
      #this provide mutation as well
      #2 children
      child1 <- parents[1, ]
      child2 <- parents[2, ]
      child1[parents[1, ] != parents[2, ]] <-
        rbinom(sum(parents[1, ] - parents[2, ] != 0), 1,
               prob = parentRank[1] / (parentRank[1] + parentRank[2]))
      child2[parents[1, ] != parents[2, ]] <-
        rbinom(sum(parents[1, ] - parents[2, ] != 0), 1,
               prob = 1 - (parentRank[1] / (parentRank[1] + parentRank[2])))
    }
    
    #check if duplicate child, and for all non-zero genes
    #if (all(!rowSums(generation_t1 == child1 |
    #                 generation_t1 == child2, na.rm = T) == C) &
    #    sum(child1) > 0 & sum(child2) > 0) {
    #not dup: add child to new generation
    generation_t1[c(i, i + 1), ] <- rbind(child1, child2)
    #update counter
    i <- i + 2
    #}
  }
  #return new new generation
  return(generation_t1)
}