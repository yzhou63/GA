select <- function(Y, X, iter, objFun = c("AIC", "BIC", "logLik", "user"),
                  family = "gaussian",
                  crossMeth = c("method1", "method2", "method3", "user"),
                  user_cross_fn = NULL,
                  numChromosomes = NULL,
                  pCrossover = 0.8,
                  mutation_rate = NULL,
                  user_objFun = NULL,
                  converge = TRUE, minimize = TRUE, parallel = FALSE,
                  userfun = NULL) {

    #input objects
    if(missing(iter)) {iter <- 100}
    if(missing(objFun)) {objFun <- "AIC"}
    if(missing(family)) {family <- "gaussian"}
    if(missing(converge)) {converge <- FALSE}
    if(missing(minimize)) {minimize <- TRUE}
    if(missing(parallel)) {parallel <- FALSE}
    if(missing(user_objFun)) {user_objFun <- NULL}
    if(missing(mutation_rate)) {mutation_rate <- "1 / (P * sqrt(C))"}



    require(abind)
    require(parallel)

    #error checking

    ################
    # functions
    ################


    # Selection of parents function ----------------
    select_parents <- function(generation_t0) {

        rank <- get("rank", mode = "any", envir = parent.frame())
        P <- get("P", mode = "any", envir = parent.frame())

        # probability of selection ----------------
        phi <- (2 * rank) / (P * (P + 1))

        #method 1: both parents selected by rank
        #parentInd <- sample(1:P, 2, prob=phi, replace = F) #this is better than method 2

        #method 2: first parent by rank, second random
        parentInd <- c(sample(1:P, 1, prob=phi, replace = F),
                       sample(1:P, 1, replace = F))


        return(parentInd)
    }

    # crossover function ----------------
    crossover_parents <- function(generation_t0, parentInd, crossMeth, pCrossover)  {

        # get parent info
        parent1 <- generation_t0[parentInd[1], ]
        parent2 <- generation_t0[parentInd[2], ]
        rank <- get("rank", mode = "any", envir = parent.frame())
        C <- length(parent1)
        parent1r <- rank[parentInd[1]]
        parent2r <- rank[parentInd[2]]

        if (rbinom(1, 1, pCrossover) == 1 ) {
            if (crossMeth == "method1") {

                #METHOD 1 ----------------
                #multipoint crossover: three crossover points
                cross <- sort(sample(seq(2,(C - 2), 2), 3, replace = F))

                child1 <- c(parent1[1:cross[1]],
                            parent2[(cross[1] + 1):cross[2]],
                            parent1[(cross[2] + 1):cross[3]],
                            parent2[(cross[3] + 1):C])
                child2 <- c(parent2[1:cross[1]],
                            parent1[(cross[1] + 1):cross[2]],
                            parent2[(cross[2] + 1):cross[3]],
                            parent2[(cross[3] + 1):C])

            } else if (crossMeth == "method2") {

                #METHOD 2 ----------------
                #method upweights parent with higher rank high
                childProb <- parent1 * parent1r[1] /
                    (parent1r + parent2r) +
                    parent2 * parent2r /
                    (parent1r + parent2r)
                child1 <- rbinom(C, 1, prob = childProb)
                child2 <- rbinom(C, 1, prob = childProb)

            } else if (crossMeth == "method3") {

                #METHOD 3 ----------------
                #randomly samples non-concordant variables between parents
                # slightly upweights parent selected by prob. proportional to rank
                child1 <- parent1
                child2 <- parent2
                child1[parent1 != parent2] <-
                    rbinom(sum(parent1 - parent2 != 0), 1,
                           prob = parent1r / (parent1r + parent2r))
                child2[parent1 != parent2] <-
                    rbinom(sum(parent1 - parent2 != 0), 1,
                           prob = parent2r / (parent1r + parent2r))

            }
            return(rbind(child1, child2))
        } else {
            child1 <- parent1
            child2 <- parent2
            return(rbind(child1, child2))
        }
    }


    # mutation function ----------------
    mutate_child <- function(mutation_rate, child) {

        C <- length(child)
        P <- get("P", mode = "any", envir = parent.frame())
        abs(round(child, 0) -
                rbinom(C, 1,
                       prob = eval(parse(text = mutation_rate))))
    }

    ##########
    # Perform genetic algorithm
    #########

    #step 1: Generate founders ----------------
        generation_t0 <- generate_founders(X)
        P <- nrow(generation_t0) #num chromosomes
        cat("1. Generate founders: ", P, "chromosomes")
        t1 <- Sys.time()

    #Step 2. Evaluate founder fitness Fitness of inital pop  ----------------
        objFunOutput_t0 <- evaluate_Fitness(generation_t0, Y, X)
        cat("\n2. Evaluate founder fitness")
        t1 <- c(t1, Sys.time())

        #create array to store fitness data for each iteration
        convergeData <- array(dim = c(P, 2, 1)) #P x 2 x iter

        #save founder fitness evaluation data
        convergeData[, , 1] <- objFunOutput_t0[
                                order(objFunOutput_t0[, 2], decreasing = T),
                                c(1, 3)]

    #Step 3. loop through successive generations until either:
        tol = 1e12 * .Machine$double.eps #tolerance for converage check
        converged <- 0

        cat("\n3. Begin breeding \n Generations: ")

        for (i in 1:iter) {

                # 1. create next generation ----------------
                if (i == 1) {
                    generation_t1 <- generation_t0
                    objFunOutput_t1 <- objFunOutput_t0
                }

                generation_t1 <- create_next_generation(generation_t1,
                                                        objFunOutput_t1,
                                                        select_parents,
                                                        crossover_parents,
                                                        mutate_child)

                # 2. evaluate children fitness ----------------
                objFunOutput_t1 <- evaluate_Fitness(generation_t1, Y, X)

                # store fitness data
                convergeData <- abind(convergeData,
                                        objFunOutput_t1[order(objFunOutput_t1[, 2],
                                                            decreasing = T), c(1, 3)])

                # cat generation and save timing
                cat(i,"-", sep = "")
                t1 <- c(t1, Sys.time())

                # 3. check convergence ----------------
                if (i > 10 & converge == TRUE) {
                    if(isTRUE(all.equal(mean(convergeData[1:(P * 0.25), 2, i]),
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

    #Step 4: process output ----------------
    bestModel <- generation_t1[convergeData[1, 1, i], ]
    value <- convergeData[1, 2, dim(convergeData)[3]]
    if(dim(convergeData)[3] < iter) {converged <- "Yes"
    } else {converged <- "No"}

    output <- list("BestModel" = colnames(X)[bestModel == 1],
                   objFun = list("Obj_fn" = objFun,
                                 value = as.numeric(round(value, 4))),
                   iter = dim(convergeData)[3],
                   converged = converged,
                   convergeData = convergeData,
                   timing = t1)
    return(output)
}




