
rm(list=ls())

##############
#parent function for package
##############


#X <- x
#Y <- y
#objFun <- "AIC"
#family <- "gaussian"
#crossover_method <- "method1"
#pCrossover = 0.8
#iter <- 100
#parallel <- FALSE
#minimize <- T
#converge <- T
#calc_objective_function = calc_objective_function
#start_chrom = NULL
#mutation_rate = NULL

GAfun <- function(Y, X, iter, family = "gaussian",
                  objFun = c("AIC", "BIC", "logLik", "user"),
                  calc_objective_function = calc_objective_function,
                  crossover_parents = crossover_parents,
                  crossover_method = c("method1", "method2", "method3"),
                  start_chrom = NULL,
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

    #calc_objective_function <- match.fun(calc_objective_function)
    
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
    crossover_parents <- function(generation_t0, parentInd, 
                                  crossover_method, pCrossover, rank)  {
        
        # get parent info
        parent1 <- generation_t0[parentInd[1], ]
        parent2 <- generation_t0[parentInd[2], ]
        C <- length(parent1)
        parent1r <- rank[parentInd[1]]
        parent2r <- rank[parentInd[2]]
        
        if (rbinom(1, 1, pCrossover) == 1 ) {
            if (crossover_method == "method1") {
                
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

            } else if (crossover_method == "method2") {
                
                #METHOD 2 ----------------
                #method upweights parent with higher rank high
                childProb <- parent1 * parent1r[1] /
                    (parent1r + parent2r) +
                    parent2 * parent2r /
                    (parent1r + parent2r)
                child1 <- rbinom(C, 1, prob = childProb)
                child2 <- rbinom(C, 1, prob = childProb)

            } else if (crossover_method == "method3") {
                
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
    mutate_child <- function(mutation_rate, child, P) {
        
        C <- length(child)
        abs(round(child, 0) -
                rbinom(C, 1, 
                       prob = eval(parse(text = mutation_rate))))
    }
    
    # calculate obj function ----------------
    calc_objective_function <- function(mod) {
        
        objFun <- get("objFun", mode = "any", envir = parent.frame())
        #user_func <-match.fun(user_func)
        
        if (objFun == "AIC") {
            extractAIC(mod)[2]
        } else if (objFun == "BIC") {
            BIC(mod)
        } else if (objFun == "logLik") {
            logLik(mod)
        } else {
            ## need to add user function
        }
    }
    
    # rank obj function output ----------------
    rank_objective_function <- function(objFunOutput) {
        
        objFun <- get("objFun", mode = "any", envir = parent.frame())
        minimize <- get("minimize", mode = "any", envir = parent.frame())
        
        if (objFun %in% c("AIC", "BIC") | minimize == TRUE) {
            r <- rank(-objFunOutput, na.last = TRUE, ties.method = "first")
        } else if (objFun == "logLike" | minimize == FALSE) {
            r <- rank(objFunOutput, na.last = TRUE, ties.method = "first")
        }
        return(cbind(chr = 1:P, rank = r, objFunOutput))
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
        objFunOutput_t0 <- evaluate_fitness(generation_t0, Y, X,
                                            family, objFun,
                                            parallel, minimize, 
                                            calc_objective_function,
                                            rank_objective_function)
        
        
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
                                                        select_parents,
                                                        crossover_method,
                                                        crossover_parents,
                                                        pCrossover,
                                                        mutate_child,
                                                        mutation_rate)
   
                

                # 2. evaluate children fitness ----------------
                objFunOutput_t1 <- evaluate_fitness(generation_t1, Y, X,
                                                    family, objFun,
                                                    parallel, minimize, 
                                                    calc_objective_function,
                                                    rank_objective_function)

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

        

##############
# subfunctions
##############

#initiative founding chromosomes
generate_founders <- function(X, start_chrom = NULL) {

    #d: #number of genes/variables in design matrix
    #P: #number of parent chromosomes
    #geneSample: Random sample of genes/vars for initial parent dataet
    #firstgeneration:
    #generation_t:

    # number of predictors ---------------
    C <- dim(X)[2]

    # number of founders ----------------
    if (is.null(start_chrom)) {
        P <- 2 * C
        if (P > 200) {P <- 200
        } else if (P < 20) {
        P <- 10 * C}}
    else {P <- start_chrom}

    #randomly generate founders ----------------
    geneSample <- sample(c(0, 1),
                         replace = TRUE,
                         size = ceiling(1.2 * C * P))

    #update geneSample to make sure that each gene/variable
    #will exist in at least one chrome, but not all.
    #geneSample <- c(rep(0, C - 1), rep(1, C), 0, geneSample)

    #create a first generation
    x <- seq_along(geneSample)
    firstGen <- split(geneSample, ceiling(x / C))

    generation_t0 <- matrix(unlist(unique(firstGen)[1:P]),
                       ncol = C, byrow = TRUE)

    #generation_t0 <- generation_t0[apply(generation_t0[, -1], 1,
    #                                 function(x) !all(x == 0)), ]
    generation_t0 <- generation_t0[apply(generation_t0, 1,
                                         function(x) !all(x == 0)), ]
       

    return(generation_t0)
}

evaluate_fitness <- function(generation_t0, Y, X, family, objFun,
                             parallel, minimize, 
                             calc_objective_function, 
                             rank_objective_function) {

    # get objects ----------------
    #family <- get("family", mode = "any", envir = parent.frame())
    #objFun <- get("objFun", mode = "any", envir= parent.frame())
    #parallel <- get("parallel", mode = "any", envir= parent.frame())
    #minimize <- get("minimize", mode = "any", envir= parent.frame())

    #number parent chromosomes  
    P <- dim(generation_t0)[1]

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
        if(dim(X)[1] < 1000) {preschedule <- FALSE
        } else {preschedule <- TRUE}

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

    # return rankings ----------------
    return(rank)
}

create_next_generation <- function(generation_t0, objFunOutput_t0,
                                   select_parents,
                                   crossover_method,
                                   crossover_parents,
                                   pCrossover,
                                   mutate_child,
                                   mutation_rate) {

    # set variables
    P <- dim(generation_t0)[1]
    C <- dim(generation_t0)[2]
    rank <- objFunOutput_t0[, 2]
    
    # inherit variables
    #match.fun(user_cross_function)
    #crossover_method <- get("crossover_method", mode = "any", envir = parent.frame())
    #pCrossover <- get("pCrossover", mode = "any", envir = parent.frame())
    #mutation_rate <- get("mutation_rate", mode = "any", envir = parent.frame())
    
    select_parents <- match.fun(select_parents)
    crossover_parents <- match.fun(crossover_parents)
    mutate_child <- match.fun(mutate_child)

    #Create matrix for next generation
    generation_t1 <- matrix(NA, dim(generation_t0)[1], dim(generation_t0)[2])

    #########
    #Selection, Crossover, and Mutation
    #########

    i <- 1 #initialize while loop
    while(i <= dim(generation_t1)[1]) {

        
        # Selection ----------------
        parentInd <- select_parents(generation_t0)
        print("select")
        
        # Crossover ---------------- 
        children <- crossover_parents(generation_t0, parentInd, 
                                      crossover_method, pCrossover, rank)
        print("crossover")
        # Mutation ----------------
        child1 <- mutate_child(mutation_rate, children[1, ], P)
        child2 <- mutate_child(mutation_rate, children[2, ], P)
        print("mutate")
        # Check child dups and all zeros -----------------
        #if (all(!rowSums(generation_t1 == child1, na.rm = T) == C,
        #        !rowSums(generation_t1 == child2, na.rm = T) == C) &
        #        sum(child1) > 0 & sum(child2) > 0) {
            generation_t1[c(i, i + 1), ] <- rbind(child1, child2)
            #update counter
        i <- i + 2
        #}
        cat(i)
    }
    
    #return new new generation
    return(generation_t1)
}



##############
# simulation
##############

#simulate toy dataset
library(simrel)
n <- 100
p <- 10
m <- 5
q <- 5
gamma <- 2
R2 <- 0.8
relpos <- sample(1:p, m, replace = F)
dat <- simrel(n, p, m, q, relpos, gamma, R2)
x <- dat$X
y <- dat$Y

##############
#test function with simulated data
##############

output.meth1 <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                    crossover_method = "method1", converge = T, family = "gaussian",
                      pCrossover = 0.9)
output.meth1.all <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                          crossover_method = "method1", converge = F, family = "gaussian",
                      pCrossover = 0.9)
output.meth2 <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                          crossover_method = "method2", converge = T, family = "gaussian",
                          pCrossover = 0.9)
output.meth2.all <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                      crossover_method = "method2", converge = F, family = "gaussian",
                      pCrossover = 0.9)
output.meth3 <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                          crossover_method = "method3", converge = T, family = "gaussian",
                          pCrossover = 0.9)
output.meth3.all <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                      crossover_method = "method3", converge = F, family = "gaussian",
                      pCrossover = 0.9)


test <- rbind(c(output.meth1$objFun, output.meth1$iter),
        c(output.meth1.all$objFun, output.meth1.all$iter),
        c(output.meth2$objFun,     output.meth2$iter),
        c(output.meth2.all$objFun, output.meth2.all$iter),
        c(output.meth3$objFun,     output.meth3$iter),
        c(output.meth3$objFun,     output.meth3.all$iter))

colnames(test) <- c("ObjFun", "value", "iterations")
rownames(test) <- grep("output.meth", ls(), value = T)
test

#check variables selected
dat$relpred
output.meth1$BestModel
output.meth1.all$BestModel
output.meth2$BestModel
output.meth2.all$BestModel
output.meth3$BestModel
output.meth3.all$BestModel

#plot results
plots <- grep("output.meth", ls(), value = T)
for (i in 1:length(plots)) {

    output <- eval(parse(text = plots[i]))
    convergeData <- output$convergeData
    require(scales)

    #png(paste0("~/repos/STAT243/project/GAplots_method_0", i, ".png"),
    #    width = 640, heigh = 480)

    par(mfrow = (c(1, 1)))

        #all AIC across interations
        plot(jitter(rep(1, nrow(convergeData))), convergeData[, 2, 1], type = "p",
               pch = 19, col = alpha("blue", 0.1),
               ylim = c(min(convergeData[, 2, ], na.rm = T),
                        max(convergeData[, 2, ], na.rm = T)),
               xlim = c(1, output$iter), xlab = "Generations", ylab = "AIC",
               main = paste("GA performance using \n objFun = AIC, lm \n, ",
                            plots[i]))
        for (i in 2:output$iter) {
            points(jitter(rep(i, nrow(convergeData))), convergeData[, 2, i], type = "p",
                   pch = 19, col = alpha("blue", 0.25))
        }

        #best AIC per iteration
        lines(1:output$iter, sapply(1:output$iter,
                                   function(x) convergeData[1, 2, x]), type = "l",
             xlab = "iterations", ylab = "AIC", col = "red", lwd = 2,
             main = paste("Best AIC per iteration \n objFun = AIC, lm \n",
                          plots[i]))
        lines(1:output$iter, sapply(1:output$iter,
                                   function(x) mean(convergeData[, 2, x])),
                                   lty = 2, lwd = 2, col = "red")


    #dev.off()
}
stop()
#tryCatch((tryCatch(match.fun(he))),
#warning = function(w) {print(paste("The user specified function does not exist")); }, error = function(e) {print(paste("The user specified function does not exist")); stop(e) })

