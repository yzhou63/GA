
rm(list=ls())

##############
#parent function for package
##############

select <- function(Y, X, iter, objFun = c("AIC", "BIC", "logLik", "user"),
                  family = "gaussian",
                  crossMeth = c("method1", "method2", "method3"),
                  numChromosomes = NULL,
                  pCrossover = 1,
                  user_objFun = NULL,
                  converge = TRUE, minimize = TRUE, parallel = FALSE) {

    #input objects
    if(missing(iter)) {iter <- 100}
    if(missing(objFun)) {objFun <- "AIC"}
    if(missing(family)) {family <- "gaussian"}
    if(missing(converge)) {converge <- FALSE}
    if(missing(minimize)) {minimize <- TRUE}
    if(missing(parallel)) {parallel <- FALSE}
    if(missing(user_objFun)) {user_objFun <- NULL}
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
        #1. finish default or user specified # of iterations
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
                objFunOutput_t1 <- evaluate_Fitness(generation_t1, Y, X)

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

##############
# subfunctions
##############

#initiative founding chromosomes
generate_founders <- function(X) {

    #d: #number of genes/variables in design matrix
    #P: #number of parent chromosomes
    #geneSample: Random sample of genes/vars for initial parent dataet
    #firstgeneration:
    #generation_t:

    # number of predictors ---------------
    C <- dim(X)[2]

    # number of founders ----------------
    P <- 2 * C
    if (P > 200) {P <- 200
    } else if (P < 20) {
        P <- 20}

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

    generation_t0 <- generation_t0[apply(generation_t0[, -1], 1,
                                     function(x) !all(x == 0)), ]


    return(generation_t0)
}

evaluate_Fitness <- function(generation_t0, Y, X) {

    # get objects ----------------
    family <- get("family", mode = "any", envir = parent.frame())
    objFun <- get("objFun", mode = "any", envir= parent.frame())
    parallel <- get("parallel", mode = "any", envir= parent.frame())
    minimize <- get("minimize", mode = "any", envir= parent.frame())

    #test parms
    #family <- "gaussian"
    ##objFun <- "AIC"
    #minimize <- TRUE
    #parallel <- FALSE

    # number parent chromosomes  ----------------
    P <- dim(generation_t0)[1]

    # calculate obj function ----------------
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

    #return rankings ----------------
    return(rank)
}

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
                #crossover: method upweights parent with higher rank high
                #create one child
                childProb <- parents[1, ] * parentRank[1] /
                        (parentRank[1] + parentRank[2]) +
                                parents[2, ] * parentRank[2]  /
                        (parentRank[1] + parentRank[2])
                #mutation: this method has A LOT of mutation
                child1 <- rbinom(C, 1, prob = childProb)
                child2 <- rbinom(C, 1, prob = childProb)

        } else if (crossMeth == "method2") {

            #METHOD 2:
                #crossover: two crossover points, non-overlapping
                #creates two children
                cross1 <- sample(2:(C/2-1), 1)
                cross2 <- sample((C/2+1):(C-1), 1)
                child1 <- c(parents[1, 1:cross1],
                            parents[2, (cross1+1):(C/2)],
                            parents[1, (C/2+1):cross2],
                            parents[2, (cross2+1):C])

                child2 <- c(parents[2, 1:cross1],
                            parents[1, (cross1+1):(C/2)],
                            parents[2, (C/2+1):cross2],
                            parents[1, (cross2+1):C])

                #mutation:
                child1 <- abs(round(child1, 0) -
                                rbinom(C, 1, prob = 1 / (P * sqrt(C))))
                child2 <- abs(round(child2, 0) -
                                  rbinom(C, 1, prob = 1 / (P * sqrt(C))))

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



##############
# simulation
##############

#simulate toy dataset
library(mlbench)
n <- 100
p <- 40
sigma <- 1
set.seed(30)
sim <- mlbench.friedman1(n, sd = sigma)
colnames(sim$x) <- c(paste("real", 1:5, sep = ""),
                     paste("bogus", 1:5, sep = ""))
bogus <- matrix(rnorm(n * p), nrow = n)
colnames(bogus) <- paste("bogus", 5 + (1:ncol(bogus)), sep = "")
x <- cbind(sim$x, bogus)
y <- sim$y

#binomial
#y[y < quantile(y, probs = 0.5)] <- 0
#y[y >= quantile(y, probs = 0.5)] <- 1

#poisson
sort(runif(10, min = 0, max = 10))/10
dim(x);length(y)
dim(x)

##############
#test function with simulated data
##############


a1 <- array(NA, dim = c(5, 4, 5))
a2 <- array(rnorm(100), dim = c(5, 4, 4))
a1[, , 2:5] <- a2

#GAfun(Y = iris$Sepal.Length, X = iris[ , - 1], objFun = "AIC")
#run against three different crossover/mutation methods
output.meth1 <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                crossMeth = "method1", converge = T, family = "gaussian")
output.meth1.all <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                      crossMeth = "method1", converge = F) #with defaults
output.meth2 <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                crossMeth = "method2", converge = T) #with defaults
output.meth2.all <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                      crossMeth = "method2", converge = F) #with defaults
output.meth3 <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                crossMeth = "method3", converge = T) #with defaults
output.meth3.all <- GAfun(y, x, iter = 100, objFun = "AIC", parallel = F,
                      crossMeth = "method3", converge = F) #with defaults

test <- rbind(c(output.meth1$objFun, output.meth1$iter),
        c(output.meth1.all$objFun, output.meth1.all$iter),
        c(output.meth2$objFun,     output.meth2$iter),
        c(output.meth2.all$objFun, output.meth2.all$iter),
        c(output.meth3$objFun,     output.meth3$iter),
        c(output.meth3$objFun,     output.meth3.all$iter))

colnames(test) <- c("ObjFun", "value", "iterations")
rownames(test) <- grep("output.meth", ls(), value = T)
test

#timing
plot(diff(output.meth1$timing), type = "l", pch = 20, col = "red")
lines(diff(output.meth2.all$timing), pch = 20, col = "blue")
lines(diff(output.meth3.all$timing), pch = 20, col = "cyan")
mean(diff(output.meth1.all$timing))
mean(diff(output.meth2.all$timing))
mean(diff(output.meth3.all$timing))
sum(diff(output.meth1.all$timing))
sum(diff(output.meth2.all$timing))
sum(diff(output.meth3.all$timing))

#check variables selected
colnames(x)[output.meth1$BestModel == 1]
colnames(x)[output.meth2$BestModel == 1]
colnames(x)[output.meth3$BestModel == 1]



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
