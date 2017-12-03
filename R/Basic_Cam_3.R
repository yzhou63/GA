
rm(list=ls())

##############
#parent function for package
##############
GAfun <- function(Y, X, iter, objFun = c("AIC", "BIC", "logLike", "user"), 
                  family = "gaussian", 
                  crossMeth = c("method1", "method2", "method3"), 
                  minimize = TRUE, parallel = NULL) {
    
    #input objects
    if(missing(iter)) {iter <- 100}
    if(missing(objFun)) {objFun <- "AIC"}
    if(missing(family)) {family <- "gaussian"}
    if(missing(minimize)) {minimize <- TRUE}
    if(missing(parallel)) {parallel <- FALSE
    } else if (parallel == TRUE) {
        require(parallel)
        require(foreach)
        require(doParallel)
        
        parallel <- TRUE
    }

    #error checking
    
    
    #function objects
    
    
    ##########
    # Perform genetic algorithm
    #########
    
    #step 1: inital generation
        #Generate initial population
        cat("1. Create founders")
        
        generation_t0 <- initialGA(X)
    
    #Step 2. Evaluate Fitness of inital pop
        cat("\n2. Eval founder fitness")
    
        objFunOutput_t0 <- evalGeneration(generation_t0, Y, X)
    
        #create array to store fitness data for each iteration
        convergeData <- array(dim = c(nrow(generation_t0), 2, iter))
    
        #save evaluation data
        convergeData[, , 1] <- objFunOutput_t0[
                                order(objFunOutput_t0[, 2], decreasing = T), 
                                c(1, 3)] #save top 5 chromosomes
    
        
        
    #Step 3. loop through successive generations until either:
        #1. finish default or user specified # of iterations
        #2. convergence to specific optimum
        tol = 0.0001#sqrt(.Machine$double.eps) #tolerance for converage check
        converged <- 0
        
        cat("\n3. Begin breeding \n Generations: ")
        
        for (i in 1:iter) {
            
            #breeding: create children
            if (i == 1) { generation_t1 <- generation_t0 
                            objFunOutput_t1 <- objFunOutput_t0 }
            
            generation_t1 <- createNextGen(generation_t1, objFunOutput_t1, iter)
            
            #eval fitness of children
            objFunOutput_t1 <-evalGeneration(generation_t1, 
                                                   Y, X)
            
            #store fitness data
            convergeData[, , i] <- 
                objFunOutput_t1[order(objFunOutput_t1[, 2], 
                                            decreasing = T), c(1, 3)]
            
            #check convergence
                # after 10 iterations, check if new value is different from moving 
                #average of prev 10 iterations
                #this isn't working great right now
            if (i > 20) {
                if(isTRUE(all.equal(convergeData[1:5, 2, i-1], 
                                convergeData[1:5, 2, i],
                                tolerance = tol))) {
                    converged <- converged + 1 
                }
            }
            if (converged >= 10) { 
                cat("\n#### Converged! ####") 
                break }
            
            cat(i,"-", sep = "")
        }
        

    #process output
    bestModel <- generation_t1[convergeData[1, 1, i], ]
    value <- round(convergeData[1, 2, i], 4)
    if(i < iter) {converged <- "Yes"} else {converged <- "No"}
    
    output <- list("BestModel" = bestModel, 
                   objFun = c(objFun, value), 
                   iter = i,
                   converged = converged,
                   convergeData = convergeData)
    return(output)
}


##############
# subfunctions
##############

#initiative founding chromosomes
initialGA <- function(X) {
    
    #d: #number of genes/variables in design matrix
    #P: #number of parent chromosomes
    #geneSample: Random sample of genes/vars for initial parent dataet
    #firstgeneration: 
    #generation_t:
    
    #get number of genes/variables in dataset
    C <- dim(X)[2]

    # I choose initial parent chromosomes, the most recent 
    # literacture said it should be 10 times the number of genes/vars
    # Cam: reading says choose P to satisfy C ≤ P ≤ 2C
    #P <- 10 * C 
    P <- 2 * C
    if (P > 100) {P <- 100}
    
    #randomly generate parent generate
    geneSample <- sample(c(0, 1), replace = TRUE, size = ceiling(1.2 * C * P))
        #question: why 1.2?
    
    #update geneSample to make sure that each gene/variable 
    #will exist in at least one chrome, but not all.
    geneSample <- c(rep(0, C - 1), rep(1, C), 0, geneSample)
  
    #create a first generation
    x <- seq_along(geneSample)
    firstGen <- split(geneSample, ceiling(x / C))
    
    generation_t0 <- matrix(unlist(unique(firstGen)[1:P]), 
                       ncol = C, byrow = TRUE)
    
    generation_t0 <- generation_t0[apply(generation_t0[, -1], 1, 
                                     function(x) !all(x == 0)), ]
    

    return(generation_t0)
}

#eval geneartion using lm/glm and objective function
#default is AIC
evalGeneration <- function(generation_t0, Y, X) {
    
    family <- get("family", mode = "any", envir = parent.frame())
    objFun <- get("objFun", mode = "any", envir= parent.frame())
    parallel <- get("parallel", mode = "any", envir= parent.frame())
    minimize <- get("minimize", mode = "any", envir= parent.frame())
    
    #family <- "gaussian"
    #objFun <- "AIC"
    #minimize <- TRUE
    #parallel <- FALSE
    
    #input objects
        #generation_t: chromosomes to be evaluated
        #Y: dependant var
        #X: covariates
        #family: regression distribution
        ##optim: min or max obj func
        #objFun: objective fun
    
    #function objects
        #mod: glm or lm object
        #P: number of chromosomes
    
    #number of parent chromosomes
    P <- dim(generation_t0)[1]
    
    ######
    #evaluate each chromosome with selected objective function
    ######
    
    
    #foreach options
    mcoptions <- list(preschedule=T)
    
    if(parallel == TRUE) {
        nCores <- detectCores() - 1
        registerDoParallel(nCores)
    } else {nCores <- 1}
    
    #initialize cores
    objFunOutput <- foreach(i = 1:P,
                      .options.multicore = mcoptions,
                      .combine = c,             # how to combine results
                      .errorhandling=c("pass"),
                      .verbose = F) %dopar% {
              
            if (family == "gaussian") {
                  mod <- lm(Y ~ X[, generation_t0[i, ] == 1])
                } else {mod <- glm(Y ~ ., X[generation_t0[i, ] == 1], 
                                 family = family)}
              
              #evaluate with objective function
              if (objFun == "AIC") {                 #AIC
                  return(extractAIC(mod)[2])
                } else if (objFun == "BIC") {      #BIC
                  return(BIC(mod))
                    } else if (objFun == "logLike") {  #loglikelihood
                  return(logLik(mod)[1])
                    } else {
                    return(objFun(mod))}    #user defined function
                      }
    
    #rank obj function output
    if (objFun %in% c("AIC", "BIC") & !objFun %in% c("loglike")) {
        r <- rank(-objFunOutput, na.last = TRUE, ties.method = "first")
    } else if (objFun == "logLike" | minimize == FALSE) {
        r <- rank(objFunOutput, na.last = TRUE, ties.method = "first")
        }
    return(cbind(chr = 1:P, rank = r, objFunOutput))
}

#create next generation function
createNextGen <- function(generation_t0, objFunOutput_t0, iter) {
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
    require(scales)
    #plot(phi, score)
    #legend("bottomleft", c("parents", "sampled parents"), pch = c(21, 19))
    
    #Create matrix for next generation
    generation_t1 <- matrix(NA, dim(generation_t0)[1], dim(generation_t0)[2])
    
    #keep high rank parents and put in new generation
    goodGeneration_t0 <- generation_t0[rank > quantile(rank, probs = 0.95), ]
    generation_t1[1:dim(goodGeneration_t0)[1], ] <- goodGeneration_t0
    
    #########
    #Selection, Crossover, and Mutation
    #########
    
    i <- 1#dim(goodGeneration_t0)[1] + 1 #initialize while loop
    while(i <= dim(generation_t1)[1]) {
            
        #SELECTION: select Parents
        #method 1: both parents selected by rank
            parentInd <- sample(1:P, 2, prob=phi, replace = F) #this is better than method 2
                       
        #method 2: first parent by rank, second random
            #parentInd <- c(sample(1:P, 1, prob=phi, replace = F),
            #               sample(1:P, 1, replace = F))
        
        parents <- generation_t0[parentInd, ]
        parentRank <- rank[parentInd]
        parentScore <- score[parentInd]
        
        #update plot with sampled parents
        #points(phi[parentInd], score[parentInd], pch = 19,
               #col = alpha("black", 0.5))
        
        #CROSSOVER and MUTATION parent to  create child
         
        if (crossMeth == "method1") {
            
            #METHOD 1:
                #crossover: method upweights parent with higher rank high
                child <- parents[1, ] * parentRank[1] / 
                        (parentRank[1] + parentRank[2]) + 
                                parents[2, ] * parentRank[2]  / 
                        (parentRank[1] + parentRank[2])
                #mutation: this method has A LOT of mutation
                child <- rbinom(C, 1, prob = child)
        
        } else if (crossMeth == "method2") {
            
            #METHOD 2: 
                #crossover: two crossover points, non-overlapping
                cross1 <- sample(2:(C/2-1), 1)
                cross2 <- sample((C/2+1):(C-1), 1)
                child <- c(parents[1, 1:cross1], 
                            parents[2, (cross1+1):(C/2)], 
                            parents[1, (C/2+1):cross2],
                            parents[2, (cross2+1):C])
                #mutation:
                child <- abs(round(child, 0) - 
                                rbinom(C, 1, prob = 1 / (P * sqrt(C))))
                
        } else if (crossMeth == "method3") {    
            
            #METHOD 3
                #randomly samples non-concordant variables
                #between parents, slightly upweights parent selected
                #by prob. proportional to rank
                #this provide mutation as well
                child <- parents[1, ]
                child[parents[1, ] != parents[2, ]] <- 
                    rbinom(sum(parents[1, ] - parents[2, ] != 0), 1, 
                         prob = parentRank[1] / (parentRank[1] + parentRank[2]))
        }
        
        #check if duplicate child, and for all non-zero genes
        if (all(!rowSums(generation_t1 == child, na.rm = T) == C) & 
            (sum(child) > 0)) {
            #not dup: add child to new generation    
            generation_t1[i, ] <- child 
            #update counter
            i <- i + 1
        } 
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
colnames(bogus) <- paste("bogus", 5+(1:ncol(bogus)), sep = "")
x <- cbind(sim$x, bogus)
y <- sim$y
dim(x);length(y)

##############
#test function with simulated data
##############

#run against three different crossover/mutation methods
output.meth1 <- GAfun(y, x, iter = 500, objFun = "AIC", parallel = T,
                crossMeth = "method1") #with defaults
output.meth2 <- GAfun(y, x, iter = 500, objFun = "AIC", parallel = T,
                crossMeth = "method2") #with defaults
output.meth3 <- GAfun(y, x, iter = 500, objFun = "AIC", parallel = T,
                crossMeth = "method3") #with defaults


test <- rbind(c(output.meth1$objFun, output.meth1$iter),
      c(output.meth2$objFun, output.meth2$iter),
      c(output.meth3$objFun, output.meth3$iter))

colnames(test) <- c("ObjFun", "value", "iterations")
rownames(test) <- c("cross.meth1", "cross.meth2", "cross.meth3")
test

#parent selection random by p(given rank)
    #ObjFun value      iterations
    #cross.meth1 "AIC"  "154.4266" "35"      
    #cross.meth2 "AIC"  "153.0102" "41"      
    #cross.meth3 "AIC"  "152.2202" "39" 

#first parent selected by rank, second random
    #ObjFun value      iterations
    #cross.meth1 "AIC"  "152.1328" "39"      
    #cross.meth2 "AIC"  "159.7817" "122"     
    #cross.meth3 "AIC"  "164.4339" "95" 

#check variables selected
colnames(x)[output$BestModel == 1]
colnames(x)[!output$BestModel == 1]

#plot results    
plots <- paste0("output.meth", 1:3)
for (i in 1:3) {
   
    output <- eval(parse(text = plots[i]))
    convergeData <- output$convergeData
    require(scales)
    
    png(paste0("~/repos/STAT243/project/GAplots_method_0", i, ".png"), 
        width = 640, heigh = 480)
    
    par(mfrow = (c(1, 2)))
    
        #best AIC per iteration
        plot(1:output$iter, sapply(1:output$iter, 
                                   function(x) convergeData[1, 2, x]), type = "l",
             xlab = "iterations", ylab = "AIC", 
             main = paste("Best AIC per iteration \n objFun = AIC, lm \n",
                          plots[i]))
        
        #all AIC across interations
        plot(jitter(rep(1, 100)), convergeData[, 2, 1], type = "p", 
             pch = 19, col = alpha("blue", 0.1), 
             ylim = c(min(convergeData[, 2, ], na.rm = T), 
                      max(convergeData[, 2, ], na.rm = T)), 
             xlim = c(1, output$iter), xlab = "Generations", ylab = "AIC",
             main = paste("GA performance using \n objFun = AIC, lm \n, ", 
                          plots[i]))
        for (i in 2:output$iter) {
            points(jitter(rep(i, 100)), convergeData[, 2, i], type = "p",
                   pch = 19, col = alpha("blue", 0.25))
        }
        mtext(paste(text, collapse = ""), side = 3)
        
    dev.off()  
}
