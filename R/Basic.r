



# X is the dataset. Row is prediction variables and columns are prediction variables
initialGA <- function(X) {
  d<-dim(X)[2]
# In the following step, I choose initial number, the most recent literacture said it should be 10 times the number of variable
  if(d<20){
    P <- 10*d
  }else{
    P <- sample(d:d*2,1)
  }
  #randomly sample genes
  sample <- sample(c(0,1), replace=TRUE, size=ceiling(1.2*d*P))
  #To make sure that each variable will exist in at least one chrome, but not all.
  sample <- c(rep(0, d-1), rep(1, d), 0, sample)
  x <- seq_along(sample)
  #create a first generation
  firstgeneration<-split(sample, ceiling(x/d))
  fgmatrix<-matrix(unlist(unique(firstgeneration)[1:P]), ncol = dim(X)[2], byrow = TRUE)
  fgmatrix <- fgmatrix[apply(fgmatrix[,-1], 1, function(x) !all(x==0)),]
  return(fgmatrix)
}






#Finding the best 5% AIC and Storing all AIC in a matrix
AIC<-function(fgmatrix, Y){
dim<-dim(fgmatrix)[1]
run.aic <- vector()
#Find AIC
for(i in 1:as.integer(dim)){
run.vars = X[, fgmatrix[i,]==1]
g = lm(Y~., X[fgmatrix[i,]==1])
run.aic[i] = extractAIC(g)[2]


}

# rank function gives the rank in ascending order. 
#The number with rank 1 is the smallest in the list. 
#Here we want the smallest AIC to be the higest ranked number in the list.
return(run.aic)
}  



ChildrenCreating <- function(fgmatrix, runs.aic){
# rank function gives the rank in ascending order. 
#The number with rank 1 is the smallest in the list. 
#Here we want the smallest AIC to be the higest ranked number in the list.
dim = dim(fgmatrix)[1]
r = rank(-runs.aic)
#Assign Probability according to rank
phi = 2*r/(dim*(dim+1))
#Create A matrix for next generations
NextGeneration <- matrix(0, dim(fgmatrix)[1], dim(fgmatrix)[2])
i = 1
while(i<=dim(fgmatrix)[1]){
#Crossover and Mutation
#Step1: Create Parents
#Parent1
P1<-sample(1:dim,1,prob=phi)
parent.1 = fgmatrix[P1,]
P1r<-r[P1]
#Parent2
P2<-sample(1:dim,1,prob=phi)
parent.2 = fgmatrix[P2,]
P2r<-r[P2]
#Do cross over.
P<-parent.1 *P1r/(P1r+P2r) + parent.2 *P2r/(P1r+P2r)
#Change it to dimension 1
#To do mutation
#Store children
NextGeneration[i,]<-rbinom(10, 1, P)
i = i + 1
}
return(NextGeneration)
}




#In the following function, we select the next generations
selection <- function(fgmatrix, child){
  #old dimension
  oldim<-dim(child)[1]
  #Before Selection
fgmatrixO <- rbind(child,fgmatrix)
fgmatrixO <- fgmatrixO[apply(fgmatrixO[,-1], 1, function(x) !all(x==0)),]
fgmatrixO <- unique(fgmatrixO)
#Selection

runs.aic<-AIC(fgmatrixO, Y)
#Give the thereshold
threshold<- sort(runs.aic)[oldim]
newfgmatrix<-fgmatrixO[runs.aic>threshold,]
return(newfgmatrix)
}

#Toy example
X<-mtcars[2:dim(mtcars)[2]]
Y<-mtcars$mpg

#Y is the variable u wanna predict
#X is the variable u can use to predict
ga<- function(Y, X){
  #Step 1 initialization
  fgmatrix<-initialGA(X)
  #consider the scope of top 5 percent
  #create a matrix to store top 5 percent AIC
  best.aic <- matrix(0,1,floor(dim(fgmatrix)[1]/20))
  #Step 2 evaluate aic
  runs.aic<-AIC(fgmatrix, Y)
  best.aic <- rbind(best.aic,sort(runs.aic)[1:floor(dim(fgmatrix)[1]/20)])
  i = 1
  tol = .Machine$double.eps ^ 0.5
  while(!isTRUE(all.equal(best.aic[i+1,], best.aic[i,], tolerance = tol))){
  #Step 3 create children
  child<-ChildrenCreating(fgmatrix,  runs.aic)
  #Step 4 select next generation
  nextgeneration<-selection(fgmatrix, child)
  #evaluate aic
  runs.aic<-AIC(nextgeneration, Y)
  best.aic<-rbind(best.aic, sort(runs.aic)[1:floor(dim(fgmatrix)[1]/20)])
  i = i + 1
  #print(best.aic)
  }
  return(best.aic)
}

