
# STARTING UP:
#
# 1. Make sure you have installed the fastICA package as in:
#    > install.packages('fastICA') # if not already done
#
# 2. Make sure you have 'sourced' all files in this directory and
#    in ../common/
#
# 3. Now you can run the demo using e.g.
#    > rlingamdemo(2)

# this gets called by rlingamdemo() from prompt
rlingamdemo <- function( task=1 ) {
    
  # Use the fastICA R package
  library('fastICA')

  #------------------------------------------------------------------------
  # Task 1: simple three-variable test
  #------------------------------------------------------------------------

  if (task==1) {

    x <- matrix(rnorm(10000)^3,1,10000)
    y <- 0.3*x + matrix(rnorm(10000)^3,1,10000)
    z <- -2*y + matrix(rnorm(10000)^3,1,10000)
    D <- rbind(z,y,x)
    res <- lingam(D)    
    print(res)
    plotgraph(res$Bpruned)
    
  } 

  #------------------------------------------------------------------------
  # Task 2: randomly choose DAG and params, simulate data, and estimate
  #------------------------------------------------------------------------

  if (task==2) {

    # Creates a random graph, generates data, then estimates the graph 
    simpletestlingam()
	
  } 

}

