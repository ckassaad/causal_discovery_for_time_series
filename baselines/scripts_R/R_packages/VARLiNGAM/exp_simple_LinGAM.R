set.seed(1)
n <- 1000
w <- rep(0,n)
x <- rep(0,n)
y <- rep(0,n)
epsw <- rnorm(n)^3
epsx <- rnorm(n)^3
epsy <- rnorm(n)^3

for(i in 3:n)
{
  x[i] <- 0.3*x[i-1]+0.5*epsx[i]
  y[i] <- 0.8*y[i-1]+0.5*epsy[i] #0.8*x[i-1]+
  w[i] <- -0.6*w[i-1]+0.8*y[i-1]+0.8*x[i-2]+0.5*epsw[i]
}


X <- cbind(x,y,w)




setwd('/home/kassaad/Documents/Codes/R - codes/VARLiNGAM')
source("sourcedir.R")
source("main1.R")

nlags = 2


X_can <- tsdata2canonicalform(X,nlags) # put data into canonical form
result <- VARLiNGAM(X_can,"ols", pruning=FALSE)

