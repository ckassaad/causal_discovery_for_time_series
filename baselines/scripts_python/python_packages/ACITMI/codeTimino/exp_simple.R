# Copyright (c) 2010-2012  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

start_time <- Sys.time()

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

# c1 <- 0.3 
# c2 <- 1
# 
# x <- rnorm(n) * sqrt(c1) * sqrt(1 / 252.)
# y <- rnorm(n) * exp(c2)
# w <- x+y

# "fork", "v_structure", "cycle", "diamond", "hidden", "complex"
struct_name <- "complex"
setwd('/home/kassaad/Documents/Codes/R - codes/simulated_ts_data/')
filenames <- list.files(struct_name, pattern="*.csv")
for (i in 1:length(filenames)){
  data = read.csv(paste(struct_name,filenames[i],sep="/"))
  data$X <- NULL
  # traints_linear, traints_gam or traints_gp
  # indtestts_hsic or indtestts_crosscov
  unit_graph <- timino_dag(data, alpha = 0.05, max_lag = 5, model = traints_linear, indtest = indtestts_crosscov, output = TRUE)
  unit_graph[is.na(unit_graph)] <- 3
  path<-paste(struct_name,filenames[i],sep="/results/res_")
  write.csv(unit_graph,path)
}


# g <- granger_dag_pairwise_nl(x,y,w)
# print(g)

show("====")
show("DONE")
show("====")

show("true summary time graph:")
show(cbind(c(0,0,0),c(0,0,0),c(1,1,0)))
#show(cbind(c(0,0,0),c(1,0,0),c(1,0,0)))

show("estimated summary time graph:")
show(d)

end_time <- Sys.time()

print(start_time)
print(end_time)
print(end_time - start_time)


setwd('/home/kassaad/Documents/Codes/Causality-time-series/data/dairymarkets')
b = read.csv('Butter.csv')
c = read.csv('Cheese.csv')
m = read.csv('Milk.csv')
start_b = which(b$Date=='2008-01-31')
end_b = which(b$Date=='2018-12-31')
start_c = which(c$Date=='2008-01-31')
end_c = which(c$Date=='2018-12-31')
start_m = which(m$Date=='2008-01-31')
end_m = which(m$Date=='2018-12-31')
b = b$Values[start_b:end_b]
c = c$Values[start_c:end_c]
m = m$Values[start_m:end_m]
length(c) ==length(b)
length(c) ==length(m)

data <- cbind(b, c)
data <- cbind(data, m)

start_time <- Sys.time()
d <- timino_dag(data, alpha = 0.05, max_lag = 10, model = traints_linear, indtest = indtestts_crosscov, output = TRUE)
end_time <- Sys.time()

print(start_time)
print(end_time)
print(end_time - start_time)

# length(x)
# length(x*y)
# z <- x*y
# a <- acf(cbind(x,y))
# c1 <- rev(a$acf[,,1][,2])
# c2 <- a$acf[,,2][,1]
# length(c(c1,c2[-1]))
# c <- ccf(x,y)
# c$acf
# c$type
# c$snames
# 
# 
# mean(x)
# mean(y)
# var(x)
# var(y)
# 
# round(10*log10(1000000/2))