# Copyright (c) 2010-2013 Jonas Peters  [jonas.peters@tuebingen.mpg.de]
# All rights reserved.  See the file COPYING for license terms. 
#
#
# This file contains independence tests for time series and one wrapper
#


# wrapper  
indtestts <- function(f,z,r,alpha,max_lag,plotit=FALSE)
{
    if(plotit == TRUE & 1 == 0)
    {
        par(mfrow = c(2, max_lag + 1))
        len <- dim(as.matrix(r))[1]
        for(i in (-max_lag):(max_lag))
        {
            zz<-as.matrix(z)
            zz<-zz[(max(0,i)+1):(len+min(i,0)),]
            resres<-r[(max(0,i)+1-i):(len+min(i,0)-i)]
            plot(zz,resres)  
        }
        readline(prompt = "Plots show the scatter plots between residuals and time series. Press <Enter> to continue...")        
        dev.off()
    }
    
    if(plotit)
    {
        plot.new()
        # readline(prompt = "Plot shows the cross correlation. Press <Enter> to continue...")        
        #par(mfrow = c(2,1))
        ccf(z,r,lag.max=max_lag,type = "correlation",plot = FALSE)
        acf(r)
        readline(prompt = "Plots show cross and auto correlation of the residuals. Press <Enter> to continue...")        
        #dev.off()
    }
    
    result <- f(z,r,alpha,max_lag,plotit)
    if(alpha == 0)
    {
        result$crit.value <- Inf
    }
    return(result)
}







# works only if eps is white noise!! (i.e. acf(k) = 0 for k!=0)
# todo: make it a function from max_lag_min and max_lag_max  
cross_cov_test <- function(x,eps,alpha,max_lag,plotit=FALSE)
{
    corr1 <- ccf(x,eps,lag.max=max_lag,type = "correlation",plot = FALSE)
    T <- max(abs(corr1$acf))
    sigma <- matrix(0,2*max_lag,2*max_lag)
    # estimate cov of x-cor.
    # Theorem 11.2.3 in brockwell and davis: "bartletts formula"
    # H_0: rho_{12} == 0 => non-zero summands only for j = k - h 
    acr <- acf(x,lag.max=2*max_lag,type = "correlation",plot=FALSE)
    for(i in 1:(2*max_lag))
    {
        for(j in 1:i)
        {
            sigma[i,j]=acr$acf[1+(i-j)]
            sigma[j,i]=acr$acf[1+(i-j)]
        }
    }
    sigma <- sigma/length(x)
    R <- chol(sigma)
    num_simulations<-20000
    z <- matrix(rnorm(2*max_lag*num_simulations),num_simulations,2*max_lag)
    z <- z%*%R
    maxz <- apply(abs(z),1,max)
    maxzorder <- sort(maxz)
    quan <- maxzorder[ceiling(num_simulations-alpha*num_simulations)]
    pval <- sum(maxzorder>T)/num_simulations
    resu <- list(statistic = T, crit.value = quan, p.value = pval)
    return(resu)
}


# bonferroni correction included
indtestts_crosscov <- function(z,r1,alpha,max_lag,plotit=FALSE)
{
    z<-as.matrix(z)
    pdim2<-dim(z)
    Tvecquanvec<-c()
    #    Tvecquanvec<-matrix(0,pdim2[2],3)
    Tvec<-matrix(0,pdim2[2],1)
    quanvec<-matrix(0,pdim2[2],1)
    for(i in 1:pdim2[2])
    {
        Tvecquanvec[[i]]<-cross_cov_test(z[,i],r1,alpha/pdim2[2],max_lag,plotit)
        Tvec[i,1]<-Tvecquanvec[[i]]$statistic
        quanvec[i,1]<-Tvecquanvec[[i]]$crit.value
    }
    bb <- which.max(Tvec-quanvec)
    T <- Tvec[bb]
    quan <- quanvec[bb]
    pval <- Tvecquanvec[[bb]]$p.value*pdim2[2]
    resu <- list(statistic = T, crit.value = quan, p.value = pval)
    return(resu)
    
}









#tests all z against each residuals (bonferroni with 2*max_lag+1)
#make it a function from max_lag_min and max_lag_max  
indtestts_hsic <- function(z,res,alpha,max_lag,plotit=FALSE)
{
    #ifelse(is.matrix(res)==TRUE,
    #    len<-dim(res)[1],
    #    len<-length(res))
    ifelse(is.matrix(z)==TRUE,
           len<-dim(z)[1],
           len<-length(z))
    jj<-0
    T<-c(0)
    quan<-c(0)
    pval<-c(0)
    for(i in (-max_lag):(max_lag))
    {
        jj<-jj+1
        ifelse(is.matrix(z),
               zz<-z,
               zz<-matrix(z))
        zz<-zz[(max(0,i)+1):(len+min(i,0)),]
        # is this correct?
        # zz<-z[(max(0,i)+1):(len+min(i,0))]
        resres<-res[(max(0,i)+1-i):(len+min(i,0)-i)]
        Tquani<-indtestHsic(zz,resres,alpha/(2*max_lag+1), pars = list(method = "ExactFastTrace"))
        T[jj] <- Tquani$statistic
        quan[jj] <- Tquani$crit.value
        pval[jj] <- Tquani$p.value
    }
    #bb<-which.max(T-quan)
    bb<-which.min(pval)
    T_final<-T[bb]
    quan_final<-quan[bb]
    pval <- pval[bb] * (2*max_lag+1)
    if(plotit)    
    {
        par(mfrow = c(1,1))
        i <- -max_lag -1 + bb
        ifelse(is.matrix(z),
               zz<-z,
               zz<-matrix(z))
        zz<-zz[(max(0,i)+1):(len+min(i,0)),]
        resres<-res[(max(0,i)+1-i):(len+min(i,0)-i)]
        plot(resres,zz)
        
        show(sprintf("shift of residuals: %i. p-val of cor=%.3e. p-val of HSIC=%.3e", i, cor.test(zz,resres)$p.value, indtestHsic(zz,resres,0.9)[3]))
        readline(prompt = "This is the scatter plot between the residuals and the shifted time series.")
    }
    resu <- list(statistic = T_final, crit.value = quan_final, p.value = pval)
    return(resu)
    
}


#tests all z against a matrix of residuals (no bonferroni)
indtestts_hsic2 <- function(z,res,alpha,max_lag,plotit=FALSE)
{
    ifelse(is.matrix(res)==TRUE,
           len<-dim(res)[1],
           len<-length(res))
    quan<-c(0)
    resres<-matrix(0,len-2*max_lag,2*max_lag+1)	
    for(i in 1:(2*max_lag+1))
    {
        # show(c(i,min(len-2*max_lag-1+i,len)))
        resres[,i]<-res[i:min(len-2*max_lag-1+i,len)]
    }
    ifelse(is.matrix(z),
           zz<-z[(max_lag+1):(len-max_lag),],
           zz<-z[(max_lag+1):(len-max_lag)])
    Tquan1 <- indtestHsic(zz,resres,alpha)
    pval <- -1
    resu <- list(statistic = Tquan1$statistic, crit.value = Tquan1$crit.value, p.value = pval)
    return(resu)
    
}


#tests each z against a matrix of residuals (bonferroni with dim(z))
indtestts_hsic3 <- function(z,res,alpha,max_lag,plotit=FALSE)
{
    ifelse(is.matrix(res)==TRUE,
           len<-dim(res)[1],
           len<-length(res))
    T<-c(0)
    quan<-c(0)
    resres<-matrix(0,len-2*max_lag,2*max_lag+1)	
    if(is.matrix(z))
    {
        zz<-z[(max_lag+1):(len-max_lag),]
        dimz<-dim(z)[2]
    }
    else
    {
        zz<-matrix(z[(max_lag+1):(len-max_lag)])
        dimz<-1
    }
    
    for(i in 1:(2*max_lag+1))
    {
        # show(c(i,min(len-2*max_lag-1+i,len)))
        resres[,i]<-res[i:min(len-2*max_lag-1+i,len)]
    }
    for(jj in 1:dimz)
    {
        Tquani<-indtestHsic(zz[,jj],resres,alpha/dimz)
        T[jj]<-Tquani$statistic
        quan[jj]<-Tquani$crit.value
    }
    bb<-which.max(T-quan)
    pval <- -1
    resu <- list(statistic = T[bb], crit.value = quan[bb], p.value = pval)
    return(resu)
    
}

