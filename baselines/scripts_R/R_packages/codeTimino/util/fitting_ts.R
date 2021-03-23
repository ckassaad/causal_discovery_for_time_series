# Copyright (c) 2010-2013 Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
#
# This file contains two kinds of functions
#
# =========
# 1. Many instances of traints_***
# =========
# traints_*** <- function(x, Y, pars)
#
# Fits a time series model for x using Y (optional)
# 
# INPUT:
#   x         vectors of target time series (N data points, 1 dimension)
#   Y         Nxd matrix of explanatory time series (N data points, d dimensions)
#   pars      list containing parameters of the regression method
#	$maxOrder	
#
#
# OUTPUT:
#   result    list with the result of the regression
#      $model        list of learned model, e.g.
#          $order        order of the fitted process
#          $degrf        degree of freedom
#      $Xfit         fitted outputs for training inputs according to the learned model
#      $residuals    noise values (e.g., residuals in the additive noise case)
#
#
# =========
# 2. One instance of traints_model
# =========
# traints_model <- function(f, x, Y, pars)
# 
# does the same as above, but takes the time series fitting method as input.


omitGP <- FALSE # some users have problems with the gptk package 


# =========
# 1. Many instances of traints_***
# =========

# ===
# gams 
# ===
# todo: gams choose maxOrder as the fixedOrder 
traints_gam <- function(x, Y, pars,plot = FALSE)
{
    if(dim(as.matrix(Y))[1] == 0)
    {
        mod <- gamts.fit2(x,pars$maxOrder,plot)
    }
    else
    {
        mod <- gamts.fit(x,Y,pars$maxOrder,plot)
    }
    
    if(length(x) == 0)
    {
        show("Warning. Is the time series in the correct format? (length == 0)")
    }
    model <- list(order = mod$order, degrf = mod$degrf)
    result <- list(residuals = mod$resid, model = model)
}

# fits X_t using X_t-k 
gamts.fit2<-function(x, max_order, plot = FALSE, order.fixed=FALSE, output = FALSE)
{
    bestaic <- 0
    bestorder <- 0
    for(order in 1:max_order)
    {
        if(!order.fixed | order == max_order)
        {
            datamatrix<-matrix(0,length(x)-order,1+order)
            for(j in 0:order)
            {
                datamatrix[1:(length(x)-order),j+1]<-x[(-j+1+order):(length(x)-j)]
            }
            dat<-data.frame(datamatrix)
            labs<-"X1~"
            if(order>1)
            {
                for(i in 2:(order))
                {
                    labs<-paste(labs,"s(X",i,")+",sep="")
                }
            }
            labs<-paste(labs,"s(X",order+1,")",sep="")
            mod_gam<-gam(formula=formula(labs),data=dat)
            if(bestaic == 0 || mod_gam$aic < bestaic)
            {
                if(output)
                {
                    show(sprintf("best aic so far: %.3f. For order %i, the aic is %.3f", bestaic, order, mod_gam$aic))
                }
                bestaic <- mod_gam$aic
                bestmodel <- mod_gam
                bestorder <- order
            }
        }
    }
    res<-rep(0,length(x))
    res[(length(x)-length(bestmodel$residuals)+1):length(x)]<-bestmodel$residuals
    output<-list(resid=res, order=bestorder, degrf=bestmodel$df.residual)
}

#fits X_t using X_t-k and Y_t-k
gamts.fit<-function(x, y, max_order, plot = FALSE, order.fixed=FALSE)
{
    bestaic <- 0
    bestorder <- 0
    ifelse(length(y)==length(x),
           y<-matrix(y,length(y),1),aaa<-0)
    pdim<-dim(y)
    for(order in 1:max_order)
    {
        if(!order.fixed | order == max_order)
        {
            datamatrix<-matrix(0,length(x)-order,1+order+order*pdim[2])
            for(j in 0:order)
            {
                datamatrix[1:(length(x)-order),j+1]<-x[(-j+1+order):(length(x)-j)]
            }
            for(i in 1:pdim[2])
            {
                for(j in 1:order)
                {
                    datamatrix[1:(length(x)-order),(i-1)*order+order+1+j]<-y[(-j+1+order):(length(x)-j),i]
                }
            }
            dat<-data.frame(datamatrix)
            labs<-"X1~"
            for(i in 2:(pdim[2]*order+order))
            {
                labs<-paste(labs,"s(X",i,")+",sep="")
            }
            labs<-paste(labs,"s(X",pdim[2]*order+order+1,")",sep="")
            mod_gam<-gam(formula=formula(labs),data=dat)
            if(bestaic == 0 || mod_gam$aic < bestaic)
            {
                bestaic <- mod_gam$aic
                bestmodel <- mod_gam
                bestorder <- order
            }
        }
    }
    res<-rep(0,length(x))
    res[(length(x)-length(bestmodel$residuals)+1):length(x)]<-bestmodel$residuals
    output<-list(resid=res, order=bestorder, degrf=bestmodel$df.residual)
}




# ===
# gp 
# ===
if(!omitGP)
{
    library(kernlab)
    library("gptk")
    
    gp_regression <- function(X,y, pars=list())
    {
        options=gpOptions("ftc")
        options$kern$comp=list("rbf","white")
        #options$learnScales=TRUE
        model<-gpCreate(dim(X)[2],1,X,y,options)
        y2<-gpOut(model,X)
        model$Yfit<-y2
        model$residuals<-y-y2
        return(model)
    }
    
    # todo: gp choose maxOrder as the fixedOrder 
    traints_gp <- function(x, Y, pars,plot = FALSE)
    {
        if(dim(as.matrix(Y))[1] == 0)
        {
            mod <- gpts.fit2(x,pars$fixedOrder)
        }
        else
        {
            mod <- gpts.fit(x,Y,pars$fixedOrder)
        }
        model <- list(order = mod$order)
        result <- list(residuals = mod$resid, model = model)
    }
    
    #fits X_t using X_t-k 
    gpts.fit2<-function(x,max_order,plot = FALSE)
    {
        datamatrix<-matrix(0,length(x)-max_order,max_order)
        for(j in 1:max_order)
        {
            datamatrix[1:(length(x)-max_order),j]<-x[(-j+1+max_order):(length(x)-j)]
        }
        target<-x[(1+max_order):length(x)]
        mod_gp<-gp_regression(datamatrix,matrix(target))
        res<-rep(0,length(x))
        res[(length(x)-length(mod_gp$residuals)+1):length(x)]<-mod_gp$residuals
        output<-list(resid=res, order=max_order)
    }
    
    #fits X_t using X_t-k and Y_t-k
    gpts.fit<-function(x,y,max_order,plot= FALSE)
    {
        if(length(y)==length(x))
        {
            y <- matrix(y,length(y),1)
        }
        pdim <- dim(y)
        datamatrix <- matrix(0,length(x)-max_order,1+max_order+max_order*pdim[2])
        for(j in 1:max_order)
        {
            datamatrix[1:(length(x)-max_order),j]<-x[(-j+1+max_order):(length(x)-j)]
        }
        for(i in 1:pdim[2])
        {
            for(j in 1:max_order)
            {
                datamatrix[1:(length(x)-max_order),(i-1)*max_order+max_order+j]<-y[(-j+1+max_order):(length(x)-j),i]
            }
        }
        target<-x[(1+max_order):length(x)]
        mod_gp<-gp_regression(datamatrix,matrix(target))
        res<-rep(0,length(x))
        res[(length(x)-length(mod_gp$residuals)+1):length(x)]<-mod_gp$residuals
        output<-list(resid=res, order=max_order)
    }
}


# ===
# linear 
# ===
traints_linear <- function(x, Y, pars,plot = FALSE)
{
    if(dim(as.matrix(Y))[1] == 0)
    {
        mod <- ar(x,aic=TRUE,order.max=pars$maxOrder)
        resi <- mod$resid
    }
    else
    {
        mod<-ar(cbind(x,Y),aic=TRUE,order.max=pars$maxOrder)
        resi <- mod$resid[,1] 
    }
    model <- list(order = max(mod$order,1))
    result <- list(residuals = resi, model = model)
}





# =========
# 2. train_model
# =========

traints_model <- function(f,x,Y,pars = list(),plot = FALSE)
{
    result <- f(x,Y,pars,plot)
    # sanity check (can be removed)
    if (length(result$residuals) != length(x))
    {
        stop("The length of time series and residuals do not coincide. Check!")
    }
    
    return(result)
}
