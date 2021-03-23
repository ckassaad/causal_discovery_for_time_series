indtestHsic <- function(x,y,alpha=0.05, pars = list(method = "IncChol"))
    # Copyright (c) 2010-2013 Jonas Peters [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms.
{
    # outputs the test statistic (N*HSIC) and the critical value (according to alpha). If the test statistic is 
    # larger than the critical value, 
    # H_0 (X and Y are independent) is rejected.
    # requires
    # library(kernlab)

    
    if(is.matrix(x)==FALSE){
        x<-as.matrix(x)}
    if(is.matrix(y)==FALSE){
        y<-as.matrix(y)}
    len <- dim(x)[1]
    
    # compute distance matrices
    xnorm<-as.matrix(dist(x,method="euclidean",diag=TRUE,upper=TRUE))
    xnorm<-xnorm^2
    ynorm<-as.matrix(dist(y,method="euclidean",diag=TRUE,upper=TRUE))
    ynorm<-ynorm^2
    
    # choose median heuristic for bandwidth
    if(len>1000)
    {
        sam <- sample(1:len,1000)
        xhilf<-xnorm[sam,sam]
        yhilf<-ynorm[sam,sam]
    }
    else
    {
        xhilf<-xnorm
        yhilf<-ynorm
    }
    sigmax<-sqrt(0.5*median(xhilf[lower.tri(xhilf,diag=FALSE)]))
    sigmay<-sqrt(0.5*median(yhilf[lower.tri(yhilf,diag=FALSE)]))
    
    
    if(pars$method == "Exact" || pars$method == "ExactFastTrace")
    {
        ###
        # Compute GramMat
        ###
        ptm <- proc.time()
        KX <- exp(-xnorm/(2*sigmax^2))
        KY <- exp(-ynorm/(2*sigmay^2))
        timeGramMat <- (proc.time() - ptm)[1]
        
        ###
        # Compute HSIC
        ###
        if(pars$method == "Exact")
        {
            ptm <- proc.time()
            H<-diag(1,len)-1/len*matrix(1,len,len)
            HSIC <- 1/(len^2)*sum(diag(KX%*%H%*%KY%*%H))
            timeHSIC <- (proc.time() - ptm)[1]
        }
        if(pars$method == "ExactFastTrace")
        {
            ptm <- proc.time()
            H<-diag(1,len)-1/len*matrix(1,len,len)
            HSIC <- 1/(len^2) * sum((KX - 1/len*(KX%*%rep(1,len))%*%t(rep(1,len)))*t(KY - 1/len*(KY%*%rep(1,len))%*%t(rep(1,len))))
            timeHSIC <- (proc.time() - ptm)[1]
        }
        
        ###
        # Compute Gamma Approximation
        ###
        ptm <- proc.time()
        mux <- (sum(KX)-len)/(len*(len-1))
        muy <- (sum(KY)-len)/(len*(len-1))
        
        mean_h0 <- 1/len*(1+mux*muy-mux-muy)
        var_h0 <- (2*(len-4)*(len-5))/(len*(len-1)*(len-2)*(len-3)) * 1/((len-1)^2)*sum((KX - 1/len*(KX%*%rep(1,len))%*%t(rep(1,len)))*t(KX - 1/len*(KX%*%rep(1,len))%*%t(rep(1,len)))) * 1/((len-1)^2)*sum((KY - 1/len*(KY%*%rep(1,len))%*%t(rep(1,len)))*t(KY - 1/len*(KY%*%rep(1,len))%*%t(rep(1,len))))
        timeGamma <- (proc.time() - ptm)[1]
        
    }
    
    if(pars$method == "IncChol" || pars$method == "IncCholFastTrace")
    {
        ###
        # Compute GramMat
        ###
        ## incomplete cholesky decomposition
        ptm <- proc.time()
        LX <- inchol(x, kernel="rbfdot", kpar=list(sigma=1/(2*sigmax^2)), tol = 0.0001, maxiter = 300)
        LX <- matrix(LX,nrow=dim(LX)[1], ncol=dim(LX)[2])
        LY <- inchol(y, kernel="rbfdot", kpar=list(sigma=1/(2*sigmay^2)), tol = 0.0001, maxiter = 300)
        LY <- matrix(LY,nrow=dim(LY)[1], ncol=dim(LY)[2])
        LXc <- LX-1/len*(as.matrix(rep(1,len))%*%colSums(LX))
        LYc <- LY-1/len*(as.matrix(rep(1,len))%*%colSums(LY))
        timeGramMat <- (proc.time() - ptm)[1]
        
        #  tr( H*LX*LX'*H*LY*LY')
        # =tr( LXc*LX'*LYc*LY')
        # =tr( (LX'*LYc) * (LY'*LXc) )
        
        ###
        # Compute HSIC
        ###
        if(pars$method == "IncChol")
        {
            ptm <- proc.time()
            HSIC <- 1/(len^2)*sum(diag((t(LX)%*%LYc)%*%(t(LY)%*%LXc)))
            timeHSIC <- (proc.time() - ptm)[1]
        }
        if(pars$method == "IncCholFastTrace") # doesn't make a difference
        {
            ptm <- proc.time()
            HSIC <- 1/(len^2)*sum( (t(LX)%*%LYc) * t((t(LY)%*%LXc)))
            timeHSIC <- (proc.time() - ptm)[1]
        }
        
        ###
        # Compute Gamma Approximation
        ###
        ptm <- proc.time()
        mux <- (crossprod(colSums(LX))-len)/(len*(len-1))
        muy <- (crossprod(colSums(LY))-len)/(len*(len-1))
        
        mean_h0 <- 1/len*(1+mux*muy-mux-muy)
        var_h0 <- (2*(len-4)*(len-5))/(len*(len-1)*(len-2)*(len-3))*1/((len-1)^2)*sum(diag((t(LX)%*%LXc)%*%(t(LX)%*%LXc)))*1/((len-1)^2)*sum(diag((t(LY)%*%LYc)%*%(t(LY)%*%LYc)))
        timeGamma <- (proc.time() - ptm)[1]
    }        
    
    
    
    a <- (mean_h0^2)/var_h0
    b <- len*var_h0/mean_h0
    critical_value <- qgamma(1-alpha,shape=a,scale=b)
    p_value <- 1-pgamma(len*HSIC,shape=a,scale=b)
    resu <- list(statistic = len*HSIC, crit.value = critical_value, p.value = p_value, time = c(timeGramMat,timeHSIC,timeGamma))
    return(resu)
}
