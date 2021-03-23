# Copyright (c) 2010-2013 Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
#
# This file contains several implementations of Granger causality:
#
# biv_granger_causality
# biv_nl_granger_causality
# partial_granger_causality
#
# pairwise_granger_causality
# pairwise_nl_granger_causality
# multiv_granger_causality


granger_pairwise<-function(x,y,alpha,max_lag,output = FALSE)
{
    #tests, whether y is causing x.
    #mod2 = mod_full
    #mod1 = mod_restr
    #choose the order of mod_full to be the order of mod_restr 
    
    mod2_fit<-ar(cbind(x,y),aic=TRUE,order.max=max_lag)
    order2<-max(mod2_fit$order,1)
    mod_fit<-ar(x,aic=FALSE,order.max=order2)
    order1<-max(mod_fit$order,1)
    
    if(output)
    {
        show(c(order1,order2))
    }
    r1<-mod_fit$resid[-c(1:order1)]
    r2<-mod2_fit$resid[-c(1:order2),1]
    r2b<-mod2_fit$resid[-c(1:order2),2]
    RSS1<-t(r1)%*%r1
    RSS2<-t(r2)%*%r2
    
    T_num <- (RSS1-RSS2)/(2*order2-order1)
    T_den <- RSS2/(length(x)-(order2+order2+1))
    T <- T_num/T_den
    if(2*order2-order1==0)
    {
        stop("The choice of time orders leads to a problem.")
    }
    quan <- qf(1-alpha,2*order2-order1,length(x)-(order2+order2+1))
    pval <- 1-pf(T,2*order2-order1,length(x)-(order2+order2+1))
    Tquan<-c(T,quan)   
    
    if(output)
    {
        show("test statistic, critical value and p-val are:")
        show(c(Tquan, pval))
        if(Tquan[1]<Tquan[2])
        {	
            show("Thus, H0 (no reduction) is not rejected: 2nd var DOES NOT CAUSE 1st var.")
        }
        else
        {
            show("Thus, H0 (no reduction) is rejected: 2nd var CAUSES 1st var.")
        }
    }
    Tquan<-Tquan
}




granger_pairwise_nl<-function(x,y,alpha,max_lag, output = FALSE)
{
    #tests, whether y is causing x.
    #choose the order of mod_full to be the order of mod_restr 
    mod_fit<-gamts.fit2(x,max_order=max_lag,order.fixed=FALSE)
    mod2_fit<-gamts.fit(x,y,max_order=mod_fit$order,order.fixed=TRUE)
    df1<-length(x)-mod_fit$degrf
    df2<-length(x)-mod2_fit$degrf
    
    order1<-mod_fit$order
    order2<-mod2_fit$order
    r1<-mod_fit$resid
    r2<-mod2_fit$resid
    RSS1<-t(r1)%*%r1
    RSS2<-t(r2)%*%r2
    
    T_num <- (RSS1-RSS2)/(df2-df1)
    T_den <- RSS2/(length(x)-df2)
    T <- T_num/T_den
    show(c(df2, df1, RSS1, RSS2)) 
    quan <- qf(1-alpha,df2-df1,length(x)-df2)
    pval <- 1 - pf(T,df2-df1,length(x)-df2)
    Tquan<-c(T,quan)  
    
    if(output)
    {
        show("test statistic, critical value and p-val are:")
        show(c(Tquan, pval))
        if(Tquan[1]<Tquan[2])
        {	
            show("Thus, H0 (no reduction) is not rejected: 2nd var DOES NOT CAUSE 1st var.")
        }
        else
        {
            show("Thus, H0 (no reduction) is rejected: 2nd var CAUSES 1st var.")
        }
    }
    Tquan<-Tquan 
}


granger_partial<-function(x,y,z,alpha,max_lag,output = FALSE)
{
    #tests, whether y is causing x given the information from z.
    #mod_2 = mod_full
    #mod_1 = mod_restr
    #choose the order of mod_full to be the order of mod_restr 
    
    mod2_fit<-ar(cbind(x,y,z),aic=TRUE,order.max=max_lag)
    order2<-max(mod2_fit$order,1)
    mod_fit<-ar(cbind(x,z),aic=FALSE,order.max=order2)
    order1<-max(mod_fit$order,1)
    
    r1<-mod_fit$resid[-c(1:order1),1]
    r2<-mod2_fit$resid[-c(1:order2),1]
    RSS1<-t(r1)%*%r1
    RSS2<-t(r2)%*%r2
    
    dim_cond_set<-dim(z)[2]
    if(length(z)==length(x)){dim_cond_set=1}
    
    T_num <- (RSS1-RSS2)/order2
    T_den <- RSS2/(length(x)-((dim_cond_set+1)*order2+1))
    T <- T_num/T_den
    
    quan <- qf(1-alpha,order2,length(x)-((dim_cond_set+1)*order2+1))
    Tquan <- c(T,quan)
    if(output)
    {
        show("test statistic and critical value are:")
        show(Tquan)
        if(Tquan[1]<Tquan[2])
        {	
            show("Thus, H0 (no reduction) is not rejected: 2nd var DOES NOT CAUSE 1st var.")
        }
        else
        {
            show("Thus, H0 (no reduction) is rejected: 2nd var CAUSES 1st var.")
        }
    }
    Tquan <- Tquan
    
}

granger_dag_pairwise<-function(M,alpha,max_lag)
{
    #M (nxp) should consists of p cols, each of which is the realization of one
    #variable, that has n data points.
    #(i,j)==1 in the answer matrix means causal link i->j
    #(i,j)==0 in the answer matrix means no causal link i->j
    np<-dim(M)
    CC<-matrix(88,np[2],np[2])
    for(i in 1:(np[2]-1))
    {
        for(j in (i+1):np[2])
        {
            Fc<-granger_pairwise(M[,j],M[,i],alpha,max_lag)
            ifelse(Fc[1]<Fc[2],
                   CC[i,j]<-0,
                   CC[i,j]<-1)
            Fc<-granger_pairwise(M[,i],M[,j],alpha,max_lag)
            ifelse(Fc[1]<Fc[2],
                   CC[j,i]<-0,
                   CC[j,i]<-1)
        }
    }
    return(CC)
}





granger_dag_pairwise_nl<-function(M,alpha,max_lag)
{
    #M (nxp) should consists of p cols, each of which is the realization of one
    #variable, that has n data points.
    #(i,j)==1 in the answer matrix means causal link i->j
    #(i,j)==0 in the answer matrix means no causal link i->j
    np<-dim(M)
    CC<-matrix(88,np[2],np[2])
    for(i in 1:(np[2]-1))
    {
        for(j in (i+1):np[2])
        {
            Fc<-granger_pairwise_nl(M[,j],M[,i],alpha,max_lag)
            ifelse(Fc[1]<Fc[2],
                   CC[i,j]<-0,
                   CC[i,j]<-1)
            Fc<-granger_pairwise_nl(M[,i],M[,j],alpha,max_lag)
            ifelse(Fc[1]<Fc[2],
                   CC[j,i]<-0,
                   CC[j,i]<-1)
        }
    }
    return(CC)
}





granger_dag_partial<-function(M,alpha,max_lag,output = TRUE)
{
    #M (nxp) should consists of p cols, each of which is the realization of one
    #variable, that has n data points.
    #(i,j)==1 in the answer matrix means causal link i->j
    #(i,j)==0 in the answer matrix means no causal link i->j
    np<-dim(M)
    CC<-matrix(88,np[2],np[2])
    for(i in 1:(np[2]-1))
    {
        for(j in (i+1):np[2])
        {
            Fc<-granger_partial(M[,j],M[,i],M[,-c(i,j)],alpha,max_lag, output)
            ifelse(Fc[1]<Fc[2],
                   CC[i,j]<-0,
                   CC[i,j]<-1)
            Fc<-granger_partial(M[,i],M[,j],M[,-c(i,j)],alpha,max_lag, output)
            ifelse(Fc[1]<Fc[2],
                   CC[j,i]<-0,
                   CC[j,i]<-1)
        }
    }
    return(CC)
}

