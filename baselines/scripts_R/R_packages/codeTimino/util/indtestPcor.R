# Copyright (c) 2010-2013 Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.


partial.cor <- function(x,y,z)
{
    a <- pcor.test(x,y,z)
    result <- a$p.value
}


nonlin.partial.cor <- function(x,y,z, dimZ = 1)
{
    if(dimZ == 1)
    {
        epsx <- gam(x~s(z))$residuals
        epsy <- gam(y~s(z))$residuals
    }
    if(dimZ == 2)
    {
        epsx <- gam(x~s(z[,1]+z[,2]))$residuals
        epsy <- gam(y~s(z[,1]+z[,2]))$residuals
    }
    if(dimZ == 4)
    {
        epsx <- gam(x~s(z[,1]+z[,2]+z[,3]+z[,4]))$residuals
        epsy <- gam(y~s(z[,1]+z[,2]+z[,3]+z[,4]))$residuals
    }
    result <- cor.test(epsx,epsy)$p.value
}

