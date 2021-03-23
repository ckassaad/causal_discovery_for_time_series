# Copyright (c) 2012-2013 Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.

indtestAll <- function(f,x,y,alpha,pars = list())
{
    result <- f(x,y,alpha,pars)
}

indtestMutualAll <- function(f,x,alpha,pars = list())
{
    result <- f(x,alpha,pars)
}
