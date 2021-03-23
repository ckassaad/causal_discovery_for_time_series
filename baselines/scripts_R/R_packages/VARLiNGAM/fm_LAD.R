fm_LAD<-function(Data, u_res){

  # this function finds the correction to be added to the LAD estimation in
  # order to get a FM-LAD estimation (fully modified LAD estimation see
  # Phillips ET 1995)
  #
  # INPUT
  # Data: original data (not in canonical form)
  # u_res: residuals of LAD estimation
  # 
  # OUTPUT
  # correction terms: Mhat_fmLAD = Mhat_LAD + correction term
  #
  # e.g. Mhat_fmLAD = res$Mhat + fm_LAD(Data, res$residuals)

  require(sm)

  T <- nrow(u_res)
  Td <- nrow(Data)
  k <- ncol(u_res)
  nl = Td-T # number of time lags
  p <- (k*nl)
  FMterm <- matrix(nrow=k, ncol=p)

  Xt <- as.data.frame(matrix(nrow=T, ncol=p))
  for(j in 1:nl){
    Xt[,((j*k-k+1):(j*k))] <- Data[(nl-j+1):(Td-j),]
  }
  Xt <- as.matrix(Xt)
  XX <- matrix(0,nrow=p,ncol=p)
  for (j in 1:T){
    mm <- Xt[j,] %*% t(Xt[j,])
    XX <- XX+mm
  }
  DXt <- Xt[-1,]-Xt[-T,]
  XDX <- matrix(0,nrow=p,ncol=p)
  for (j in 1:(T-1)){
    mm <- Xt[(j+1),] %*% t(DXt[j,])
    XDX <- XDX+mm
  }
  t1 <- T-1
  U <- matrix(nrow=t1, ncol=(p+1))

  ## start loop
  for(i in 1:k){
    u0 <- u_res[,i]
    #h(0): probability density of u0 at 0
    h0 <- sm.density(u0, eval.points=0, display="none")$estimate
    U[,1] <- u0[2:T]
    U[,2:(p+1)] <- DXt

    ###### function Gammaj
    Gammaj <- function(u, tauu){
      ta <- abs(tauu)
      N <- nrow(u)
      G1 <- matrix(0, nrow=ncol(u), ncol=ncol(u))
      for (jj in (ta +1):N){
	mm <- u[jj,] %*%  t(u[jj-ta,])
	G1 <- mm+G1
      }
      rr <- G1/N
      if(tauu<0){ rr <- t(G1/N)}
      rr
    }
    
    ###### function plug_in_bw
    # find the right bandwidth with the plug-in Andrews (Econometrica 1991)
    # method
    plug_in_bw <- function(UU){
      aa1 <- 1:ncol(UU)
      aa2 <- 1:ncol(UU)
      for (j in 1:ncol(UU)){
	yt <- UU[-1,j]
	yt_1 <- UU[-nrow(UU),j]
	ar1 <- lm(yt~yt_1-1)
	cc <- ar1$coefficients
	vr <- var(ar1$residuals)
	aa1[j] <- (4*cc^2*vr^2)/(1-cc)^8
	aa2[j] <- (vr^2)/(1-cc)^4
      }
      alpha2 <- sum(aa1, na.rm=TRUE)/sum(aa2, na.rm=TRUE)
      bw <- 2.6614*(alpha2*t1)^(1/5)
      bw
    }

    ###### function K_parzen
    # Parzen Kernel
    K_parzen <- function(x){
      rr <- 0
      if (abs(x) <= 0.5){
	rr <- 1-6*x^2+ 6*abs(x)^3
      }
      if (abs(x) > 0.5 & abs(x) <= 1 ){
	rr <- 2*(1-abs(x))^3
      }
      rr
    }

    Ot <- matrix(0, ncol=(p+1), nrow=(p+1))
    for (j in (-t1+1):(t1-1)){
      mm <- Gammaj(U, abs(j)) *K_parzen(j/plug_in_bw(U))
      Ot <- mm+Ot
    }
    Duu <- matrix(0, ncol=(p+1), nrow=(p+1))
    for (j in 0:(t1-1)){
      mm <- Gammaj(U, abs(j)) *K_parzen(j/plug_in_bw(U))
      Duu <- mm+Duu
    }
    
    Oxx <- Ot[-1,-1]
    Dxx <- Duu[-1,-1]
    
    nu <- sign(u0[2:T])
    Wt <- U
    Wt[,1] <- nu

    Oww <- matrix(0, ncol=(p+1), nrow=(p+1))
    for (j in (-t1+1):(t1-1)){
      mm <- Gammaj(Wt, abs(j)) *K_parzen(j/plug_in_bw(Wt))
      Oww <- mm+Oww
    }
    Dww <- matrix(0, ncol=(p+1), nrow=(p+1))
    for (j in 0:(t1-1)){
      mm <- Gammaj(Wt, abs(j)) *K_parzen(j/plug_in_bw(Wt))
      Dww <- mm+Dww
    }

    Oxv <- Oww[2:(p+1),1]
    Dxv <- Dww[2:(p+1),1]
    
    Dxvpp <- (XDX/t1-Dxx)%*% solve(Oxx)%*%Oxv  + Dxv
    FMterm[i,] <- - solve(2*h0*XX)%*%(t1*Dxvpp)
  }
  ## end loop

  rr <- as.list(1:nl)
  for(j in 1:nl){
    rr[[j]] <- FMterm[,(k*j-k+1):(k*j)]
  }

  rr
}