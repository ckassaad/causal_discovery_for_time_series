canonicalform_USfirmsdata <- function(Growth, nlags=1) {

  # INPUT
  # Growth: firm growth data
  # nlags: time lags in analysis
  # 
  # OUTPUT
  # firm growth data in canonical form

  ## canonical format - columns:
  # 1st col: 'id number', each firm gets an id number
  # 2nd col: 'time-point t of current value, starting from 1'
  # following: 'past values', nvar*nlags colums:
  #            nvar variables at time t-1
  #            nvar variables at time t-2
  #            nvar variables at time t-nlags
  # following: 'current values', nvar columns:
  #            nvar variables at time t
  # where nvar = number of variables

  # only use observation at time t if:
  # 1) obs. at time t-1 ... t-nlags have the same index (same firm) as obs.
  #    at time t AND
  # 2) all obs. from time t to t-nlags are in in consecutive years AND
  # 3) no observations are missing (NA)

  dims <- dim(Growth)

  # get id-number (change company name in number) and time points
  idvec <- array(0,dim=c(dims[1],1))
  timevec <- array(0,dim=c(dims[1],1))
  cnt <- 1
  temp1 <- Growth$coname[1]
  idvec[1] <- cnt
  timevec[1] <- Growth$ayear[1]-1972 
  for (i in 2:dims[1]) {
    timevec[i] <- Growth$ayear[i]-1972 # earliest year = 1973 = timepoint 1
    temp2 <- Growth$coname[i]
    if (temp1 == temp2) idvec[i] <- cnt
    else {
      cnt <- cnt+1
      idvec[i] <- cnt
      temp1 <- temp2
    }
  }
  Growth$coname <- idvec
  colnames(Growth)[1] <- "id"
  Growth$ayear <- timevec
  colnames(Growth)[2] <- "time"
  colnames(Growth)[3:6] <- c("empl.gr","sales.gr","rnd.gr","opinc.gr")

  cat("Put data in canonical format. This can take some time ...")

  data <- array(0,dim=c(dims[1]-nlags,2+4*(nlags+1))) # 4 = nvar
  cnt <- 1

  # get canonical form, this is slow
  for (i in (nlags+1):dims[1]) {
    if ( all(Growth$id[i] == Growth$id[(i-1):(i-nlags)]) &
         all(Growth$time[i] == Growth$time[(i-1):(i-nlags)]+1:nlags) &
         all(!is.na(Growth[i:(i-nlags),3])) ) {
      temp <- t(Growth[c((i-1):(i-nlags),i),3:6])
      dim(temp) <- c(1,4*nlags+4)
      data[cnt,] <- c(Growth$id[i],Growth$time[i],temp)
      cnt <- cnt+1
    } # end if
  } # end for i

  data <- data[1:(cnt-1),]

  # get labels for colnames
  temp1 <- array(0,dim=c(1,4)); temp2 <- NULL
  for (i in 1:nlags) {
    temp1 <- paste(sep="",'pastval', 1:4, i)
    temp2 <- c(temp2,temp1)
  }
  temp1 <- paste(sep="",'curval', 1:4)
  temp2 <- c(temp2,temp1)

  colnames(data) <- c("id","curtime",temp2)

  cat("Done!\n")

  data

}
