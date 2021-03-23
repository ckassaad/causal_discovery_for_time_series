
#' Estimate response function for each person using smoothed Finite Impulse Response.
#' @param data The data to be used to estimate response function 
#' @param stimuli A vector containing '0' when the stimuli of interest is not present 
#' and '1' otherwise. Number of observations across time must equal the length of data. 
#' @param interval Time between observations; for fMRI this is the repetition time. 
#' Defaults to 1.
#' @return Shape of response function and convolved time series vector. 
#' @keywords internal 
sFIR <- function(data, 
                 stimuli, 
                 response_length = 16, 
                 interval = 1){
  
  run_length <- length(stimuli)
  t = seq(from = 1, to = response_length, by =interval)
  
  X_fir <- matrix(0, run_length, (response_length/interval))
  
  # set up basis vectors
  c_onsets <- which(stimuli == 1)
  id_cols <- seq(from=1, to = response_length/interval)
  for (j in 1: length(c_onsets)){
    id_rows <- seq(from =c_onsets[j], to= (c_onsets[j]+response_length/interval-1))
    for (k in 1:min(length(id_rows),(length(X_fir[,1])-c_onsets[j]))){
      X_fir[id_rows[k], id_cols[k]]<- 1
    }
  }
  
  ### estimate beta 
  # get R2s to select which vector to use
  if (dim(as.matrix(data))[2]> 1){
  R2 <- matrix(,length(data[1,]), 1)
  for (p in 1:length(data[1,]))
    R2[p]<- summary(lm(data[,p]~X_fir))$r.squared
  best <- which(R2 == max(R2))
  } else 
    best <- 1
  
  ### For smoothing
  C <- seq(1:length(t))%*%matrix(1, 1, length(t))
  h <- sqrt(1/(7/interval))
  
  C2 <- apply((C-t(C)), c(1,2), function(i) i^2)
  RI <- solve(.1*exp(-h/2*(C2)))
  # MRI <- matrix(0,1, n_cond*length(t)+1)
  MRI <- matrix(0,length(t),length(t))
  MRI[1:length(t),1:length(t)] = RI
  ## end smooth 
  
  if (dim(as.matrix(data))[2]> 1){
  est_hrf <- solve(t(X_fir)%*%X_fir + 1^2*MRI)%*%t(X_fir)%*%data[,best]
  } else
    est_hrf <- solve(t(X_fir)%*%X_fir + 1^2*MRI)%*%t(X_fir)%*%data
#plot(ts(est_hrf))
  conv_onsets <- stats::convolve(as.numeric(stimuli), rev(est_hrf), type = c("open"))
  
  res <- list(est_rf = est_hrf, 
              conv_stim_onsets = conv_onsets)
  return(res)
}