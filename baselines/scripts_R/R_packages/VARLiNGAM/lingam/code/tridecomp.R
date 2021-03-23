tridecomp <- function( W, choice='qr' ) {

  # SYNTAX:
  # res <- tridecomp( W, choice )
  # QR, RQ, QL, or LQ decomposition specified by
  # choice = 'qr', 'rq', 'ql', or 'lq' respectively
  #
  # res$A is the first matrix and res$B is the second
  #
  # Based on MATLAB code kindly provided by Matthias Bethge
  # Adapted for R by Patrik Hoyer

  m <- nrow(W)
  n <- ncol(W)
  Jm <- matrix(0,m,m)
  Jm[m:1,] <- diag(m)
  Jn <- matrix(0,n,n)
  Jn[n:1,] <- diag(n)

  res <- list()

  if (choice == 'qr') {
    r <- qr(W)
    res$A <- qr.Q(r)
    res$B <- qr.R(r)
  } else if (choice == 'lq') {
    r <- qr(t(W))
    res$A <- t(qr.R(r))
    res$B <- t(qr.Q(r))
  } else if (choice == 'ql') {
    r <- qr(Jm %*% W %*% Jn)
    res$A <- Jm %*% qr.Q(r) %*% Jm
    res$B <- Jm %*% qr.R(r) %*% Jn
  } else if (choice == 'rq') {
    r <- qr(Jn %*% t(W) %*% Jm)
    res$A <- t(Jn %*% qr.R(r) %*% Jm)
    res$B <- t(Jn %*% qr.Q(r) %*% Jn)
  }

  res
  
}
