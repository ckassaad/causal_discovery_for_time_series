all.perm <- function(n) {
  p <- matrix(1, ncol = 1)
  for (i in 2:n) {
    p <- pp <- cbind(p, i)
    v <- c(1:i, 1:(i - 1))
    for (j in 2:i) {
      v <- v[-1]
      p <- rbind(p, pp[, v[1:i]])
    }
  }
  p
}

