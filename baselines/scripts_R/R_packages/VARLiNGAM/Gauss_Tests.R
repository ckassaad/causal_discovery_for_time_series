Gauss_Tests <- function(X) {

  # INPUT
  # X: data vector or data matrix (assume more observations than variables)
  #
  # OUTPUT
  # p-values of normality tests: Shapiro-Wilk, Shapiro-Francia and Jarque-Bera

  require(nortest)
  require(tseries)

  X <- as.matrix(X)
  dims <- dim(X)

  if (dims[1]<dims[2]) {
    X <- t(X)
    dims <- dim(X)
  }

  m <- dims[1]
  if (m > 5000) {
    cat("get subsample of residuals to perform normality-test\n")
    subsamp <- sample(m,5000)
    X1 <- X[subsamp,]
  }
  else X1 <- X

  shapiro <- list(); shapfranc <- list(); jarque <- list();
  pval1 <- NULL; pval2 <- NULL; pval3 <- NULL

  for (i in 1:dims[2]) {
    # Shapiro-Wilk test
    shapiro[[i]] <- shapiro.test(X1[,i])
    pval1 <- c(pval1, unclass(shapiro[[i]])$p.value)
    # Shapiro-Francia test
    shapfranc[[i]] <- sf.test(X1[,i])
    pval2 <- c(pval2, unclass(shapfranc[[i]])$p.value)
    # Jarque-Bera test
    jarque[[i]] <- jarque.bera.test(X1[,i])
    pval3 <- c(pval3, unclass(jarque[[i]])$p.value)
  }

  ls <- list()
  ls$ShapiroWilk <- pval1
  ls$ShapiroFrancia <- pval2
  ls$JarqueBera <- pval3

  ls

}
