ols_est <- function(X,y,regstats=FALSE) {

  # ordinary least squares estimation of a regression model

  # INPUT
  # X, y: Data such that of y=c+Xb (with c, b the parameters to be estimated)
  # regstats: boolean, if TRUE print coefficients, standard errors, t-values,
  #       p-values, R^2 (and some more)
  # 
  # OUTPUT
  # b: regression coefficients s.t. y=c+Xb

  fit <- lm(y~X)

  if (regstats) {
    print(summary(fit))
  }

  b <- t(fit$coef)
  b

}