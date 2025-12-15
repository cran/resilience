#' Identifying Predictors of Resilience to Stressors in Single-Arm Studies of Pre-Post Change
#'
#' Studies of resilience in older adults are typically conducted with a single-arm 
#' where everyone experiences the stressor. The simplistic approach of regressing 
#' change versus baseline yields biased estimates due to mathematical coupling 
#' and regression-to-the-mean. This function provides a method to correct the bias.
#'
#' @param formula Formula object where LHS is Post response variable and RHS is 
#'   pre-response variable plus the covariates. Note: the first variable of the RHS 
#'   of the formula must be the pre-response variable. For example, `y2 ~ y1 + x1 + x2`.
#' @param data Data frame containing all the variables. Only complete cases are used 
#'   in the analysis, i.e. rows of dataframe with missing values in any of the 
#'   predictors are automatically deleted.
#' @param change A logical variable. If `TRUE` (default) the dependent variable 
#'   of regression is the pre-post change. If `FALSE`, the post response is 
#'   used as the dependent variable.
#' @param k A sensitivity analysis parameter. Typically, it is greater than or equal 
#'   to 1.0. It is recommended that the user provide at least two values to examine 
#'   how the effects vary with `k`. Default setting allows three values: k = 1.0, 1.5, 
#'   and 2.0. For more details about this parameter refer to the manuscript.
#' @param m Another sensitivity analysis parameter. It is set equal to 1.0. 
#'   It is recommended that the user not change this unless there is information 
#'   from an external study to justify a different value.
#' @param nboot Number of bootstrap samples for calculating the confidence intervals 
#'   of corrected regression coefficients. Default is 1000.
#' @param ci.level Confidence coefficient for confidence interval. Default is 0.95 
#'   (95\% confidence intervals).
#' @param boot.method The bootstrap method for confidence interval. Four options 
#'   are provided: "perc" (percentile), "norm" (normal approximation), "basic", 
#'   and "bca" (bias-corrected accelerated bootstrap). Default is "perc".
#' @param ncores Number of cores available for parallel computing. Default is set 
#'   to 2 due to CRAN requirements. If more cores are available, the user can 
#'   utilize all available cores with the command: 
#'   `ncores = parallel::detectCores()`
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{naive.beta}}{Unadjusted, naive estimates of regression coefficients}
#'   \item{\code{corrected.beta}}{The corrected coefficients of the variables. A matrix 
#'         with one column of parameter values for each value of sensitivity parameter \code{k}}
#'   \item{\code{CI}}{A list of length equal to the number of sensitivity values. Each 
#'         component of the list is a matrix with two columns of lower and upper 
#'         confidence interval for each parameter}
#'}
#'
#' @section Missing Data:
#' This function uses complete-case analysis only. If your data contains missing 
#' values, the function will issue a warning and exclude cases with any missing 
#' values. For more robust inference with missing data, consider using multiple 
#' imputation followed by \code{prepost_resilience_mi_list()}, which pools 
#' results across multiply imputed datasets using Rubin's rules.
#'
#' @section Parallel Processing:
#' The function uses the \code{parallel} and \code{foreach} packages to perform 
#' parallel computations of bootstrap confidence intervals for different values 
#' of the sensitivity parameter, \code{k}.
#'
#' @references 
#' Varadhan, R., Zhu, J., and Bandeen-Roche, K (2024). Identifying Predictors of 
#' Resilience to Stressors in Single-Arm Studies of Pre-Post Change. 
#' \emph{Biostatistics}. 25(4): 1094-1111.
#'
#' @seealso 
#' \code{\link{prepost_mi}} for analysis with multiply imputed data.
#'
#' @examples
#' \dontrun{
#' data(tkr)
#' 
#' names(tkr.dat)
#' dim(tkr.dat)
#' 
#' # pre-post change regression
#' ans1 <- prepost(post.Y ~ pre.Y + I(age-mean(age)) + I((age - mean(age))^2) + 
#'   bmi + gender + as.factor(smoker), data=tkr.dat, k=c(1.2, 1.5), nboot=200)
#' print(ans1)
#' 
#' # Post regression
#' ans2 <- prepost(post.Y ~ pre.Y + I(age-mean(age)) + I((age - mean(age))^2) + 
#'   bmi + gender + as.factor(smoker), data=tkr.dat, 
#'   k=c(1.2, 1.5), change=FALSE, nboot=200, boot.method="norm")
#' print(ans2)
#' 
#' # without any covariates
#' ans3 <- prepost(post.Y ~ pre.Y, data=tkr.dat, k=c(1.2, 2.0), nboot=200)
#' print(ans3)
#' 
#' # Bootstrapping using "bca" - relatively slow
#' # Not run
#' # ans4 <- prepost(post.Y ~ pre.Y, data=tkr.dat, k=c(1.2, 2.0), change=FALSE,
#' #   boot.method = "bca")
#' }
#'
#' @author Ravi Varadhan
#' @export
#' @importFrom stats lm model.frame cor coef terms
#' @importFrom nptest np.boot
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom graphics abline barplot legend matplot mtext par
#' @importFrom stats as.formula complete.cases qnorm qt sd var
#'
prepost <- function(formula, data, change=TRUE, k=c(1.0,1.5,2.0), m=1, nboot=1000, 
                    ci.level=0.95, boot.method=c("perc", "norm", "basic", "bca"),
                    ncores=2){
  
  statfun <- function(x, bdata, kb){
    datax <- bdata[x, ]
    mf <- model.frame(formula, data = datax)
    rho <- cor(mf[,1], mf[,2])
    mod.naive <- lm(formula.new, data = datax)
    beta.naive <- coef(mod.naive)
    if (npar > 2) {
      mod1 <- lm(formula.1, data = datax)
      S <- summary(mod1)$r.squared } else S <- 0
    krho <- kb*rho
    correction.1b <- (m - krho) / (1-S) - (m - 1*change)
    if (npar > 2) {
      correction.xb <- (krho - m)/(1-S)*coef(mod1)[-1]
      correction.b <- c(correction.1b, correction.xb)
    } else correction.b <- correction.1b
    beta <- beta.naive + c(0, correction.b)
    return(beta)
  }
  
  # Check for missing data and warn
  formula_vars <- all.vars(formula)
  available_vars <- formula_vars[formula_vars %in% names(data)]
  
  if (length(available_vars) > 0) {
    missing_count <- sum(is.na(data[, available_vars, drop = FALSE]))
    if (missing_count > 0) {
      warning(sprintf(
        "Data contains %d missing values in model variables. ",
        missing_count
      ), "Results are based on complete cases only.", "\n", 
      "Consider multiply imputated analysis with prepost_mi() ",
      "for more robust inference.")
    }
  }
  
  if (length(showConnections()) > 0) closeAllConnections()
  cl <- makeCluster(ncores)  
  registerDoParallel(cl)  
  nbc <- max(2, round(ncores/2))
  
  varbs <- labels(terms(formula)) # independent variables
  y1 <- varbs[1]
  npar <- length(varbs) + 1 
  
  intcpt <- attr(terms(formula), "intercept")
  formula.text <- if(intcpt==1) paste(formula[[2]], "-", y1, "~", paste(varbs, collapse = " + ")) else paste(formula[[2]], "-", y1, "~", paste(varbs,  collapse = " + "), "-1") 
  formula.text1 <- formula.text
  formula.new <- if(change) as.formula(formula.text1) else formula
  if (npar > 2) formula.1 <- if(intcpt==1) paste(y1, "~", paste(varbs[-1], collapse = " + ")) else paste( y1, "~", paste(varbs[-1],  collapse = " + "), "-1") 
  
  mf <- model.frame(formula, data = data)
  rho <- cor(mf[,1], mf[,2])
  mod.naive <- lm(formula.new, data = data)
  beta.naive <- coef(mod.naive)
  
  if (npar > 2) {
    mod1 <- lm(formula.1, data = data)
    S <- summary(mod1)$r.squared} else S <- 0
  
  beta <- matrix(NA, length(coef(mod.naive)), length(k))
  
  for (j in 1:length(k)){
    krho <- k[j]*rho
    correction.1b <- (m - krho) / (1-S) - (m - 1*change)
    if (npar > 2) {
      correction.xb <- (krho - m)/(1-S)*coef(mod1)[-1]
      correction.b <- c(correction.1b, correction.xb)
    } else correction.b <- correction.1b
    beta[,j] <-   beta.naive + c(0, correction.b)
  }
  
  boot.method <- match.arg(boot.method)
  
  ci <- foreach (j = 1:length(k), combine=`cbind`) %dopar% {
    npbs <- nptest::np.boot(x=1:nrow(data), statistic = statfun, kb=k[j], bdata = data, 
                            R=nboot, level=ci.level, 
                            method=boot.method, boot.dist=FALSE, parallel = TRUE, 
                            cl=parallel::makeCluster(nbc))
    t(eval(parse(text=paste("npbs", boot.method, sep="$"))))
  }
  
  closeAllConnections()
  
  names(ci) <- paste("k", k, sep="=")
  rownames(beta) <- names(beta.naive)
  colnames(beta) <- paste("k", k, sep="=")
  list(naive.beta = coef(summary(mod.naive)), corrected.beta = beta, CI = ci)
}

.onAttach <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    requireNamespace("nptest", quietly = TRUE)
  })
}

