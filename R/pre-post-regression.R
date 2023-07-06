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
    correction.1b <- (1 - krho) / (1-S) - (1 - 1*change)
    if (npar > 2) {
      correction.xb <- (krho - 1)/(1-S)*coef(mod1)[-1]
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


