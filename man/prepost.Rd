\name{prepost}
\alias{prepost}
\title{Identifying Predictors of Resilience to Stressors in Single-Arm Studies of Pre-Post Change}
\description{
Studies of resilience in older adults are typically conducted with a single-arm where everyone experiences the stressor. The simplistic
approach of regressing change versus baseline yields biased estimates due to mathematical coupling and regression-to-the-mean. This function provides a method to correct the bias.
}
\usage{
prepost(formula, data, change=TRUE, k=c(1.0,1.5, 2.0), m=1, nboot=1000, 
    ci.level=0.95, boot.method=c("perc", "norm", "basic", "bca"),
    ncores=2)
}
\arguments{
  \item{formula}{
formula object where LHS is Post response variable and RHS is pre-response variable plus the covariates; note: the first variable of the RHS of the formula must be the pre-response variable. For example, \code{y2 ~ y1 + x1 + x2}.
}
  \item{data}{
data-frame containing all the variables. Only complete cases are used in the analysis, i.e. rows of dataframe with missing values in any of the predictors are automatically deleted.
}
  \item{k}{
A sensitivity analysis parameter.  Typically, it is greater than or equal to 1.0.  It is recommended that the user provide at least two values to examine how the effects vary with `k'.  Default setting allows three values:  k = 1.0, 1.5, and 2.0.  For more details about this parameter refer to the manuscript.
}
  \item{m}{
Another sensitivity analysis parameter.  It is set equal to 1.0.  It is recommended that the user not change this unless there is information from an external study to justify a diferent value.  For more details about this parameter refer to the manuscript.
}
  \item{change}{
A logical variable.  If \code{TRUE} the dependent variable of regression is the pre-post change. If \code{FALSE}, the post response is used as the dependent variable.
}
  \item{nboot}{
Number of boot strap samples for calculating the confidence intervals of corrected regression coefficients.  Default is 1000. 
}
  \item{ci.level}{
  Confidence coefficient for confidence interval.  Default is 95\% confidence intervals.
}
  \item{boot.method}{
The bootstrap method for confidence interval.  Four options are provided: percentile, normal approximation, basic, and bias-corrected accelerated bootstrap. Default is ``percentile.''}
  \item{ncores}{
  Number of cores avaialble for parallel computing.  Default is set to 2 due to CRAN requirements.  If more cores are available, the user can utilize all available cores with the command: \cr
  \code{ncores = parallel::detectCores()}
}
}
\details{
The function uses the \code{parallel} and \code{foreach} packages to perform parallel computations of bootstrap confidence intervals for different values of the sensitivity parameter, `k.'
}
\value{
  \item{beta.naive}{Unadjusted, naive estimates of regression coefficients}
  \item{corrected.beta}{The corrected coefficients of the variables. A matrix with one column of parameter values for each value of sensitivity parameter `k.'}
  \item{CI}{A list of length equal to the number of sensitivity values. Each component of the list is a matrix with two columns of lower and upper confidence interval for each parameter.}
}
\references{
Varadhan, R., Zhu, J., and Bandeen-Roche, K (2023). Identifying Predictors of Resilience to Stressors
in Single-Arm Studies of Pre-Post Change. \emph{Biostatistics}. To appear. 
}
\author{
Ravi Varadhan
}

\examples{
data(tkr)

names(tkr.dat)
dim(tkr.dat)

# pre-post change regression
ans1 <- prepost(post.Y ~ pre.Y + I(age-mean(age)) + I((age - mean(age))^2) + 
bmi + gender + as.factor(smoker), data=tkr.dat, k=c(1.2, 1.5), nboot=200)
print(ans1)

# Post regression
ans2 <- prepost(post.Y ~ pre.Y + I(age-mean(age)) + I((age - mean(age))^2) + 
bmi + gender + as.factor(smoker), data=tkr.dat, 
k=c(1.2, 1.5), change=FALSE, nboot=200, boot.method="norm")
print(ans2)

# without any covariates
ans3 <- prepost(post.Y ~ pre.Y, data=tkr.dat, k=c(1.2, 2.0), nboot=200)
print(ans3)

# Bootstrapping using "bca" - realtively slow
# Not run
#ans4 <- prepost(post.Y ~ pre.Y, data=tkr.dat, k=c(1.2, 2.0), change=FALSE,
#boot.method = "bca")

}
\keyword{pre-post}
\keyword{resilience}
