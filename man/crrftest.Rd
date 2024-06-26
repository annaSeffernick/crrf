% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crrftest.R
\name{crrftest}
\alias{crrftest}
\title{Penalized Likelihood Ratio Test in Fine-Gray Competing Risks Regression}
\usage{
crrftest(form, etype, dset, test = ~1, firth = TRUE, alpha = 0.05, eps = 1e-10)
}
\arguments{
\item{form}{A formula with the survival model, with usual survival package syntax for the time to event data: Surv(time, event) notation. Values must be column names in dset.}

\item{etype}{A character to indicate event type of interest, as defined by the event indicator in dset.}

\item{dset}{A data.frame including the survival variables defined in the formula and any additional covariates to include in the model.}

\item{test}{A right-hand formula of parameters to test (e.g. \code{~ B + D}). As default the null hypothesis that all parameters are 0 is tested.}

\item{firth}{A Boolean of whether to fit with Firth penalty. Default is TRUE.}

\item{alpha}{A numeric of significance level. Default is 0.05.}

\item{eps}{A numeric of threshold to begin confidence interval calculations. Default is eps=1e-10. Not recommended to decrease, as this may lead to errors in calculating the Firth Penalty.}
}
\value{
A crrftest object, contianing a list with the LRT statistic, p-value, degrees of freedom, log-likelihoods, description of method, and function call.
}
\description{
Penalized Likelihood Ratio Test in Fine-Gray Competing Risks Regression
}
\examples{
set.seed(12345)
n <- 100
tm <- rexp(n)
ev <- factor(rbinom(n,2,0.3))
x <- rnorm(n)
y <- rbinom(n,1,0.5)
u <- runif(n)
dset <- cbind.data.frame(tm=tm,ev=ev,x=x,y=y,u=u)
fgm.res <- crrf(Surv(tm,ev)~x+y+u,etype="1",dset,firth=TRUE,CI=FALSE)
fgm.test <- crrftest(Surv(tm,ev)~x+y+u,etype="1",dset,test=~x+y,firth=TRUE)
}
