#' Penalized Likelihood Ratio Test in Fine-Gray Competing Risks Regression
#'
#' @param form  A formula with the survival model, with usual survival package syntax for the time to event data: Surv(time, event) notation. Values must be column names in dset.
#' @param etype A character to indicate event type of interest, as defined by the event indicator in dset.
#' @param dset  A data.frame including the survival variables defined in the formula and any additional covariates to include in the model.
#' @param test  A right-hand formula of parameters to test (e.g. \code{~ B + D}). As default the null hypothesis that all parameters are 0 is tested.
#' @param firth A Boolean of whether to fit with Firth penalty. Default is TRUE.
#' @param alpha A numeric of significance level. Default is 0.05.
#' @param eps   A numeric of threshold to begin confidence interval calculations. Default is eps=1e-10. Not recommended to decrease, as this may lead to errors in calculating the Firth Penalty.
#'
#' @return A crrftest object, contianing a list with the LRT statistic, p-value, degrees of freedom, log-likelihoods, description of method, and function call.
#' @import stats
#' @import survival
#' @export
#'
#' @examples
#' set.seed(12345)
#' n <- 100
#' tm <- rexp(n)
#' ev <- factor(rbinom(n,2,0.3))
#' x <- rnorm(n)
#' y <- rbinom(n,1,0.5)
#' u <- runif(n)
#' dset <- cbind.data.frame(tm=tm,ev=ev,x=x,y=y,u=u)
#' fgm.res <- crrf(Surv(tm,ev)~x+y+u,etype="1",dset,firth=TRUE,CI=FALSE)
#' fgm.test <- crrftest(Surv(tm,ev)~x+y+u,etype="1",dset,test=~x+y,firth=TRUE)
crrftest <- function(form, etype, dset, test=~1,
                     firth=TRUE, alpha=0.05, eps=1e-10)
{
  # Extract full model covs
  cov1.name <- unlist(strsplit(as.character(form)[3], "[+]"))
  cov1.name <- trimws(cov1.name)
  # check models are nested
  form2 <- as.formula(paste(as.character(form)[2],
                            as.character(test)[2], sep="~"))
  cov2.name <- unlist(strsplit(as.character(test)[2], "[+]"))
  cov2.name <- trimws(cov2.name)
  if(!all(cov2.name %in% cov1.name)&&(cov2.name !=1))
    stop("test formula is not a subset of whole formula!")
  # If nested, fit full model
  full <- crrf(form, etype, dset, firth, CI=FALSE, alpha, eps, CI.print = F)
  # If intercept only model, create constant intercept variable, fit with regular fine gray
  # because penalty does not affect
  if(length(cov2.name)==1 && cov2.name==1){
    intcpt <- rep(1, times=nrow(dset))
    dset2 <- cbind.data.frame(dset, intcpt)
    form2.m <- as.formula(paste(as.character(form2)[2], "intcpt", sep="~"))
    fgdat <- survival::finegray(form2.m, data=dset2, etype=etype)
    red <- survival::coxph(Surv(fgstart, fgstop, fgstatus)~intcpt, weight=fgdat$fgwt, data=fgdat)
    loglik.vec <- c(full$loglik, red$loglik[2])
    df <- length(cov1.name)
    method <- ifelse(firth, "Standard ML for Intercept-Only Model and Penalized ML for Full Model",
                     paste0("Standard ML for Intercept-Only and Full Models"))
  }
  else{
    # Fit the reduced model
    red <- crrf(form2, etype, dset, firth, CI=FALSE, alpha, eps, CI.print = F)
    # Save the full and reduced log-likelihoods
    loglik.vec <- c(full$loglik, red$loglik)
    # Degrees of freedom = difference in model size
    df <- abs(length(cov1.name)-length(cov2.name))
    method <- ifelse(firth, "Penalized LRT", "Standard LRT")
  }
  names(loglik.vec) <- c("full", "reduced")
  lr.stat <- 2*abs(loglik.vec[2] - loglik.vec[1])
  lr.stat <- unname(lr.stat)
  p <- stats::pchisq(lr.stat, df, lower.tail = F)
  res <- list(statistic=lr.stat, pvalue=p, df=df, loglik=loglik.vec,
              method=method, call=match.call())
  class(res) <- "crrftest"
  return(res)
}
