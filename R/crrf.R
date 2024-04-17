#' Firth-penalized Fine Gray model
#'
#' @param form     A formula with the survival model, with usual survival package syntax for the time to event data: Surv(time, event) notation. Values must be column names in dset.
#' @param etype    A character to indicate event type of interest, as defined by the event indicator in dset.
#' @param dset     A data.frame including the survival variables defined in the formula and any additional covariates to include in the model.
#' @param firth    A Boolean of whether to fit with Firth penalty. Default is TRUE.
#' @param CI       A Boolean of whether or not to compute the (penalized) partial profile likelihood confidence intervals. Default is FALSE.
#' @param alpha    A numeric of significance level. Default is 0.05.
#' @param eps      A numeric of threshold to begin confidence interval calculations. Default is eps=1e-10. Not recommended to decrease, as this may lead to errors in calculating the Firth Penalty.
#' @param CI.print A Boolean of whether to print the CI calculations for each variable while calculating to monitor progress. Default is FALSE.
#'
#' @return A crrf object, with a list of output, including a table with estimates and CI intervals if calculated and the output from nlminb optimization.
#' @import stats
#' @import survival
#' @export
#'
#' @examples
#' set.seed(12345)
#' n=100
#' tm=rexp(n)
#' ev=factor(rbinom(n,2,0.3))
#' x=rnorm(n)
#' y=rbinom(n,1,0.5)
#' u=runif(n)
#' dset=cbind.data.frame(tm=tm,ev=ev,x=x,y=y,u=u)
#' fgm.res=crrf(Surv(tm,ev)~x+y+u,etype="1",dset,firth=TRUE,CI=FALSE)
crrf <- function(form,etype,dset,firth=TRUE,CI=FALSE,alpha=0.05,eps=1e-10,
                 CI.print=F){
  fg.res <- maxlogL.fg(form,etype,dset,firth=firth)
  CI.tbl <- NULL
  if(CI){
    beta <- fg.res$beta
    k <- length(beta)
    CI.tbl <- matrix(NA,k,3)
    CI.clm <- c("ml.beta","lower","upper")
    colnames(CI.tbl) <- CI.clm
    rownames(CI.tbl) <- names(beta)
    for(i in 1:k){
      ci.res <- one.cb.fg(fg.res,i,alpha,eps)
      if(CI.print){print(ci.res)}
      CI.tbl[i,1:3] <- ci.res[1:3]
    }
  }

  fit <- list(coefficients=fg.res$beta,
              alpha=alpha, eps=eps, firth=fg.res$firth, etype=fg.res$etype,
              loglik=-1*fg.res$nlogL, iter=fg.res$iterations,
              conv=ifelse(fg.res$convergence==0, TRUE, FALSE),
              formula=fg.res$form, call=match.call(),
              CI.tbl=CI.tbl,fg.res=fg.res)
  class(fit) <- "crrf"
  return(fit)

}



##################################
# find confidence bound for one parameter

one.cb.fg <- function(ml.res,      # result of maxlogL.fg
                   index,       # index of parameter to compute bounds for
                   alpha=0.05,  # alpha for confidence bound
                   eps=1e-10,   # small value used to set lb and ub
                   n.iter=10)   # number of iterations to try modifying eps when produces error computing firth penalty
{
  beta <- ml.res$beta[index]      # extract beta of interest
  theta <- exp(beta)/(1+exp(beta)) # transform to unit scale
  delta <- stats::qchisq(1-alpha,1)/2

  # uniroot for lower bound on unit scale
  lb0 <- eps  # lower limit to test for lower CI bound
  for(i in 1:n.iter){
    nlogL.lb0 <- try(profile.nlogL.fg(lb0,index,ml.res,T), silent=TRUE) # test if profile penalized likelihood can be calculated
    if(inherits(nlogL.lb0, "try-error")){
      lb0 <- 10*lb0 # if profile penalized likelihood cannot be calculated, increase lower limit
    }
  }
  lb1 <- theta # MLE - upper limit

  root.lb <- stats::uniroot(profile.dif, c(lb0, lb1), index=index,
                     ml.res=ml.res, alpha=alpha)$root # find the root
  cilb <- log(root.lb)-log(1-root.lb) # convert back to original scale
  cilb.nlogL <- profile.nlogL.fg(root.lb,index,ml.res,T) # save profile penalized negative log likelihood at confidence bound

  # uniroot for upper bound on unit scale
  ub0 <- theta # MLE - lower limit
  ub1 <- 1-eps # upper limit to test for upper CI bound
  for(i in 1:n.iter){
    nlogL.ub1 <- try(profile.nlogL.fg(ub1,index,ml.res,T), silent=TRUE) # test if profile penalized likelihood can be caluclated
    if(inherits(nlogL.ub1, "try-error")){
      ub1 <- 1-10*(1-ub1) # if profile penalized likelihood cannot be caluclated, decrease upper limit
    }
  }

  root.ub <- stats::uniroot(profile.dif, c(ub0, ub1), index=index,
                     ml.res=ml.res, alpha=alpha)$root # find root
  ciub <- log(root.ub)-log(1-root.ub) # convert back to original scale
  ciub.nlogL <- profile.nlogL.fg(root.ub,index,ml.res,T) # save profile penalized negative log likelihood at confidence bound

  # organize output
  res <- c(ml.beta=beta,
           lower=cilb,
           upper=ciub,
           ml.nlogL=ml.res$nlogL,
           lb.nlogL=cilb.nlogL,
           ub.nlogL=ciub.nlogL,
           delta=delta,
           alpha=alpha)

  return(res)


}

#################################
# Profile likelihood CI function to find root
profile.dif <- function(beta,    # a scalar
                        index,   # the index of the beta in the result
                        ml.res,  # result of maxlogL.fg
                        alpha)   # confidence limit
{
  beta.ml <- ml.res$beta[index]      # extract beta of interest
  theta.ml <- exp(beta)/(1+exp(beta)) # transform to unit scale
  delta <- stats::qchisq(1-alpha,1)/2     # delta for CI calculation
  best.nlogL <- ml.res$nlogL
  nlogL.line <- best.nlogL+delta
  nlogL.temp <- profile.nlogL.fg(beta,index,ml.res,T)
  res <- nlogL.temp - nlogL.line
  return(res)
}


##################################
# Profile negative log-likelihood for Fine-Gray model

profile.nlogL.fg <- function(beta,         # a scalar
                             index,        # the index of the beta in the result
                             ml.res,       # result of maxlogL.fg
                             unit.scale=F) # indicates whether beta has been transformed to be on the unit scale
{
  if (unit.scale){
    beta <- log(beta)-log(1-beta)
  }

  all.beta <- ml.res$beta
  oth.beta <- all.beta[-index]

  ml.res <- maxlogL.fg(ml.res$form,
                       ml.res$etype,
                       ml.res$dset,
                       beta0=oth.beta,
                       firth=ml.res$firth,
                       fix.beta=beta,
                       fix.index=index)

  return(ml.res$nlogL)
}


###################################
# maximum likelihood for Fine-Gray model with optional Firth penalty

maxlogL.fg <- function(form,
                       etype,
                       dset,
                       beta0=NULL,
                       firth=F,
                       fix.beta=NULL,
                       fix.index=NULL)
{
  init.beta <- beta0
  notes <- NULL
  if (is.null(beta0)){
    fg.res <- fg.coxph(form,etype,dset)
    beta0 <- stats::coef(fg.res)

    if (!is.null(fix.index)){
      fix.beta <- beta0[fix.index]
      beta0 <- beta0[-fix.index]
    }
    init.beta <- beta0
    notes <- "beta0 computed by fg.coxph"
  }

  res <- stats::nlminb(beta0,
                       nlogL.fg,
                       form=form,
                       etype=etype,
                       dset=dset,
                       firth=firth,
                       fix.beta=fix.beta,
                       fix.index=fix.index)

  res$beta0 <- init.beta
  res$firth <- firth
  res$fix.beta <- fix.beta
  res$fix.index <- fix.index
  res$beta <- construct.beta(res$par,fix.beta,fix.index)
  res$nlogL <- res$objective
  res$notes <- notes
  res$dset <- dset
  res$form <- form
  res$etype <- etype

  return(res)
}


#################################
# negative log-likelihood of Fine-Gray model with optional Firth penalty

nlogL.fg <- function(beta,            # vector of coefficients allowed to vary in likelihood optimization
                     dset,            # data set
                     form,            # model formula for finegray function of survival package
                     etype,           # event of interest
                     firth=F,         # indicates whether to add a Firth penalty term
                     fix.beta=NULL,   # scalar or vector for a coefficient to fix for confidence interval calculation
                     fix.index=NULL)  # indices of fixed coefficients
{
  # construct the beta vector
  new.beta <- construct.beta(beta,fix.beta,fix.index)                              # initialize as beta

  # use coxph to compute likelihood for Fine-Gray model given new.beta

  fg.res0 <- fg.coxph(form,etype=etype,data=dset,
                      init=new.beta,                     # specified beta values
                      control=coxph.control(iter.max=0)) # no iterative likelihood optimization

  fpen <- firth.penalty(fg.res0,firth)

  logL <- fg.res0$loglik[1]              # this is the likelihood at new.beta


  res <- -logL-fpen

  return(res)
}


###############################
# compute Firth penalty

firth.penalty <- function(fg.res,firth=F)

{
  if(!firth) return(0)                # return 0 if no penalty

  fpen <- NA
  fg.details <- survival::coxph.detail(fg.res)    # extract coxph details
  info.mat <- fg.details$imat           # get the information matrix (for each time point)

  if(is.vector(info.mat)){
    fpen <- 0.5*log(sum(info.mat))
  }

  if(is.matrix(info.mat)){
    tot.imat <- info.mat
    fpen <- 0.5*log(det(tot.imat))
  }

  if(is.array(info.mat))
  {
    dim.imat <- dim(info.mat)
    tot.imat <- matrix(0,dim.imat[1],      # initialize the total information matrix as zero
                       dim.imat[2])
    for (i in 1:dim.imat[3]){
      tot.imat <- tot.imat+info.mat[,,i]
    }            # compute total information matrix as sum of risk-set specific information matrix

    fpen <- 0.5*log(det(tot.imat))
  }

  if (is.na(fpen))
    stop("Technical error computing Firth penalty.")

  return(fpen)
}


################################
# wrapper to use coxph to fit Fine-Gray model

fg.coxph <- function(form,etype,data,init,control=survival::coxph.control())
{
  # prepare data to use coxph to fit Fine-Gray model
  split.form <- unlist(strsplit(as.character(form), split="~"))
  fg.form <- stats::as.formula(paste0(split.form[2],"~."))
  fgdat <- survival::finegray(fg.form,etype=etype,data=data)

  # prepare formula for Fine-Gray model
  split.form <- unlist(strsplit(as.character(form), split="~"))
  new.form <- stats::as.formula(paste0("Surv(fgstart, fgstop, fgstatus)~",split.form[3]))

  fg.res <- survival::coxph(new.form,        # formula in Fine-Gray format for coxph
                            data=fgdat,      # data in Fine-Gray format for coxph
                            weight=fgdat$fgwt,    # Fine-Gray weights for coxph
                            init=init,       # specified initial values for beta
                            control=control) # controls for iterative optimization

  return(fg.res)
}


################################
# construct beta vector

construct.beta <- function(beta,
                           fix.beta=NULL,
                           fix.index=NULL)
{
  new.beta <- beta
  if(!is.null(fix.beta))                     # if fix.beta is provided, insert it into new.beta
  {
    if (length(fix.index)!=length(fix.beta))
      stop("length(fix.beta) must equal length(fix.index)")

    fix.index <- fix.index-0.5                   # subtract from fixed index so it will insert appropriately later
    new.beta <- c(beta,fix.beta)                 # append fix.beta to beta
    new.index <- c(1:length(beta),fix.index)     # append fix.index to indices of beta
    new.order <- order(new.index)                # order by new.index
    new.beta <- new.beta[new.order]              # this is the new beta
  }
  return(new.beta)
}

