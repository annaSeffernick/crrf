#' @exportS3Method print crrf
print.crrf <- function(x,...)
{
  firth.res <- unclass(x)
  firth.res <- as.list(firth.res)
  method <- ifelse(x$firth, "Firth penalized Fine-Gray", "Fine-Gray")

  cat("Model fitted by", method, "\n\n")
  if(is.null(x$CI.tbl)){
    cat("\n No Confidence Intervals Calculated \n")
    out <- cbind(x$coefficients, exp(x$coefficients))
    dimnames(out) <- list(names(x$coefficients), c("coef", "exp(coef)"))
    print(out)
  }
  else{
    cat("\n Confidence Intervals Calcualted by Profile Likelihood \n\n")
    out <- cbind(x$coefficients,
                 x$CI.tbl[,2], x$CI.tbl[,3],
                 exp(x$coefficients), exp(x$CI.tbl[,2]), exp(x$CI.tbl[,3]))
    dimnames(out) <- list(names(x$coefficients), c("coef",
                                                   paste(c("lower", "upper"), 1 - x$alpha),
                                                   "exp(coef)",
                                                   paste(c("exp lower", "exp upper"), 1 - x$alpha)))
    print(out)
  }
  invisible(x)

}
