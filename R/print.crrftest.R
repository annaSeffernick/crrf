#' @exportS3Method print crrftest
print.crrftest <- function(x, ...)
{
  # model call
  print(x$call)
  # print method
  cat("Model Fitted By", x$method, "\n\n")
  # print likelihoods
  out <- c(x$loglik[2], x$loglik[1], x$statistic/2)
  names(out) <- c("Reduced Model", "Full Model", "Difference")
  cat("\n Likelihoods: \n")
  print(out)
  # result summary
  cat("\n Likelihood ratio test=", x$statistic, " on ", x$df,
      " df, p=", x$p, "\n", sep="")
  invisible(x)
}
