toString.gconsensus <- function(x, ...) {
  ss1 <- x$subset
#  if (x$config$MC_use.HKSJ.adjustment) {
#    if (!is.null(x$w.i)) {
#      mm <- sum(ss1)
#	  wsse <- sum(x$w.i*(x$reported.data$mean[ss1] - x$mu)^2)
#      qq <- max(1, sqrt(1/(mm - 1)*wsse))
#    } else {
#      qq <- 1
#    }
#  } else {
#    qq <- 1
#  }

  str <- paste0("Study: ", x$study, "\n")
  str <- paste0(str, "Measurand: ", x$measurand, "\n")
  str <- paste0(str, "Number of included sources: ", sum(ss1), "\n")
  str <- paste0(str, "Consensus method: ", x$method, "\n\n")
  str2 <- capture.output(x$fit, file = NULL)
  for (i in 1:length(str2)) str <- paste0(str, (str2[i]), "\n")
  
  return(str)
}
