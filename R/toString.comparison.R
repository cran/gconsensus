toString.comparison <- function(x, ...) {
  ss1 <- x$gconsensus$subset

  str <- paste0("Consensus Comparison by Statistical Method\n")
  str <- paste0(str, "Exercise: ", x$gconsensus$exercise, "\n")
  str <- paste0(str, "Measurand: ", x$gconsensus$measurand, "\n")
  str <- paste0(str, "Number of included sources: ", sum(ss1), "\n\n")
  
  str2 <- capture.output(x$fit, file = NULL)
  for (i in 1:length(str2)) str <- paste0(str, (str2[i]), "\n")

  return(str)
}
