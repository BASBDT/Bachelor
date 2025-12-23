# R/rcpp_numeric.R

#' Internal: number of threads for parallel ops
#' @keywords internal
.tbc_threads <- function() {
  n <- getOption("TBC.threads", NA_integer_)
  if (is.na(n)) {
    if (requireNamespace("RcppParallel", quietly = TRUE)) {
      return(max(1L, getNamespace("RcppParallel")$defaultNumThreads()))
    }
    return(1L)
  }
  as.integer(max(1L, n))
}

#' Internal: set thread pool size for RcppParallel (optional)
#' @keywords internal
.tbc_set_threads <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0L) return(invisible(NULL))
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    RcppParallel::setThreadOptions(numThreads = n)
  }
  invisible(NULL)
}
