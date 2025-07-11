# msea_ora_range : wrapper function of over-representation analysis using one of the following options:
#          "ora_full","bino_naive","bino_weighted", "bino_shrink"

#' Wrapper function for Over-Representation Analysis with p-value range estimation
#'
#' This function performs over-representation analysis (ORA) to assess metabolite set enrichment 
#' while considering uncertainty due to undetected metabolites. It wraps different methods for 
#' estimating a p-value range, including full enumeration and binomial resampling.
#'
#' @param SIG A character vector of statistically significant metabolite IDs.
#' @param DET A character vector of all detected metabolite IDs (background). Required for all 
#'   methods except \code{ora_full}.
#' @param M A named list of metabolite sets, each containing a character vector of metabolite IDs.
#' @param option Method to use for estimating the p-value range. One of \code{"ora_full"}, 
#'   \code{"bino_naive"}, \code{"bino_weighted"}, or \code{"bino_shrink"}.
#' @param probs Numeric vector of quantile probabilities for binomial simulation 
#'   (e.g., \code{c(0.025, 0.975)}). Ignored for \code{ora_full}.
#' @param nsim Number of simulations for binomial-based estimation. Ignored for \code{ora_full}.
#' @param lambda Shrinkage parameter for \code{"bino_shrink"} option.
#'
#' @return A list containing a matrix with p-value range results for each metabolite set. 
#' Columns include lower, median, and upper p-values.
#'
#' @details
#' This wrapper function allows switching between multiple ORA implementations that estimate the 
#' uncertainty due to undetected metabolites. The \code{ora_full} method uses exhaustive enumeration 
#' of all possible detection patterns, while the other methods use binomial resampling with different 
#' estimation strategies (naive, weighted, or shrinkage-based).
#'
#' @author Hiroyuki Yamamoto
#'
#' @examples
#' ## Example
#' M <- list(
#'   set1 = c("A", "B", "C"),
#'   set2 = c("C", "D", "E", "F"),
#'   set3 = c("G", "H")
#' )
#' DET <- c("A", "B", "C", "D", "E")  # Detected
#' SIG <- c("A", "C")                # Significant
#' 
#' msea_ora_range(SIG, DET, M, option = "bino_naive")
#' 
#' @keywords MSEA ORA enrichment
#' @export
msea_ora_range <- function(SIG, DET = NULL, M, option = "ora_full", probs = NULL, nsim = NULL, lambda = NULL) {
  if (option == "ora_full") {
    B <- ora_full(SIG, DET, M)
  } else if (option == "bino_naive") {
    B <- ora_bino(SIG, DET, M, method="naive", probs = c(0.025, 0.975), nsim = 1000)
  } else if (option == "bino_weighted") {
    B <- ora_bino(SIG, DET, M, method="weighted", probs = c(0.025, 0.975), nsim = 1000)
  } else if (option == "bino_shrink") {
    B <- ora_bino(SIG, DET, M, method="shrink", probs = c(0.025, 0.975), nsim = 1000, lambda = lambda)
  } else {
    stop("Invalid option. Use 'ora_full', 'bino_naive', 'bino_weighted', or 'bino_shrink'.")
  }
  return(B)
}
