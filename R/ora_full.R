#' Over-representation analysis with full enumeration of undetected metabolite patterns
#'
#' This function performs over-representation analysis (ORA) by enumerating all possible patterns of significant and non-significant assignments among undetected metabolites for each metabolite set. It returns the minimum, median, and maximum p-values from Fisher's exact tests across these patterns, thereby estimating the full uncertainty range due to undetected metabolites.
#'
#' @param SIG A character vector of statistically significant metabolite IDs.
#' @param DET A character vector of all detected metabolite IDs (the background).
#' @param M A named list of metabolite sets, where each element is a character vector of metabolite IDs.
#'
#' @return A list containing:
#' \item{\code{Range of p-values}}{
#'   A matrix with rows corresponding to metabolite sets and three columns:
#'   \code{lower p-value}, \code{p-value(median)}, and \code{upper p-value}.
#' }
#'
#' @details
#' For each metabolite set, the number of undetected metabolites is calculated.
#' The function then considers all possible numbers of significant metabolites (from 0 to the total number of undetected ones) among those undetected.
#' For each case, a 2x2 contingency table is constructed and Fisher's exact test is applied.
#' The resulting p-values are aggregated to report the minimum, median, and maximum values.
#'
#' @author Hiroyuki Yamamoto
#'
#' @references
#' Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.\cr
#' Global functional profiling of gene expression. \emph{Genomics}. 2003;81(2):98–104.
#'
#' @examples
#' ## Example: small set
#' M <- list(
#'   set1 = c("A", "B", "C"),
#'   set2 = c("C", "D", "E", "F"),
#'   set3 = c("G", "H")
#' )
#' DET <- c("A", "B", "C", "D", "E")  # Detected
#' SIG <- c("A", "C")                # Significant
#'
#' ora_full(SIG, DET, M)$`Range of p-values`
#'
#' @keywords MSEA ORA enrichment
#' @export
ora_full <- function(SIG, DET, M) {

  ALL <- unique(as.character(unlist(M)))

  # Step 1: Label assignment for each metabolite group
  L1 <- setlabel(ALL, M) # All metabolites
  L2 <- setlabel(DET, M) # Detected metabolites
  L3 <- setlabel(SIG, M) # Significant metabolites

  # Step 2: Perform base ORA on detected metabolites
  B <- ora_det(SIG, DET, M)

  # Initialize vectors to store results
  P <- NULL
  P_range <- NULL

  # Step 3: Loop through each pathway
  for (i in 1:length(M)) {
    # Counts for each pathway
    l1 <- colSums(L1)[i]  # Total metabolites
    l2 <- colSums(L2)[i]  # Detected metabolites
    l3 <- colSums(L3)[i]  # Significant metabolites

    # Proportion of significant among detected, and estimate for undetected
    r <- l3 / l2
    n <- l1 - l2  # Count of undetected metabolites

    # Construct adjusted 2x2 table
    a <- round(B$TAB[[i]][1,1] + n * r)
    b <- round(B$TAB[[i]][1,2] + n * (1 - r))

    # Baseline significance proportion from detected metabolites
    p <- length(SIG) / length(DET)

    # Count of non-significant detected and total substances
    c <- max(0, round(length(ALL) * p - a))
    d <- max(0, round(length(ALL) * (1 - p) - b))

    # Conditional branching based on option
    # Calculate p-value range
    possible_values <- 0:n  # Full range of undetected metabolites
    p_values <- sapply(possible_values, function(x) {
      # Construct 2x2 table for all patterns
      a_var <- l3 + x
      b_var <- l2 - l3 + (n - x)
      c_var <- max(0, round(length(ALL) * p - a_var))
      d_var <- max(0, round(length(ALL) * (1 - p) - b_var))

      tab_var <- matrix(c(a_var, b_var, c_var, d_var), nrow = 2)

    # Fisher's exact test for each pattern if valid table
      if (all(tab_var >= 0) && all(is.finite(tab_var))) {
        fisher.test(tab_var, alternative = "greater")$p.value
      } else {
        NA  # Invalid table entries
      }
    })

    # Obtain the range of p-values (minimum and maximum)
    p_values <- na.omit(p_values)  # Remove NAs from invalid tables
    p_min <- min(p_values)
    p_max <- max(p_values)

    # Store the range result
    P_range <- rbind(P_range, c(p_min, median(p_values), p_max))
  }

  # Output the range of lower and upper p-values
  rownames(P_range) <- names(M)
  colnames(P_range) <- c("lower p-value", "p-value(median)", "upper p-value")

  result <- list(P_range)
  names(result) <- c("Range of p-values")

  return(result)
}
