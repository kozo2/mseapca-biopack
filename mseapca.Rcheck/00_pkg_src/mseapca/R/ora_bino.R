# ora_bino : ORA that adjusts for undetected metabolites using binomial resampling.

#source("C:/Users/yamamoto/Documents/R/msea/mseapca/dev/mseapca/R/ora_det.R")

#' Over-representation analysis with binomial resampling adjustment
#'
#' Performs ORA while adjusting for undetected metabolites by binomial resampling.
#'
#' @param SIG Character vector of significant metabolites.
#' @param DET Character vector of detected metabolites (background).
#' @param M Named list of metabolite sets.
#' @param method \code{"naive"}, \code{"weighted"}, or \code{"shrink"}.
#' @param probs Quantiles for the confidence interval (default 95\%).
#' @param nsim Number of binomial simulations (default 1000).
#' @param lambda Shrinkage constant used when \code{method = "shrink"}.
#'
#' @return A list containing one matrix.
#' Rows = metabolite sets; columns = lower, median, and upper p-values.
#'
#' @author Hiroyuki Yamamoto
#'
#' @examples
#' M   <- list(set1 = c("A","B","C"), set2 = c("C","D","E"))
#' DET <- c("A","B","C","D")
#' SIG <- c("A","C")
#' 
#' # Binomial-adjusted ORA (shrink method)
#' ora_bino(SIG, DET, M, method = "shrink", nsim = 200)
#' 
#' @keywords MSEA ORA simulation
#' @export
ora_bino <- function(SIG, DET, M, method="naive", probs = c(0.025, 0.975), nsim = 1000, lambda = 5) {

  # Set a random seed for reproducibility
  set.seed(1)

  ALL <- unique(as.character(unlist(M)))

  # Perform base ORA on detected metabolites
  B <- ora_det(SIG, DET, M)

  # Initialize vectors to store results
  P <- NULL
  P_range <- NULL

  p <- length(SIG) / length(DET)

  # Perform calculations for each pathway
  for (i in 1:length(M)) {

    # Get metabolite counts for each pathway
    l1 <- sum(ALL %in% M[[i]])  # Total metabolites in the pathway
    l2 <- sum(DET %in% M[[i]])   # Detected metabolites in the pathway
    l3 <- sum(SIG %in% M[[i]])   # Significant metabolites in the pathway

    # Directly use the proportion of significant metabolites among detected metabolites
    n <- l1 - l2  # Number of undetected metabolites in the pathway
    if (method == "naive") {
      r <- if (l2 > 0) l3 / l2 else p
    } else if (method == "weighted") { # weighted
      r <- (l3 + n * p) / l1
    } else if (method == "shrink"){
      r <- if (l2 + lambda > 0) (l3 + lambda * p) / (l2 + lambda) else p
    }

    # Initial 2x2 table based on undetected metabolites
    a <- B$TAB[[i]][1,1]
    b <- B$TAB[[i]][1,2]

    # Overall count of non-significant metabolites
    c <- round(length(ALL) * p - a)
    d <- round(length(ALL) * (1 - p) - b)

    # Resampling based on binomial distribution and p-value range calculation
    simulated_p_values <- numeric(nsim)

    for (j in 1:nsim) {
      # Resample significant metabolites among undetected ones using a binomial distribution
      sampled_significant <- rbinom(1, size = n, prob = r)

      # Reconstruct the 2x2 table based on resampling
      a_var <- a + sampled_significant
      b_var <- b + (n - sampled_significant)
      c_var <- round(length(ALL) * p - a_var)
      d_var <- round(length(ALL) * (1 - p) - b_var)

      tab_var <- t(matrix(c(a_var, b_var, c_var, d_var), nrow = 2))

      # Calculate p-value using Fisher's exact test
      simulated_p_values[j] <- fisher.test(tab_var, alternative = "greater")$p.value
    }

    # Obtain the range of p-values from simulations
    p_min <- quantile(simulated_p_values, probs = probs[1])
    p_median <- median(simulated_p_values)
    p_max <- quantile(simulated_p_values, probs = probs[2])

    # Store the range of p-values
    P_range <- rbind(P_range, c(p_min, p_median, p_max))
  }

  # Set row and column names for p-value range output
  rownames(P_range) <- names(M)
  colnames(P_range) <- c("lower p-value", "p-value(median)","upper p-value")

  # Display results
  list("Range of p-values for each pathway (95% confidence interval)" = P_range)

}
