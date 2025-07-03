#' Over-representation analysis adjusting for undetected metabolites
#'
#' This function performs metabolite set enrichment analysis using over-representation analysis (ORA), incorporating point estimates to adjust for potentially undetected metabolites.
#' It supports three estimation methods: naive, weighted, and shrinkage-based adjustment.
#'
#' @param SIG Character vector of metabolite IDs considered statistically significant.
#' @param DET Character vector of all detected metabolite IDs (background set).
#' @param M A named list, where each element is a metabolite set (e.g., pathway) containing character vectors of metabolite IDs.
#' @param method A character string specifying the estimation method to use. One of \code{"naive"}, \code{"weighted"}, or \code{"shrink"}. Default is \code{"naive"}.
#' @param lambda A numeric value used in the \code{"shrink"} method as a shrinkage parameter. Default is \code{5}.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{Result of MSEA (ORA with adjustment)}: A matrix with raw p-values and adjusted q-values (BH correction) for each metabolite set.}
#' }
#'
#' @details
#' This function estimates the impact of undetected metabolites on enrichment results. It builds upon the ORA results from detected metabolites, then adjusts the contingency tables by estimating how many undetected metabolites might be significant, based on a specified method.
#'
#' @author Hiroyuki Yamamoto
#'
#' @references
#' Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.\cr
#' Global functional profiling of gene expression.\cr
#' \emph{Genomics}. 2003 Feb;81(2):98-104.
#'
#' @examples
#' # Example using metabolomics data
#' data(fasting)
#' data(pathway)
#'
#' # Compute PCA loadings
#' pca <- prcomp(fasting$X, scale=TRUE)
#' pca <- pca_loading(pca)
#'
#' # Detected and significant metabolites
#' metabolites <- colnames(fasting$X)
#' SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]
#' DET <- metabolites
#' M <- pathway$fasting
#'
#' # Run adjusted ORA (e.g., shrinkage method)
#' B <- ora_est(SIG, DET, M, method = "shrink", lambda = 5)
#' B$`Result of MSEA (ORA with adjustment)`
#'
#' @keywords msea
#' @export
ora_est <- function(SIG, DET, M, method="naive", lambda = 5) {

  ALL0 <- unique(unlist(M)) # 対象物質の情報
  ALL <- unique(c(ALL0, DET))

  # Step 1: Label assignment for each metabolite group
  L1 <- setlabel(ALL, M) # All metabolites
  L2 <- setlabel(DET, M) # Detected metabolites
  L3 <- setlabel(SIG, M) # Significant metabolites

  # Step 2: Perform base ORA on detected metabolites
  B <- ora_det(SIG, DET, M)

  # Initialize vectors to store results
  P <- NULL
  P_range <- NULL

  # Baseline significance proportion from detected metabolites
  p <- length(SIG) / length(DET)

  # Step 3: Loop through each pathway
  for (i in 1:length(M)) {
    # Counts for each pathway
    l1 <- colSums(L1)[i]  # Total metabolites
    l2 <- colSums(L2)[i]  # Detected metabolites
    l3 <- colSums(L3)[i]  # Significant metabolites

    n <- l1 - l2  # Count of undetected metabolites

    # 推定r: naive or weighted or shrink
    if (method == "naive") {
      r <- if (l2 > 0) l3 / l2 else p
    } else if (method == "weighted") {
      r <- (l3 + n * p) / l1
    } else if (method == "shrink"){
      r <- if (l2 + lambda > 0) (l3 + lambda * p) / (l2 + lambda) else p
    }

    # Construct adjusted 2x2 table
    a <- round(B$TAB[[i]][1,1] + n * r)
    b <- round(B$TAB[[i]][1,2] + n * (1 - r))

    # Count of non-significant detected and total substances
    c <- max(0, round(length(ALL) * p - a))
    d <- max(0, round(length(ALL) * (1 - p) - b))

    # Perform Fisher's test in the default case
    tab <- matrix(c(a, b, c, d), nrow = 2)
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- resfish$p.value

    print(i)
    print(tab)

  }

  # Adjust p-values for multiple testing (default)
  Q <- p.adjust(P, method = "BH")
  PQ <- cbind(P, Q)
  rownames(PQ) <- names(M)
  colnames(PQ) <- c("p.value", "q.value")
  result <- list(PQ)
  names(result) <- c("Result of MSEA (ORA with adjustment)")

  return(result)
}
