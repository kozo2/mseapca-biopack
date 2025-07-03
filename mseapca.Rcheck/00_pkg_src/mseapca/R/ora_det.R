# ora_det : Over-representation analysis (ORA) using detected metabolites only.

#source("C:/Users/yamamoto/Documents/R/msea/mseapca/dev/mseapca/R/setlabel.R")

#' Over-representation analysis using detected metabolites only
#'
#' This function performs metabolite set enrichment analysis using over-representation analysis (ORA) 
#' under the assumption that only detected metabolites are used as the background.
#' A one-sided Fisher's exact test is applied to each metabolite set.
#'
#' @param SIG Character vector of metabolite IDs considered statistically significant.
#' @param DET Character vector of all detected metabolite IDs (background set).
#' @param M A named list, where each element is a metabolite set (e.g., pathway) containing 
#'   character vectors of metabolite IDs.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\code{Result of MSEA(ORA)}: A matrix with raw p-values and adjusted q-values 
#'       (BH correction) for each metabolite set.}
#'     \item{\code{significant metabolites}: A list of significant metabolites overlapping 
#'       with each metabolite set.}
#'     \item{\code{TAB}: A list of 2x2 contingency tables used for Fisher's exact tests.}
#'   }
#'
#' @author Hiroyuki Yamamoto
#'
#' @references
#' Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.
#' Global functional profiling of gene expression.
#' \emph{Genomics}. 2003 Feb;81(2):98-104.
#'
#' @examples
#' # Example: Simple ORA with metabolomics data
#' data(MetaboAnalyst)
#' 
#' # Use pre-computed significant and detected metabolites
#' SIG <- MetaboAnalyst$fasting$data$SIG
#' DET <- MetaboAnalyst$fasting$data$DET
#' M <- MetaboAnalyst$fasting$pathway
#' 
#' # Perform ORA using detected metabolites only
#' B <- ora_det(SIG, DET, M)
#' B$`Result of MSEA(ORA)`
#' 
#' @keywords msea
#' @export
ora_det <- function (SIG, DET, M)
{
  DET <- as.character(as.matrix(DET))
  SIG <- as.character(as.matrix(SIG))
  num_all <- length(DET)
  num_sig <- length(SIG)
  Lall0 <- setlabel(DET, M)
  #Lall <- Lall0[, colSums(Lall0) != 0]
  Lall <- Lall0
  if (ncol(Lall) < 2) {
    stop("more than two metabolite set are necessary")
    return()
  }
  Lsig <- setlabel(SIG, M)
  #Lsig <- Lsig[, colSums(Lall0) != 0]
  #l <- colSums(Lall0) != 0
  P <- NaN
  TAB <- NULL
  for (i in 1:length(M)) {
    a1 <- sum(Lsig[, i])
    a2 <- sum(Lall[, i]) - sum(Lsig[, i])
    a3 <- length(SIG) - a1
    a4 <- (length(DET) - length(SIG)) - a2
    tab <- t(matrix(c(a1, a2, a3, a4), 2))
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- resfish$p.value
    TAB[i] <- list(tab)
  }
  Q <- p.adjust(P, method = "BH")
  LES <- NaN
  for (i in 1:ncol(Lsig)) {
    les <- SIG[Lsig[, i] == 1]
    LES[i] <- list(les)
  }
  names(LES) <- colnames(Lsig)
  PQ <- cbind(P, Q)
  rownames(PQ) <- colnames(Lsig)
  colnames(PQ) <- c("p.value", "q.value")
  RES <- list(PQ, LES,TAB)
  names(RES) <- c("Result of MSEA(ORA)", "significant metabolites", "TAB")
  return(RES)
}
