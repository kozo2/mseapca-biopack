#' ORA using all metabolites in the pathway
#'
#' This function performs over-representation analysis (ORA) using all metabolites present in the given metabolite set list as the background, without specifying a reference metabolome. This corresponds to the behavior of MetaboAnalyst when no reference metabolome is uploaded.
#'
#' @param SIG Character vector of significant metabolites
#' @param M Named list of metabolite sets
#'
#' @return A list with:
#' \item{\code{Result of MSEA(ORA)}}{Matrix of p-values and q-values}
#' \item{\code{significant metabolites}}{List of significant IDs per set}
#' \item{\code{Contingency tables}}{A list of 2×2 contingency tables used in Fisher's exact tests.}
#'
#' @author Hiroyuki Yamamoto
#'
#' @references
#' Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.\cr
#' Global functional profiling of gene expression.\cr
#' \emph{Genomics}. 2003 Feb;81(2):98-104.\cr\cr
#' Yamamoto H, Fujimori T, Sato H, Ishikawa G, Kami K, Ohashi Y,\cr
#' Statistical hypothesis testing of factor loading in principal component analysis and its application to metabolite set enrichment analysis.\cr
#' BMC Bioinformatics, (2014) 15(1):51.\cr\cr
#' Yamamoto H.\cr
#' Probabilistic Over-Representation Analysis for Metabolite Set Enrichment Analysis Considering Undetected Metabolites", Jxiv, (2024).
#'
#' @examples
#' # Example1 : Metabolome data
#' data(fasting)
#' data(pathway)
#'
#' # pca and pca loading
#' pca <- prcomp(fasting$X, scale=TRUE)
#' pca <- pca_loading(pca)
#'
#' # all detected metabolites
#' metabolites <- colnames(fasting$X)
#'
#' # statistically significant negatively correlated metabolites in PC1 loading
#' SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]
#'
#' # metabolite set list
#' M <- pathway$fasting
#'
#' # MSEA by over representation analysis
#' B <- ora_all(SIG, M)
#' B$`Result of MSEA(ORA)`
#'
#' @export
ora_all <- function (SIG, M)
{
  ALL <- unique(as.character(unlist(M)))
  SIG <- as.character(as.matrix(SIG))
  SIG <- SIG[SIG %in% ALL]
  num_all <- length(ALL)
  num_sig <- length(SIG)
  Lall0 <- setlabel(ALL, M)
  Lall <- Lall0[, colSums(Lall0) != 0]
  if (ncol(Lall) < 2) {
    stop("more than two metabolite set are necessary")
    return()
  }
  Lsig <- setlabel(SIG, M)
  Lsig <- Lsig[, colSums(Lall0) != 0]
  l <- colSums(Lall0) != 0
  P <- NaN
  TAB <- NULL
  for (i in 1:sum(l)) {
    a1 <- sum(Lsig[, i]) # ここはok
    a2 <- sum(Lall[, i]) - sum(Lsig[, i]) # ここもok
    a3 <- length(SIG) - a1 # ここはok
    a4 <- (length(ALL) - length(SIG)) - a2 # ここも駄目
    tab <- t(matrix(c(a1, a2, a3, a4), 2))
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- resfish$p.value
    TAB[i] <- list(tab)
  }
  names(TAB) <- colnames(Lsig)
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
  RES <- list(PQ, LES, TAB)
  names(RES) <- c("Result of MSEA(ORA)", "significant metabolites", "Contingency tables")
  return(RES)
}
