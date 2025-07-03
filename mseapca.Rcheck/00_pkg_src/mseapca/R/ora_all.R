# ora_all: over-representation analysis under the same conditions as MetaboAnalyst,
# when no reference metabolome is specified.

#source("C:/Users/yamamoto/Documents/R/msea/mseapca/dev/mseapca/R/setlabel.R")

#' ORA using all metabolites in the pathway
#'
#' This function performs over-representation analysis (ORA) using all metabolites present in the 
#' given metabolite set list as the background, without specifying a reference metabolome. 
#' This corresponds to the behavior of MetaboAnalyst when no reference metabolome is uploaded.
#'
#' @param SIG Character vector of significant metabolites
#' @param M Named list of metabolite sets
#'
#' @return A list with:
#'   \item{\code{Result of MSEA(ORA)}}{Matrix of p-values and q-values}
#'   \item{\code{significant metabolites}}{List of significant IDs per set}
#'   \item{\code{Contingency tables}}{A list of 2×2 contingency tables used in Fisher's exact tests.}
#'
#' @references
#' Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.
#' Global functional profiling of gene expression.
#' Genomics. 2003 Feb;81(2):98-104.
#' 
#' Yamamoto H, Fujimori T, Sato H, Ishikawa G, Kami K, Ohashi Y,
#' Statistical hypothesis testing of factor loading in principal component analysis and its 
#' application to metabolite set enrichment analysis.
#' BMC Bioinformatics, (2014) 15(1):51.
#' 
#' Yamamoto H.
#' Probabilistic Over-Representation Analysis for Metabolite Set Enrichment Analysis 
#' Considering Undetected Metabolites", Jxiv, (2024).
#'
#' @author Hiroyuki Yamamoto
#'
#' @examples
#' # Example1 : Metabolome data
#' data(MetaboAnalyst)
#' 
#' # Use pre-computed significant metabolites
#' SIG <- MetaboAnalyst$fasting$data$SIG
#' 
#' # metabolite set list
#' M <- MetaboAnalyst$fasting$pathway
#' 
#' # MSEA by over representation analysis
#' B <- ora_all(SIG, M)
#' B$`Result of MSEA(ORA)`
#' 
#' ## Example2 : Proteome data (commented out - data not available)
#' # data(covid19)
#' # data(pathway)
#' # 
#' # X <- covid19$X$proteomics
#' # Y <- covid19$Y
#' # D <- covid19$D
#' # tau <- covid19$tau
#' # 
#' # protein_name <- colnames(X)
#' # 
#' # # pls-rog and pls-rog loading
#' # plsrog <- pls_rog(X,Y,D)
#' # plsrog <- plsrog_loading(plsrog)
#' # 
#' # # statistically significant proteins
#' # index_prot <- which(plsrog$loading$R[,1]>0 & plsrog$loading$p.value[,1]<0.05)
#' # sig_prot <- protein_name[index_prot]
#' # 
#' # # protein set list
#' # M <- pathway$covid19$proteomics
#' # 
#' # # MSEA by over representation analysis
#' # B <- ora_all(sig_prot, M)
#' # B$`Result of MSEA(ORA)`
#' 
#' # Example3: Metabolome data
#' data(MetaboAnalyst)
#' 
#' SIG <- MetaboAnalyst$fasting$data$SIG
#' M <- MetaboAnalyst$fasting$pathway
#' 
#' # Perform ORA using detected metabolites only
#' B <- ora_all(SIG, M)
#' B$`Result of MSEA(ORA)`
#' 
#' @keywords htest
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
