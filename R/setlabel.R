#' Binary labeling of metabolite IDs across metabolite sets
#'
#' Assigns binary labels (1 or 0) to indicate whether each metabolite ID is included in each metabolite set. This function is used for over-representation analysis (ORA) in metabolomics or other omics data.
#'
#' @param M_ID A character vector of metabolite IDs (e.g., detected or significant metabolites).
#' @param M A named list of metabolite sets, where each element is a character vector of metabolite IDs.
#' @param option Optional argument. If \code{"anno"}, input IDs are assumed to be annotation strings (e.g., "C00031,C00022") and are split by semicolon after replacing commas. If \code{NULL} (default), IDs are used as-is.
#'
#' @return A binary matrix with rows representing input metabolite IDs and columns representing metabolite sets. Each cell contains 1 if the metabolite ID is found in the set, and 0 otherwise.
#'
#' @details
#' When \code{option = "anno"}, each element of \code{M_ID} may contain multiple IDs separated by commas (e.g., from an annotation field). These are split and matched against the metabolite sets.
#'
#' This function is used internally in various ORA methods (e.g., \code{ora_det}, \code{ora_all}, \code{ora_est}) to compute contingency tables for enrichment analysis.
#'
#' @author Hiroyuki Yamamoto
#'
#' @examples
#' # Example usage
#' M <- list(
#'   set1 = c("A", "B", "C"),
#'   set2 = c("C", "D"),
#'   set3 = c("E", "F")
#' )
#'
#' M_ID <- c("A", "C", "E", "G")
#'
#' setlabel(M_ID, M)
#'
#' @keywords MSEA labeling ORA
#' @export
setlabel <- function(M_ID, M, option = NULL) {
  n <- length(M_ID)
  p <- length(M)

  L <- matrix(0, nrow = n, ncol = p)
  colnames(L) <- names(M)

  for (i in seq_len(p)) {
    m <- unique(as.character(unlist(M[[i]])))

    for (j in seq_len(n)) {
      if (is.null(option)){
        b <- M_ID[j]
      }
      else if (option == "anno") {
        a <- chartr(",", ";", M_ID[j])
        b <- unique(unlist(strsplit(a, ";")))
      } else {
        stop("Unknown option: ", option)
      }

      if (any(b %in% m)) {
        L[j, i] <- 1
      }
    }
  }

  return(L)
}
