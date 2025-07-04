% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ora_det.R
\name{ora_det}
\alias{ora_det}
\title{Over-representation analysis using detected metabolites only}
\usage{
ora_det(SIG, DET, M)
}
\arguments{
\item{SIG}{Character vector of metabolite IDs considered statistically significant.}

\item{DET}{Character vector of all detected metabolite IDs (background set).}

\item{M}{A named list, where each element is a metabolite set (e.g., pathway) containing 
character vectors of metabolite IDs.}
}
\value{
A list containing:
  \itemize{
    \item{\code{Result of MSEA(ORA)}: A matrix with raw p-values and adjusted q-values 
      (BH correction) for each metabolite set.}
    \item{\code{significant metabolites}: A list of significant metabolites overlapping 
      with each metabolite set.}
    \item{\code{TAB}: A list of 2x2 contingency tables used for Fisher's exact tests.}
  }
}
\description{
This function performs metabolite set enrichment analysis using over-representation analysis (ORA) 
under the assumption that only detected metabolites are used as the background.
A one-sided Fisher's exact test is applied to each metabolite set.
}
\examples{
# Example: Simple ORA with metabolomics data
data(MetaboAnalyst)

# Use pre-computed significant and detected metabolites
SIG <- MetaboAnalyst$fasting$data$SIG
DET <- MetaboAnalyst$fasting$data$DET
M <- MetaboAnalyst$fasting$pathway

# Perform ORA using detected metabolites only
B <- ora_det(SIG, DET, M)
B$`Result of MSEA(ORA)`

}
\references{
Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.
Global functional profiling of gene expression.
\emph{Genomics}. 2003 Feb;81(2):98-104.
}
\author{
Hiroyuki Yamamoto
}
\keyword{msea}
