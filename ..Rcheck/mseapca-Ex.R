pkgname <- "mseapca"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mseapca')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MetaboAnalyst")
### * MetaboAnalyst

flush(stderr()); flush(stdout())

### Name: MetaboAnalyst
### Title: Example metabolomics data extracted from MetaboAnalyst (fasting
###   study)
### Aliases: MetaboAnalyst
### Keywords: datasets metabolomics

### ** Examples

## Load data
data(MetaboAnalyst)

SIG <- MetaboAnalyst$fasting$data$SIG
DET <- MetaboAnalyst$fasting$data$DET
M   <- MetaboAnalyst$fasting$pathway

## Simple ORA on detected metabolites
res <- ora_det(SIG, DET, M)
head(res$`Result of MSEA(ORA)`)




cleanEx()
nameEx("msea_ora")
### * msea_ora

flush(stderr()); flush(stdout())

### Name: msea_ora
### Title: Wrapper for several over-representation analysis (ORA) flavours
### Aliases: msea_ora
### Keywords: MSEA ORA enrichment

### ** Examples

## Toy data ----------------------------------------------------
M   <- list(set1 = c("A","B","C"),
            set2 = c("C","D","E"))
DET <- c("A","B","C","D")
SIG <- c("A","C")

## Plain ORA with detected metabolites
msea_ora(SIG, DET, M, option = "det")

## ORA with shrinkage adjustment
msea_ora(SIG, DET, M, option = "est_shrink", lambda = 5)




cleanEx()
nameEx("msea_ora_range")
### * msea_ora_range

flush(stderr()); flush(stdout())

### Name: msea_ora_range
### Title: Wrapper function for Over-Representation Analysis with p-value
###   range estimation
### Aliases: msea_ora_range
### Keywords: MSEA ORA enrichment

### ** Examples

## Example
M <- list(
  set1 = c("A", "B", "C"),
  set2 = c("C", "D", "E", "F"),
  set3 = c("G", "H")
)
DET <- c("A", "B", "C", "D", "E")  # Detected
SIG <- c("A", "C")                # Significant

msea_ora_range(SIG, DET, M, option = "bino_naive")




cleanEx()
nameEx("ora_all")
### * ora_all

flush(stderr()); flush(stdout())

### Name: ora_all
### Title: ORA using all metabolites in the pathway
### Aliases: ora_all
### Keywords: htest

### ** Examples

# Example1 : Metabolome data
data(MetaboAnalyst)

# Use pre-computed significant metabolites
SIG <- MetaboAnalyst$fasting$data$SIG

# metabolite set list
M <- MetaboAnalyst$fasting$pathway

# MSEA by over representation analysis
B <- ora_all(SIG, M)
B$`Result of MSEA(ORA)`

## Example2 : Proteome data (commented out - data not available)
# data(covid19)
# data(pathway)
# 
# X <- covid19$X$proteomics
# Y <- covid19$Y
# D <- covid19$D
# tau <- covid19$tau
# 
# protein_name <- colnames(X)
# 
# # pls-rog and pls-rog loading
# plsrog <- pls_rog(X,Y,D)
# plsrog <- plsrog_loading(plsrog)
# 
# # statistically significant proteins
# index_prot <- which(plsrog$loading$R[,1]>0 & plsrog$loading$p.value[,1]<0.05)
# sig_prot <- protein_name[index_prot]
# 
# # protein set list
# M <- pathway$covid19$proteomics
# 
# # MSEA by over representation analysis
# B <- ora_all(sig_prot, M)
# B$`Result of MSEA(ORA)`

# Example3: Metabolome data
data(MetaboAnalyst)

SIG <- MetaboAnalyst$fasting$data$SIG
M <- MetaboAnalyst$fasting$pathway

# Perform ORA using detected metabolites only
B <- ora_all(SIG, M)
B$`Result of MSEA(ORA)`




cleanEx()
nameEx("ora_bino")
### * ora_bino

flush(stderr()); flush(stdout())

### Name: ora_bino
### Title: Over-representation analysis with binomial resampling adjustment
### Aliases: ora_bino
### Keywords: MSEA ORA simulation

### ** Examples

M   <- list(set1 = c("A","B","C"), set2 = c("C","D","E"))
DET <- c("A","B","C","D")
SIG <- c("A","C")

# Binomial-adjusted ORA (shrink method)
ora_bino(SIG, DET, M, method = "shrink", nsim = 200)




cleanEx()
nameEx("ora_det")
### * ora_det

flush(stderr()); flush(stdout())

### Name: ora_det
### Title: Over-representation analysis using detected metabolites only
### Aliases: ora_det
### Keywords: msea

### ** Examples

# Example: Simple ORA with metabolomics data
data(MetaboAnalyst)

# Use pre-computed significant and detected metabolites
SIG <- MetaboAnalyst$fasting$data$SIG
DET <- MetaboAnalyst$fasting$data$DET
M <- MetaboAnalyst$fasting$pathway

# Perform ORA using detected metabolites only
B <- ora_det(SIG, DET, M)
B$`Result of MSEA(ORA)`




cleanEx()
nameEx("ora_est")
### * ora_est

flush(stderr()); flush(stdout())

### Name: ora_est
### Title: Over-representation analysis adjusting for undetected
###   metabolites
### Aliases: ora_est
### Keywords: msea

### ** Examples

# Example using metabolomics data
data(MetaboAnalyst)

# Use pre-computed significant and detected metabolites
SIG <- MetaboAnalyst$fasting$data$SIG
DET <- MetaboAnalyst$fasting$data$DET
M <- MetaboAnalyst$fasting$pathway

# Run adjusted ORA (e.g., shrinkage method)
B <- ora_est(SIG, DET, M, method = "shrink", lambda = 5)
B$`Result of MSEA (ORA with adjustment)`




cleanEx()
nameEx("ora_full")
### * ora_full

flush(stderr()); flush(stdout())

### Name: ora_full
### Title: Over-representation analysis with full enumeration of undetected
###   metabolite patterns
### Aliases: ora_full
### Keywords: MSEA ORA enrichment

### ** Examples

## Example: small set
M <- list(
  set1 = c("A", "B", "C"),
  set2 = c("C", "D", "E", "F"),
  set3 = c("G", "H")
)
DET <- c("A", "B", "C", "D", "E")  # Detected
SIG <- c("A", "C")                # Significant

ora_full(SIG, DET, M)$`Range of p-values`




cleanEx()
nameEx("setlabel")
### * setlabel

flush(stderr()); flush(stdout())

### Name: setlabel
### Title: Binary labeling of metabolite IDs across metabolite sets
### Aliases: setlabel
### Keywords: MSEA ORA labeling

### ** Examples

# Example usage
M <- list(
  set1 = c("A", "B", "C"),
  set2 = c("C", "D"),
  set3 = c("E", "F")
)

M_ID <- c("A", "C", "E", "G")

setlabel(M_ID, M)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
