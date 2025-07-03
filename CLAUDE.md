# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**mseapca** is an R package for Metabolite Set Enrichment Analysis (MSEA) and Over-Representation Analysis (ORA) in metabolomics data. The package provides multiple statistical methods for pathway enrichment analysis using metabolite sets.

## Common Commands

### Package Development
```bash
# Build the package
R CMD build .

# Check the package
R CMD check mseapca_*.tar.gz

# Install the package locally
R CMD INSTALL .

# Load package in R for testing
R -e "devtools::load_all()"

# Document the package (requires devtools)
R -e "devtools::document()"

# Run R interactively to test functions
R
```

### Testing Functions
```r
# Example usage after loading the package
library(mseapca)
data(MetaboAnalyst)
result <- msea_ora(loadings, lib, type="det")
```

## Code Architecture

### Core Components

1. **Main Entry Point**: `R/msea_ora.R`
   - Wrapper function that dispatches to different ORA methods based on `type` parameter
   - Supported types: "det", "all", "est_naive", "est_weighted", "est_shrink"

2. **ORA Implementations**:
   - `ora_det.R`: Fisher's exact test on detected metabolites only
   - `ora_all.R`: ORA using all metabolites in the library
   - `ora_est.R`: Estimation-based methods with three variants (naive, weighted, shrink)
   - `ora_bino.R`: Binomial test-based approach
   - `ora_full.R`: Complete ORA implementation with all features

3. **Utility Functions**:
   - `setlabel.R`: Helper for labeling metabolite sets
   - `msea_ora_range.R`: Range-based analysis for multiple loadings

### Key Design Patterns

- All ORA functions follow similar input/output patterns:
  - Input: metabolite loadings, metabolite library, optional parameters
  - Output: data frame with pathway names, p-values, adjusted p-values, and metabolite counts
  
- P-value adjustment using Benjamini-Hochberg method is standard across all methods

- The metabolite library format is consistent: named list where names are pathway identifiers and values are metabolite vectors

### Function Signatures

Main functions expect:
- `x`: Numeric vector of metabolite loadings or scores
- `lib`: List of metabolite sets (pathways)
- `minsize`: Minimum pathway size (default: 1)
- Additional method-specific parameters

## Important Implementation Notes

- The package depends on the `loadings` package - ensure it's available
- XML package is imported for potential KEGG/MetaboAnalyst data parsing
- Fisher's exact test is implemented using R's built-in `fisher.test()`
- All p-value adjustments use `p.adjust()` with method="BH"
- The MetaboAnalyst.RData contains example data for testing (fasting study)

## Development Considerations

When modifying the package:
1. Maintain consistency in function return formats (data.frame with standard columns)
2. Preserve the wrapper pattern in `msea_ora()` when adding new methods
3. Update documentation in corresponding .Rd files when changing function signatures
4. Test with the included MetaboAnalyst dataset to ensure compatibility