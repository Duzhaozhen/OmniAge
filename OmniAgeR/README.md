# OmniAgeR

Provides a comprehensive suite of tools for calculating and evaluating various aging-related clocks from DNA methylation and transcriptomic data.

## Installation

You can install the development version of OmniAgeR from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Duzhaozhen/OmniAgeR")
```

## 📖 Quick Start
```r
library(OmniAgeR)

# 1. Load example data (requires OmniAgeRData package installed)
# Using the Hannum example dataset as a demonstration
data_list <- loadOmniAgeRdata("omniager_hannum_example")
beta_matrix <- data_list[[1]]

# 2. Call the core function to calculate aging scores
# You can specify multiple clocks of interest in the clockNames argument
results <- epiMarker(
    betaM = beta_matrix, 
    clockNames = c("Horvath2013", "Hannum", "PhenoAge")
)

# 3. View the calculation results
head(results)
```

## 📖 Tutorials
For comprehensive details on function usage, parameter specifications, and benchmarking case studies, please refer to the package vignette:
* [OmniAgeR: User Guide and Tutorials](vignettes/OmniAgeR.Rmd) - Comprehensive guide for the R-based workflow.
