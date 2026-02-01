
# OmniAge 🧬 

**OmniAge** is a cross-platform computational suite designed for the robust estimation and analysis of biological aging biomarkers. It provides a unified framework supporting both **R** and **Python**, integrating a vast suite of aging clocks and biomarkers.

---

## 📂 Repository Structure (Monorepo)

This repository is organized as a monorepo to ensure version consistency across platforms:

* **`OmniAgeR/`**: An R package providing the core implementation of aging clocks and biomarkers.
* **`OmniAgePy/`**: A Python package (`omniage`) optimized for high-throughput aging omic biomarker prediction.
---

## 🚀 Installation

### 1. R Version
Install the development version directly from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("Duzhaozhen/OmniAge", subdir = "OmniAgeR")
```

### 2. Python Version
Install via pip with the subdirectory flag:

```r
pip install git+[https://github.com/Duzhaozhen/OmniAge.git#subdirectory=OmniAgePy](https://github.com/Duzhaozhen/OmniAge.git#subdirectory=OmniAgePy)
```
---

## 📖 Quick Start
### Python API
```r
import omniage
import pandas as pd
import numpy as np

# Initialize the clock (automatically handles model decompression)
clock = DNAmCTFClock()

print("1. Loading Lung Pre-cancerous Lesions Dataset...")
beta_df = pd.read_csv("../example/data/LungInv_beta.csv", index_col=0)
meta_df = pd.read_csv("../example/data/LungInv_pheno.csv", index_col=0)
print(meta_df['Group'].value_counts())

print("2. Calculating Mitotic Ages...")
# Batch calculate all mitotic clocks using the "mitotic" group alias.
mitotic_ages = omniage.cal_epimarker(beta_df, clocks="Mitotic",ages=meta_df['Age'],return_dict=False)
```

### R API
```r
library(OmniAgeR)

download_OmniAgeR_example("LungInv")
load_OmniAgeR_example("LungInv")
my_comparisons <- list(c("N\nN=21", "LCIS\nN=13"), c("LCIS\nN=13", "LCIS->LC\nN=22"))
EpiAge.o<-EpiAge(data.m = bmiq.m,clock_names = "mitotic",ages.v = df$Age)
```

---

## 📖 Tutorials & Documentation

We provide step-by-step guides to help you get started with OmniAge:

### Python (omniage)
* [Python Package Tutorial](tutorials/OmniAgePy_tutorial.pdf) - Comprehensive guide for the Python-based workflow

### R (OmniAgeR)
* [R Package Tutorial](tutorials/OmniAgeR_tutorial.pdf) - Comprehensive guide for the R-based workflow.

---


