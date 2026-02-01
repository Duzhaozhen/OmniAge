# OmniAge

**A comprehensive Python framework for aging research and biological clock calculation**

**Authors:**
* **Zhaozhen Du** (Shanghai Institute of Nutrition and Health, CAS) - duzhaozhen2022@sinh.ac.cn
* **Andrew E. Teschendorff** (Shanghai Institute of Nutrition and Health, CAS) - andrew@sinh.ac.cn

---

## 📖 Overview

**OmniAge** offers high-level interfaces for calculating a vast suite of aging-related clocks and biomarkers. It is designed to integrate seamlessly with standard bioinformatics workflows (e.g., `scanpy`, `pandas`).

The package covers three main categories:
1.  **Epigenetic (DNAm) Clocks**: Including Chronological, Biological (e.g., GrimAge, DunedinPACE), Cellular (Mitotic), and Gestational age clocks.
2.  **Transcriptomic (RNA) Clocks**: Including sc-ImmuAging, Brain_CT_clock, and Pasta.
3.  **Surrogate Biomarkers**: Proxies for proteins (CRP, IL6) and lifestyle traits (Smoking, BMI).

> **✨ For a complete list of supported clocks and detailed usage examples, please refer to our [Tutorial Notebook](notebooks/tutorial.ipynb).**

## 🛠 Installation

You can install `omniage` via pip:

```bash
pip install omniage