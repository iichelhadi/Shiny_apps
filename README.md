Here's an updated **top-level `README.md`** for your `Shiny_apps` repository that introduces both apps (RNA-seq and scRNA-seq), while preserving and integrating the detailed information you provided for the RNA-seq app.

---

````markdown
# Shiny Apps for RNA and scRNA-seq Analysis

This repository contains two interactive Shiny applications developed for RNA-seq and single-cell RNA-seq (scRNA-seq) data analysis. Each app is built to help researchers explore, analyze, and visualize transcriptomic datasets through a user-friendly browser interface without requiring advanced coding skills.

## Contents

- [`RNA-seq_analysis_app/`](./RNA-seq_analysis_app) – Differential expression analysis and visualization of bulk RNA-seq data
- [`scRNA-seq_analysis_app/`](./scRNA-seq_analysis_app) – End-to-end single-cell RNA-seq analysis using Seurat and related tools

---

## 📦 RNA-seq Analysis Shiny App

An interactive Shiny application for exploring and analyzing RNA-seq count data, supporting human, mouse, and rat datasets.

### Key Features

- **Data Upload**: Accepts raw count data and metadata (CSV/TSV/space/semicolon-separated)
- **Species Selection**: Analyze datasets from human, mouse, or rat
- **Exploratory Data Analysis**: Visualize expression distributions, PCA plots, and sample clustering
- **Differential Expression**: Identify differentially expressed genes using DESeq2
- **Custom Visualizations**: Generate heatmaps, volcano plots, and boxplots

### Installation

#### Requirements

- R (v4.4.1 or later)
- Required packages:
  ```r
  install.packages(c("shiny", "ggplot2", "dplyr", "plotly"))
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("DESeq2")
````

### Running the App

```bash
git clone https://github.com/iichelhadi/Shiny_apps.git
cd Shiny_apps/RNA-seq_analysis_app
```

```r
shiny::runApp("app.R")
```

---

## 🔬 scRNA-seq Analysis Shiny App

An interactive Shiny app for analyzing single-cell RNA-seq datasets, using the Seurat framework, scDblFinder, SingleR, and other tools.

Key features include:

* 10X Genomics data import
* Quality control filtering
* Doublet detection with `scDblFinder`
* Dimensionality reduction (PCA, UMAP)
* Clustering and cell type annotation
* Marker gene identification and visualization

For full documentation, visit:
[`scRNA-seq_analysis_app/readme.md`](./scRNA-seq_analysis_app/readme.md)

---

## 💬 Contact

For questions, suggestions, or contributions:

**Elhadi Iich**
📧 [iichelhadi@gmail.com](mailto:iichelhadi@gmail.com)
🌐 [GitHub Profile](https://github.com/iichelhadi)

---

## 📄 License

This repository is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

```

---

This version gives a clean overview of both apps, with the **bulk RNA-seq app explained in detail**, and the **scRNA-seq app linked to its own full README** for clarity and maintainability. Let me know if you'd like to embed app screenshots or badges.
```
