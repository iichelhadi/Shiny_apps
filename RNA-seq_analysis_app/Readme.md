# RNA-seq Analysis Shiny App

This repository contains a Shiny application designed for analyzing RNA-seq data. The app provides an interactive interface for users to upload RNA-seq datasets, perform various analyses, and visualize the results directly in their browser.

## Overview

The RNA-seq Analysis Shiny App includes the following features:
- **Data Upload**: Upload RNA-seq count data and associated metadata.
- **Exploratory Data Analysis**: Visualize data distributions, PCA plots, and more.
- **Differential Expression Analysis**: Identify differentially expressed genes using methods like DESeq2.
- **Data Visualization**: Generate customizable plots such as heatmaps, volcano plots, and boxplots.

## Installation

### Requirements

Ensure that you have the following installed:
- R (version 4.4.1 or later)
- R packages: `shiny`, `DESeq2`, `ggplot2`, `dplyr`, `plotly`, and others as required.

You can install the required packages by running:
```r
install.packages(c("shiny", "ggplot2", "dplyr", "plotly"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Running the App
# Clone the repository:
```bash
git clone https://github.com/iichelhadi/Shiny_apps.git
cd Shiny_apps/RNA-seq_analysis_app
```

# Run the Shiny app: Open R or RStudio and run the following command:

```r
shiny::runApp("app.R")
```

Interact with the App: The app should open in your default web browser, where you can start uploading your RNA-seq data and performing analyses.

# Features

- Data Upload: Easily upload RNA-seq count data in CSV/TSV/semicolon/space separated format.
- Species: RNA-seq data can be from analyzed from Human, Mouse or Rat.
- Exploratory Data Analysis: Interactive plots for data exploration, including PCA and expression distributions.
- Differential Expression Analysis: Use DESeq2 to identify differentially expressed genes.
- Customizable Visualizations: Generate publication-quality plots that can be customized directly in the app.

### Contact

For any questions or suggestions, please reach out:
- Name: Elhadi Iich
- Email: e.iich@hotmail.nl
