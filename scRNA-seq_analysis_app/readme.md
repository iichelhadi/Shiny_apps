# scRNA-seq Explorer

A comprehensive Shiny application for analyzing single-cell RNA sequencing data using the Seurat framework. This interactive tool provides a user-friendly interface for performing common scRNA-seq analysis tasks without extensive coding knowledge.

## Features

- **Data Import**: Load 10X Genomics data formats
- **Quality Control**: Interactive filtering based on metrics like mitochondrial percentage, number of genes, and UMI counts
- **Doublet Detection**: Identify and remove cell doublets using scDblFinder
- **Dimensionality Reduction**: PCA and UMAP visualizations
- **Clustering**: Flexible clustering options with multiple algorithms and resolution settings
- **Gene Expression Visualization**: Feature plots and violin plots for gene expression patterns
- **Marker Gene Identification**: Find differentially expressed genes between clusters
- **Cell Type Annotation**: Automated annotation using SingleR or manual annotation

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/scRNA-seq-explorer.git
cd scRNA-seq-explorer
```

2. Install the required R packages:

```r
# Install CRAN packages
install.packages(c("shiny", "shinythemes", "shinyFiles", "DT", "tidyverse", 
                  "ggplot2", "patchwork", "cowplot", "parallelly", "reticulate"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(c("scDblFinder", "BiocParallel", "MAST"))

# Install Seurat and scCustomize
install.packages("Seurat")
devtools::install_github("samuel-marsh/scCustomize")
```

3. For Leiden clustering, install the leidenalg Python package:
```bash
pip install leidenalg
```

4. Run the application:
```r
shiny::runApp("path/to/scRNA-seq-explorer")
```

## Usage

### Data Input
1. Select a 10X Genomics output directory containing the matrix, barcodes, and features files
2. Set QC parameters to filter cells based on quality metrics
3. Optionally detect and filter out doublets

### Preprocessing
1. Set parameters for variable feature selection and dimensionality reduction
2. Run the dimension reduction to visualize cells in UMAP space
3. Run clustering with desired resolution and algorithm

### Visualization
1. Select genes of interest to visualize their expression on UMAP or as violin plots
2. Explore the relationship between gene expression and cell clusters

### Marker Genes
1. Find marker genes for all clusters or between specific clusters
2. Visualize top marker genes in a heatmap
3. Download marker gene lists for further analysis

### Cell Annotation
1. Run automated cell type annotation using reference datasets
2. Visualize cell types on UMAP plot
3. Export annotated Seurat object for further analysis

## Example Data

Example 10X Genomics datasets can be downloaded from:
- [10X Genomics Datasets](https://www.10xgenomics.com/resources/datasets)
- [Satija Lab](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

## File Structure

```
scRNA-seq-explorer/
├── app.R             # Main application file
├── global.R          # Global variables and functions
├── www/              # Static assets
│   └── styles.css    # Custom CSS styles
├── data/             # Example datasets (optional)
└── README.md         # Documentation
```

## Customization

The app can be customized by modifying the following files:
- `app.R`: UI layout and server logic
- `global.R`: Functions for data processing and analysis
- `www/styles.css`: Custom styling for the application

## Performance Considerations

- The app performs memory-intensive operations. We recommend at least 16GB of RAM for analyzing datasets with more than 5,000 cells.
- For larger datasets (>20,000 cells), consider subsampling or using a server with more RAM.
- Some operations (dimensionality reduction, clustering, marker gene identification) may take several minutes to complete depending on dataset size.

## Troubleshooting

Common issues:
- **Memory errors**: Increase the memory allocation for R and/or use a smaller dataset
- **Missing leidenalg**: Ensure Python is properly configured with reticulate and leidenalg is installed
- **Slow performance**: Adjust parameters (reduce number of variable features, use fewer PCs for analysis)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this application in your research, please cite this repo:
```
https://github.com/iichelhadi/Shiny_apps/tree/main/scRNA-seq_analysis_app
```
And the core software packages:
- Seurat: Stuart et al. (2019). Comprehensive Integration of Single-Cell Data. Cell, 177(7), 1888-1902.
- scDblFinder: Germain et al. (2021). scDblFinder: a pipeline for filtering cells with artifactual doublet-like expression profiles. Bioinformatics, 37(19), 3333-3335.
- MAST: Finak et al. (2015). MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology, 16, 278.
- SingleR: Aran et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology, 20(2), 163-172.

Contact

For any questions or suggestions, please reach out:

    Name: Elhadi Iich
    Email: iichelhadi@gmail.com

