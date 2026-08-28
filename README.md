# gRaphiaExtra

'Graphia' (https://graphia.app/) is a powerful open source visual analytics application developed to aid the interpretation of large and complex datasets. For more details, see Freeman et al. (2022) (doi:10.1371/journal.pcbi.1010310). 'gRaphia' (unreleased) is an extension of the 'Graphia' application  within the R environment, providing tools for network analysis and visualisation. 'gRaphiaExtra' provides additional functionality specifically designed for single-cell RNA-sequencing data, enabling users to seamlessly integrate and utilise existing 'Seurat' analysis outputs in 'gRaphia'. The package also provides supplementary functions to support and enhance the 'gRaphia' analysis framework.

Current version: 0.26.8
Currently under review by CRAN. Expected release in September 2026.

### Installation
To install the package:
```
# check for remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# install gRaphiaExtra package
remotes::install_github("mskgrg/gRaphiaExtra")
```

### Functions
Functions included in this version of the package:
- testSeuratClusterEnrichment()
