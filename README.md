<div align="center">
<h1>gRaphiaExtra: An R Package for Integrating ‘Seurat’ Objects into ‘gRaphia’</h1>

**Nilabhra R. Das<sup>1</sup>, 
Kaitlyn A. Flynn<sup>1</sup>,
John P. Kemp<sup>1</sup>**

<sup>1</sup> 
Mater Research Institute, The University of Queensland, Brisbane, Australia
</div>

### Introduction
Utilising graph-based network analysis frameworks, 'Graphia' <https://graphia.app/> is a powerful open source visual analytics application developed to aid the interpretation of large and complex datasets. For more details, see article by Freeman et al. (2022) <doi:10.1371/journal.pcbi.1010310>. 'gRaphia' is an extension of the 'Graphia' application  within the R environment, providing tools for network analysis and visualisation. 'gRaphiaExtra' provides additional functionality specifically designed for single-cell RNA-sequencing data, enabling users to seamlessly integrate and utilise existing 'Seurat' analysis outputs in 'gRaphia'. The package also provides supplementary functions to support and enhance the 'gRaphia' analysis framework.

Current CRAN version: 0.26.9

Current GitHub version: 0.26.8

### Installation
To install the package from CRAN:
```
# install gRaphiaExtra package
install.packages("gRaphiaExtra")
```

To install the package from GitHub:
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
- **testSeuratClusterEnrichment():** Test enrichment for disease genes (or, other gene sets of interest) within clusters in the marker gene list created using the FindAllMarkers() from the Seurat R package.
