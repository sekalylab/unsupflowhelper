# unsupflowhelper v1.0

Unsupflowhelper is a R suite of tools for unsupervised analysis and visualization of flow cytometry data. 

Unsupflowhelper is compatible with other flow cytometry libraries, like flowCore, flowWorkspace, MetaCyto and CytoExploreR.
These tools excel at compensating, setting up gating strategies, or plotting more traditional flow cytometry plots, but are not built for modern unsupervised workflows.

Unsupflowhelper will help you perform dimension reduction and clustering, and vizualise or annotate your results in a user-friendly fashion.

## Installation
You can install the package in R using 

```
library(devtools)
install_github("sekalylab/unsupflowhelper")
```

A few additional optional packages are also recommended for more advanced analysis.
```
devtools::install_github("cvarrichio/Matrix.utils")
devtools::install_github("cole-trapnell-lab/monocle3")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SummarizedExperiment")
```

## Vignettes
A Step-by-step guide to using unsuperflowhelper is available in vignettes.

[Basic vignette](https://sekalylab.github.io/unsupflowhelper/guides/unsupervised_flow_vignette)


For a history/changelog, please see the [NEWS](https://github.com/sekalylab/unsupflowhelpr/NEWS.md).

