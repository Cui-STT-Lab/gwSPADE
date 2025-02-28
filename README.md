# gwSPADE: Gene Frequency-weighted Reference-free Deconvolution in Spatial Transcriptomics

<img width="1467" alt="figure_overview" src="https://github.com/user-attachments/assets/f5e33750-773a-46c5-8fa0-2f3a51e900dd" />

gwSPADE is a gene frequency-weighted reference-free SPAtial DEconvolution method for spatial transcriptomics data. gwSPADE requires only the gene count matrix and utilizes appropriate weighting schemes within a topic model to accurately recover cell-type transcriptional profiles and their proportions at each spatial location, without relying on external single-cell reference information.

## Installation

You can install the development version of gwSPADE from [GitHub](https://github.com/Cui-STT-Lab/gwSPADE) with:

```r
# install.packages("devtools")
devtools::install_github("Cui-STT-Lab/gwSPADE")
