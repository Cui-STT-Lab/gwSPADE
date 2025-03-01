# gwSPADE: Gene Frequency-weighted Reference-free Deconvolution in Spatial Transcriptomics

<img width="1467" alt="figure_overview" src="https://github.com/user-attachments/assets/f5e33750-773a-46c5-8fa0-2f3a51e900dd" />

gwSPADE is a gene frequency-weighted reference-free SPAtial DEconvolution method for spatial transcriptomics data. gwSPADE requires only the gene count matrix and utilizes appropriate weighting schemes within a topic model to accurately recover cell-type transcriptional profiles and their proportions at each spatial location, without relying on external single-cell reference information.

## Installation

You can install the development version of gwSPADE from [GitHub](https://github.com/Cui-STT-Lab/gwSPADE) with:

```r
# install.packages("devtools")
devtools::install_github("Cui-STT-Lab/gwSPADE")
```

#### **Run gwSPADE with Example Data**
## Run gwSPADE with Example Data

The function WLDA takes the spatial transcriptomics data matrix `corpus` (spots x genes) as inputs.

The `corpus` needs to be non-transformed counts (or round the data to the nearest integers).

```r
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot
## remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 100)
odGenes <- getOverdispersedGenes(as.matrix(counts),
                                 gam.k=5,
                                 alpha=0.05,
                                 plot=FALSE,
                                 use.unadjusted.pvals=FALSE,
                                 do.par=TRUE,
                                 max.adjusted.variance=1e3,
                                 min.adjusted.variance=1e-3,
                                 verbose=FALSE, details=TRUE)
genes <- odGenes$ods
length(genes)

mobCorpus <- preprocess(t(cd),
                         selected.genes = genes,
                         # can then proceed to filter this list, if desired
                         # min.reads = 1, 
                         min.lib.size = 1, # can still filter pixels
                         min.detected = 1, # can still filter to make sure the selected genes are present in at least 1 pixel
                         ODgenes = FALSE, # don't select the over dispersed genes
                         verbose = TRUE)

corpus = as.matrix(mobCorpus$corpus)

Ks = seq(2,20,1)
ldas <- fitLDA(as.matrix(corpus), Ks = Ks, plot=FALSE, verbose=FALSE)

ncores = 10
wldas = parallel::mclapply(Ks, function(k){
  ldamodel = optimalModel(ldas, k)
  wlda_bdc = WLDA(corpus, k, type = 'bdc', ldamodel = ldamodel)
  return(wlda_bdc)
}, mc.cores = ncores)
names(wldas) = paste0('k=',Ks)
plt = PerplexityPlot(wldas, corpus = corpus)
print(plt) #Find the optimal number of deconvolved cell type
```
![plt](https://github.com/user-attachments/assets/89e05d74-e529-4c70-a2b2-277edd46b308)


To make a scatterpie plot showing the predicted proportions for all cell types:
```r
topicCols = c('#d42728', '#F9BD3F', '#2c9f2c', '#1e77b4', '#6D1A9C', "#f4f1de", "#f4a99f")
names(topicCols) = paste0('Topic_', c(5,2,1,4,7,3,6))
order = c(5,2,1,4,7,3,6)

STdeconvolve::vizAllTopics(wldas$`k=7`$theta[, order], pos, 
                           topicCols = topicCols,
                           groups = annot, 
                           group_cols = rainbow(length(levels(annot))),
                           r=0.4)
```
![mOB](https://github.com/user-attachments/assets/7bd22da0-c96b-4fd5-b544-cc3e5b5e780c)


