setwd('~/Desktop/2022-2023/wLDA/')

library(Matrix)
library(Seurat)
library(STdeconvolve)
library(patchwork)
library(dplyr)
library(hdf5r)
library(ggplot2)
library(parallel)
source('code/WLDA.R')

output_path =  'WLDA_simulation'
# create output folder if doesn't exist
if(!file.exists(output_path)){
  dir.create(output_path)
}

out_appdix = 'mOB'
output_file = paste0('./',output_path,'/', out_appdix,'.rds')
output_pdf = paste0('./',output_path,'/', out_appdix,".pdf")

ncores = 19

data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot
## remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 100)# supplementary
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
corpus[1:5,1:5]
dim(corpus)

pdf(output_pdf)

Ks = seq(2,20,1)

ldas <- fitLDA(as.matrix(corpus),
               Ks = Ks,
               ncores = ncores, # number of cores to fit LDA models in parallel
               plot=TRUE, verbose=FALSE)

dev.off()

wldas = parallel::mclapply(Ks, function(k){
  ldamodel = optimalModel(ldas, k)
  wlda_bdc = WLDA(corpus, k, type = 'bdc', ldamodel = ldamodel)
  return(wlda_bdc)
}, mc.cores = ncores)
names(wldas) = paste0('k=',Ks)
plt = PerplexityPlot(wldas, corpus = corpus)

results = list(ldas = ldas, wldas = wldas)
saveRDS(results, output_file)

optLDA <- optimalModel(models = model$lda, opt = 12)
results <- getBetaTheta(optLDA,
                        perc.filt = 0,
                        betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

STdeconvolve::vizAllTopics(deconProp, pos, 
                           groups = annot, 
                           group_cols = rainbow(length(levels(annot))),
                           r=0.4)

mobProxyTheta <- model.matrix(~ 0 + annot)
rownames(mobProxyTheta) <- names(annot)
# fix names
colnames(mobProxyTheta) <- unlist(lapply(colnames(mobProxyTheta), function(x) {
  unlist(strsplit(x, "annot"))[2]
}))

mobProxyGexp <- counts %*% mobProxyTheta

rownames(deconGexp) = colnames(deconProp) = paste0('Topic_', seq(ncol(deconProp)))

corPlot_matchall(deconGexp, t(mobProxyGexp), deconProp, mobProxyTheta)
corPlot_matchall(wldas$`k=7`$phi, t(mobProxyGexp), wldas$`k=7`$theta, mobProxyTheta)

topicCols = c('#d42728', '#F9BD3F', '#2c9f2c', '#1e77b4', '#6D1A9C', "#f4f1de", "#f4a99f")
names(topicCols) = paste0('Topic_', c(5,2,1,4,7,3,6))
order = c(5,2,1,4,7,3,6)

STdeconvolve::vizAllTopics(wldas$`k=7`$theta[, order], pos, 
                           topicCols = topicCols,
                           groups = annot, 
                           group_cols = rainbow(length(levels(annot))),
                           r=0.4)

pdf('mOB_simulation/mOB_bdc_filter.pdf')
theta = filterTheta(wldas$`k=7`$theta, 0.0)
for (i in seq(ncol((wldas$`k=7`$theta)))) {
  m <- as.data.frame(theta[,i])
  p <- pos
  other <- 1 - rowSums(m)
  m <- cbind(m, other)
  colnames(m) <- c(colnames(model$wlda$`k=12`$theta)[i], "Other")
  p = STdeconvolve::vizAllTopics(theta = m,
                                 pos = p,
                                 topicOrder=seq(ncol(m)),
                                 topicCols=c("black", "white"), # colors for cell-type 7, and "other"
                                 groups = annot, 
                                 group_cols = rainbow(length(levels(annot))), # make scatterpie borders white to only show the cell-type proportions.
                                 r = 0.4,
                                 lwd = 0.1,
                                 showLegend = TRUE,
                                 # BONUS: plot the scatterpies on top of a raster image of the H&E tissue, if this argument is equal to the rgb matrix representing the image
                                 overlay = NA) 
  print(p)
}

dev.off()

expr = t(corpus)
spatial_location = pos[rownames(corpus),]
vizGenesCounts(expr, spatial_location, genes = 'Beta-s',
               size = 7, stroke = 0.1,
               alpha = 1,
               plotTitle = 'Beta-s',
               showLegend = TRUE)





