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

ncores = 19

data = readRDS('./Data/MPOA.rds')

out_appdix = 'MPOA'
output_file = paste0('./',output_path,'/', out_appdix, '.rds')
output_pdf = paste0('./',output_path,'/', out_appdix,".pdf")

corpus = as.matrix(data$sim)
dim(corpus)
corpus[1:4,1:4]

pdf(output_pdf)

Ks = seq(2,20,1)
ldas <- fitLDA(counts = corpus,
               Ks = Ks,
               perc.rare.thresh = 0.01,
               seed = 0,
               ncores = ncores,
               plot = TRUE)

dev.off()

wldas = parallel::mclapply(Ks, function(k){
  ldamodel = optimalModel(ldas, k)
  wlda_bdc = WLDA(corpus, k, type = 'bdc', ldamodel = ldamodel)
  return(wlda_bdc)
}, mc.cores = ncores)

names(wldas) = paste0('k=',Ks)

results = list(ldas = ldas, wldas = wldas)
saveRDS(results, output_file)

beta = data[[1]]$gtCtGenes
theta = data[[1]]$gtSpotTopics

lda_beta = exp(models_A2$ldas$models$`9`@beta)
colnames(lda_beta) = models$ldas$models$`9`@terms
rownames(lda_beta) = paste0('Topic_', seq(nrow(lda_beta)))

lda_theta = ldas$models$`9`@gamma
colnames(lda_theta) = paste0('Topic_', seq(ncol(lda_theta)))
rownames(lda_theta) = ldas$models$`9`@documents
#apply(lda_theta,2,mean)
#lda_theta = filterTheta(lda_theta, 0.05)

corPlot_matchall(lda_beta, beta, lda_theta, theta)
corPlot_matchall(wldas$`k=9`$phi, beta, wldas$`k=9`$theta, theta)


cell_type_colors <- c(
  "Inhibitory" = "red",    # Red
  "Excitatory" = "#3ab9f9",    # Blue
  "Endothelial" = "yellow",   # Yellow
  "Astrocyte" = "#a07623e9",     # Brown
  "OD Mature" = "#4DAF4A",     # Green
  "OD Immature" = "#a33883",   # Purple
  "Ependymal" = "blue",     # Blue (same as Excitatory, assuming this is correct)
  "Microglia" = "#f7399e",     # Pink
  "Pericytes" = "orange"      # Orange
)

cell_type_colors = cell_type_colors[rownames(beta)]
decon_colors = cell_type_colors

beta_match = function(beta, beta_true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                       type = "b") # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  pairs <- lsatPairs(t(corMtx))
  
  beta = beta[pairs$colsix,]
  #rownames(beta) = rownames(beta_true)
  return(beta)
}

beta_bdc = beta_match(wldas$`k=9`$phi, beta)
names(decon_colors) = rownames(beta_bdc)

data_list = readRDS('Data/MPOA_simulation_data.rds')
position_list = lapply(data_list, function(data){
  pos = data$cellCounts[,1:2]
})

pdf('./MPOA/wLDA_bdc_proportion_plot.pdf', width = 8)
for(i in seq(length(data_list))){
  spots = rownames(position_list[[i]])
  theta_bregma = wldas$`k=9`$phi[spots, ]
  p = vizAllTopics(theta_bregma, position_list[[i]],
                   topicCols = decon_colors, r = 42) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank()
    )
  print(p)
}
dev.off()

