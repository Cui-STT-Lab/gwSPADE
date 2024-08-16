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

out_appdix = '10X_Visium'
output_file = paste0('./',output_path,'/', out_appdix,'.rds')
output_pdf = paste0('./',output_path,'/', out_appdix,'.pdf')

ncores = 19

se = readRDS('./Data/se.rds')
cd = se@assays$Spatial@layers$counts
cd[1:5,1:5]

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
dim(counts)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.01, alpha = 0.01, nTopOD = 1000)
corpus[1:4,1:4]
corpus = t(corpus)
dim(corpus)

pdf(output_pdf)

Ks = seq(8,20,1)
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

pos = se@images$slice1$centroids@coords
rownames(pos) = se@images$slice1$centroids@cells
head(pos)
colnames(pos) <- c("y", "x")

flip = function(location){
  colnames(location) <- c("y", "x")
  location[,1] = max(location[,1]) - location[,1]
  return(as.data.frame(location))
}
spatial_location = flip(pos)

topicCols = c("#00FFFF", "#FF0000", "#FF006D", "#FF6D00", "#0092FF", "#B6FF00","#FFEB00",
              '#2c9f2c', "#FEB915", "#6D1A9C", "#0024FF", "#B600FF", "#FF00DB")

theta_bdc = wldas$`k=13`$theta
colnames(theta_bdc) = paste0('Topic', seq(ncol(theta_bdc)))

pdf('Plot/10x_proportion.pdf')
plt <- STdeconvolve::vizAllTopics(theta = theta_bdc,
                                  pos = spatial_location,
                                  topicCols = topicCols,
                                  r = 50,
                                  lwd = 0,
                                  showLegend = TRUE,
                                  plotTitle = NA)
print(plt)

for (i in seq(ncol(wldas$`k=13`$theta))) {
  p = VizTopic(wldas$`k=13`$theta, i, spatial_location, colors = topicCols[i], r = 50)
  print(p)
}

dev.off()

marker_list = markergene_list(wldas$`k=13`$theta*1000, 5)

markers = c("Hpca", "Plp1")
df <- merge(as.data.frame(flip(pos)), 
            as.data.frame(t(as.matrix(t(corpus)[markers,]))), 
            by = 0)

ps <- lapply(markers, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 1, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.00,
                showLegend = TRUE) +
    
    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +
    
    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})

ps[[1]]+ps[[2]]


