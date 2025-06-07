setwd('~/Desktop/2022-2023/wLDA/')
library(ggplot2)
library(RcppEigen)
library(Rcpp)
# library(keyATM)
library(quanteda)
library(magrittr)
library(gtools)
source('code/functions.R')
source('code/WLDA.R')

set.seed(2012)

n_topics <- k <- 4
n_docs <- 1000
n_vocab <- V <- 100

b1 = gtools::rdirichlet(1, c(rep(30,3),rep(1, 97)))
#plot(1:100, b1)
b2 = gtools::rdirichlet(1, c(rep(2,5),rep(20,15),rep(2,80)))
b3 = gtools::rdirichlet(1, c(rep(2,5),rep(15,10),rep(20,10),rep(2,75)))
b4 = gtools::rdirichlet(1, c(rep(5,10),rep(20,20),rep(5,70)))

beta = rbind(b1,b2,b3,b4)
beta[,1:5]
beta.true = beta
beta.true[,1:5]
rowSums(beta.true)
rownames(beta.true) = paste0('Cell', seq(nrow(beta.true)))
colnames(beta.true) = paste0('gene', seq(ncol(beta.true)))

df_long <- reshape2::melt(beta.true)
colnames(df_long) <- c("cell", "gene", "expr")

ResultView = function(mat, type = 'b', xLabs = '', yLabs = '', title = ''){
  if(type == 'b'){
    mat = t(apply(mat, 1, function(row){
      (row-min(row))/(max(row)-min(row))
    }))
  }
  
  df_long <- reshape2::melt(mat)
  colors = c("#f1faee","#FEB915","red")
  
  p = ggplot(df_long, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    #scale_fill_gradient(low = "blue", high = "red") +
    scale_fill_gradientn(colours = colors, limits = c(0, 1)) + 
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Hide x-axis text
      axis.ticks.x = element_blank(), # Hide x-axis ticks
      axis.text.y = element_text(angle = 0, hjust = 0, size = 12),
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 15)
    ) +
    #ggtitle('Gene Expression Profile') +
    guides(fill = ggplot2::guide_colorbar(title = "Level",
                                          title.position = "top",
                                          title.hjust = 0.5,
                                          ticks.colour = "black",
                                          frame.colour= "black",
                                          label.hjust = 0
    )) +
    labs(x = xLabs, y = yLabs, title = title)
  return(p)
}


p_sim_beta = ResultView(beta.true, type = 'b', xLabs = 'Genes', title = 'Gene Expression Profile')
p_sim_beta
#theta.true <- gtools::rdirichlet(260, rep(1, 6)) 
theta.true <- gtools::rdirichlet(n_docs, rep(1/n_topics, n_topics)) 
colnames(theta.true) = rownames(beta.true)
rownames(theta.true) = paste0('Doc', seq(nrow(theta.true)))

p_sim_theta = ResultView(t(theta.true), type = 't', xLabs = 'Spots', title = 'Cell Type Distribution')
p_sim_theta

combined_plot <- p_sim_beta + p_sim_theta + plot_layout(nrow = 2) +
  patchwork::plot_annotation('Profile of True Composition')

docs <- matrix(0, nrow = n_docs, ncol = n_vocab) 
ksai <- 1000 # average words per doc
for (i in 1:n_docs) {
  # draw topics for each word
  tops <- rmultinom(1, rpois(1, ksai), theta.true[i, ]) # draw words
  for (j in 1:n_topics) {
    docs[i, ] <- docs[i, ] + rmultinom(1, tops[j], beta.true[j, ]) 
  }
}
rownames(docs) = paste0('Doc', seq(nrow(docs)))
colnames(docs) = paste0('gene', seq(ncol(docs)))
docs[1:5,1:5]
corpus = docs

lda = STdeconvolve::fitLDA(corpus, Ks = c(4), ncores = 8,
                                plot = F,
                                verbose = F)
beta_lda = exp(lda$models$`4`@beta)
theta_lda = lda$models$`4`@gamma
colnames(beta_lda) = colnames(corpus)
rownames(beta_lda) = colnames(theta_lda) = paste0('Topic_',seq(nrow(beta_lda)))
rownames(theta_lda) = rownames(corpus)

types = c('infor_bdc', 'inverse_bdc', 'bdc', 'infor', 'inverse', 'mutual_point', 'tf_idf')
wldas = parallel::mclapply(types, function(type){
  wlda = WLDA(corpus, k=4, type = type, ldamodel = lda$models$`4`)
  return(wlda)
}, mc.cores = 9)
names(wldas) = paste0('wlda_', types)

Decon_betas = c(list(lda = beta_lda),
                lapply(wldas, function(wlda){beta = wlda$phi}))
Decon_thetas = c(list(lda = theta_lda),
                 lapply(wldas, function(wlda){theta = wlda$theta}))

beta = beta.true
theta = theta.true

for (i in seq(Decon_betas)) {
  corPlot_matchall(Decon_betas[[i]], beta.true, Decon_thetas[[i]], theta.true,
                   title = names(Decon_betas)[i])
}

