require(keyATM)
require(quanteda)
require(magrittr)
require(ggplot2)
require(dplyr)
require(STdeconvolve)
require(gridExtra)
require(patchwork)

filterTheta <- function(theta, perc.filt=0.05, verbose=TRUE){
  ## remove rare cell-types in pixels
  theta[theta < perc.filt] <- 0
  ## re-adjust pixel proportions to 1
  theta <- theta/rowSums(theta)
  ## if NAs because all cts are 0 in a spot, replace with 0
  theta[is.na(theta)] <- 0
  ## drop any cts that are 0 for all pixels
  dropped_cts <- names(which(colSums(theta) == 0))
  if(length(dropped_cts) > 0){
    if(verbose){
      message("Cell-types with no proportions in pixels after filtering dropped:\n",
              dropped_cts, "\n")
    }
  }
  theta <- theta[,which(colSums(theta) > 0)]
  
  empty_pixels <- names(which(rowSums(theta) == 0))
  if(length(empty_pixels) > 0){
    if(verbose){
      message(length(empty_pixels), " pixels with no cell-types after filtering.", "\n")
    }
  }
  
  return(theta)
}


corPlot = function(compare_matrix, ground_matrix, type = 'b', annotation=TRUE){
  corMtx <- getCorrMtx(m1 = as.matrix(compare_matrix), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(ground_matrix), # the reference `beta` (celltypes x genes)
                       type = type) # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  #rownames(corMtx) <- paste0(substitute(compare_matrix),'_', seq(nrow(corMtx)))
  pairs <- lsatPairs(t(corMtx))
  m <- t(corMtx)[pairs$rowix, pairs$colsix]
  
  p = correlationPlot(mat = t(m), # transpose back
                  colLabs = 'Deconvolved cell-types', # aka x-axis, and rows of matrix
                  rowLabs = 'Ground truth cell-types', # aka y-axis, and columns of matrix
                  title = "Transcriptional correlation", annotation = annotation) +
    
    ## this function returns a `ggplot2` object, so can add additional aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
  return(p)
  
}


vizAllTopics_adjust <- function(theta, pos,
                                topicOrder=seq(ncol(theta)),
                                topicCols=rainbow(ncol(theta)),
                                groups = NA,
                                group_cols = NA,
                                r = max(0.4, max(pos)/nrow(pos)*4),
                                lwd = 0.01,
                                showLegend = TRUE,
                                plotTitle = NA,
                                overlay = NA) {
  
  ## check that theta and pos are either data.frames or matrices
  if( !is.matrix(theta) & !is.data.frame(theta) ){
    stop("`theta` must be a matrix or data.frame.")
  }
  if( !is.matrix(pos) & !is.data.frame(pos) ){
    stop("`pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }
  
  if( (any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] != 2) ){
    stop("`pos` must have exactly 2 columns named `x` and `y`.")
  }
  
  # pixel cell-type distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  colnames(theta_ordered) <- paste0("X", colnames(theta_ordered))
  
  # ensure that `theta` and `pos` pixel rownames maintain same order
  # after the merge so as to not mess up the order of `groups`
  # if provided
  # make sure only using the shared pixels
  pixels <- intersect(rownames(theta_ordered), rownames(pos))
  pixels <- rownames(theta_ordered)[which(rownames(theta_ordered) %in% pixels)]
  
  # add columns "x", "y" with document positions from `pos`
  theta_ordered_pos <- merge(data.frame(theta_ordered),
                             data.frame(pos), by=0)
  rownames(theta_ordered_pos) <- theta_ordered_pos[,"Row.names"]
  ## make sure pixels in the original order before the merge
  theta_ordered_pos <- theta_ordered_pos[pixels,]
  
  # first column after merge is "Row.names", last two are "x" and "y"
  # problem is that data frame will replace "-" and " " with "."
  topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]
  
  # color of piechart groups (lines of piechart):
  if (is.na(groups[1])) {
    groups <- rep("0", dim(theta_ordered_pos)[1])
    theta_ordered_pos$Pixel.Groups <- groups
  } else {
    theta_ordered_pos$Pixel.Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c("0" = "gray")
  }
  
  message("Plotting scatterpies for ", dim(theta_ordered_pos)[1], " pixels with ", length(topicColumns),
          " cell-types...this could take a while if the dataset is large.", "\n")
  a = which(apply(theta_ordered_pos[,topicColumns], 2, function(col) all(col == 0)))
  
  if(length(a) == 0){
    topicCols = topicCols
  }else{
    topicCols = topicCols[-a]
  }
  
  if (!is.na(overlay[1])){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], y = 0:dim(overlay)[1])) +
      ggplot2::coord_equal(xlim = c(0,dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
      ) +
      ggplot2::annotation_raster(overlay, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "CellTypes") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
      ) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "CellTypes") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  }
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
}


scale0_1 <- function(x) {
  xAdj <- (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
  return(xAdj)
}


corPlot_matchall = function(beta, beta_true, theta, theta_true, title=NULL,annotation=TRUE){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                       type = "b") # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  pairs <- lsatPairs(t(corMtx))
  m <- t(corMtx)[pairs$rowix, pairs$colsix]
  
  p1 = correlationPlot(mat = t(m), # transpose back
                       colLabs = 'Deconvolved cell-types', # aka x-axis, and rows of matrix
                       rowLabs = 'Ground truth cell-types', # aka y-axis, and columns of matrix
                  title = "Gene Expression Profile", annotation = annotation) +
    
    ## this function returns a `ggplot2` object, so can add additional aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
 
  theta_test = theta[,pairs$colsix]

  corMtx_theta <- getCorrMtx(m1 = as.matrix(theta_test), # the deconvolved cell-type `beta` (celltypes x genes)
                             m2 = as.matrix(theta_true), # the reference `beta` (celltypes x genes)
                             type = "t")

  p2 = correlationPlot(mat = corMtx_theta, # transpose back
                       colLabs = 'Deconvolved cell-types', # aka x-axis, and rows of matrix
                       rowLabs = 'Ground truth cell-types', # aka y-axis, and columns of matrix
                  title = "Cell Distribution", annotation = annotation) +

    ## this function returns a `ggplot2` object, so can add additional aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
  
  if(is.null(title)) title = 'Deconvolved Results Compare'
  #grid.arrange(p1, p2, ncol = 2)
  #combined_plot <- plot_grid(p1, p2, ncol = 2)
  combined_plot <- p1 + p2
  p = combined_plot + plot_layout(ncol = 2) + patchwork::plot_annotation(title)
  return(p)
}


Terms <- function(phi, n){
  names = colnames(phi)
  terms_by_topic = apply(phi, 1, function(celltype){
    names[order(celltype, decreasing = TRUE)[1:n]]
  } )
  return(terms_by_topic)
}


Perplexity = function(beta, theta, corpus){
  beta = beta[, colnames(corpus)]
  theta = theta[rownames(corpus), ]
  x = slam::as.simple_triplet_matrix(corpus)
  perplexity = exp(-sum(log(colSums(beta[,x$j]*t(theta)[,x$i]))*x$v)/sum(x$v))
  # inner = log(theta%*%beta)
  # perplexity_check = exp(-sum(corpus*inner)/sum(corpus))
  return(perplexity)
}


PerplexityPlot <- function(models, corpus=NULL, perc.rare.thresh = 0.05){
  
  if(!is.null(models$models)){
    fitted_models <- models$models
    Ks <- as.vector(unlist(sapply(fitted_models, slot, "k")))
    pScores <- models$perplexities
    
    out <- lapply(1:length(Ks), function(i) {
      apply(getBetaTheta(fitted_models[[i]], corpus = corpus, verbose = FALSE, perc.filt = 0)$theta, 2, mean)
    })
    # original function, no perc.filt = 0
    ## number of cell-types present at fewer than `perc.rare.thresh` on average across pixels
    numrare <- unlist(lapply(out, function(x) sum(x < perc.rare.thresh)))
    alphas = unlist(sapply(fitted_models, slot, "alpha"))
  } else {
    Ks <- sapply(models, function(wlda){k = wlda$number_of_topics})
    pScores <- sapply(models, function(wlda){
      if(is.null(wlda$perplexity)){
        pscore = Perplexity(wlda$phi, wlda$theta, corpus)
        }else{
          pscores = wlda$perplexity
      }
      return(pscores)
    })
    
    out <- lapply(models, function(wlda) {
      apply(wlda$theta, 2, mean)
    })
    ## number of cell-types present at fewer than `perc.rare.thresh` on average across pixels
    numrare <- unlist(lapply(out, function(x) sum(x < perc.rare.thresh)))
    
    alphas <- sapply(models, function(wlda) {
      wlda$priors$alpha[1]
    }) 
  }
  
  dat <- data.frame(K = as.double(Ks),
                    rareCts = numrare,
                    perplexity = pScores,
                    rareCtsAdj = scale0_1(numrare),
                    perplexAdj = scale0_1(pScores),
                    alphas = alphas)
  dat[["alpha < 1"]] <- ifelse(dat$alphas < 1, 'gray90', 'gray50')
  dat$alphaBool <- ifelse(dat$alphas < 1, 0, 1)
  
  prim_ax_labs <- seq(min(dat$rareCts), max(dat$rareCts))
  prim_ax_breaks <- scale0_1(prim_ax_labs)
  ## if number rareCts stays constant, then only one break. scale0_1(prim_ax_labs) would be NaN so change to 0
  if(length(prim_ax_labs) == 1){
    prim_ax_breaks <- 0
    ## also the rareCtsAdj <- scale0_1(rareCts) would be NaN, so set to 0, so at same position as the tick,
    ## and its label will still be set to the constant value of rareCts
    dat$rareCtsAdj <- 0
  }
  
  if(max(dat$rareCts) < 1){
    sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), (max(dat$perplexity)-min(dat$perplexity))/1)
  } else {
    sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), (max(dat$perplexity)-min(dat$perplexity))/max(dat$rareCts))
  }
  sec_ax_breaks <- scale0_1(sec_ax_labs)
  
  if(length(sec_ax_labs) == 1){
    sec_ax_breaks <- 0
    dat$perplexAdj <- 0
  }
  
  plt <- ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(y=rareCtsAdj, x=K), col="blue", lwd = 2) +
    ggplot2::geom_point(ggplot2::aes(y=perplexAdj, x=K), col="red", lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y=rareCtsAdj, x=K), col="blue", lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y=perplexAdj, x=K), col="red", lwd = 2) +
    ggplot2::geom_bar(ggplot2::aes(x = K, y = alphaBool), fill = dat$`alpha < 1`, stat = "identity", width = 1, alpha = 0.5) +
    ggplot2::scale_y_continuous(name=paste0("# cell-types with mean proportion < ", round(perc.rare.thresh*100, 2), "%"), breaks = prim_ax_breaks, labels = prim_ax_labs,
                                sec.axis= ggplot2::sec_axis(~ ., name="perplexity", breaks = sec_ax_breaks, labels = round(sec_ax_labs, 2))) +
    ggplot2::scale_x_continuous(breaks = min(dat$K):max(dat$K)) +
    ggplot2::labs(title = "Fitted model K's vs deconvolved cell-types and perplexity",
                  subtitle = "LDA models with \u03b1 > 1 shaded") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size=15, face=NULL),
      # plot.subtitle = ggplot2::element_text(size=13, face=NULL),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "black", size = 0.1),
      panel.ontop = TRUE,
      axis.title.y.left = ggplot2::element_text(color="blue", size = 13),
      axis.text.y.left = ggplot2::element_text(color="blue", size = 13),
      axis.title.y.right = ggplot2::element_text(color="red", size = 15, vjust = 1.5),
      axis.text.y.right = ggplot2::element_text(color="red", size = 13),
      axis.text.x = ggplot2::element_text(angle = 0, size = 13),
      axis.title.x = ggplot2::element_text(size=13)
    )
  return(plt)
}

vizGenesCounts <- function(counts, spatial_location, genes,
                          groups = NA,
                          group_cols = NA,
                          size = 7, stroke = 0.5,
                          alpha = 1,
                          plotTitle = NA,
                          showLegend = TRUE) {
  
  #counts = sweep(counts,2,colSums(counts),"/")
  if(length(genes) == 1){
    df <- merge(as.data.frame(spatial_location), 
                as.data.frame(as.matrix(counts[genes,])), 
                by = 0)
    colnames(df) = c(colnames(df)[1:3], genes)
    counts <- df[,genes]
  }else{
    df <- merge(as.data.frame(spatial_location), 
                as.data.frame(t(as.matrix(counts[genes,]))), 
                by = 0)
    counts = rowSums(df[,genes])
    }
  
  
  # color spots by group:
  if (is.na(groups[1])) {
    groups <- " "
  } else {
    groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c(" " = "white")
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df, ggplot2::aes(x=x, y=y, fill=counts, color = groups),
                        shape = 21,
                        stroke = stroke, size = size, 
                        alpha = alpha) +
    viridis::scale_fill_viridis(option = "A", direction = -1) +
    ggplot2::scale_color_manual(values = group_cols)
  
  p <- p +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()) + 
    
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    )) + ## remove the pixel "groups", which is the color aesthetic for the pixel borders
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

   
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
}


VizTopic = function(theta_match, celltype, spatial_location, colors = NULL,
                    groups = NA, group_cols = NA,
                    seed = 0, r = 0.4,lwd = 0.2, plotTitle = NA){
  m <- as.data.frame(theta_match[,celltype])
  other <- 1 - rowSums(m)
  m <- cbind(other, m)
  colnames(m) <- c("Other",colnames(theta_match)[celltype])
  
  colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
                     "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
                     "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
                     "#dcf0d0","#f4a99f","#c8bad8",
                     "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
                     "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
                     "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
                     "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785",
                     "#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080")
  
  if(is.null(colors)){
    set.seed(seed)
    colors = colorCandidate[sample(1:length(colorCandidate),length(celltype))]
  }
  
  topicCols = c('white',colors)
  vizTopic = STdeconvolve::vizAllTopics(theta = m,
                          pos = spatial_location,
                          topicOrder=seq(ncol(m)),
                          topicCols=topicCols, # colors for cell-type 7, and "other"
                          groups = groups,
                          group_cols = group_cols,
                          r = r,
                          lwd = lwd,
                          plotTitle = plotTitle,
                          showLegend = TRUE,
                          overlay = NA) 
  return(vizTopic)
}


CARD.viz.gene <- function(spatial_expression,spatial_location,gene.visualize,NumCols,
                          size = 7, stroke = 0.3,
                          alpha = 1){
  expression = as.data.frame(as.matrix(spatial_expression))
  expression = sweep(expression,2,colSums(expression),"/")
  location = as.data.frame(spatial_location)
  if(sum(colnames(expression)==rownames(location))!= nrow(location)){
    stop("The colnames of expression data does not match with the rownames of spatial location data")
  }
  gene.select = gene.visualize
  if(sum(toupper(gene.select) %in% toupper(rownames(spatial_expression))) != length(gene.select)){
    stop("There existing selected genes that are not in the expression data!")
  }
  Data = NULL
  for(i in 1:length(gene.select)){
    #### load spatial dataset
    gene = gene.select[i]
    ind = which(toupper(rownames(expression)) == toupper(gene))
    df = as.numeric(expression[ind,])
    names(df) = colnames(expression)
    df = (df - min(df)) / (max(df) - min(df))
    d = data.frame(value = df,
                   x=as.numeric(location$x),
                   y = as.numeric(location$y))
    d$gene = gene
    Data = rbind(Data,d)
  }
  Data$gene = factor(Data$gene,levels = gene.select)
  p = ggplot() + 
    geom_point(data = Data, aes(x=x, y=y, fill = value), shape = 21#15
               , position = position_dodge2(width = 0.05),
               stroke = stroke, size = size, 
               alpha = alpha)+
    viridis::scale_fill_viridis(option = "A", direction = -1) +
    #scale_fill_gradientn(colours = c('#f1faee',"#b4e0ea",'darkblue'), limits = c(0, 1)) + 
    scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
    coord_equal()+
    facet_wrap(~gene,ncol = NumCols)+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          legend.position="bottom",
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.title =element_blank(),
          legend.title=element_text(size = 15,face="bold"),
          legend.text=element_text(size = 14),
          strip.text = element_text(size = 18,face="bold"),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(1.0, 'cm'))
    #guides(color = guide_colourbar(title = "Expression"))

  return(p)
}


gene_list = function(deconGexp, threshold){
  gene_list = list()
  gene = colnames(deconGexp)
  names = rownames(deconGexp)
  for (i in 1:nrow(deconGexp)) {
    gene_list[[names[i]]] = names(sort(deconGexp[i,gene[deconGexp[i,]>=threshold]], decreasing = TRUE))
  }
  return(gene_list)
}

markergene_list = function(deconGexp, threshold){
  ProxyLayerMarkers <- list()
  names = rownames(deconGexp)
  
  for (i in seq(nrow(deconGexp))){
    celltype <- i
    ## log2FC relative to other cell-types
    ## highly expressed in cell-type of interest
    highgexp <- names(which(deconGexp[celltype,] >= threshold))
    ## high log2(fold-change) compared to other deconvolved cell-types and 
    if(nrow(deconGexp)>2){
      log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
    } else {
      log2fc <- sort(log2(deconGexp[celltype,highgexp]/deconGexp[-celltype,highgexp]), decreasing=TRUE)
    }
    ## for gene set of the ground truth cell-type, get the genes
    ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
    markers <- names(log2fc[log2fc > 1])
    ProxyLayerMarkers[[names[i]]] <- markers
  }
  return(ProxyLayerMarkers)
}

