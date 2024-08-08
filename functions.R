#' Spatial transcriptomics deconvolution with the base model of SMART
#'
#' @param stData a spot(rows)-by-gene(columns) spatial transcriptomics matrix.
#' @param markerGs a list with marker genes for each cell type.
#' @param noMarkerCts the number of cell types without marker genes.
#' @param outDir the path to the output directory.
#' @param seed the starting value.
#' @param iterations the number of iterations (default: 2000).
#' @param priors do not change the priors if you don't know what they are.
#'
#' @return a model object
#' @return a cell type proportion matrix
#' @return a cell type-specific gene expression matrix
#' @import keyATM
#' @import quanteda
#' @export
#'
#' @examples
#' \dontrun{
#' library(SMART)
#' data(MPOA)
#' res <- SMART_base(stData=MPOA$stMat, markerGs=MPOA$markers, noMarkerCts=1,
#'                   outDir='SMART_results', seed=1, iterations=2000, prior=NULL)
#' }
require(keyATM)
require(quanteda)
require(magrittr)
require(ggplot2)
# wLDA_base <- function(stData, number_of_topics, seed=1, weight = 'inv-freq',
#                       iterations=1500, priors=NULL){
#   
#   # read data
#   stData <- keyATM_read(texts = as.dfm(stData), keep_docnames = TRUE)
#   
#   if(is.null(priors)){
#     # options
#     priors = list()
#     my_options <- list(seed          = seed,
#                        iterations    = iterations,
#                        verbose       = FALSE,
#                        use_weights   = TRUE,
#                        weights_type  = weight)
#   }else{
#     my_options <- list(seed          = seed,
#                        iterations    = iterations,
#                        verbose       = FALSE,
#                        use_weights   = TRUE,
#                        weights_type  = weight,
#                        estimate_alpha = 0)
#   }
#   
#   
#   # modelling
#   out <- weightedLDA(docs             = stData,
#                      number_of_topics = number_of_topics,
#                      model            = "base",
#                      priors = priors,
#                      options          = my_options)
#   
#   
#   colnames(out$theta)=gsub('[0-9]{1,2}_(.*)','\\1',colnames(out$theta))
#   rownames(out$phi)=gsub('[0-9]{1,2}_(.*)','\\1',rownames(out$phi))
#   
#   return(list(model=out, ct_proportions = out$theta, ct_spec_gexp = out$phi))
# }


#' Build hash table of simulated spots for each bregma layer of a given
#' MERFISH experiment.
#' 
#' @description uses `hash` to create a hash table of each bregma for a given
#'     input MERFISH experiment. Goal is to generate simulated spots for each bregma
#'     by partitioning cells into spots based on their Centroid coordinates. The
#'     size of the simulated spots is chosen as well but default is 100um x 100um.
#'     Note that for a given experiment/dataset there are multiple bregma layers. 
#'     Important: "Cell_class" column is where the cell labels will be pulled
#'     from. These will be used to construct the: "cellTypeTable", where downstream
#'     in `buildBregmaCorpus`, will be used to make "gtSpotTopics" and "gtCtGenes".
#' 
#' @param cellCentroidsAndClass an "annot.table" data.frame from
#'     "mpoa_merfish_clean.RData" that has been filtered for cells (the rows)
#'     of a given merfish experiment and has the following columns:
#'     c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")
#' @param counts cell x gene count matrix for the individual cells in all the
#'     bregmas for which spots will be simulated
#' @param patch_size dimension in um to make the simulated spots for a MERFISH
#'     experiment. The Centroid coords are already in um. (default: 100)
#' 
#' @return a hash table where each key is a bregma ID
#'     and the returned object is a list that contains
#' \itemize{
#' \item bregmaFullDf: the "annot.table" data.frame for the specific bregma with
#'     a new column "patch_id" indicating which patch a cell is assigned to. Written
#'     as "xcoord_ycoord". Simulated spots on the edges are dropped so some cells
#'     are not assigned to a patch and their patch ID is ""
#' \item cellTypeTable: table of counts of different cell types in each simulated patch
#' \item patchTotalCells: vector of total cells in each simulated patch
#' \item cellTypeCount: vector of counts of unique cell types in each simulated patch
#' \item cellGexp: individual cell x gene counts matrix for the specific cells
#'     in the bregma
#' \item patchGexp: patch x gene count matrix; ie the simulated gene counts for
#'     each simulated patch
#' }
#' 
#' @noRd

simulateBregmaSpots <- function (cellCentroidsAndClass, counts, patch_size = 100) {
  
  # dictionary hash table
  h <- hash::hash()
  
  data <- cellCentroidsAndClass
  
  for (bregma in unique(data$Bregma)) {
    
    bregma_key <- as.character(bregma)
    print(bregma_key)
    
    selected_bregma <- data[which(data$Bregma == bregma),]
    
    # 1. Get patch edge coordinates:
    
    # Sequence of X-coord positions for left edge of each patch:
    x_edges <- seq(min(selected_bregma$Centroid_X), max(selected_bregma$Centroid_X), patch_size)
    # drop first and last to avoid any issues with the edges of the whole region
    inner_x_edges <- x_edges[2:(length(x_edges)-1)]
    # Sequence of Y-coord positions for bottom edge of each patch:
    y_edges <- seq(min(selected_bregma$Centroid_Y), max(selected_bregma$Centroid_Y), patch_size)
    inner_y_edges <- y_edges[2:(length(y_edges)-1)]
    
    selected_bregma$patch_id <- character(length(rownames(selected_bregma)))
    
    # 2. add patch IDs to bregma cells, for the patch they belong to:
    
    for (x in inner_x_edges) {
      for (y in inner_y_edges) {
        patch_id <- paste0(as.character(x), "_", as.character(y))
        patch <- selected_bregma[which( (selected_bregma$Centroid_X > x) &
                                          (selected_bregma$Centroid_X < x+patch_size) &
                                          (selected_bregma$Centroid_Y > y) &
                                          (selected_bregma$Centroid_Y < y+patch_size) ),]
        
        if (length(rownames(patch)) > 0) {
          selected_bregma[rownames(patch),]$patch_id <- patch_id
        }
      }
    }
    
    # 3. get table of counts of cell types for each patch in bregma
    patch_cell_types <- selected_bregma[which(selected_bregma$patch_id != ""),
                                        c("patch_id", "Cell_class")]
    cellTypeTable <- table(patch_cell_types[])
    
    # 4. total cell counts in each patch
    patchTotalCells <- rowSums(cellTypeTable)
    
    # 5. counts of unique cell types in each patch
    cellTypeCount <- c()
    for (i in seq_len(length(rownames(cellTypeTable)))) {
      patch_num_cell_types <- length(which(cellTypeTable[i,] != 0))
      cellTypeCount <- append(cellTypeCount, patch_num_cell_types)
    }
    
    # 6. collapse gene counts for cells in same spot to make simulation
    patches <- unique(selected_bregma$patch_id[which(!selected_bregma$patch_id == "")])
    patchGexp <- do.call(rbind, lapply(patches, function(patch){
      cells <- rownames(selected_bregma[which(selected_bregma$patch_id == patch),])
      mat <- as.matrix(counts[cells,])
      if (length(cells) == 1){
        patch_counts <- as.vector(mat)
      } else if (length(cells) > 1){
        patch_counts <- colSums(mat)
      } else if (length(cells) == 0){
        cat("WARNING:", bregma, "patch", patch, "had no cells in `counts` to use for simulated gene expression", "\n")
      }
      patch_counts
    }))
    rownames(patchGexp) <- patches
    
    # 7. gene count matrix for individual cells in the bregma
    # this also includes cells not assigned to patches
    bregma_cells <- rownames(selected_bregma)
    cellGexp <- counts[bregma_cells,]
    
    # 8. combine data objects and append to hash table
    h[[bregma_key]] <- list(bregmaFullDf = selected_bregma,
                            cellTypeTable = cellTypeTable,
                            patchTotalCells = patchTotalCells,
                            cellTypeCount = cellTypeCount,
                            cellGexp = cellGexp,
                            patchGexp = patchGexp)
  }
  return(h)
}


#' Generate simulated corpus as input into `topicmodels` and the matched ground
#' truth spot-topic proportion and topic-gene proportion matrices for a given
#' simulated bregma spot dataset (built with `simulateBregmaSpots`)
#' 
#' @description Important: "Cell_class" column in "bregmaFullDf" as input for `simulateBregmaSpots`
#'     is where the cell labels will be pulled from. These will be used to construct
#'     the: "cellTypeTable", where downstream in `buildBregmaCorpus`, will be used
#'     to make "gtSpotTopics" and "gtCtGenes".
#' 
#' @param hashTable output hash table of simulated bregma spots from
#'     `simulateBregmaSpots`
#' @param bregmaID ID to reference a given bregma in the hashTable. ex: "-0.04"
#' 
##' @return a list that contains
#' \itemize{
#' \item sim: slam::as.simple_triplet_matrix(corpus); required format for topicmodels::LDA input
#' \item gtSpotTopics: the ground truth of major cell class proportions in each simulated patch.
#'     Use as ground truth "theta" (document x topix proportions)
#' \item gtCtGenes: the normalized gene count proportions of each of the major cell classes.
#'     Use as ground truth "beta" (topics x word frequencies)
#' \item cellCounts: data.frame that has patch names, "x", and "y" centroid coordinates,
#'     and counts of total cells in each patch
#' \item classColors: vector of colors for each cell class
#' \item annotDf: data.frame of the individual cells of the bregma,
#'     their position, cell type label, and assigned patch
#' }
#' 
#' @noRd

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}


buildBregmaCorpus <- function (hashTable, bregmaID) {
  
  bregmaID <- as.character(bregmaID)
  
  # corpus in slam format
  sim <- hashTable[[bregmaID]][["patchGexp"]]
  sim <- sim[order(rownames(sim)),]
  # print(sim)
  # remove "Blanks" from data:
  sim <- sim[,!grepl("Blank", colnames(sim))]
  sim <- slam::as.simple_triplet_matrix(sim)
  
  # ground truth spot - cell type proportions
  gtSpotTopics <- hashTable[[bregmaID]][["cellTypeTable"]]/rowSums(hashTable[[bregmaID]][["cellTypeTable"]])
  gtSpotTopics <- as.data.frame.matrix(gtSpotTopics[order(rownames(sim)),])
  
  # reformat `gtDocTopic` proportions into data frame with spot coordinates
  # tmp_positions <- do.call(rbind, lapply(rownames(gtSpotTopics), function(x){
  #   coords <- strsplit(x, "_")[[1]]
  #   as.numeric(coords)
  # }))
  # colnames(tmp_positions) <- c("x", "y")
  # rownames(tmp_positions) <- rownames(gtSpotTopics)
  # tmp_proportions <- lapply(colnames(gtSpotTopics), function(i) {
  #   gtSpotTopics[,i]
  # })
  # names(tmp_proportions) <- colnames(gtSpotTopics)
  # gtSpotTopics <- merge(tmp_positions, as.data.frame(tmp_proportions), by="row.names")
  
  # the individual cell annotation df but only for cells assigned to spots
  df <- hashTable[[bregmaID]][["bregmaFullDf"]]
  df <- df[which(df$patch_id != ""),]
  
  # cells and their assigned class (only cells in simulated patches)
  cellTypes <- df[,c("Cell_class")]
  cells <- rownames(df)
  
  # ground truth cell type - gene expression frequencies 
  mat <- hashTable[[bregmaID]][["cellGexp"]][cells,]
  mm <- stats::model.matrix(~ 0 + factor(cellTypes))
  colnames(mm) <- levels(factor(cellTypes))
  gtCtGenes <- t(t(as.matrix(mat)) %*% mm)
  # remove "Blanks" from data:
  gtCtGenes <- gtCtGenes[,!grepl("Blank", colnames(gtCtGenes))]
  gtCtGenes <- gtCtGenes/rowSums(gtCtGenes)
  
  # number of total cells in each spot
  cell_counts <- hashTable[[bregmaID]]$patchTotalCells
  count_df <- do.call(rbind, lapply(names(cell_counts), function(x){
    coords <- strsplit(x, "_")[[1]]
    as.numeric(coords)
  }))
  colnames(count_df) <- c("x", "y")
  rownames(count_df) <- names(cell_counts)
  count_df <- as.data.frame(count_df)
  count_df$counts <- cell_counts
  # same order as simulated spots
  count_df <- count_df[order(rownames(sim)),]
  
  # colors for each cell class
  classColors <- gg_color_hue(length(unique(df$Cell_class)))
  names(classColors) <- names(unique(df$Cell_class))
  
  bregma <- list(sim = sim,
                 gtSpotTopics = gtSpotTopics,
                 gtCtGenes = gtCtGenes,
                 cellCounts = count_df,
                 classColors = classColors,
                 annotDf = df)
  
  return(bregma)
}

require(dplyr)
simulateSpots <- function (cellCentroidsAndClass, counts, patch_size = 1000) {

  data <- cellCentroidsAndClass

    # 1. Get patch edge coordinates:

    # Sequence of X-coord positions for left edge of each patch:
    x_edges <- seq(min(data$Centroid_X), max(data$Centroid_X), patch_size)
    # drop first and last to avoid any issues with the edges of the whole region
    inner_x_edges <- x_edges[2:length(x_edges)-1]
    # Sequence of Y-coord positions for bottom edge of each patch:
    y_edges <- seq(min(data$Centroid_Y), max(data$Centroid_Y), patch_size)
    inner_y_edges <- y_edges[2:length(y_edges)-1]

    data$patch_id <- character(length(rownames(data)))
    patch_data = expand.grid(xmin = x_edges, ymin = y_edges) %>%
      mutate(xmax = xmin + 1000, ymax = xmin + 1000)
                            
    # 2. add patch IDs to cells, for the patch they belong to:

    for (x in inner_x_edges) {
      for (y in inner_y_edges) {
        patch_id <- paste0(as.character(x), "_", as.character(y))
        patch <- data[which( (data$Centroid_X > x) &
                                          (data$Centroid_X < x+patch_size) &
                                          (data$Centroid_Y > y) &
                                          (data$Centroid_Y < y+patch_size) ),]

        if (length(rownames(patch)) > 0) {
          data[rownames(patch),]$patch_id <- patch_id
        }
      }
    }

    # 3. get table of counts of cell types for each patch
    patch_cell_types <- data[which(data$patch_id != ""),
                                        c("patch_id", "Cell_class")]
    cellTypeTable <- table(patch_cell_types[])
    cellTypeTable = as.data.frame.matrix(cellTypeTable)

    # 4. total cell counts in each patch
    patchTotalCells <- rowSums(cellTypeTable)
    # patchTotalCells_move = names(patchTotalCells[patchTotalCells < 5])
    # patchTotalCells_filter = names(patchTotalCells[patchTotalCells >= 5])
    # patchTotalCells = patchTotalCells_filter
    # cellTypeTable = cellTypeTable[patchTotalCells_filter,]
    # data[data$patch_id %in% patchTotalCells_move,]$patch_id
    
    # 5. counts of unique cell types in each patch
    cellTypeCount <- c()
    for (i in seq_len(length(rownames(cellTypeTable)))) {
      patch_num_cell_types <- length(which(cellTypeTable[i,] != 0))
      cellTypeCount <- append(cellTypeCount, patch_num_cell_types)
    }

    # 6. collapse gene counts for cells in same spot to make simulation
    patches <- unique(data$patch_id[which(!data$patch_id == "")])
    patchGexp <- do.call(rbind, lapply(patches, function(patch){
      cells <- rownames(data[which(data$patch_id == patch),])
      mat <- as.matrix(counts[cells,])
      if (length(cells) == 1){
        patch_counts <- as.vector(mat)
      } else if (length(cells) > 1){
        patch_counts <- colSums(mat)
      } else if (length(cells) == 0){
        cat("WARNING:", "patch", patch, "had no cells in `counts` to use for simulated gene expression", "\n")
      }
      patch_counts
    }))
    rownames(patchGexp) <- patches

    # 7. gene count matrix for individual cells in the tissue
    # this also includes cells not assigned to patches
    cells <- rownames(data)
    cellGexp <- counts[cells,]

    # 8. combine data objects and append to hash table
    h <- list(DataFullDf = data,
              cellTypeTable = cellTypeTable,
              patchTotalCells = patchTotalCells,
              cellTypeCount = cellTypeCount,
              cellGexp = cellGexp,
              patchGexp = patchGexp,
              patch_data = patch_data)

  return(h)
}
#' 
#' 
#' #' Generate simulated corpus as input into `topicmodels` and the matched ground
#' #' truth spot-topic proportion and topic-gene proportion matrices for a given
#' #' simulated bregma spot dataset (built with `simulateBregmaSpots`)
#' #' 
#' #' @description Important: "Cell_class" column in "bregmaFullDf" as input for `simulateBregmaSpots`
#' #'     is where the cell labels will be pulled from. These will be used to construct
#' #'     the: "cellTypeTable", where downstream in `buildBregmaCorpus`, will be used
#' #'     to make "gtSpotTopics" and "gtCtGenes".
#' #' 
#' #' @param hashTable output hash table of simulated bregma spots from
#' #'     `simulateBregmaSpots`
#' #' @param bregmaID ID to reference a given bregma in the hashTable. ex: "-0.04"
#' #' 
#' ##' @return a list that contains
#' #' \itemize{
#' #' \item sim: slam::as.simple_triplet_matrix(corpus); required format for topicmodels::LDA input
#' #' \item gtSpotTopics: the ground truth of major cell class proportions in each simulated patch.
#' #'     Use as ground truth "theta" (document x topix proportions)
#' #' \item gtCtGenes: the normalized gene count proportions of each of the major cell classes.
#' #'     Use as ground truth "beta" (topics x word frequencies)
#' #' \item cellCounts: data.frame that has patch names, "x", and "y" centroid coordinates,
#' #'     and counts of total cells in each patch
#' #' \item classColors: vector of colors for each cell class
#' #' \item annotDf: data.frame of the individual cells of the bregma,
#' #'     their position, cell type label, and assigned patch
#' #' }
#' #' 
#' #' @noRd
buildCorpus <- function (hashTable) {

  # corpus in slam format
  sim <- hashTable[["patchGexp"]]
  sim = sim[rownames(hashTable[["cellTypeTable"]]), ]
  sim <- sim[order(rownames(sim)),]
  # print(sim)
  # remove "Blanks" from data:
  sim <- sim[,!grepl("Blank", colnames(sim))]
  sim <- slam::as.simple_triplet_matrix(sim)

  # ground truth spot - cell type proportions
  cellTypeTable = hashTable[["cellTypeTable"]]
  gtSpotTopics <- hashTable[["cellTypeTable"]]/rowSums(hashTable[["cellTypeTable"]])
  # gtSpotTopics <- as.data.frame.matrix(gtSpotTopics[order(rownames(sim)),])

  # reformat `gtDocTopic` proportions into data frame with spot coordinates
  # tmp_positions <- do.call(rbind, lapply(rownames(gtSpotTopics), function(x){
  #   coords <- strsplit(x, "_")[[1]]
  #   as.numeric(coords)
  # }))
  # colnames(tmp_positions) <- c("x", "y")
  # rownames(tmp_positions) <- rownames(gtSpotTopics)
  # tmp_proportions <- lapply(colnames(gtSpotTopics), function(i) {
  #   gtSpotTopics[,i]
  # })
  # names(tmp_proportions) <- colnames(gtSpotTopics)
  # gtSpotTopics <- merge(tmp_positions, as.data.frame(tmp_proportions), by="row.names")

  # the individual cell annotation df but only for cells assigned to spots
  df <- hashTable[["DataFullDf"]]
  df <- df[which(df$patch_id != ""),]

  # cells and their assigned class (only cells in simulated patches)
  cellTypes <- df[,c("Cell_class")]
  cells <- rownames(df)

  # ground truth cell type - gene expression frequencies
  mat <- hashTable[["cellGexp"]][cells,]
  mm <- stats::model.matrix(~ 0 + factor(cellTypes))
  colnames(mm) <- levels(factor(cellTypes))
  gtCtGenes <- t(t(as.matrix(mat)) %*% mm)
  # remove "Blanks" from data:
  gtCtGenes <- gtCtGenes[,!grepl("Blank", colnames(gtCtGenes))]
  gtCtGenes <- gtCtGenes/rowSums(gtCtGenes)

  # number of total cells in each spot
  cell_counts <- hashTable$patchTotalCells
  count_df <- do.call(rbind, lapply(names(cell_counts), function(x){
    coords <- strsplit(x, "_")[[1]]
    as.numeric(coords)
  }))
  colnames(count_df) <- c("x", "y")
  rownames(count_df) <- names(cell_counts)
  count_df <- as.data.frame(count_df)
  count_df$counts <- cell_counts
  # same order as simulated spots
  count_df <- count_df[order(rownames(sim)),]

  # colors for each cell class
  classColors <- gg_color_hue(length(unique(df$Cell_class)))
  names(classColors) <- names(unique(df$Cell_class))

  Data <- list(sim = sim,
                 gtSpotTopics = gtSpotTopics,
                 gtCtGenes = gtCtGenes,
                 cellCounts = count_df,
                 classColors = classColors,
                 annotDf = df)

  return(Data)
}



#' Pull out cell-type proportions across pixels (theta) and
#' cell-type gene probabilities (beta) matrices from fitted LDA models from fitLDA
#'
#' @param lda an LDA model from "topicmodels" R package. From list of models returned by
#'     fitLDA
#' @param perc.filt proportion threshold to remove cell-types in pixels (default: 0.05)
#' @param betaScale factor to scale the predicted cell-type gene expression profiles (default: 1)
#' @param verbose Boolean for verbosity (default: TRUE)
#'
#' @return A list that contains
#' \itemize{
#' \item beta: cell-type (rows) by gene (columns) distribution matrix.
#'     Each row is a probability distribution of a cell-type expressing each gene
#'     in the corpus
#' \item theta: pixel (rows) by cell-types (columns) distribution matrix. Each row
#'     is the cell-type composition for a given pixel
#' }
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3, ncores=2)
#' optLDA <- optimalModel(models = ldas, opt = 3)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' head(results$theta)
#' head(results$beta)
#' 
# getBetaTheta_wlda <- function(wlda, perc.filt=0.05, betaScale=1, verbose=TRUE) {
#   result = wlda
#   
#   theta <- result$ct_proportions
#   beta <- result$ct_spec_gexp
#   
#   ## filter out cell-types with low proportions in pixels
#   if(verbose){
#     message("Filtering out cell-types in pixels that contribute less than ",
#             perc.filt, " of the pixel proportion.", "\n")
#   }
#   theta <- filterTheta(theta, perc.filt=perc.filt, verbose=verbose)
#   
#   ## scale the beta
#   beta <- beta * betaScale
#   
#   return(list(beta=beta,
#               theta=theta))
# }


#' Function to filter out cell-types in pixels below a certain proportion
#' 
#' @description Sets cell-types in each pixel to 0 that are below a given proportion.
#'     Then renormalizes the pixel proportions to sum to 1.
#'     Cell-types that result in 0 in all pixels after this filtering step are removed.
#'
#' @param theta pixel (rows) by cell-types (columns) distribution matrix. Each row
#'     is the cell-type composition for a given pixel
#' @param perc.filt proportion threshold to remove cell-types in pixels (default: 0.05)
#' @param verbose Boolean for verbosity (default: TRUE)
#' 
#' @return A filtered pixel (rows) by cell-types (columns) distribution matrix.
#' 
#' @noRd
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


# ## Estimate model parameters for data generation
# getparams <- function(x, pval_cutoff = 0.05,maxiter=100){
#   p <- nrow(x)
#   n <- ncol(x)
#   ## parameter estimation
#   ## parallel possible here, which will increase the efficiency
#   params <- t(apply(x, 1, function(gene){
#     m <- mean(gene)
#     v <- stats::var(gene)
#     mle_NB <- fit_nb_optim(gene,maxiter=maxiter)
#     c(mle_NB[1], mle_NB[2])
#   }))
#   names(params) = c('theta', 'mu')
#   # return(list(params = params))
#   return(params)
# }
# 
# 
# fit_nb_optim <- function(x,maxiter= 500){
#   m = mean(x)
#   v = stats::var(x)
#   
#   if(v>m){
#     size = m^2/(v-m)
#   }else{
#     size = 100
#   }
#   ## fix the mu through the moment matching to facilitate the estimation
#   fitres <- tryCatch({
#     opt_out <- optim(size,nb_loglik_mom,mu=m,y=x,control = list(maxit = maxiter), method="Brent",lower=0,upper=1000)
#     c(opt_out$par,m,-opt_out$value,1-opt_out$convergence)
#   },
#   error = function(cond){
#     # library(MASS)
#     glmfit <- glm.nb(x~1)
#     c(glmfit$theta, exp(glmfit$coefficients),as.numeric(logLik(glmfit)),min(glmfit$converged,is.null(glmfit$th.warn))) 
#   })
#   
#   
#   names(fitres) <- c("theta","mu","llk","convergence")
#   return(fitres)
# }
# 
# # Function to find vectors containing the character
# find_vectors_with_char <- function(char, list_of_vectors) {
#   names(which(sapply(list_of_vectors, function(vec) char %in% vec)))
# }
# 
# dominate_gene = function(gene_list){
#   unique_chars <- unique(unlist(gene_list))
#   
#   # Create a list to hold the results
#   result_list <- list()
#   
#   # Iterate over each character
#   for (char in unique_chars) {
#     vectors_with_char <- find_vectors_with_char(char, gene_list)
#     
#     # Check if the character appears in at least 2 vectors
#     if (length(vectors_with_char) >= 2) {
#       result_list[[char]] <- vectors_with_char
#     }
#   }
#   return(result_list)
# }

require(STdeconvolve)
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


RMSE_spot <- function(beta, beta_true, theta, theta_true){
  mtx = getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                   m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                   type = "b")
  
  # mtx must not more rows than columns
  # values in matrix converted to 0-1 scale relative to all values in mtx
  pairing <- clue::solve_LSAP(scale0_1(t(mtx)), maximum=TRUE)
  # clue::lsat returns vector where for each position the first element is a row
  # and the second is the paired column
  # rowsix <- seq_along(pairing)
  colsix <- as.numeric(pairing)
  
  theta_test = theta[,colsix]
  RMSE = sqrt(rowSums((theta_test - theta_true)^2)/ncol(theta_true))
  
  return(RMSE)
}

Spearman = function(beta, beta.true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta.true), # the reference `beta` (celltypes x genes)
                       type = 'b') # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  rownames(corMtx) <- paste0(substitute(beta),'_', seq(nrow(corMtx)))
  pairing <- clue::solve_LSAP(scale0_1(t(corMtx)), maximum=TRUE)
  colsix <- as.numeric(pairing)
  beta <- beta[colsix,]
  sharedGenes <- intersect(colnames(beta), colnames(beta.true))
  
  Spearman = sapply(seq(nrow(beta)), function(row){
    cor(beta[row,sharedGenes], beta.true[row,sharedGenes], method = 'spearman')
  })
  names(Spearman) = rownames(beta.true)
  return(Spearman)
}

PCC_beta = function(beta, beta.true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta.true), # the reference `beta` (celltypes x genes)
                       type = 'b') # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  rownames(corMtx) <- paste0(substitute(beta),'_', seq(nrow(corMtx)))
  pairing <- clue::solve_LSAP(scale0_1(t(corMtx)), maximum=TRUE)
  colsix <- as.numeric(pairing)
  m <- corMtx[colsix,]
  PCC = diag(m)
  names(PCC) = colnames(m)
  return(PCC)
}


PCC_theta = function(beta, beta_true, theta, theta_true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                       type = 'b') # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  rownames(corMtx) <- paste0(substitute(beta),'_', seq(nrow(corMtx)))
  pairing <- clue::solve_LSAP(scale0_1(t(corMtx)), maximum=TRUE)
  colsix <- as.numeric(pairing)
  #m <- corMtx[colsix,]
  theta_test = theta[,colsix]
  corMtx_theta <- getCorrMtx(m1 = as.matrix(theta_test), # the deconvolved cell-type `beta` (celltypes x genes)
                             m2 = as.matrix(theta_true), # the reference `beta` (celltypes x genes)
                             type = "t")
  
  PCC = diag(corMtx_theta)
  names(PCC) = colnames(corMtx_theta)
  return(PCC)
}

require(gridExtra)
require(patchwork)
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


JS_Divergence = function(beta, beta_true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                       type = "b") # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  pairs <- lsatPairs(t(corMtx))
  
  beta_match = beta[pairs$colsix,]
  
  JSD = lapply(1:nrow(beta), function(i){
    P = beta_match[i,]
    Q = beta[i,]
    M = 1/2*(P+Q)
    JSD = 1/2*sum(P*log(P/M)) + 1/2*sum(Q*log(Q/M))
    return(JSD)
  })
  JSD = sum(unlist(JSD))
  
  return(JSD)
}

KL_Divergence = function(beta, beta_true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                       type = "b") # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  pairs <- lsatPairs(t(corMtx))
  
  beta_match = beta[pairs$colsix,]
  
  KL = lapply(1:nrow(beta), function(i){
    P = beta_match[i,]
    Q = beta[i,]
    KL = sum(P*log(P/Q))
    return(KL)
  })
  KL = sum(unlist(KL))
  
  return(KL)
}


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

RankGenePlot = function(beta, beta_true, size=5, title = NULL){
  
  summary_df = RankGenes(beta, beta_true)
  if(is.null(title)) title = "Cell-type transcriptional profiles"
  
  p <- ggplot2::ggplot(data = summary_df) +
    ggplot2::geom_hex(ggplot2::aes(x = gexp_rank, y = beta_rank), bins=10) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::labs(title = title,
                  x = "Ground truth gene expression level",
                  y = "Deconvolved gene expression level") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=size, color = "black"),
                   axis.text.y = ggplot2::element_text(size=size, color = "black"),
                   axis.title.y = ggplot2::element_text(size=size),
                   axis.title.x = ggplot2::element_text(size=size),
                   plot.title = ggplot2::element_text(size=size),
                   legend.text = ggplot2::element_text(size = size, colour = "black"),
                   legend.title = ggplot2::element_text(size = size, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank()
                   # legend.position="none"
    ) +
    ggplot2::coord_equal() +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Count",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 0.5,
                                                   frame.colour= "black",
                                                   frame.linewidth = 0.5,
                                                   label.hjust = 0
    ))
  
  return(p)
}


RankGenes = function(beta, beta_true){
  corMtx <- getCorrMtx(m1 = as.matrix(beta), # the deconvolved cell-type `beta` (celltypes x genes)
                       m2 = as.matrix(beta_true), # the reference `beta` (celltypes x genes)
                       type = "b") # "b" = comparing beta matrices, "t" for thetas
  ## row and column names need to be characters
  pairs <- lsatPairs(t(corMtx))
  #beta = beta[pairs$colsix,]
  matches <- rownames(beta_true)
  names(matches) = rownames(beta)[pairs$colsix]
  sharedGenes <- intersect(colnames(beta), colnames(beta_true))
  
  summary_df <- do.call(rbind, lapply(names(matches), function(i){
    # get paired ground truth ct
    ct <- as.vector(matches[i])
    ## the deconvolved cell-type
    topic <- i
    topic_betas <- beta[topic, sharedGenes]
    ct_betas <- beta_true[ct, sharedGenes]
    ## values for a given topic=ground truth pair.
    ## append these together into final summary df
    df <- data.frame(gexp = ct_betas,
                     beta = topic_betas,
                     gexp_rank = rank(ct_betas),
                     beta_rank = rank(topic_betas),
                     pair = paste0(topic, " vs ", ct))
  }))
  return(summary_df)
}


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


Theta_domain = function(theta){
  binary_matrix = t(apply(theta, 1, function(row) {
    row == max(row) # Logical vector: TRUE where the element is the max in the row
  }))
  binary_matrix <- matrix(as.numeric(binary_matrix), nrow = nrow(theta))
  rownames(binary_matrix) = rownames(theta)
  colnames(binary_matrix) = colnames(theta)
  
  return(binary_matrix)
}



Terms <- function(phi, n){
  names = colnames(phi)
  terms_by_topic = apply(phi, 1, function(celltype){
    names[order(celltype, decreasing = TRUE)[1:n]]
  } )
  return(terms_by_topic)
}


Coherence_score = function(counts, phi, top_n_tokens = 10,
                           smoothing_beta = 1, type='PMI',
                           genelist = NULL, threshold=2/1000, verbose = TRUE){
  if(is.null(genelist)){
    #top_terms = Terms(phi, top_n_tokens)
    genelist = markergene_list(phi, threshold)
    }
  
  top_terms = do.call(cbind, lapply(genelist, function(genes){genes[seq(top_n_tokens)]}))
  
  dtm_data = slam::as.simple_triplet_matrix(counts)
  n = min(colSums(!is.na(top_terms)))
  if(verbose) cat('The number of top marker genes is', n, '\n')
  
  top_terms = top_terms[1:n,]
  a = apply(top_terms, 2, coherence, dtm_data = dtm_data, smoothing_beta = smoothing_beta,
            type = type)
  return(a)
}

coherence = function(dtm_data, top_terms, smoothing_beta, type){
  # Get the relevant entries of the document-term matrix
  rel_dtm <- dtm_data[,top_terms]
  # Turn it into a logical representing co-occurences
  df_dtm <- rel_dtm > 0
  D = nrow(dtm_data)
  # Calculate document frequencies for each term and all of its co-occurences
  cooc_mat <- slam::tcrossprod_simple_triplet_matrix(t(df_dtm))
  # Quickly get the number of top terms for the for-loop below
  top_n_tokens <- length(top_terms)
  
  if(type == 'UMass'){
    # Using the syntax from the paper, calculate coherence
    c_l <- 0
    for (m in 2:top_n_tokens) {
      for (l in 1:(m - 1)) {
        df_ml <- cooc_mat[m,l]
        df_l <- cooc_mat[l,l]
        c_l <- c_l + log((df_ml + smoothing_beta) / df_l)
      }
    }
  }else if(type == 'PMI'){
    c_l <- 0
    for (m in 2:top_n_tokens) {
      for (l in 1:(m-1)) {
        df_ml <- cooc_mat[m,l]
        df_l <- cooc_mat[l,l]
        df_m = cooc_mat[m,m]
        c_l <- c_l + log((df_ml+smoothing_beta) / (df_l*df_m/D))
      }
    }
  }
  return(c_l)
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

# CoherenceScore_Plot = function(models, corpus=NULL, threshold = 2/1000, perc.rare.thresh = 0.05,
#                                check='mean', type = 'PMI', verbose = TRUE){
#   if(!is.null(models$models)){
#     fitted_models <- models$models
#     Ks <- as.vector(unlist(sapply(fitted_models, slot, "k")))
#     pScores <- models$perplexities
#     out <- lapply(1:length(Ks), function(i) {
#       beta = getBetaTheta(fitted_models[[i]], verbose = verbose, perc.filt = 0)$beta
#       rownames(beta) = paste0('Topic_', seq(nrow(beta)))
#       return(beta)
#     })
# 
#     out_theta <- lapply(1:length(Ks), function(i) {
#       apply(getBetaTheta(fitted_models[[i]], verbose = verbose, perc.filt = 0)$theta, 2, mean)
#     })
#     numrare <- unlist(lapply(out_theta, function(x) sum(x < perc.rare.thresh)))
#     
#     # original function, no perc.filt = 0
#     alphas = unlist(sapply(fitted_models, slot, "alpha"))
#   } else {
#     Ks <- sapply(models, function(wlda){k = wlda$number_of_topics})
#     pScores <- sapply(models, function(wlda){
#       if(is.null(wlda$perplexity)){
#         pscore = Perplexity(wlda$phi, wlda$theta, corpus)
#       }else{
#         pscores = wlda$perplexity
#       }
#       return(pscores)
#     })
#     
#     out <- lapply(models, function(wlda) {
#       wlda$phi
#     })
#     
#     out_theta <- lapply(models, function(wlda) {
#       apply(wlda$theta, 2, mean)
#     })
#     numrare <- unlist(lapply(out_theta, function(x) sum(x < perc.rare.thresh)))
#     
#     alphas <- sapply(models, function(wlda) {
#       wlda$priors$alpha[1]
#     }) 
#   }
#   
#   out_marker = lapply(out, function(beta){
#     markergene_list(beta, threshold)
#   })
#   n_markers = min(unlist(lapply(out_marker, function(genelist){lengths(genelist)})))
#   cs = lapply(seq(out), function(i){
#     Coherence_score(corpus, out[[i]], top_n_tokens = n_markers, type = type, genelist = out_marker[[i]])
#   })
#   
#   if(check == 'median'){
#     cScores = sapply(cs, function(cs_score){median(cs_score)})
#   }else if(check == 'mean'){
#     cScores = sapply(cs, function(cs_score){mean(cs_score)})
#   }
#   
#   dat <- data.frame(K = as.double(Ks),
#                     rareCts = numrare,
#                     perplexity = pScores,
#                     cScores = cScores,
#                     rareCtsAdj = scale0_1(numrare),
#                     perplexAdj = scale0_1(pScores),
#                     cScoresAdj = scale0_1(cScores),
#                     alphas = alphas)
#   
#   dat[["alpha < 1"]] <- ifelse(dat$alphas < 1, 'gray90', 'gray50')
#   dat$alphaBool <- ifelse(dat$alphas < 1, 0, 1)
#   
#   prim_ax_labs <- seq(min(dat$rareCts), max(dat$rareCts))
#   prim_ax_breaks <- scale0_1(prim_ax_labs)
#   ## if number rareCts stays constant, then only one break. scale0_1(prim_ax_labs) would be NaN so change to 0
#   if(length(prim_ax_labs) == 1){
#     prim_ax_breaks <- 0
#     ## also the rareCtsAdj <- scale0_1(rareCts) would be NaN, so set to 0, so at same position as the tick,
#     ## and its label will still be set to the constant value of rareCts
#     dat$rareCtsAdj <- 0
#   }
#   
#   sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), length.out = max(dat$rareCts))
#   sec_ax_breaks <- scale0_1(sec_ax_labs)
#   sec_cs_ax_labs <- seq(min(dat$cScores), max(dat$cScores), length.out = max(dat$rareCts))
#   
#   # Create custom labels combining pScores and cScores with different colors
#   sec_ax_labels <- sapply(1:length(sec_ax_labs), function(i) {
#     paste0(
#       '<span style="color:red;">', round(sec_ax_labs[i], 2), '</span>',
#       ' (<span style="color:green;">', round(sec_cs_ax_labs[i], 3), '</span>)'
#     )
#   })
#   
#   # Create the plot
#   plt <- ggplot2::ggplot(dat) +
#     ggplot2::geom_point(ggplot2::aes(y = rareCtsAdj, x = K), col = "blue", size = 2) +
#     ggplot2::geom_point(ggplot2::aes(y = perplexAdj, x = K), col = "red", size = 2) +
#     ggplot2::geom_point(ggplot2::aes(y = cScoresAdj, x = K), col = "green", size = 2) +
#     ggplot2::geom_line(ggplot2::aes(y = rareCtsAdj, x = K), col = "blue", size = 2) +
#     ggplot2::geom_line(ggplot2::aes(y = perplexAdj, x = K), col = "red", size = 2) +
#     ggplot2::geom_line(ggplot2::aes(y = cScoresAdj, x = K), col = "green", size = 2) +
#     ggplot2::geom_bar(ggplot2::aes(x = K, y = alphaBool), fill = dat$`alpha < 1`, stat = "identity", width = 1, alpha = 0.5) +
#     ggplot2::scale_y_continuous(
#       name = paste0("# cell-types with mean proportion < ", round(perc.rare.thresh*100, 2), "%"), 
#       breaks = prim_ax_breaks, 
#       labels = prim_ax_labs,
#       sec.axis = ggplot2::sec_axis(
#         ~ ., 
#         name = "<span style='color:red;'>Perplexity</span> (<span style='color:green;'>Coherence Score</span>)", 
#         breaks = sec_ax_breaks, 
#         labels = sec_ax_labels
#       )
#     ) +
#     ggplot2::scale_x_continuous(breaks = min(dat$K):max(dat$K)) +
#     ggplot2::labs(
#       title = "Fitted model K's vs deconvolved cell-types and criteria scores",
#       subtitle = "LDA models with α > 1 shaded"
#     ) +
#     ggplot2::theme_classic() +
#     ggplot2::theme(
#       panel.background = ggplot2::element_blank(),
#       plot.title = ggplot2::element_text(size = 15),
#       panel.grid.minor = ggplot2::element_blank(),
#       panel.grid.major = ggplot2::element_line(color = "black", size = 0.1),
#       panel.ontop = TRUE,
#       axis.title.y.left = ggplot2::element_text(color = "blue", size = 13),
#       axis.text.y.left = ggplot2::element_text(color = "blue", size = 13),
#       axis.title.y.right = ggtext::element_markdown(size = 15, vjust = 1.5),
#       axis.text.y.right = ggtext::element_markdown(size = 13),
#       axis.text.x = ggplot2::element_text(angle = 0, size = 13),
#       axis.title.x = ggplot2::element_text(size = 13)
#     )
#   
#   return(plt)
#   #return(cs)
# }

colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
                   "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
                   "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
                   "#dcf0d0","#f4a99f","#c8bad8",
                   "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
                   "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
                   "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
                   "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785",
                   "#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080")
#colors = colorRampPalette(colorCandidate)(8)
#"#1E77B4" "#D571B0" "#AB6531" "#DCC8B0" "#64ACA7" "#C9596F" "#439691" "#F08080"

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


GeneRank.ProfilePlot = function(deconGexp, celltype, threshold = 5, markers = NULL){
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > threshold))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  if(is.null(markers)){markers <- names(log2fc)[1:4]}
  
  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[match(markers, dat$genes)] <- markers
  
  plt <- ggplot2::ggplot(data = dat) +
    # Set bars to dark blue
    ggplot2::geom_col(aes(x = order, y = values), fill = "darkblue", color = "darkblue", width = 1) +
    
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +
    
    ggplot2::labs(title = paste0('Deconvolved cell-type ', celltype,' transcriptional profile'),
                  x = "Gene expression rank",
                  y = "log2(FC)") +
    
    ggrepel::geom_text_repel(aes(x = order, y = values, label = selectedLabels, size = 1),
                             color = "red",min.segment.length = 0, box.padding = 1) +
    
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.x = ggplot2::element_text(size=15, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  return(plt)
}


DeconPieplot <- function(DeconData,
                         TumorST,
                         spatial_location,
                         img_path,
                         topicCols = NULL,
                         pie_scale = 0.4,
                         scatterpie_alpha = 1,
                         border_color = 'darkgrey') {
  plot_col = colnames(DeconData)
  if(is.null(topicCols)){
    topicCols = hue_pal(l = 80)(length(plot_col))
  }
  DeconData = as.data.frame(DeconData)
  DeconData$cell_ID = rownames(DeconData)
  DeconData <- DeconData[rownames(DeconData) %in% rownames(TumorST@meta.data), ]
  
  ## Preprocess data
  slice <- names(TumorST@images)[1]
  
  spatial_coord <- as.data.frame(spatial_location) %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::mutate(
      imagerow_scaled =
        x * TumorST@images[[slice]]@scale.factors$lowres,
      imagecol_scaled =
        y * TumorST@images[[slice]]@scale.factors$lowres
    ) %>%
    dplyr::inner_join(DeconData, by = "cell_ID")
  
  ### Load histological image into R
  img_path <- img_path # lowers image png(input dir)
  img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))
  
  if (img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  } else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }
  
  # Convert image to grob object
  img_grob <- grid::rasterGrob(img,
                               interpolate = FALSE,
                               width = grid::unit(1, "npc"),
                               height = grid::unit(1, "npc")
  )
  
  a = which(apply(DeconData[,plot_col], 2, function(col) all(col == 0)))
  
  if(length(a) == 0){
    topicCols = topicCols
  }else{
    topicCols = topicCols[-a]
  }
  
  
  ## Plot spatial scatterpie plot
  scatterpie_plt <- ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = 0,
      ymax = -nrow(img)
    ) +
    scatterpie::geom_scatterpie(
      data = spatial_coord,
      ggplot2::aes(
        x = imagecol_scaled,
        y = imagerow_scaled
      ),
      cols = plot_col,
      color = border_color,
      alpha = scatterpie_alpha,
      pie_scale = pie_scale,
      lwd=0.1
    ) +
    ggplot2::scale_fill_manual(values = as.vector(topicCols))+
    ggplot2::scale_y_reverse() +
    ggplot2::ylim(nrow(img), 0) +
    ggplot2::xlim(0, ncol(img)) +
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_void() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = NULL,
      ylim = NULL,
      expand = TRUE,
      clip = "on"
    ) + ggplot2::theme(legend.text=element_text(size = 12))
  return(scatterpie_plt)
}




