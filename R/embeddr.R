#
# R functions for Laplacian Eigenmaps
# High resolution pseudotemporal of single-cell RNA-seq data using
# laplacian eigenmaps and principal curves.
#
# kieran.campbell@sjc.ox.ac.uk
#
# (c) 2015

#' Laplacian eigenmaps embedding of single-cell RNA-seq data.
#'
#' @param sce The SCESet object
#' @param genes_for_embedding A vector of gene indices or names to subset the sce for the embedding. The returned
#' object contains the full original gene set found in sce.
#' @param kernel The choice of kernel. 'nn' will give nearest neighbours, 'dist' gives minimum distance and
#' 'heat' gives a heat kernel. Discussed in detail in 'Laplacian Eigenmaps and Spectral Techniques for Embedding and Clustering',
#' Belkin & Niyogi
#' @param metric The metric with which to assess 'closeness' for nearest neighbour selection, one of
#' "correlation" (pearson) or "euclidean". Default is "correlation".
#' @param nn Number of nearest neighbours if kernel == 'nn'
#' @param eps Maximum distance parameter if kernel == 'dist'
#' @param t 'time' for heat kernel if kernel == 'heat'
#' @param symmetrize How to make the adjacency matrix symmetric. Note that slightly
#' counterintuitively, node i having node j as a nearest neighbour doesn't guarantee node
#' j has node i. There are several ways to get round this:
#' \itemize{
#' \item \code{mean} If the above case occurs make the link weight 0.5 so the adjacency matrix becomes \eqn{0.5(A + A')}
#' \item \code{ceil} If the above case occurs set the link weight to 1 (ie take the ceiling of the mean case)
#' \item \code{floor} If the above case occurs set the link weight to 0 (ie take the floor of the mean case)
#' }
#' @param measure_type Type of laplacian eigenmap, which corresponds to the constraint on the eigenvalue problem. If
#' type is 'unorm' (default), then the graph measure used is the identity matrix, while if type is 'norm' then the measure
#' used is the degree matrix.
#' @param p Dimension of the embedded space, default is 2
#'
#' @import scater
#' @importFrom Biobase exprs
#'
#' @export
#' @return An object of class SCESet
#' 
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
embeddr <- function(sce, genes_for_embedding = NULL,
                    kernel=c('nn','dist','heat'),
                    metric=c('correlation', 'euclidean','cosine'),
                    nn = round(log(ncol(sce))), eps = NULL, t = NULL,
                    symmetrize = c('mean','ceil','floor'),
                    measure_type = c('unorm','norm'), p = 2) {

  sce <- remove_pseudotime(sce)

  if(is.null(genes_for_embedding)) genes_for_embedding <- 1:dim(sce)[1]
  W <- weighted_graph(exprs(sce[genes_for_embedding,]), kernel = kernel, metric = metric,
                      nn = nn, eps = eps, t = t, symmetrize = symmetrize)

  cellDist(sce) <- W

  LE <- laplacian_eigenmap(W, measure_type = measure_type, p = p)

  redDim(sce) <- as.matrix(LE$embedding)
  pData(sce)$connected_component <- LE$connected_component

  validObject( sce )
  return( sce )
}

#' Construct a weighted graph adjacency matrix
#'
#' @param x The feature-by-sample (e.g. genes are rows, cells are columns) data matrix
#' @param kernel The choice of kernel. 'nn' will give nearest neighbours, 'dist' gives minimum distance and
#' 'heat' gives a heat kernel. Discussed in detail in 'Laplacian Eigenmaps and Spectral Techniques for Embedding and Clustering',
#' Belkin & Niyogi
#' @param metric The metric with which to assess 'closeness' for nearest neighbour selection, one of
#' "correlation" (pearson) or "euclidean". Default is "correlation".
#' @param nn Number of nearest neighbours if kernel == 'nn'
#' @param eps Maximum distance parameter if kernel == 'dist'
#' @param t 'time' for heat kernel if kernel == 'heat'
#' @param symmetrize How to make the adjacency matrix symmetric. Note that slightly
#' counterintuitively, node i having node j as a nearest neighbour doesn't guarantee node
#' j has node i. There are several ways to get round this:
#' \itemize{
#' \item \code{mean} If the above case occurs make the link weight 0.5 so the adjacency matrix becomes \eqn{0.5(A + A')}
#' \item \code{ceil} If the above case occurs set the link weight to 1 (ie take the ceiling of the mean case)
#' \item \code{floor} If the above case occurs set the link weight to 0 (ie take the floor of the mean case)
#' }
#'
#' @importFrom lsa cosine
#'
#' @export
#' @return An n by n adjacency matrix
#' @examples 
#' x <- matrix(rnorm(50), ncol = 10) # synthetic data matrix
#' wg <- weighted_graph(x)
weighted_graph <- function(x, kernel = c('nn','dist','heat'),
                           metric=c('correlation', 'euclidean', 'cosine'),
                           nn = round(log(ncol(x))), eps = NULL, t = NULL,
                           symmetrize = c('mean','ceil','floor')) {
  ## sanity checking
  kernel <- match.arg(kernel)
  symmetrize <- match.arg(symmetrize)
  metric <- match.arg(metric)

  if(kernel == 'heat' && is.null(t)) stop('Please specify t for using heat kernel')
  if(kernel == 'dist' && is.null(eps)) stop('Please specify eps for using distance kernel')

  ## compute distance matrix
  dm <- NULL
  if(metric == 'euclidean') {
    dm <- as.matrix(dist(t(x)))
  } else if(metric == 'correlation') {
    dm <- cor(x, method='pearson')
  } else if(metric == 'cosine') {
    dm <- cosine(x)
  }

  W <- NULL
  if(kernel == 'nn' || kernel == 'heat') {
    if(metric == 'euclidean') {
      W <- apply(dm, 2, function(r) {
        ind <- order(r)[1:nn]
        rep_vec <- rep(0, length(r))
        rep_vec[ind] <- 1
        names(rep_vec) <- names(r)
        return( rep_vec )
      })
    } else if(metric == 'correlation' | metric == 'cosine') {
      W <- apply(dm, 2, function(r) {
        ind <- order(r, decreasing = TRUE)[1:nn]
        rep_vec <- rep(0, length(r))
        rep_vec[ind] <- 1
        names(rep_vec) <- names(r)
        return( rep_vec )
      })
    }
    ## symmetrise as A -> B => B -> A
    W <- 0.5 * (W + t(W))
    if(symmetrize == 'ceil') {
      W <- ceiling(W)
    } else if(symmetrize == 'floor') {
      W <- floor(W)
    }
  }

  if(kernel == 'heat') {
    if(metric == 'correlation') stop('Heat hernel not supported for correlation metric')
    We <- exp(-dm * dm / t)
    W <- W * We
  }

  if(kernel == 'dist') {
    W <- apply(dm, 1, function(r) as.numeric(r < eps))
  }
  return( W )
}

#' Laplacian eigenmaps
#'
#' Construct a laplacian eigenmap embedding
#' @param W The weighted graph adjacency matrix
#' @param measure_type Type of laplacian eigenmap (norm for normalised, unorm otherwise)
#' @param p Dimension of the embedded space, default is 2
#' @param eig_tol The value below which an eigenvalue is deemed to be 0. Ideally these would all
#' be \code{.Machine$double.eps} but numerical instabilities mean they are often larger. A good
#' cutoff seems to be around 1e-10. If more than one zero eigenvalue is found but this is inconsistent
#' with the number of connected clusters then the user will be warned to adjust eig_tol.
#'
#' @import igraph
#'
#' @export
#'
#' @return A list containing two entries:
#' \itemize{
#' \item{\code{embedding}}{A \code{p} by \code{nrow(W)} matrix with the laplacian eigenmaps embedding}
#' \item{\code{connected_components}}{A vector of length \code{nrow(W)} of integers specifying to which
#' connected component each cell belongs}
#' }
#' @examples 
#' ## create synthetic weight matrix
#' W <- matrix(sample(0:1, 25, replace = TRUE), nrow=5)
#' diag(W) <- 0 ; W <- 0.5(W + t(W))
#' le <- laplacian_eigenmap(W) # returns a list with entries embedding and connected_components
laplacian_eigenmap <- function(W, measure_type = c('unorm','norm'), p = 2, eig_tol = 1e-10) {
  measure_type <- match.arg(measure_type)
  if(nrow(W) != ncol(W)) {
    print('Input weight matrix W must be symmetric')
    return( NULL )
  }

  L  <- diag(rowSums(W)) - W
  l <- nrow(L)
  M <- NULL # object to be returned
  eig <- NULL # eigenvectors/values

  if(measure_type == 'norm') {
    Ds <- diag(1/sqrt(rowSums(W)))
    L_sym <- Ds %*% L %*% Ds
    eig <- eigen(L_sym, symmetric = TRUE) # values ordered in decreasing order
    M <- diag(Ds) * eig$vectors[,(l-1):(l-p)]
  }
  if(measure_type == 'unorm') {
    eig <- eigen(L, symmetric = TRUE) # values ordered in decreasing order
    M <- eig$vectors[,(l-1):(l-p)]
  }

  ## deal with multiple connected components
  cc <- NULL

  if(sum(eig$values < eig_tol) > 1) {
    warning(paste('More than one non-zero eigenvalue - disjoint clusters.'))
    g <- graph.adjacency(ceiling(W), mode = 'undirected')
    if(no.clusters(g) != sum(eig$values < eig_tol)) {
      stop("Number of connected components doesn't match number of eigenvalues of the laplacian identified as zero. Consider changing eig_tol")
    }
    cc <- clusters(g)$membership
  } else {
    cc <- rep(1, l)
  }

  colnames(M) <- paste0('component_', 1:p)
  M <- data.frame(M)

  return( list(embedding = M, connected_components = cc) )
}

#' removes the pseudotime and low dimensional trajectory
#' @param sce The \code{SCESet}
#'
#' @return An \code{SCESet} with the pseudotime & trajectory removed
remove_pseudotime <- function(sce) {
  sce$pseudotime <- NULL
  sce$trajectory_1 <- NULL
  sce$trajectory_2 <- NULL
  return( sce )
}



#' Cluster the resulting embedding
#'
#' Cluster the embedded representation using either kmeans or mixture models
#' from the mclust package
#'
#' @param sce The SCESet object
#' @param k The number of clusters to find in the data
#' @param method Either 'kmeans' or 'mm' to use \code{mclust}
#' @export
#' @import mclust
#' @importFrom dplyr select
#' @importFrom Biobase phenoData
#' @importFrom Biobase phenoData<-
#'
#' @return The dataframe M with a new numeric variable `cluster` containing the assigned cluster
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- cluster_embedding(sce)
cluster_embedding <- function(sce, method=c('mm', 'kmeans'), k = NULL) {
  M_xy <- redDim(sce)
  method <- match.arg(method)
  if(method == 'kmeans') {
    km <- kmeans(M_xy, k)
    phenoData(sce)$cluster <- km$cluster
  } else if(method == 'mm') {
    mc <- Mclust(M_xy, G=k)
    phenoData(sce)$cluster <- as.factor(mc$classification)
  }
  return( sce )
}

# fit_pseudotime_thinning <- function(M, clusters=NULL, ...) {
#   load_all("/net/isi-scratch/kieran/embeddr/curver/")
#
#   Mp <- M
#   if(!is.null(clusters)) Mp <- filter(Mp, cluster %in% clusters)
#   Y <- curver::reconstruct(select(Mp, component_1, component_2), niter=1, h=0.2)
#   #Y <- curver::reconstruct(select(Mp, component_1, component_2), ...)
#
#   D <- as.matrix(dist(Y))
#   g <- graph.adjacency(D, weighted=TRUE, mode='undirected')
#   g_mst <- minimum.spanning.tree(g)
#   A <- as.matrix(get.adjacency(g_mst))
#   endpoints <- which(rowSums(A) == 1)
#   ordering <- get.shortest.paths(g_mst, from=endpoints)
#   paths <- ordering$vpath
#   best_path <- paths[[ which.max(sapply(paths, length)) ]]
#
#   Z <- Y[ordering$vpath[[1]],]
#
#   ## now we have the ordering want to work out the arc-length
#   n <- dim(Z)[1]
#   Z_start <- Z[1:n-1,]
#   Z_end <- Z[2:n,]
#   Z_diff <- Z_end - Z_start
#   pst <- sqrt(rowSums(Z_diff^2))
#   pseudotime <- c(0, cumsum(pst))
#   Mp$pseudotime <- pseudotime[invPerm(as.integer(ordering$vpath[[1]]))]
#   Mp$trajectory_1 <- Y[,1]
#   Mp$trajectory_2 <- Y[,2]
#   return(Mp)
# }

#' Fit the pseudotime curve
#'
#' Fits the pseudotime curve using principal curves from the princurve library
#'
#' @param sce The SCESet object
#' @param clusters The (numeric) clusters to use for the curve fitting. If NULL (default) then
#' all points are used
#' @param ... Additional arguments to be passed to \code{princurve} from \pkg{principal.curve}.
#'
#' @importFrom dplyr select
#' @importFrom princurve principal.curve
#' @import scater
#'
#' @export
#'
#' @return The dataframe \code{M} with three new variables:
#' \describe{
#' \item{pseudotime}{The pseudotime of the cell (arc-length from beginning of curve)}
#' \item{trajectory_1}{The x-coordinate of a given cell's projection onto the curve}
#' \item{trajectory_2}{The y-coordinate of a given cell's projection onto the curve}}
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
fit_pseudotime <- function(sce, clusters = NULL, ...) {
  component_1 <- component_2 <- NULL # satisfy R CMD check
  M <- as.data.frame(redDim(sce)) # select(pData(sce), component_1, component_2)
  n_cells <- dim(sce)[2]

  cells_in_cluster <- rep(TRUE, n_cells)
  if(!is.null(clusters)) cells_in_cluster <- pData(sce)$cluster %in% clusters

  Mcl <- M[cells_in_cluster, ]
  pc <- principal.curve(x = as.matrix(select(Mcl, component_1, component_2)), ...)
  pst <- pc$lambda
  pst <- (pst - min(pst)) / (max(pst) - min(pst))

  pseudotime <- trajectory_1 <- trajectory_2 <- rep(NA, n_cells)
  pseudotime[cells_in_cluster] <- pst
  trajectory_1[cells_in_cluster] <- pc$s[,1]
  trajectory_2[cells_in_cluster] <- pc$s[,2]

  phenoData(sce)$pseudotime <- pseudotime
  phenoData(sce)$trajectory_1 <- trajectory_1
  phenoData(sce)$trajectory_2 <- trajectory_2
  return( sce )
}

#' Plot the cells in the embedding
#'
#' This function takes a data frame with at least positional (component_0 & component_1) information
#' and plots the resulting embedding. If clusters are assigned it can colour by these, and if a pseudotime
#' trajectory is assigned it will plot this through the embeddding.
#'
#' @param sce The SCESet object
#' @param color_by The variable to color the embedding with (defaults to cluster)
#' @param plot_genes A vector of gene names or indices (anything with which to subset the SCESet) for
#' plotting gene expression as a size in the reduced space. This will facet wrap, so an individual plot
#' for each gene (don't choose too many).
#' @param use_short_names If fData(sce)$gene_short_names is defined and this is true then those are used
#' on the plot
#' @param plot_pseudotime If TRUE (default) the pseudotime curve (principal curve) is shown on the plot. If
#' pData(sce)$pseudotime is null then this is ignored
#'
#' @import ggplot2
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @export
#'
#' @return A \pkg{ggplot2} plot
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' plot_embedding(sce)
plot_embedding <- function(sce, color_by = 'cluster',
                           plot_genes = NULL, use_short_names = FALSE,
                           plot_pseudotime = TRUE) {
  ## satisfy R CMD check
  component_1 <- component_2 <- NULL
  trajectory_1 <- trajectory_2 <- NULL
  plot_names <- NULL

  M <- as.data.frame(redDim(sce)) # dplyr::select(pData(sce), component_1, component_2)

  if(ncol(M) < 2) stop('Please call embeddr on SCESet first')

  ## are we colouring by a factor?
  if(color_by %in% names(pData(sce))) {
    col <- match(color_by, names(pData(sce)))
    M <- cbind(M, dplyr::select(pData(sce), col))
  }

  ## are we plotting pseudotime?
  if(('pseudotime' %in% names(pData(sce))) && plot_pseudotime) {
    M <- cbind(M, select(pData(sce), pseudotime, trajectory_1, trajectory_2))
  }

  ## are we sizing by gene expression?
  if(!is.null(plot_genes)) {
    if(is.character(plot_genes)) {
      plot_names <- plot_genes
    } else {
      plot_names <- featureNames(sce)[plot_genes]
    }
    y <- exprs(sce[plot_genes,])
    y <- t(y)
    y <- data.frame(y)
    M <- cbind(M, y)
    M <- melt(M, id.vars=setdiff(names(M), plot_names),
                       variable.name='gene', value.name='count')
   if(use_short_names) M$gene <- mapvalues(M$gene, from = featureNames(sce), to = fData(sce)$gene_short_name)
  }

  if(('pseudotime' %in% names(pData(sce))) && plot_pseudotime) M <- arrange(M, pseudotime)


  plt <- ggplot(data=M)

  if(color_by %in% names(M)) { # colouring by something
    if(color_by == 'cluster') {
      mapping_str <- paste0("as.factor(", color_by, ")")
    } else {
      mapping_str <- color_by
    }
    if(is.null(plot_genes)) {
      plt <- plt + geom_point(aes_string(x = "component_1", y = "component_2",
                                         fill = mapping_str), color = 'gray20',
                              alpha = 0.65, size = 3.5, shape = 21)
    } else {
      plt <- plt + geom_point(aes_string(x = "component_1", y = "component_2",
                                         fill = mapping_str, size = "count"), color = 'gray20',
                              alpha = 0.65, shape = 21)
    }

    if(is.numeric(pData(sce)[[color_by]]) && all(pData(sce)[[color_by]] %% 1 != 0)) {
      plt <- plt + scale_fill_continuous(name = color_by)
    } else {
      plt <- plt + scale_fill_discrete(name = color_by)
    }

  } else { # not colouring by something
    warning(paste('color_by string',color_by,'not found in pData(SCESet)'))
    if(is.null(plot_genes)) {
      plt <- plt + geom_point(aes_string(x = "component_1", y = "component_2"),
                              alpha = 0.65, size = 3.5)
    } else {
      plt <- plt + geom_point(aes_string(x = "component_1", y = "component_2", size = "count"),
                              alpha = 0.65)
    }
  }

  if(("pseudotime" %in% names(M)) && plot_pseudotime) {
    ## curve has been fit so plot all
    plt <- plt + geom_path(aes_string(x = "trajectory_1", y = "trajectory_2"),
                           data=M, color='black',
                           size = 1.5, alpha = 0.8, linetype=2)
  }
  plt <- plt + xlab('Component 1') + ylab('Component 2')

  if ( library(cowplot, logical.return = TRUE) ) {
    plt <- plt + cowplot::theme_cowplot()
  } else {
    plt <- plt + theme_bw()
  }

  if(!is.null(plot_genes)) plt <- plt + facet_wrap( ~ gene )

  return( plt )
}

#' Plot cells in pseudotime
#'
#' Plot a set of genes through pseudotime using smooth LOESS curves.
#'
#'  @param sce An object of class \code{SCESet}
#'  @param nrow Number of rows of plots; passed to \code{facet_wrap}
#'  @param ncol Number of columns of plots; passed to \code{facet_wrap}
#'  @param use_short_names Logical If \code{pData(sce)} contains a \code{gene_short_name}
#'  column (such as in the \pkg{monocle} dataset \code{HSMM}) then the names in the resulting
#'  plot will be the gene short names.
#'  @param use_log_scale Logical If TRUE scale_y_log10 is added to plot
#'  @param facet_wrap_scales Passed to the scales argument in \code{facet_wrap}.
#'  Should scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x", "free_y")
#'
#'  @export
#'  @import ggplot2
#'  @importFrom reshape2 melt
#'  @importFrom Biobase fData
#'
#'  @return A \pkg{ggplot2} plot
#'  
#'  @examples 
#'  library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' plot_in_pseudotime(sce[1:4,]) # plot first four genes
plot_in_pseudotime <- function(sce, nrow = NULL, ncol = NULL,
                               use_short_names = FALSE, use_log_scale = FALSE,
                               facet_wrap_scales = "fixed") {
  xp <- data.frame(t(exprs(sce)), check.names = FALSE) # now cell-by-gene
  if(use_short_names) names(xp) <- fData(sce)$gene_short_name

  xp$pseudotime <- pData(sce)$pseudotime

  cn <- 'cluster' %in% names(pData(sce))
  if(cn) xp$cluster <- as.factor(pData(sce)$cluster)
  id_vars <- 'pseudotime'
  if(cn) id_vars <- c(id_vars, 'cluster')

  df_x <- melt(xp, id.vars=id_vars, variable.name='gene', value.name='counts')

  plt <- NULL
  if(cn) {
    plt <- ggplot(data=df_x, aes_string(x="pseudotime", y="counts", color="cluster"))
  } else {
    plt <- ggplot(data=df_x, aes_string(x="pseudotime", y="counts"))
  }
  if(use_log_scale) plt <- plt + scale_y_log10()

  return( plt + geom_point(size=1.5) +
    theme_bw() + geom_smooth(method='loess', color='firebrick') +
      facet_wrap(~ gene, nrow = nrow, ncol = ncol, scales = facet_wrap_scales) +
    ylab('Normalised log10(FPKM)') )
}

#' Reverse pseudotime
#'
#' Reverse the pseudotimes of cells via pseudotime x -> -x + max(x) + min(x)
#'
#' @param sce An object of class \code{SCESet}
#'
#' @export
#' @return \code{sce} with the pseudotime reversed.
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' sce <- reverse_pseudotime(sce)
reverse_pseudotime <- function(sce) {
  reverse <- function(x) -x + max(x) + min(x)
  phenoData(sce)$pseudotime <- reverse(pData(sce)$pseudotime)
  return( sce )
}

# plot_heatmap <- function(sce, ...) {
#   xp <- exprs(sce)[,order(pData(sce)$pseudotime)]
#
#   heatmap.2(xp, dendrogram="none", Colv=FALSE,
#             col=redblue(256), trace="none", density.info="none", scale="row", ...)
# }

#' Plot the nearest neighbour graph in the embedding
#'
#' This plots the cells in the reduced dimension and connects cells if they are
#' connected in the nearest neighbour graph W. A vertice between cells i and j is
#' coloured red if both i and j are nearest neigbours of each other and grey otherwise.
#'
#' @param sce An object of class \code{SCESet}
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr select
#' @importFrom plyr mapvalues
#' @import scater
#'
#' @return A ggplot graphic
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' plot_graph(sce)
plot_graph <- function(sce) {
  ## satisfy R CMD check
  component_1 <- component_2 <- NULL
  x <- y <- NULL

  df <- dplyr::rename(as.data.frame(redDim(sce)), x = component_1, y = component_2)
                                    #select(pData(sce), x = component_1, y = component_2)
  W <- cellDist(sce)

  diag(W) <- 0
  locs <- which((1 * lower.tri(W) * W) > 0, arr.ind = TRUE)
  from_to <- apply(locs, 1, function(xy) {
    xy_from <- df[xy[1],]
    xy_to <- df[xy[2],]
    xx <- c(xy_from, xy_to)
    names(xx) <- c('x_from','y_from','x_to','y_to')
    unlist(xx)
  })
  colnames(from_to) <- NULL
  from_to <- data.frame(t(from_to))
  from_to$connect <- mapvalues(W[locs], c(0.5, 1), c('Single','Both'))
  cols <- c('Single' = 'grey', 'Both' = 'red')

  plt <- ggplot() +
    geom_segment(data=from_to, aes_string(x="x_from", xend="x_to",
                                          y="y_from", yend="y_to", color="connect"),
                 alpha=0.5, linetype=1) +
        theme_minimal() + scale_color_manual(values = cols) +
    geom_point(data=df, aes(x=x,y=y), alpha = 0.5)
  plt + xlab('x') + ylab('y')
  plt
}

## expression~VGAM::bs(Pseudotime, df=3)

#' Fit the gene expression profile in pseudotime
#'
#' This function fits the expression profile of y as a function of
#' pseudotime using a natural cubic spline of degree three. A tobit model
#' is used to censor values less than min_expr
#'
#' @param sce An object of type SCESet
#' @param gene The gene name to fit
#' @export
#' @importFrom AER tobit
#' @importFrom splines bs
#' @importFrom survival survreg
#' @importFrom survival Surv
#'
#' @return An object of class AER::tobit
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' model <- fit_pseudotime_model(sce, 1) # fit for first gene
fit_pseudotime_model <- function(sce, gene) {
  t <- pData(sce)$pseudotime
  y <- exprs(sce)[gene ,]
  min_expr <- sce@lowerDetectionLimit # see paper
  b <- bs(t, df=3)
  fit <- AER::tobit(y ~ b, left = min_expr)
  return( fit )
}

#' Fit the null pseudotime model
#'
#' This function fits the null expression profile in y (ie y ~ 1). A tobit model is
#' used to censor values less than min_expr
#'
#' @param sce The SCESet object
#' @param gene The gene name to fit
#' @export
#'
#' @return An object returned by AER::tobit
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' model <- fit_null_model(sce, 1) # fit for first gene
fit_null_model <- function(sce, gene) {
  y <- exprs(sce)[gene,]
  min_expr <- sce@lowerDetectionLimit
  return( AER::tobit(y ~ 1, left = min_expr) )
}

#' Plot the fit in pseudotime
#'
#' @param sce An object of class SCE set
#' @param models A list of models returned by fit_pseudotime_models. If NULL (default), these will be
#' re-calculated
#' @param n_cores Number of cores to use when calculating the pseudotime models if `model` is NULL
#' @param use_log_scale Logical If TRUE scale_y_log10 is added to plot
#' @param facet_wrap_scale Passed to the scales argument in \code{facet_wrap}.
#'  Should scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x", "free_y")
#'  @param color_by The variable in pData(sce) by which to colour the points
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#'
#' @return An plot object from \code{ggplot}
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' plot_pseudotime_model(sce[1,]) # select first gene
plot_pseudotime_model <- function(sce, models = NULL, n_cores = 2,
                                  facet_wrap_scale = 'fixed',  use_log_scale = FALSE,
                                  color_by = NULL) {
  if(is.null(models)) models <- fit_pseudotime_models(sce, n_cores)
  # if(class(models) == 'list' && dim(sce)[1] != length(models)) stop('Must have a fitted model for each gene in sce')

  gene_names <- NULL
  if('gene_short_name' %in% names(fData(sce))) {
    gene_names <- fData(sce)$gene_short_name
  } else {
    gene_names <- featureNames(sce)
  }

  min_expr <- sce@lowerDetectionLimit
  y <- data.frame(t(exprs(sce)))
  names(y) <- gene_names

  y$pseudotime <- pseudotime(sce)
  y_melted <- melt(y, id.vars='pseudotime', value.name='exprs', variable.name='gene')
  pe <- data.frame(predicted_expression(sce, models))
  names(pe) <- gene_names
  pe$pseudotime <- pseudotime(sce)
  pe_melted <- melt(pe, id.vars='pseudotime', value.name='predicted', variable.name='gene')

  df <- dplyr::full_join(y_melted, pe_melted, by=c('pseudotime','gene'))
  df$predicted[df$predicted < min_expr] <- min_expr
  df$min_expr <- min_expr

  plt <- ggplot(df)  + geom_line(aes_string(x = "pseudotime", y = "predicted"), color = 'red') +
    theme_minimal() + geom_line(aes_string(x = "pseudotime", y = "min_expr"), color='grey', linetype=2) +
    scale_color_manual('', values=c('Min expr' = 'grey', 'Predicted'='red')) +
    facet_wrap(~ gene, scales = facet_wrap_scale) +
    ylab('expression')
  if(is.null(color_by)) {
    plt <- plt + geom_point(aes_string(x = "pseudotime", y = "exprs"))
  } else {
    if(color_by %in% names(pData(sce))) {
      plt <- plt + geom_point(aes_string(x = "pseudotime", y = "exprs", color = color_by))
    } else {
      plt <- plt + geom_point(aes_string(x = "pseudotime", y = "exprs"))
    }
  }
  if(use_log_scale) plt <- plt + scale_y_log10()
  return( plt )
}

# plot_pseudotime_models <- function(models, x, min_expr) {
#
# }


#' Perform likelihood ratio test
#'
#' @param model The full model y ~ pseudotime returned by AER::tobit
#' @param null_model The null model y ~ 1 returned by AER::tobit
#'
#' @export
#' @importFrom lmtest lrtest
#'
#' @return The p-value
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' full_model <- fit_pseudotime_model(sce, 1) # fit for first gene
#' null_model <- fit_null_model(sce, 1) # fit for first gene
#' pval <- compare_models(full_model, null_model)
compare_models <- function(model, null_model) {
  lrt <- lrtest(model, null_model) # VGAM <-> lmtest
  return( lrt$"Pr(>Chisq)"[2] ) # lmtest
  #lrt@Body["Pr(>Chisq)"][2, ]  # VGAM
}

#' Test a single gene as a function of pseudotime
#'
#' @param gene Name or index of the gene to be tested. 
#' @param sce An object of class SCESet
#' @param full_model If a full pseudotime model has already been calculated, use that. Otherwise (if NULL, default),
#' recalculate one on the fly.
#' 
#' @return A p-value for the given gene
#'
#' @export
#' @examples 
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' pval <- gene_pseudotime_test(sce, 1) # select first gene
gene_pseudotime_test <- function(sce, gene, full_model = NULL) {
  tryCatch({
    if(is.null(full_model)) {
      model <- embeddr::fit_pseudotime_model(sce, gene)
    } else {
      model <- full_model
    }
    null_model <- embeddr::fit_null_model(sce, gene)
    return( embeddr::compare_models(model, null_model) )
  }, error = function(e) {
    return( 1 ) # if there's an error, just return -1
  })
}

#' Pseudotime gene testing
#'
#' This function fits both smoothing B-spline tobit regression models for \strong{all}
#' genes in \code{sce} and computes p-values by comparing to the null model using
#' a likelihood ratio test.
#'
#' @param sce An object of class \code{SCESet}. If only a certain number of genes are desired to
#' be tested then \code{sce} must be subsetted before the function call.
#' @param n_cores The number of cores used in the call to \code{mcapply}
#'
#' @return A \code{data.frame} with three columns: the gene name, the unadjusted p-value and
#' the Benjamini-Hochberg adjusted q-value. Note that different multiple testing corrections
#' can be applied using the R function \code{p.adjust}.
#'
#' @export
#' @examples 
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' df_pval <- pseudotime_test(sce[1:4,]) # only for first four genes
pseudotime_test <- function(sce, n_cores = 2) {
  p_vals <- NULL
  if(n_cores == 1) {
    p_vals <- sapply(featureNames(sce),
                     function(gene_name) gene_pseudotime_test(gene_name, sce))
  } else {
    p_vals <- unlist(mclapply(featureNames(sce),
                              function(gene_name) gene_pseudotime_test(gene_name, sce)))
  }
  q_vals <- p.adjust(p_vals, method='BH')
  return( data.frame(gene=featureNames(sce), p_val = p_vals, q_val = q_vals))
}

#' Generate a list of pseudotime models corresponding to ALL genes in sce
#'
#' This will iterate over every gene in the \code{SCESet} and produce a model fit
#' for each. \strong{WARNING}: The list returned can be huge, so use only on a few genes.
#'
#' @param sce An object of class \code{SCESet}
#' @param n_cores The number of cores to use in the call to \code{mclapply}
#' @importFrom parallel mclapply
#'
#' @export
#' @return A list of models returned by AER::tobit
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' models <- fit_pseudotime_models(sce, n_cores = 1)
fit_pseudotime_models <- function(sce, n_cores = 2) {
  models <- NULL
  fpm_wrapper <- function(gene, sce) fit_pseudotime_model(sce, gene)
  if(n_cores == 1) {
    models <- lapply(featureNames(sce), fpm_wrapper, sce)
  } else {
    models <- mclapply(featureNames(sce), fpm_wrapper, sce, mc.cores = n_cores)
  }
  names(models) <- featureNames(sce)
  return( models )
}

#' Create a predicted expression matrix
#'
#' Given a list of models return a matrix corresponding to the prediction from the models.
#' Each column represents a gene and each row its expression at a given point in pseudotime.
#'
#' @param sce An object of class \code{SCESet}
#' @param models An object representing models. If of type \code{list} then for each element
#' the predicted expression is computed and a matrix returned. If of type \code{model} for
#' which a \code{predict} function is available, then a single vector corresponding to
#' \code{predict(model)} is returned. If NULL then the model is computed for all genes in
#' \code{sce} and the resulting list returned.
#' @param n_cores The number of cores to pass to \code{mclapply}.
#' 
#' @return A dataframe of predicted expression where rows are cells (pseudotime) and
#' columns are genes
#'
#' @importFrom Biobase featureNames
#' @export
#' @examples 
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' pe <- predicted_expression(sce[1:4,]) # use first four genes
predicted_expression <- function(sce, models = NULL, n_cores = 2) {
  if(!is.null(models)) {
    if(class(models) == 'list') {
      predict_list <- lapply(models, function(x) predict(x))
      return( do.call('cbind', predict_list) )
    } else {
      return(predict(models))
    }
  } else {
    ## want to calculate models one-at-a-time to make sure they don't take
    ## up too much space in memory
    names_not_null <- NULL
    predict_list <- lapply(featureNames(sce), function(gene_name) {
      fit <- fit_pseudotime_model(sce, gene_name)
      if(is.null(fit)) {
        warning(paste('Fit', gene_name, 'returned null'))
        return(NULL)
      }
      names_not_null <<- c(names_not_null, gene_name)
      return(predict(fit))
    })
    pred_mat <- do.call('cbind', predict_list)
    colnames(pred_mat) <- names_not_null
    return(pred_mat)
  }
}

#' Plot density of cells in pseudotime
#'
#' This returns a \code{ggplot} density plot of the cells across pseudotime.
#'
#' @param sce An object of class SCESet
#' @param reverse Logical If true the pseudotime will be reversed.
#' @param color_by The variable (in \code{pData(sce)}) by which to colour the density plot.
#'
#' @import ggplot2
#' @export
#' @return A `ggplot` object
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' plot_pseudotime_density(sce)
plot_pseudotime_density <- function(sce, color_by = NULL, reverse = FALSE) {
  if(reverse) sce <- reverse_pseudotime(sce)

  if(!('pseudotime' %in% names(pData(sce)))) stop('Fit pseudotime before calling plot_pseudotime_density')

  plt <- ggplot(pData(sce))
  if(!is.null(color_by)) {
    if(color_by %in% names(pData(sce))) {
    plt <- plt + geom_density(aes_string(x = "pseudotime", fill = color_by), alpha = 0.6)
   } else {
    warning(paste('Variable',color_by,'not found in names(pData(sce))'))
    plt <- plt + geom_density(aes_string(x = "pseudotime"), alpha = 0.6)
   }
    }
  else {
     plt <- plt + geom_density(aes_string(x = "pseudotime"))
  }
  plt <- plt + theme_bw() + xlab('Pseudotime') + ylab('Cellular density')
  return( plt )
}

#' Plot metrics in pseudotime
#'
#' Plot various metrics (mean, variance, CV2 and signal-to-noise ratio) using
#' a sliding window approach. Calculates metrics if not already done so.
#'
#' @param sce An object of class \code{SCESet}
#' @param ... Additional arguments to be passed to \code{calculate_metrics} concerning
#' how the metrics are calculated.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @return A ggplot graphic
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' plot_pseudotime_metrics(sce)
plot_pseudotime_metrics <- function(sce, ...) {
  df_window <- calculate_metrics(sce, ...)
  df_window <- melt(df_window, id='t')
  ggplot(df_window, aes_string(x = "t", y = "value", color = "variable")) +
    geom_line() + theme_bw() + xlab('Pseudotime')
}

#' Calculate metrics through pseudotime
#'
#' Calculate various metrics (mean, variance, CV2 and signal-to-noise ratio) using
#' a sliding window approach
#'
#' @param sce A \code{SCESet} object
#' @param gene The gene of interest on which to calculate the metrics
#' @param window_size The size of the sliding window. By default taken to be half the number of cells
#'
#' @export
#' @importFrom Biobase pData
#' @return An object of class `data.frame` where each column is a metric (window-averaged pseudotime,
#' mean, variance, CV2, signal to noise ratio)
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' metrics <- calculate_metrics(sce, 1) 
calculate_metrics <- function(sce, gene, window_size=NULL) {
  t <- pseudotime(sce)
  y <- exprs(sce)[gene,]
  y <- y[order(t)]
  t <- sort(t,decreasing = FALSE)

  if(is.null(window_size)) window_size <- length(y) / 2

  n_points <- length(y) - window_size + 1
  vt <- sapply(1:n_points, function(i) {
    interval <- i:(i+window_size - 1)
    x <- mean(t[interval])
    y_mean <- mean(y[interval])
    y_var <- var(y[interval])
    cv <- sqrt(y_var) / y_mean
    snr <- y_mean / sqrt(y_var)
    return(c(x,y_mean,y_var, cv, snr))
  })
  df_window <- data.frame(t(vt))
  names(df_window) <- c('t','Mean','Variance','CV','Signal-to-noise ratio')
  return(df_window)
}

#' Retrieve the pseudotime assignment from sce
#'
#' Equivalent to \code{pData(sce)$pseudotime}
#'
#' @param sce An object of class \code{SCESet}
#' @export
#'
#' @return A numeric vector of pseudotimes
#' @examples
#' library(scater)
#' data("sc_example_counts") ; example_sceset <- newSCESet(countData = sc_example_counts)
#' sce <- embeddr(sce)
#' sce <- fit_pseudotime(sce)
#' pseudotime(sce)
pseudotime <- function(sce) {
  if(class(sce) != 'SCESet') stop('sce must be of class SCESet')
  return(pData(sce)$pseudotime)
}

# #' Log the expression of a given SCESet
# #'
# #' @param sce An object of class SCESet
# #' @param base The base of the logarithm
# #' @param pseudocount The count to add to the expression to avoid taking
# #' the logarithms of negative numbers. If NULL this will
# #' be worked out automatically.
# #'
# #' @importFrom Biobase exprs<-
# #'
# #' #@export
# #' @return An SCESet with the expression logged
# log_expression <- function(sce, base = 10, pseudocount = NULL) {
#   if(is.null(pseudocount)) {
#       min_exprs <- min(as.vector(exprs(sce)))
#       if(min_exprs == 0) {
#         pseudocount <-  1
#       } else if(min_exprs < 0) {
#         stop('Expression contains negative values')
#       } else {
#         pseudocount <- 0
#       }
#   }
#   exprs(sce) <- log(exprs(sce) + pseudocount, base = base)
#   sce@lowerDetectionLimit <- log(sce@lowerDetectionLimit + pseudocount, base = base)
#   sce@logged <- TRUE
#   return( sce )
# }
