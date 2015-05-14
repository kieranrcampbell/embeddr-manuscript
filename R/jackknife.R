library(devtools)


## kieranrcampbell@gmail.com
## pseudotime jackknifing to estimate uncertainty

load_all('/net/isi-scratch/kieran/embeddr/embeddr/')
load('/net/isi-scratch/kieran/embeddr/embeddr/data/sce_23.Rdata')

sce <- sce_23
n_cells <- dim(sce)[2]

pseudotime_all <- pseudotime(fit_pseudotime(sce))

pseudotime_jk <- sapply(1:n_cells, function(i) {
  sce_mi <- sce[,-i]
  pj <- pseudotime(fit_pseudotime(sce_mi))
  
  pst <- rep(NA, n_cells)
  pj_inds <- match(rownames(pData(sce_mi)), rownames(pData(sce)))
  pst[pj_inds] <- pj
  pst  
})

jk_estimate <- rowMeans(pseudotime_jk, na.rm = TRUE)

