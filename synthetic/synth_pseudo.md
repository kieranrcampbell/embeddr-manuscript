## Create synthetic pseudo-temporally regulated gene expression profiles
### kieran.campbell@sjc.ox.ac.uk


```r
library(ggplot2)
library(scater)
library(embeddr)
set.seed(13)
```

We're going to use the logistic function for (in)-activation of various genes


```r
logistic <- function(x, k = 1, x_0 = 0) {
  1 / (1 + exp(-k*(x - x_0)))
}

n_cells <- 100
n_genes <- 200
pst <- runif(n_cells)

## example gating function with k=5, x_0 = .5 (to centre around 0.5)
ggplot(data.frame(x=pst), aes(x=x)) + 
  stat_function(fun=logistic, args = c(k=5, x_0 = 1)) +
  theme_bw()
```

![plot of chunk log-func](figure/log-func-1.png) 

Now we set up the pseudotemporally regulated cells. We have 70 poisson noise genes, 70
house keeping genes (roughly constant expression) and 30 genes that exhibit up-regulation
and 30 down-regulated.

For each gene there's also a constant probability (0.2) of a dropout event, in which case
the expression is drawn from a poisson distribution with rate 0.1


```r
mean_shape <- 0.39
mean_rate <- 0.008
r_shape <- 0.69
r_rate <- 1.31
```

```r
n_noise <- 70 # 70 genes exhibit poisson noise
n_hk <- 90 # 70 genes exhibit house-keeping
n_up_reg <- 20 # 30 genes upregulated
n_down_reg <- 20 # 30 genes downregulated

sample_noise <- function(n = 1, lambda_0 = 0.1) {
  rpois(n, lambda = lambda_0)
}

sample_genes <- function(n = 1, mu, r, lambda = 0.168) {
  pi <- exp(-lambda * log(mu + 1)^2) # p of dropout
  noise_component <- sample_noise(n)
  amplified_component <- rnbinom(n = n, size = r, mu = mu)
  is_amplified <- sapply(pi, function(p) sample(c(0,1), size = 1, prob = c(p, 1-p)))
  (1 - is_amplified) * noise_component + is_amplified * amplified_component
}



r_0  <- rgamma(n_up_reg + n_down_reg, shape = r_shape, rate = r_rate) #runif(n = n_up_reg + n_down_reg, 0.1, 5)
mu_0_up <- rgamma(n_up_reg, shape = mean_shape, rate = mean_rate)
mu_0_down <- rgamma(n_down_reg, shape = mean_shape, rate = mean_rate)

## one cell at a time, profile for all genes
pst_regulated_matrix <- sapply(pst, function(t) {
  mu_up <- mu_0_up * logistic(t, k = runif(1, 5, 20), x_0 = runif(1, 0, 1))
  mu_down <- mu_0_down * logistic(t, k = runif(1, -20, -5), x_0 = runif(1, 0, 1))
  mu_all <- c(mu_up, mu_down)
  sample_genes(length(mu_all), mu = mu_all, r = r_0)
})


r_0_hk <- rgamma(n_hk, shape = r_shape, rate = r_rate) # runif(n = n_hk, 0.1, 5)
mu_0_hk <- rgamma(n_hk, shape = mean_shape, rate = mean_rate)

hk_matrix <- sapply(1:n_cells, function(unused) {
  sample_genes(n_hk, mu = mu_0_hk, r = r_0_hk)
})

dropout_matrix <- matrix(sample_noise(n_noise * n_cells), ncol=n_cells)

expr <- rbind(pst_regulated_matrix, hk_matrix, dropout_matrix)
rownames(expr) <- paste0('gene', 1:nrow(expr))
colnames(expr) <- paste0('pst_cell', 1:n_cells)

gene_type <- c(rep('Pst up', n_up_reg),
                rep('Pst down', n_down_reg),
                rep('House keeping', n_hk),
                rep('Off', n_noise))
```

Next we add in an equivalent number of non-pseudotemporally regulated cells (think
contamination cells or similar). These essentially express all genes as housekeeping
with constant dropout rates as before.


```r
mu_non_pst <- rgamma(n_up_reg + n_down_reg, shape = mean_shape, rate = mean_rate)
r_non_pst <- rgamma(n_up_reg + n_down_reg, shape = r_shape, rate = r_rate)
```

```r
## want (n_up_reg + n_down_reg) cells, then n_hk then n_noise

dropout_matrix2 <- matrix(sample_noise(n_noise * n_cells), ncol=n_cells)

expr_2 <- sapply(1:(n_cells), function(unused) {
  sample_genes(n_up_reg + n_down_reg, mu = mu_non_pst, r = r_non_pst)
})

hk2 <- sapply(1:(n_cells), function(unused) {
  sample_genes(n_hk, mu = mu_0_hk, r = r_0_hk)
})

expr_2 <- rbind(expr_2, hk2, dropout_matrix2)

 # dropout_matrix_2 <- matrix(sample_noise(n_noise * n_cells), ncol=n_cells)


rownames(expr_2) <- paste0('gene', 1:nrow(expr))
colnames(expr_2) <- paste0('nonpst_cell', 1:n_cells)

X <- cbind(expr, expr_2) # 200 x 200 expression matrix

# gene_means <- rowMeans(X)
# gene_vars <- apply(X, 1, var)
# CV2 <- gene_vars / gene_means^2
# CV2_plt <- qplot(gene_means, CV2) + theme_bw() + xlab('Mean') + ylab('CV2')
```

### Calling embeddr
First we convert our count data into an `SCESet`:


```r
pd <- new('AnnotatedDataFrame', 
          data = data.frame(cell_type = c(rep('pst', n_cells), rep('non-pst', n_cells))))
rownames(pd) <- colnames(X)

sce <- newSCESet(exprsData = log10(X + 1), phenoData = pd)
```

```
## Defining 'is_exprs' using exprsData and a lower exprs threshold of 0
```

Next we reduce the dimension using laplacian eigenmaps and 10 nearest neighbours. There appears
to be good separation between the different cell types.


```r
sce  <- embeddr(sce, nn = 20)
plot_embedding(sce, color_by = 'cell_type')
```

![plot of chunk embeddr-1](figure/embeddr-1-1.png) 

```r
plot_graph(sce)
```

![plot of chunk embeddr-1](figure/embeddr-1-2.png) 

Next we cluster the embedding into the different types and check the precision and recall:


```r
sce <- cluster_embedding(sce, k = 2, method='mm') #, method='kmeans')
plot_embedding(sce)
```

![plot of chunk prec-recall](figure/prec-recall-1.png) 

```r
print(table(sce$cell_type, sce$cluster))
```

```
##          
##             1   2
##   non-pst   1  99
##   pst     100   0
```

Finally we can fit the pseudotime to cluster 1 only and check the correlation with
the synthetic pseudotime


```r
sce <- fit_pseudotime(sce, clusters = 1)
plot_embedding(sce)
```

```
## Warning: Removed 99 rows containing missing values (geom_path).
```

![plot of chunk pseudotime](figure/pseudotime-1.png) 

```r
fitted_pst <- pseudotime(sce)[1:n_cells]
qplot(pst, fitted_pst) + theme_bw() +
  ggtitle('Correlation between synthetic pseudotime and fitted')
```

![plot of chunk pseudotime](figure/pseudotime-2.png) 

```r
print(abs(cor(pst, fitted_pst, method='pearson', use = 'na.or.complete')))
```

```
## [1] 0.6802311
```

Next we want to compare the result to monocle


```r
library(monocle)
cds <- toCellDataSet(sce)
save(cds, file="/net/isi-scratch/kieran/embeddr/embeddr-manuscript/synthetic/cds.Rdata")
#load(cds)

exprs(cds) <- X
cds <- setOrderingFilter(cds, featureNames(cds))
cds <- reduceDimension(cds, use_irlba = F)
```

```
## Reducing to independent components
```

```r
cds <- orderCells(cds, num_paths = 2)

cds_pst <- cds$Pseudotime[cds$cell_type == 'pst']

pdf('/net/isi-scratch/kieran/embeddr/embeddr-manuscript/synthetic/spanning_tree.pdf')
plot_spanning_tree(cds)
plot_embedding(sce)
```

```
## Warning: Removed 99 rows containing missing values (geom_path).
```

```r
qplot(pst, cds_pst)
dev.off()
```

```
## RStudioGD 
##         2
```

```r
print(abs(cor(pst, cds_pst, method='pearson')))
```

```
## [1] 0.6815609
```

