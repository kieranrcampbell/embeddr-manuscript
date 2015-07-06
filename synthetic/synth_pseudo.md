## Create synthetic pseudo-temporally regulated gene expression profiles
### kieran.campbell@sjc.ox.ac.uk


```r
library(ggplot2)
library(scater)
library(embeddr)
set.seed(37)
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
  stat_function(fun=logistic, args = c(k=10, x_0 = .5)) +
  theme_bw()
```

![plot of chunk log-func](figure/log-func-1.png) 

Now we set up the pseudotemporally regulated cells. We have 70 poisson noise genes, 70
house keeping genes (roughly constant expression) and 30 genes that exhibit up-regulation
and 30 down-regulated.

For each gene there's also a constant probability (0.2) of a dropout event, in which case
the expression is drawn from a poisson distribution with rate 0.1


```r
n_noise <- 70 # 70 genes exhibit poisson noise
n_hk <- 70 # 70 genes exhibit house-keeping
n_up_reg <- 30 # 30 genes upregulated
n_down_reg <- 30 # 30 genes downregulated

sample_noise <- function(n = 1, lambda_0 = 0.1) {
  rpois(n, lambda = lambda_0)
}

sample_genes <- function(n = 1, mu, r, pi = 0.2) {
  noise_component <- sample_noise(n)
  amplified_component <- rnbinom(n = n, size = r, mu = mu)
  is_amplified <- sample(c(0,1), size = n, replace = TRUE, prob=c(pi, 1-pi))
  (1 - is_amplified) * noise_component + is_amplified * amplified_component
}



r_0  <- runif(n = n_up_reg + n_down_reg, 0.1, 5)
mu_0_up <- sample(50:200, size = n_up_reg, replace = TRUE)
mu_0_down <- sample(50:200, size = n_down_reg, replace = TRUE)

## one cell at a time, profile for all genes
pst_regulated_matrix <- sapply(pst, function(t) {
  mu_up <- mu_0_up * logistic(t, k = 5, x_0 = .5)
  mu_down <- mu_0_down * logistic(t, k = -5, x_0 = .5)
  mu_all <- c(mu_up, mu_down)
  sample_genes(length(mu_all), mu = mu_all, r = r_0)
})


r_0_hk <- runif(n = n_hk, 0.1, 5)
mu_0_hk <- sample(20:100, size = n_hk, replace = TRUE)

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
expr_2 <- sapply(1:n_cells, function(unused) {
  sample_genes(n_genes, mu = mu_0_hk, r = r_0_hk)
})


rownames(expr_2) <- paste0('gene', 1:nrow(expr))
colnames(expr_2) <- paste0('nonpst_cell', 1:n_cells)

X <- cbind(expr, expr_2) # 200 x 200 expression matrix
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
sce  <- embeddr(sce, nn = 10)
plot_embedding(sce, color_by = 'cell_type')
```

![plot of chunk embeddr-1](figure/embeddr-1-1.png) 

```r
plot_graph(sce)
```

![plot of chunk embeddr-1](figure/embeddr-1-2.png) 

Next we cluster the embedding into the different types and check the precision and recall:


```r
sce <- cluster_embedding(sce, k=2, method='kmeans')
plot_embedding(sce)
```

![plot of chunk prec-recall](figure/prec-recall-1.png) 

```r
print(table(sce$cell_type, sce$cluster))
```

```
##          
##             1   2
##   non-pst  98   2
##   pst       0 100
```

Finally we can fit the pseudotime to cluster 1 only and check the correlation with
the synthetic pseudotime


```r
sce <- fit_pseudotime(sce, clusters = 2)
plot_embedding(sce)
```

```
## Warning: Removed 98 rows containing missing values (geom_path).
```

![plot of chunk pseudotime](figure/pseudotime-1.png) 

```r
fitted_pst <- pseudotime(sce)[1:n_cells]
qplot(pst, fitted_pst) + theme_bw() +
  ggtitle('Correlation between synthetic pseudotime and fitted')
```

![plot of chunk pseudotime](figure/pseudotime-2.png) 

```r
print(abs(cor(pst, fitted_pst, method='spearman')))
```

```
## [1] 0.8660666
```

