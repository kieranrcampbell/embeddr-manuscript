
# Create synthetic pseudo-temporally regulated gene expression profiles
# Author: kieran.campbell@sjc.ox.ac.uk

library(ggplot2)

logistic <- function(x, k = 1, x_0 = 0) {
  1 / (1 + exp(-k*(x - x_0)))
}

n_cells <- 100
n_genes <- 200
pst <- runif(n_cells)

## example gating function with k=5, x_0 = .5 (to centre around 0.5)
ggplot(data.frame(x=pst), aes(x=x)) + stat_function(fun=logistic, args = c(k=-5, x_0 = .5))

n_noise <- 70 # 70 genes exhibit poisson noise
n_hk <- 70 # 70 genes exhibit house-keeping
n_up_reg <- 30 # 30 genes upregulated
n_down_reg <- 30 # 30 genes downregulated

sample_noise <- function(n = 1, lambda_0 = 0.1) {
  rpois(n, lambda = lambda_0)
}

sample_genes <- function(n = 1, mu, r, pi = 0.1) {
  noise_component <- sample_noise(n)
  amplified_component <- rnbinom(n = n, size = r, mu = mu)
  is_amplified <- sample(c(0,1), size = n, replace = TRUE, prob=c(pi, 1-pi))
  (1 - is_amplified) * noise_component + is_amplified * amplified_component
}



r_0  <- runif(n = n_up_reg + n_down_reg, 10, 20)
mu_0_up <- sample(50:100, size = n_up_reg, replace = TRUE)
mu_0_down <- sample(50:100, size = n_down_reg, replace = TRUE)

## one cell at a time, profile for all genes
pst_regulated_matrix <- sapply(pst, function(t) {
  mu_up <- mu_0_up * logistic(t, k = 5, x_0 = .5)
  mu_down <- mu_0_down * logistic(t, k = -5, x_0 = .5)
  mu_all <- c(mu_up, mu_down)
  sample_genes(length(mu_all), mu = mu_all, r = r_0)
})


r_0_hk <- runif(n = n_hk, 50, 100)
mu_0_hk <- sample(20:50, size = n_hk, replace = TRUE)

hk_matrix <- sapply(1:n_cells, function(unused) {
  sample_genes(n_hk, mu = mu_0_hk, r = r_0_hk)
})

dropout_matrix <- matrix(sample_noise(n_noise * n_cells), ncol=n_cells)

expr <- rbind(pst_regulated_matrix, hk_matrix, dropout_matrix)
rownames(expr) <- paste0('gene', 1:nrow(expr))
colnames(expr) <- paste0('cell', 1:ncol(expr))

gene_type <- c(rep('Pst up', n_up_reg),
                rep('Pst down', n_down_reg),
                rep('House keeping', n_hk),
                rep('Off', n_noise))
# Call in embeddr ---------------------------------------------------------

library(devtools)
library(scater)
load_all("~/oxford/embeddr//embeddr")

x <- t(log10(expr + 1))
x_mean <- colMeans(x)
x_var <- apply(x, 2, var)
genes_for_fit <- x_mean > 0.5
CV2 <- x_var[genes_for_fit] / (x_mean[genes_for_fit])^2
df_fit <- data.frame(m = x_mean[genes_for_fit], CV2 = CV2)
fit_loglin <- nls(CV2 ~ a * 10^(-k * m), data = df_fit, start=c(a=5, k=1))
ak <- coef(fit_loglin)
f <- function(x) ak[1] * 10^(-ak[2] * x)
genes_for_embedding <- (CV2 >  predict(fit_loglin))
df_fit$for_embedding <- as.factor(genes_for_embedding)
df_fit$gene_type <- as.factor(gene_type[genes_for_fit])
ggplot(df_fit, aes(x=m, y=CV2, color = gene_type)) + geom_point() +
  theme_bw() + xlab('Mean') + ylab('CV2') +
  stat_function(fun=f, color='black')


sce <- newSCESet(cellData = log10(expr + 1))

W <- weighted_graph(sce, metric='correlation')
sce <- laplacian_eigenmap(sce, W)
plot_embedding(sce)

## pseudotime doesn't work so well here
x <- pData(sce)$component_1
y <- pData(sce)$component_2
fit <- lm(y ~ poly(x, 2))
df <- data.frame(x, predict(fit))
df <- dplyr::arrange(df, x)

sce <- fit_pseudotime(sce, start = as.matrix(df))
plot_embedding(sce)
qplot(pst, pseudotime(sce))
cor(pst, pseudotime(sce))

# Compare to monocle ------------------------------------------------------

rownames(expr) <- paste0('gene', 1:nrow(expr))
colnames(expr) <- paste0('cell', 1:ncol(expr))

library(monocle)








