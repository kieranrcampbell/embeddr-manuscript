# 
# This script estimates the parameters from negative binomial distributions
# for 100 randomly selected genes (that are expressed in more than 50 cells
# at FPKM > 0). These are then used as the estimated parameters for the
# synthetic data.
# 
# kieran.campbell@sjc.ox.ac.uk


library(monocle)
library(MASS)

set.seed(12)

data(HSMM_expr_matrix)

n_cells_exprs <- rowSums(HSMM_expr_matrix > 0)
use_genes <- which((n_cells_exprs > 0.2 * ncol(HSMM_expr_matrix)) & rowMeans(HSMM_expr_matrix > 2))

X <- HSMM_expr_matrix[use_genes,]
x_nonzero <- rowSums(X > 0)

X100 <- X[sample(1:nrow(X), 400),]
X100 <- round(X100)
estimates <- apply(X100, 1, function(x) {
  fd <- fitdistr(x, densfun = 'negative binomial', start = list(mu = mean(x), size=1))
  fd$estimate
  })

means <- estimates[1,]
dispersions <- estimates[2,]

m_fit <- fitdistr(means, densfun = 'gamma')
print(m_fit)

d_fit <- fitdistr(dispersions, densfun = 'gamma')
print(d_fit)



# now infer dropout rate --------------------------------------------------

Y <- HSMM_expr_matrix[n_cells_exprs > 0.05 * ncol(HSMM_expr_matrix),]
pd <- rowSums(Y == 0) / ncol(Y)
md <- apply(Y, 1, function(x) mean(log(x[x>0] + 1)))
mean_to_12 <- md < 12

df <- data.frame(mean=md[mean_to_12], p0=pd[mean_to_12])
ggplot(df) + geom_point(aes(x=mean, y=p0))
qplot(-md^2, pd)
fit <- lm(pd ~ md^2)
lambda <- -coef(fit)[2]
fn <- function(x, lambda) exp(-lambda * x^2)
ggplot(df) + geom_point(aes(x=mean, y=p0)) + 
  stat_function(fun = fn, args = list(lambda = lambda), color='red') 




