# 
# trimmedMean <- function(x, pct = 5) {
#   cutoffs <- quantile(x, probs = c(pct / 100, 1 - pct / 100))
#   to_use <- x > cutoffs[1] & x < cutoffs[2]
#   if(length(to_use) == 0) warning('No genes in central quantile')
#   mean(x[to_use])
# }
# 
# trimmedMeanNormalise <- function(sce, pct = 5) {
#   X <- exprs(sce)
#   normalise_tm <- function(x, pct) {
#     x / trimmedMean(x, pct = pct)
#   }
#   exprs(sce) <- apply(X, 1, normalise_tm, pct)
#   sce
# }