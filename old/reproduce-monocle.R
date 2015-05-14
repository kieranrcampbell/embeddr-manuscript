

# Reproduce the differential expression lineage from monocle --------------


library(monocle)
##library(devtools)
##library(reshape2)
##library(dplyr)
##library(ggplot2)
data(HSMM)

## 'genes detectable in fewer than 50 cells at or above FPKM 1 were discarded'

FPKM_threshold <- 0.1
min_cells <- 50
is_exprs <- exprs(HSMM) > FPKM_threshold
discard_gene <- rowSums(is_exprs) < min_cells

cds <- HSMM[!discard_gene,]

## 'the remaining genes were analysed for differential expression using a Tobit-family GLM.
## A minimum FPKM value of 0.1 was used as the censoring threshold...this analysis reported genes 
## that were significantly differentially expressed between groups of cells harvested on different days.
## Only genes significant at an FDR < 0.01 (after Benjamini-Hochberg) were kept for ICA.'

diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr="expression~Media")
save(diff_test_res, file="/net/isi-scratch/kieran/embeddr/embeddr/data/diff_test_rest.Rdata")

load("/net/isi-scratch/kieran/embeddr/embeddr/data/diff_test_rest.Rdata")


# select only genes with FDR < 0.01 ---------------------------------------

alpha <- 0.01

ordering_genes <- row.names (subset(diff_test_res, qval < alpha))

cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, use_irlba = FALSE)

cds <- orderCells(cds, num_paths = 2, reverse = TRUE)

save(cds, file="/net/isi-scratch/kieran/embeddr/embeddr/data/cds.Rdata")

plot_spanning_tree(cds)

## if that doesn't work we can look at marker genes

marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5", "ANPEP", "PDGFRA",
                                                        "MYOG", "TPM1", "TPM2", "MYH2", "MYH3", "NCAM1",
                                                        "CDK1", "CDK2", "CCNB1", "CCNB2", "CCND1", "CCNA")))


cds <- setOrderingFilter(cds, marker_genes)

cds <- reduceDimension(cds, use_irlba = FALSE)

cds <- orderCells(cds, num_paths = 2, reverse = TRUE)

plot_spanning_tree(cds)

