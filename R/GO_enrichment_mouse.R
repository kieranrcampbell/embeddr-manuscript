# GO enrichment in mice using MGI symbols is surprisingly difficult. Below is code
# that does this using GO-seq and biomaRt
# kieran.campbell@sjc.ox.ac.uk

library(biomaRt)
library(org.Mm.eg.db)
library(goseq)

# BP -> AT1 transition ----------------------------------------------------


gene_names <- df_gene$gene

genes_per_class <- sapply(1:n_cuts, function(i) 1 * (df_gene$class == i))
rownames(genes_per_class) <- gene_names
# all_other_genes <- setdiff(featureNames(sce_bp_at1), gene_names)
# zero_matrix <- matrix(0, nrow=length(all_other_genes), ncol=n_cuts, dimnames=list(all_other_genes))
# genes_per_class <- rbind(genes_per_class, zero_matrix)
genes_per_class <- as.data.frame(genes_per_class)

ensembl_mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
ids <- getBM(attributes=c('mgi_symbol','ensembl_gene_id','gene_biotype'), 
             filters='mgi_symbol',values=rownames(genes_per_class), mart=ensembl_mouse)
id_lengths <- getBM(attributes=c('mgi_symbol','transcript_length'), filters='mgi_symbol',
                    values=rownames(genes_per_class),
                    mart=ensembl_mouse)

## time for some sexy dplyr
by_transcript <- dplyr::group_by(id_lengths, mgi_symbol)
gene_lengths <- dplyr::summarise(by_transcript, length=sum(transcript_length))

ids <- dplyr::left_join(ids, gene_lengths)

## found_ensembl_id_for <- gene_names %in% ensembl_ids$mgi_symbol
genes_per_class <- genes_per_class[ids$mgi_symbol, ]
rownames(genes_per_class) <- ids$ensembl_gene_id
names(genes_per_class) <- 1:n_cuts

enriched_terms <- apply(genes_per_class, 2, function(gene_set) {
  pwf <- nullp(gene_set, 'mm10', 'knownGene', bias.data=ids$length, plot.fit = FALSE)
  go <- goseq(pwf, 'mm10','ensGene', test.cats='GO:BP')
  go$log_qval <- log10(p.adjust(go$over_represented_pvalue, method='BH'))
  go <- dplyr::filter(go, log_qval < log10(0.01))
  go <- dplyr::select(go, category, log_qval, term)
  names(go) <- c('Category','log10 q-value','Term')
  return(go)
})


# BP -> AT2 transition ----------------------------------------------------

gene_names <- df_gene$gene

genes_per_class <- sapply(1:n_cuts, function(i) 1 * (df_gene$class == i))
rownames(genes_per_class) <- gene_names
# all_other_genes <- setdiff(featureNames(sce_bp_at1), gene_names)
# zero_matrix <- matrix(0, nrow=length(all_other_genes), ncol=n_cuts, dimnames=list(all_other_genes))
# genes_per_class <- rbind(genes_per_class, zero_matrix)
genes_per_class <- as.data.frame(genes_per_class)

ensembl_mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
ids <- getBM(attributes=c('mgi_symbol','ensembl_gene_id','gene_biotype'), 
             filters='mgi_symbol',values=rownames(genes_per_class), mart=ensembl_mouse)
id_lengths <- getBM(attributes=c('mgi_symbol','transcript_length'), filters='mgi_symbol',
                    values=rownames(genes_per_class),
                    mart=ensembl_mouse)

## time for some sexy dplyr
by_transcript <- dplyr::group_by(id_lengths, mgi_symbol)
gene_lengths <- dplyr::summarise(by_transcript, length=sum(transcript_length))

ids <- dplyr::left_join(ids, gene_lengths)

## found_ensembl_id_for <- gene_names %in% ensembl_ids$mgi_symbol
genes_per_class <- genes_per_class[ids$mgi_symbol, ]
rownames(genes_per_class) <- ids$ensembl_gene_id
names(genes_per_class) <- 1:n_cuts

enriched_terms <- apply(genes_per_class, 2, function(gene_set) {
  pwf <- nullp(gene_set, 'mm10', 'knownGene', bias.data=ids$length, plot.fit = FALSE)
  go <- goseq(pwf, 'mm10','ensGene')#, test.cats='GO:BP')
  go$log_qval <- log10(p.adjust(go$over_represented_pvalue, method='BH'))
  go <- dplyr::filter(go, log_qval < log10(0.01))
  go <- dplyr::select(go, category, log_qval, term)
  names(go) <- c('Category','log10 q-value','Term')
  return(go)
})

#reduced <- lapply(enriched_terms, head, 4)

