# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

#' Read in DESeq data from TSV file
#'
#' @param path Path to DESeq data
#' @export
#' @import readr magrittr
#' @return tibble of DESeq data
read_deseq <- function(path) {
  read_tsv(path) %>%
    return()
}

#' @import dplyr biomaRt
attach_biomart <- function(dat, bmart_dataset = "mmusculus_gene_ensembl", extra_attributes = c()) {
  dat <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "external_gene_name", extra_attributes),
    values = dat$id,
    mart = useMart("ensembl", dataset=bmart_dataset),
    uniqueRows = TRUE,
    verbose = FALSE
  ) %>%
    inner_join(dat, by = c("ensembl_gene_id" = "id")) %>%
    rename(EnsemblID = ensembl_gene_id, GeneName = external_gene_name)
}

#' Create a volcano plot from DESeq data set
#'
#' @param deseq_dat DESeq dataset to plot
#' @param label_top_n Annotate the top n genes of the data set
#' @export
#' @import ggplot2 ggrepel
#' @return ggplot2 object of volcano plot
volcano_plot <- function(deseq_dat, label_top_n = 20) {
  if (!("GeneName" %in% names(deseq_dat))) {
    stop("GeneName not in data frame. Please attach it with TODO")
  }
  plot <- deseq_dat %>%
    ggplot(aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(
      data = dat %>% filter(Significant == TRUE) %>% arrange(padj) %>% head(label_top_n),
      aes(label = Name)
    )
  return(plot)
}
