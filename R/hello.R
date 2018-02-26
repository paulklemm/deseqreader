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
#' @param filter_invalid Remove entry if log2Foldchange, pvalue and padj is NA
#' @export
#' @import readr magrittr
#' @return Tibble of DESeq data
read_deseq <- function(path, filter_invalid = TRUE) {
  read_tsv(path) %>%
    rename(EnsemblGeneID = id) %>%
    filter(!is.na(log2FoldChange) & !is.na(pvalue) & !is.na(padj)) %>%
    return()
}

#' Attach biomart variables to deseq data
#'
#' @import dplyr biomaRt
#' @export
#' @param dat Tibble to attach biomart results to.
#' @param join_key_dat Join key on the dataset
#' @param join_key_bm Join key on biomart
#' @param bmart_dataset biomart dataset to be used. Defaults to mmusculus_gene_ensembl
#' @param extra_attributes Any extra attributes to attach besides the gene name
#' @return Tibble with attached biomart variables
attach_biomart <- function(
    dat,
    join_key_dat = "EnsemblGeneID",
    join_key_bm = "ensembl_gene_id",
    bmart_dataset = "mmusculus_gene_ensembl",
    extra_attributes = c()
  ) {
  dat <- getBM(
    filters = join_key_bm,
    attributes = c(join_key_bm, "external_gene_name", extra_attributes),
    values = dat[join_key_dat][[1]],
    mart = useMart("ensembl", dataset=bmart_dataset),
    uniqueRows = TRUE,
    verbose = FALSE
  ) %>%
    inner_join(dat, ., by = setNames(join_key_bm, join_key_dat)) %>%
    rename("GeneName" = "external_gene_name")
}

#' Add boolean 'Significant' variable to deseq dataset
#'
#' @param deseq_dat Tibble of DESeq data
#' @param padj_min Double of minimum adjusted pvalue
#' @export
#' @import dplyr
#' @return Tibble of entries satisfying significance filter
set_significant <- function(deseq_dat, padj_min = 0.05) {
  deseq_dat %>%
    mutate(Significant = padj <= padj_min) %>%
    return()
}

#' Export deseq folder to excel with default read_deseq
#'
#' @param path Path to deseq files
#' @param outputpath Path to output Excel files
#' @param file_pattern Pattern to identify deseq files.
#' @export
#' @import writexl purrr stringr
deseq_to_excel <- function(path, outputpath, file_pattern = ".tsv$") {
  # Remove ending forward slah if it's there
  path <- path %>% str_replace(pattern = "/$", replacement = "")
  outputpath <- outputpath %>% str_replace(pattern = "/$", replacement = "")
  # Read in all files from the path
  dir(path, pattern = file_pattern) %>%
    # Walk over all path names
    iwalk(
      # Combine path with filename
      ~str_glue("{path}/{.x}") %>%
        read_deseq() %>%
        # Write excel file (replace the file ending)
        write_xlsx(path = paste0(outputpath, "/", str_replace(string = .x, pattern = "\\.[^.]*$", replacement = ".xlsx")))
    )
}

#' Create a volcano plot from DESeq data set
#'
#' @param deseq_dat DESeq dataset to plot
#' @param label_var Name of variable containing the names for plotting
#' @param label_top_n Annotate the top n genes of the data set
#' @export
#' @import ggplot2 ggrepel magrittr
#' @return ggplot2 object of volcano plot
volcano_plot <- function(deseq_dat, label_var = "GeneName", label_top_n = 20) {
  if (!(label_var %in% names(deseq_dat))) {
    stop(paste0(label_var, "not in data frame. Please attach it with attach_biomart"))
  }
  if (!("Significant") %in% names(deseq_dat)) {
    warning("Variable 'Significant' not found, added it using 'set_significant' with default values (padj <= 0.05)")
    deseq_dat <- deseq_dat %>% set_significant()
  }
  plot <- deseq_dat %>%
    ggplot(aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(
      data = deseq_dat %>% filter(Significant == TRUE) %>% arrange(padj) %>% head(label_top_n),
      aes_string(label = label_var)
    )
  return(plot)
}
