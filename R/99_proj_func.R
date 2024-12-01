plot_volcano <- function(x, title, title_plot, fc_cutoff, p_cutoff){
  
  ## Make sure title input can be used in filepath
  title <- gsub("[^A-Za-z0-9_\\-]", "_", title)
  
  ## Data wrangling
  volcano_data <- tibble(gene = rownames(x),
                         log2FoldChange = x$log2FoldChange,
                         pvalue = x$pvalue,
                         padj = x$padj)
  
  ## Annotating genes 
  volcano_data <- volcano_data |> 
    dplyr::mutate(gene_type = case_when(log2FoldChange >= fc_cutoff & 
                                          padj <= 0.05 ~ "Up regulated",
                                        log2FoldChange <= -fc_cutoff & 
                                          padj <= 0.05 ~ "Down regulated",
                                        TRUE ~ "Non significant")) 
  
  cols <- c("Up regulated" = "gold", 
            "Down regulated" = "blue", 
            "Non significant" = "darkgrey")
  
  ## Plot
  p <- volcano_data |> 
    ggplot(aes(x = log2FoldChange, 
               y = -log10(padj), 
               col=gene_type)) + 
    geom_point() +
    geom_hline(yintercept = -log10(p_cutoff), 
               linetype = "dashed") + 
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
             linetype = "dashed") +
    scale_color_manual(values=cols)+
    labs(
      title= str_c("Volcano plot - ", title_plot),
      colour = "Differencial expression"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## Create the folder to save it in
  results_dir <- "../results/DEA"
  if( !dir.exists(results_dir) ){
    dir.create(path = results_dir)
  }
  
  ## Saving and printing
  ggsave(plot = p, file.path(str_c("../results/DEA/DEA-", 
                                   title, 
                                   ".png")))
  p
}

compute_GSEA <- function(ranked_genes, gene_sets) {
  
  ## Applying GSEA
  gsea_results <- GSEA(
    geneList = ranked_genes, # Ordered ranked gene list
    minGSSize = 25, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p-value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      gene_sets,
      gs_name,
      gene_symbol
    )
  )
  
  gsea_results
}

plot_most_enriched_pathways <- function(gsea_results, 
                                        title,
                                        save_title) {
  
  ## Taking the 5 most enriched pathways
  topPathways <- gsea_results |> 
    dplyr::arrange(p.adjust) |> 
    dplyr::slice(1:5)
  
  p <- ggplot(topPathways, 
         mapping = aes(x = reorder(ID, NES), 
                       y = NES)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", 
         y = "Normalized Enrichment Score (NES)", 
         title = title) +
    theme_minimal()
  
  ## Create the folder to save it in
  results_dir <- "../results/GSEA"
  if( !dir.exists(results_dir) ){
    dir.create(path = results_dir)
  }
  
  ## Save the plot
  ggsave(plot = p, file.path(str_c("../results/GSEA/", 
                                   save_title, 
                                   ".png")))
  
  print(p)
}

# save pheatmap to pdf (https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file)
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

plot_topN_heatmap <- function(res, dds, N, annotation_data_pheno, title, title_plot){
  top_genes <- res |>
    dplyr::filter(!is.na(padj) & !is.na(gene_symbol) & !(gene_symbol=="")) |> ## Removing the genes without annotation or adjusted p-value
    dplyr::arrange(padj) |>
    dplyr::distinct(gene_symbol, .keep_all = TRUE) |> # Gene symbol can be duplicated (multiple probes for the same gene, if the same gene is multiple times in the topN gene, we skip the duplicate and take another one)
    dplyr::slice_head(n = N) |>
    dplyr::pull(transcript_ID)
  
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts <- as.data.frame(norm_counts) |>
    tibble::rownames_to_column("transcript_ID") |> 
    dplyr::filter(transcript_ID %in% top_genes) |>
    dplyr::left_join(gene_symbol.df, by = "transcript_ID")
  
  heatmap_data <- norm_counts |>
    dplyr::select(-transcript_ID) |>
    dplyr::rename(gene_symbol = gene_symbol) |>
    column_to_rownames(var = "gene_symbol")
  
  heatmap <- pheatmap(heatmap_data,
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             annotation_col = annotation_data_pheno,
                             main = str_c("Heatmap of Top ", N, " Expressed Genes \n", title_plot),
                             scale = "row")
  
  # save_pheatmap_pdf(heatmap, filename = file.path(str_c("../results/DEA/heatmap-topN-genes-", title, ".pdf")))
  save_pheatmap_png(heatmap, filename = file.path(str_c("../results/DEA/heatmap-topN-genes-", title, ".png")))
  
  print(heatmap)
}







