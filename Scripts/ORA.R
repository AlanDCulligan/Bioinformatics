run_ora <- function(
  dataset_path,
  pathway_db      = c("wikipathways", "reactome", "gobp"),
  analysis_name   = NULL,
  species         = "Homo sapiens",
  pval_cutoff     = 0.05,
  fc_cutoff       = 1,
  p_adjust_method = "BH",
  p_cutoff_ora    = 0.05,
  q_cutoff_ora    = 0.2,
  show_category   = 20,
  output_dir      = "."
) {
  # ── Libraries ────────────────────────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(EnhancedVolcano)
    library(enrichplot)
    library(msigdbr)
    library(ggplot2)
  })

  # ── Resolve pathway database ─────────────────────────────────────────────────
  pathway_db <- match.arg(pathway_db)

  db_params <- switch(pathway_db,
    wikipathways = list(
      collection    = "C2",
      subcollection = "CP:WIKIPATHWAYS",
      label         = "WikiPathways"
    ),
    reactome = list(
      collection    = "C2",
      subcollection = "CP:REACTOME",
      label         = "Reactome"
    ),
    gobp = list(
      collection    = "C5",
      subcollection = "GO:BP",
      label         = "GO_BiologicalProcess"
    )
  )

  if (is.null(analysis_name)) {
    analysis_name <- paste0("ORA_", db_params$label)
  }

  # ── Output directories ───────────────────────────────────────────────────────
  plots_dir  <- file.path(output_dir, analysis_name, "plots")
  tables_dir <- file.path(output_dir, analysis_name, "tables")

  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

  message("=== ORA | ", db_params$label, " | ", analysis_name, " ===")
  message("Output root : ", file.path(output_dir, analysis_name))
  message("Plots dir   : ", plots_dir)
  message("Tables dir  : ", tables_dir)

  # ── 1. Load & clean dataset ──────────────────────────────────────────────────
  message("\nLoading dataset...")
  dataset <- read.delim(dataset_path)
  dataset <- dataset[!is.na(dataset$ENSG_ID), ]
  dataset <- dataset[!duplicated(dataset$ENSG_ID), ]
  message("Rows after filtering: ", nrow(dataset))

  # ── 2. Volcano plot ──────────────────────────────────────────────────────────
  message("Generating volcano plot...")
  vp <- EnhancedVolcano(
    dataset,
    title    = paste0("DEG - ", analysis_name),
    lab      = dataset$external_gene_id,
    labSize  = 3,
    x        = "logFC",
    y        = "adj.P.Val",
    pCutoff  = pval_cutoff,
    FCcutoff = fc_cutoff
  )

  png(file.path(plots_dir, paste0(analysis_name, "_volcano.png")),
      width = 3200, height = 2400, res = 300)
  print(vp)
  dev.off()

  pdf(file.path(plots_dir, paste0(analysis_name, "_volcano.pdf")),
      width = 12, height = 9)
  print(vp)
  dev.off()

  message("Saved: volcano  (.png + .pdf)")

  # ── 3. Define gene lists ─────────────────────────────────────────────────────
  de.genes   <- unique(dataset[
    dataset$adj.P.Val < pval_cutoff & abs(dataset$logFC) > fc_cutoff,
    "ENSG_ID"
  ])
  bkgd.genes <- unique(dataset$ENSG_ID)

  message("DEGs:       ", length(de.genes))
  message("Background: ", length(bkgd.genes))

  if (length(de.genes) == 0) stop("No DEGs found with the current thresholds.")

  # ── 4. Retrieve gene sets ────────────────────────────────────────────────────
  message("\nFetching gene sets: ", db_params$label, "...")
  wp.sets <- msigdbr(
    species       = species,
    collection    = db_params$collection,
    subcollection = db_params$subcollection
  )
  wp.t2g <- wp.sets[, c("gs_name", "ensembl_gene")]
  message("Gene sets retrieved: ", length(unique(wp.t2g$gs_name)))

  # ── 5. ORA ───────────────────────────────────────────────────────────────────
  message("\nRunning ORA...")
  ora <- clusterProfiler::enricher(
    gene          = de.genes,
    TERM2GENE     = wp.t2g,
    universe      = bkgd.genes,
    pAdjustMethod = p_adjust_method,
    pvalueCutoff  = p_cutoff_ora,
    qvalueCutoff  = q_cutoff_ora
  )

  ora.res <- as.data.frame(ora)
  message("Significantly enriched pathways: ", nrow(ora.res))

  csv_path <- file.path(tables_dir, paste0(analysis_name, "_results.csv"))
  write.csv(ora.res, file = csv_path, row.names = FALSE)
  message("Results saved to: ", csv_path)

    # ── 7. Save text summary ─────────────────────────────────────────────────────
  txt_path <- file.path(tables_dir, paste0(analysis_name, "_summary.txt"))

  sink(txt_path)
  cat("================================================================================\n")
  cat("ORA ANALYSIS SUMMARY\n")
  cat("================================================================================\n")
  cat("Analysis name   :", analysis_name, "\n")
  cat("Pathway database:", db_params$label, "\n")
  cat("Species         :", species, "\n")
  cat("Date            :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("PARAMETERS\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Dataset path    :", dataset_path, "\n")
  cat("adj.P.Val cutoff:", pval_cutoff, "\n")
  cat("|logFC| cutoff  :", fc_cutoff, "\n")
  cat("P-adjust method :", p_adjust_method, "\n")
  cat("ORA p cutoff    :", p_cutoff_ora, "\n")
  cat("ORA q cutoff    :", q_cutoff_ora, "\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("GENE LISTS\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Background genes:", length(bkgd.genes), "\n")
  cat("DEGs (input)    :", length(de.genes), "\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("RESULTS OVERVIEW\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Enriched pathways:", nrow(ora.res), "\n")

  if (nrow(ora.res) > 0) {
    cat("\n--- TOP 10 ENRICHED PATHWAYS (by p.adjust) ---\n\n")
    top10 <- head(ora.res[order(ora.res$p.adjust), ], 10)
    for (i in seq_len(nrow(top10))) {
      cat(sprintf("[%2d] %s\n", i, top10$ID[i]))
      cat(sprintf("     GeneRatio : %s\n",   top10$GeneRatio[i]))
      cat(sprintf("     BgRatio   : %s\n",   top10$BgRatio[i]))
      cat(sprintf("     p.adjust  : %.4e\n", top10$p.adjust[i]))
      cat(sprintf("     q.value   : %.4e\n", top10$qvalue[i]))
      cat(sprintf("     Count     : %d\n",   top10$Count[i]))
      cat(sprintf("     Genes     : %s\n\n", top10$geneID[i]))
    }

    cat("--- FULL RESULTS TABLE ---\n\n")
    print(ora.res[, c("ID", "GeneRatio", "BgRatio", "pvalue",
                      "p.adjust", "qvalue", "Count")],
          row.names = FALSE)
  }

  cat("\n================================================================================\n")
  cat("Output files\n")
  cat("================================================================================\n")
  cat("Tables dir:", tables_dir, "\n")
  cat("Plots dir :", plots_dir, "\n")
  cat("CSV results:", csv_path, "\n")
  cat("This file  :", txt_path, "\n")
  sink()

  message("Saved: summary txt -> ", txt_path)

  if (nrow(ora.res) == 0) {
    warning("No enriched pathways found; skipping visualisations.")
    return(invisible(list(
      ora        = ora,
      results    = ora.res,
      pathway_db = db_params$label,
      dirs       = list(plots = plots_dir, tables = tables_dir)
    )))
  }

  # ── 6. Visualisations ────────────────────────────────────────────────────────
  save_plot <- function(filename_stem, expr,
                        png_width  = 3200, png_height = 2400, png_res = 300,
                        pdf_width  = 12,   pdf_height = 9) {

    png_file <- file.path(plots_dir, paste0(analysis_name, "_", filename_stem, ".png"))
    pdf_file <- file.path(plots_dir, paste0(analysis_name, "_", filename_stem, ".pdf"))

    png(png_file, width = png_width, height = png_height, res = png_res)
    print(expr)
    dev.off()

    pdf(pdf_file, width = pdf_width, height = pdf_height)
    print(expr)
    dev.off()

    message("Saved: ", filename_stem, "  (.png + .pdf)")
  }

  message("\nGenerating plots...")

  # Suppress known harmless deprecation warnings from enrichplot/EnhancedVolcano
  withCallingHandlers(
    {
      save_plot("barplot",
        barplot(ora, showCategory = show_category) +
          ggtitle(paste0(analysis_name, " - barplot"))
      )

      save_plot("dotplot",
        dotplot(ora, showCategory = show_category) +
          ggtitle(paste0(analysis_name, " - dotplot")) +
          theme(axis.text.y = element_text(size = 6))
      )

      edo <- pairwise_termsim(ora)

      save_plot("emapplot",
        emapplot(edo)
      )

      save_plot("treeplot",
        treeplot(edo),
        png_width = 4200, pdf_width = 16
      )

      ora.readable <- setReadable(ora, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")

      save_plot("cnetplot",
        cnetplot(
          ora.readable,
          showCategory = 5,
          cex.params   = list(gene_label = 0.6)
        )
      )
    },
    warning = function(w) {
      if (grepl("size|linewidth|by = .Count", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )

  message("\nDone! All outputs written to: ", file.path(output_dir, analysis_name))


  # ── 7. Return ────────────────────────────────────────────────────────────────
  invisible(list(
    ora        = ora,
    ora_sim    = edo,
    results    = ora.res,
    de_genes   = de.genes,
    background = bkgd.genes,
    pathway_db = db_params$label,
    dirs       = list(plots = plots_dir, tables = tables_dir)
  ))
}