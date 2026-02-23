run_gsea <- function(
  dataset_path,
  pathway_db      = c("wikipathways", "reactome", "gobp"),
  analysis_name   = NULL,
  species         = "Homo sapiens",
  ranking_metric  = c("signed_p", "logfc"),
  pval_cutoff     = 0.05,
  fc_cutoff       = 1,
  p_adjust_method = "BH",
  p_cutoff_gsea   = 0.05,
  show_category   = 20,
  top_pathways    = 3,
  seed            = 42,
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
    library(ggridges)
  })

  # ── Resolve pathway database ─────────────────────────────────────────────────
  pathway_db     <- match.arg(pathway_db)
  ranking_metric <- match.arg(ranking_metric)

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
    analysis_name <- paste0("GSEA_", db_params$label)
  }

  # ── Output directories ───────────────────────────────────────────────────────
  plots_dir  <- file.path(output_dir, analysis_name, "plots")
  tables_dir <- file.path(output_dir, analysis_name, "tables")

  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

  message("=== GSEA | ", db_params$label, " | ", analysis_name, " ===")
  message("Ranking metric : ", ranking_metric)
  message("Output root    : ", file.path(output_dir, analysis_name))
  message("Plots dir      : ", plots_dir)
  message("Tables dir     : ", tables_dir)

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

  # ── 3. Build ranked gene list ────────────────────────────────────────────────
  message("\nBuilding ranked gene list using: ", ranking_metric)

  dataset$rank <- switch(ranking_metric,
    signed_p = sign(dataset$logFC) * -log10(dataset$adj.P.Val),
    logfc    = dataset$logFC
  )

  genes.ranked <- setNames(dataset$rank, dataset$ENSG_ID)
  genes.ranked <- sort(genes.ranked, decreasing = TRUE)

  message("Ranked genes: ", length(genes.ranked))
  message("Top 3:    ", paste(head(names(genes.ranked), 3), collapse = ", "))
  message("Bottom 3: ", paste(tail(names(genes.ranked), 3), collapse = ", "))

  # ── 4. Retrieve gene sets ────────────────────────────────────────────────────
  message("\nFetching gene sets: ", db_params$label, "...")
  wp.sets <- msigdbr(
    species       = species,
    collection    = db_params$collection,
    subcollection = db_params$subcollection
  )
  wp.t2g <- wp.sets[, c("gs_name", "ensembl_gene")]
  message("Gene sets retrieved: ", length(unique(wp.t2g$gs_name)))

  # ── 5. GSEA ──────────────────────────────────────────────────────────────────
  # Note: two warnings you may see - both can be ignored:
  # 1. "There are ties in the preranked stats" - harmless, reduced by signed_p metric
  # 2. "'package:stats' may not be available when loading" - harmless parallelisation warning
  message("\nRunning GSEA (seed = ", seed, ")...")
  gsea <- suppressWarnings(
    clusterProfiler::GSEA(
      geneList      = genes.ranked,
      TERM2GENE     = wp.t2g,
      pAdjustMethod = p_adjust_method,
      pvalueCutoff  = p_cutoff_gsea,
      seed          = seed
    )
  )

  gsea.res <- as.data.frame(gsea)
  n_up   <- sum(gsea.res$NES > 0)
  n_down <- sum(gsea.res$NES < 0)
  message("Significantly enriched pathways: ", nrow(gsea.res),
          " (", n_up, " up / ", n_down, " down)")

  csv_path <- file.path(tables_dir, paste0(analysis_name, "_results.csv"))
  write.csv(gsea.res, file = csv_path, row.names = FALSE)
  message("Results saved to: ", csv_path)

  if (nrow(gsea.res) == 0) {
    warning("No enriched pathways found; skipping visualisations.")
    return(invisible(list(
      gsea       = gsea,
      results    = gsea.res,
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

  withCallingHandlers(
    {
      # Dotplot split by direction (up/down)
      save_plot("dotplot",
        dotplot(gsea, showCategory = show_category, split = ".sign") +
          facet_grid(. ~ .sign) +
          ggtitle(paste0(analysis_name, " - dotplot")) +
          theme(axis.text.y = element_text(size = 6))
      )

      edo <- pairwise_termsim(gsea)

      save_plot("emapplot",
        emapplot(edo)
      )

      save_plot("treeplot",
        treeplot(edo),
        png_width = 4200, pdf_width = 16
      )

      # Ridgeplot: expression distribution of core enriched genes per pathway
      save_plot("ridgeplot",
        ridgeplot(gsea) +
          theme(axis.text.y = element_text(size = 6)),
        png_width  = 2400, png_height = 3200,
        pdf_width  = 9,    pdf_height = 12
      )

      # Gene-concept network: version-safe cnetplot with fold change colouring
      gsea.readable <- setReadable(gsea, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")

      cnet <- tryCatch(
        cnetplot(gsea.readable, showCategory = 5, foldChange = genes.ranked,
                 cex.params = list(gene_label = 0.6)),
        error = function(e) tryCatch(
          cnetplot(gsea.readable, showCategory = 5, foldChange = genes.ranked,
                   cex_label_gene = 0.6),
          error = function(e)
            cnetplot(gsea.readable, showCategory = 5, foldChange = genes.ranked)
        )
      )
      save_plot("cnetplot", cnet)

      # Classic GSEA running enrichment score plot for top N pathways
      n_plot <- min(top_pathways, nrow(gsea.res))
      save_plot("gseaplot",
        gseaplot2(
          gsea,
          geneSetID = 1:n_plot,
          title     = paste0("Top ", n_plot, " enriched pathways")
        )
      )
    },
    warning = function(w) {
      if (grepl("size|linewidth|ties|by = .Count|package:stats", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )

  message("\nDone! All outputs written to: ", file.path(output_dir, analysis_name))

    txt_path <- file.path(tables_dir, paste0(analysis_name, "_summary.txt"))

  sink(txt_path)
  cat("================================================================================\n")
  cat("GSEA ANALYSIS SUMMARY\n")
  cat("================================================================================\n")
  cat("Analysis name   :", analysis_name, "\n")
  cat("Pathway database:", db_params$label, "\n")
  cat("Species         :", species, "\n")
  cat("Date            :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("PARAMETERS\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Dataset path    :", dataset_path, "\n")
  cat("Ranking metric  :", ranking_metric, "\n")
  cat("adj.P.Val cutoff:", pval_cutoff, "\n")
  cat("|logFC| cutoff  :", fc_cutoff, "(volcano plot only)\n")
  cat("P-adjust method :", p_adjust_method, "\n")
  cat("GSEA p cutoff   :", p_cutoff_gsea, "\n")
  cat("Seed            :", seed, "\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("GENE LIST\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Total ranked genes:", length(genes.ranked), "\n")
  cat("Top gene (rank)   :", names(genes.ranked)[1],
      sprintf("(%.4f)\n", genes.ranked[1]))
  cat("Bottom gene (rank):", names(genes.ranked)[length(genes.ranked)],
      sprintf("(%.4f)\n", genes.ranked[length(genes.ranked)]))
  cat("--------------------------------------------------------------------------------\n")
  cat("RESULTS OVERVIEW\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Enriched pathways:", nrow(gsea.res), "\n")
  cat("  Up-regulated (NES > 0)  :", n_up, "\n")
  cat("  Down-regulated (NES < 0):", n_down, "\n")

  if (nrow(gsea.res) > 0) {
    cat("\n--- TOP 10 UP-REGULATED PATHWAYS (NES > 0, by p.adjust) ---\n\n")
    up_res <- gsea.res[gsea.res$NES > 0, ]
    up_res <- head(up_res[order(up_res$p.adjust), ], 10)
    if (nrow(up_res) > 0) {
      for (i in seq_len(nrow(up_res))) {
        cat(sprintf("[%2d] %s\n", i, up_res$ID[i]))
        cat(sprintf("     NES       : %+.4f\n", up_res$NES[i]))
        cat(sprintf("     p.adjust  : %.4e\n",  up_res$p.adjust[i]))
        cat(sprintf("     q.value   : %.4e\n",  up_res$qvalue[i]))
        cat(sprintf("     Set size  : %d\n",    up_res$setSize[i]))
        cat(sprintf("     Core genes: %s\n\n",  up_res$core_enrichment[i]))
      }
    } else {
      cat("None.\n\n")
    }

    cat("--- TOP 10 DOWN-REGULATED PATHWAYS (NES < 0, by p.adjust) ---\n\n")
    dn_res <- gsea.res[gsea.res$NES < 0, ]
    dn_res <- head(dn_res[order(dn_res$p.adjust), ], 10)
    if (nrow(dn_res) > 0) {
      for (i in seq_len(nrow(dn_res))) {
        cat(sprintf("[%2d] %s\n", i, dn_res$ID[i]))
        cat(sprintf("     NES       : %+.4f\n", dn_res$NES[i]))
        cat(sprintf("     p.adjust  : %.4e\n",  dn_res$p.adjust[i]))
        cat(sprintf("     q.value   : %.4e\n",  dn_res$qvalue[i]))
        cat(sprintf("     Set size  : %d\n",    dn_res$setSize[i]))
        cat(sprintf("     Core genes: %s\n\n",  dn_res$core_enrichment[i]))
      }
    } else {
      cat("None.\n\n")
    }

    cat("--- FULL RESULTS TABLE ---\n\n")
    print(gsea.res[order(gsea.res$p.adjust),
                   c("ID", "setSize", "NES", "pvalue", "p.adjust", "qvalue")],
          row.names = FALSE)
  }

  cat("\n================================================================================\n")
  cat("Output files\n")
  cat("================================================================================\n")
  cat("Tables dir :", tables_dir, "\n")
  cat("Plots dir  :", plots_dir, "\n")
  cat("CSV results:", csv_path, "\n")
  cat("This file  :", txt_path, "\n")
  sink()

  message("Saved: summary txt -> ", txt_path)

  # ── 7. Return ────────────────────────────────────────────────────────────────
  invisible(list(
    gsea         = gsea,
    gsea_sim     = edo,
    results      = gsea.res,
    genes_ranked = genes.ranked,
    pathway_db   = db_params$label,
    dirs         = list(plots = plots_dir, tables = tables_dir)
  ))
}
