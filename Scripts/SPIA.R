run_spia <- function(
  dataset_path,
  analysis_name   = NULL,
  organism        = "hsa",
  species         = "Homo sapiens",
  pval_cutoff     = 0.05,
  fc_cutoff       = 1,
  n_bootstrap     = 2000,
  combine_method  = "fisher",
  p_threshold     = 0.05,
  output_dir      = "."
) {
  # ── Libraries ────────────────────────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(SPIA)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(EnhancedVolcano)
    library(ggplot2)
  })

  if (is.null(analysis_name)) {
    analysis_name <- "SPIA_KEGG"
  }

  # ── Output directories ───────────────────────────────────────────────────────
  plots_dir  <- file.path(output_dir, analysis_name, "plots")
  tables_dir <- file.path(output_dir, analysis_name, "tables")

  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

  message("=== SPIA | KEGG | ", analysis_name, " ===")
  message("Organism    : ", organism)
  message("Output root : ", file.path(output_dir, analysis_name))
  message("Plots dir   : ", plots_dir)
  message("Tables dir  : ", tables_dir)

  # ── 1. Load & clean dataset ──────────────────────────────────────────────────
  message("\nLoading dataset...")
  dataset <- read.delim(dataset_path)
  dataset <- dataset[!is.na(dataset$ENSG_ID), ]
  dataset <- dataset[!duplicated(dataset$ENSG_ID), ]
  message("Rows after filtering: ", nrow(dataset))

  # ── 2. Identifier mapping: Ensembl -> Entrez ─────────────────────────────────
  message("Mapping Ensembl IDs to Entrez IDs...")
  ids <- suppressMessages(
    clusterProfiler::bitr(
      dataset$ENSG_ID,
      fromType = "ENSEMBL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  )

  dataset <- merge(dataset, ids, by.x = "ENSG_ID", by.y = "ENSEMBL")
  dataset <- dataset[!duplicated(dataset$ENTREZID), ]
  message("Genes after Entrez mapping: ", nrow(dataset))

  if (nrow(dataset) == 0) stop("No genes remaining after Entrez ID mapping.")

  # ── 3. Volcano plot ──────────────────────────────────────────────────────────
  message("Generating volcano plot...")
  vp <- suppressWarnings(
    EnhancedVolcano(
      dataset,
      title    = paste0("DEG - ", analysis_name),
      lab      = dataset$external_gene_id,
      labSize  = 3,
      x        = "logFC",
      y        = "adj.P.Val",
      pCutoff  = pval_cutoff,
      FCcutoff = fc_cutoff
    )
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

  # ── 4. Build SPIA input vectors ──────────────────────────────────────────────
  de      <- dataset[dataset$adj.P.Val < pval_cutoff & abs(dataset$logFC) > fc_cutoff, ]
  deg.fc  <- setNames(de$logFC, de$ENTREZID)
  all.genes <- dataset$ENTREZID

  message("DEGs:       ", length(deg.fc))
  message("Background: ", length(all.genes))

  if (length(deg.fc) == 0) stop("No DEGs found with the current thresholds.")

  # ── 5. SPIA ──────────────────────────────────────────────────────────────────
  # Note: downloads KEGG pathway data on first run — may take a few minutes.
  message("\nRunning SPIA (nB = ", n_bootstrap, ")...")
  message("Note: KEGG data will be downloaded on first run.")

  spia.res <- tryCatch(
    spia(
      de       = deg.fc,
      all      = all.genes,
      organism = organism,
      nB       = n_bootstrap,
      combine  = combine_method,
      plots    = FALSE
    ),
    error = function(e) stop("SPIA failed: ", conditionMessage(e))
  )

  # Remove incomplete cases before saving/plotting
  spia.res.clean <- spia.res[complete.cases(spia.res), ]

  n_sig    <- sum(spia.res.clean$pG    < p_threshold, na.rm = TRUE)
  n_sig_both <- sum(spia.res.clean$pG  < p_threshold &
                    spia.res.clean$pNDE < p_threshold &
                    spia.res.clean$pPERT < p_threshold, na.rm = TRUE)
  n_act    <- sum(spia.res.clean$tA > 0 & spia.res.clean$pG < p_threshold, na.rm = TRUE)
  n_inh    <- sum(spia.res.clean$tA < 0 & spia.res.clean$pG < p_threshold, na.rm = TRUE)

  message("Significant pathways (pG < ", p_threshold, "): ", n_sig)
  message("  Activated (tA > 0): ", n_act)
  message("  Inhibited (tA < 0): ", n_inh)
  message("  Significant for both enrichment AND perturbation: ", n_sig_both)

  # ── 6. Save results ──────────────────────────────────────────────────────────
  csv_path <- file.path(tables_dir, paste0(analysis_name, "_results.csv"))
  write.csv(spia.res, file = csv_path, row.names = FALSE)
  message("Results saved to: ", csv_path)

  # ── 7. Save text summary ─────────────────────────────────────────────────────
  txt_path <- file.path(tables_dir, paste0(analysis_name, "_summary.txt"))

  write_summary <- function() {
    con <- file(txt_path, open = "wt")
    on.exit(close(con))

    cat("================================================================================\n", file = con)
    cat("SPIA ANALYSIS SUMMARY\n", file = con)
    cat("================================================================================\n", file = con)
    cat("Analysis name   :", analysis_name, "\n", file = con)
    cat("Database        : KEGG signaling pathways\n", file = con)
    cat("Organism        :", organism, "\n", file = con)
    cat("Species         :", species, "\n", file = con)
    cat("Date            :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = con)
    cat("--------------------------------------------------------------------------------\n", file = con)
    cat("PARAMETERS\n", file = con)
    cat("--------------------------------------------------------------------------------\n", file = con)
    cat("Dataset path    :", dataset_path, "\n", file = con)
    cat("adj.P.Val cutoff:", pval_cutoff, "\n", file = con)
    cat("|logFC| cutoff  :", fc_cutoff, "\n", file = con)
    cat("Bootstrap iter  :", n_bootstrap, "\n", file = con)
    cat("Combine method  :", combine_method, "\n", file = con)
    cat("Sig. threshold  :", p_threshold, "\n", file = con)
    cat("--------------------------------------------------------------------------------\n", file = con)
    cat("GENE LISTS\n", file = con)
    cat("--------------------------------------------------------------------------------\n", file = con)
    cat("Background genes:", length(all.genes), "\n", file = con)
    cat("DEGs (input)    :", length(deg.fc), "\n", file = con)
    cat("--------------------------------------------------------------------------------\n", file = con)
    cat("RESULTS OVERVIEW\n", file = con)
    cat("--------------------------------------------------------------------------------\n", file = con)
    cat("Total pathways tested          :", nrow(spia.res.clean), "\n", file = con)
    cat("Significant (pG <", p_threshold, ")       :", n_sig, "\n", file = con)
    cat("  Activated (tA > 0)           :", n_act, "\n", file = con)
    cat("  Inhibited (tA < 0)           :", n_inh, "\n", file = con)
    cat("  Both enriched AND perturbed  :", n_sig_both, "\n", file = con)

    if (n_sig > 0) {
      sig <- spia.res.clean[spia.res.clean$pG < p_threshold, ]
      sig <- sig[order(sig$pG), ]

      cat("\n--- ACTIVATED PATHWAYS (tA > 0, pG < ", p_threshold, ") ---\n\n", file = con)
      act <- sig[sig$tA > 0, ]
      if (nrow(act) > 0) {
        for (i in seq_len(nrow(act))) {
          cat(sprintf("[%2d] %s\n",            i, act$Name[i]),   file = con)
          cat(sprintf("     KEGG ID   : %s\n", act$ID[i]),        file = con)
          cat(sprintf("     tA        : %+.4f\n", act$tA[i]),     file = con)
          cat(sprintf("     pNDE      : %.4e\n",  act$pNDE[i]),   file = con)
          cat(sprintf("     pPERT     : %.4e\n",  act$pPERT[i]),  file = con)
          cat(sprintf("     pG        : %.4e\n",  act$pG[i]),     file = con)
          cat(sprintf("     pGFdr     : %.4e\n\n",act$pGFdr[i]),  file = con)
        }
      } else {
        cat("None.\n\n", file = con)
      }

      cat("--- INHIBITED PATHWAYS (tA < 0, pG < ", p_threshold, ") ---\n\n", file = con)
      inh <- sig[sig$tA < 0, ]
      if (nrow(inh) > 0) {
        for (i in seq_len(nrow(inh))) {
          cat(sprintf("[%2d] %s\n",            i, inh$Name[i]),   file = con)
          cat(sprintf("     KEGG ID   : %s\n", inh$ID[i]),        file = con)
          cat(sprintf("     tA        : %+.4f\n", inh$tA[i]),     file = con)
          cat(sprintf("     pNDE      : %.4e\n",  inh$pNDE[i]),   file = con)
          cat(sprintf("     pPERT     : %.4e\n",  inh$pPERT[i]),  file = con)
          cat(sprintf("     pG        : %.4e\n",  inh$pG[i]),     file = con)
          cat(sprintf("     pGFdr     : %.4e\n\n",inh$pGFdr[i]),  file = con)
        }
      } else {
        cat("None.\n\n", file = con)
      }

      cat("--- FULL SIGNIFICANT RESULTS TABLE ---\n\n", file = con)
      capture.output(
        print(sig[, c("Name", "ID", "tA", "pNDE", "pPERT", "pG", "pGFdr", "pGFWER")],
              row.names = FALSE),
        file = con
      )
    }

    cat("\n================================================================================\n", file = con)
    cat("Output files\n", file = con)
    cat("================================================================================\n", file = con)
    cat("Tables dir :", tables_dir, "\n", file = con)
    cat("Plots dir  :", plots_dir,  "\n", file = con)
    cat("CSV results:", csv_path,   "\n", file = con)
    cat("This file  :", txt_path,   "\n", file = con)
  }

  write_summary()
  message("Saved: summary txt -> ", txt_path)

  if (n_sig == 0) {
    warning("No significant pathways found (pG < ", p_threshold, "); skipping visualisations.")
    return(invisible(list(
      results    = spia.res,
      deg_fc     = deg.fc,
      background = all.genes,
      dirs       = list(plots = plots_dir, tables = tables_dir)
    )))
  }

  # ── 8. Visualisations ────────────────────────────────────────────────────────
  save_plot <- function(filename_stem, expr,
                        png_width  = 2400, png_height = 2400, png_res = 300,
                        pdf_width  = 10,   pdf_height = 10) {

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

  # Two-evidence plot: enrichment (pNDE) vs perturbation (pPERT)
  # Pathways in bottom-left are significant for BOTH
  png(file.path(plots_dir, paste0(analysis_name, "_evidence_plot.png")),
      width = 2400, height = 2400, res = 300)
  plotP(spia.res.clean, threshold = p_threshold)
  dev.off()

  pdf(file.path(plots_dir, paste0(analysis_name, "_evidence_plot.pdf")),
      width = 10, height = 10)
  plotP(spia.res.clean, threshold = p_threshold)
  dev.off()

  message("Saved: evidence_plot  (.png + .pdf)")

  # Bar plot of significant pathways by combined p-value (pG)
  sig <- spia.res.clean[spia.res.clean$pG < p_threshold, ]
  sig <- sig[order(sig$tA), ]
  sig$direction <- ifelse(sig$tA > 0, "Activated", "Inhibited")
  sig$Name <- factor(sig$Name, levels = sig$Name)

  bar <- ggplot(sig, aes(x = tA, y = Name, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("Activated" = "#d73027", "Inhibited" = "#4575b4")) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = paste0(analysis_name, " - Pathway perturbation (tA)"),
      x     = "Total perturbation accumulation (tA)",
      y     = NULL,
      fill  = NULL
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7))

  save_plot("barplot", bar, png_width = 3200, png_height = 2400,
            pdf_width = 14, pdf_height = 9)

  # Scatter plot: -log10(pNDE) vs -log10(pPERT), sized by |tA|
  scatter <- ggplot(spia.res.clean, aes(
    x     = -log10(pNDE),
    y     = -log10(pPERT),
    size  = abs(tA),
    color = ifelse(pG < p_threshold,
                   ifelse(tA > 0, "Activated", "Inhibited"),
                   "Not significant")
  )) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c(
      "Activated"       = "#d73027",
      "Inhibited"       = "#4575b4",
      "Not significant" = "grey70"
    )) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = -log10(p_threshold), linetype = "dashed", colour = "grey40") +
    labs(
      title  = paste0(analysis_name, " - Enrichment vs Perturbation"),
      x      = "-log10(pNDE)  [enrichment]",
      y      = "-log10(pPERT) [perturbation]",
      color  = NULL,
      size   = "|tA|"
    ) +
    theme_bw()

  save_plot("scatter", scatter, png_width = 3200, png_height = 2400,
            pdf_width = 12, pdf_height = 9)

  message("\nDone! All outputs written to: ", file.path(output_dir, analysis_name))

  # ── 9. Return ────────────────────────────────────────────────────────────────
  invisible(list(
    results    = spia.res,
    results_clean = spia.res.clean,
    deg_fc     = deg.fc,
    background = all.genes,
    dirs       = list(plots = plots_dir, tables = tables_dir)
  ))
}