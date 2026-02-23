# --- Enrichment Analysis Runner --- #

# 1. ORA (Over-Representation Analysis)
# 2. GSEA (Gene Set Enrichment Analysis)
# 3. SPIA (Signaling Pathway Impact Analysis)

# Parameters for fold chagne, p-value and other thresholds can be set according to results from DE methods

# Set working directory to the location of the project dir (assuming this script is in the project root)
# Copy the functions to the same directory as this script for easier sourcing
setwd("/Users/alan/Desktop/Genetic_Analysis/5_Enrichment")

# Path to dataset statistics file (logFC + p-values from DE analysis)
data_path <- "dataset-statistics.txt"

# ── Analysis switches ─────────────────────────────────────────────────────────
run_ORA  <- TRUE
run_GSEA <- TRUE
run_SPIA <- TRUE

# ── 1. ORA ───────────────────────────────────────────────────────────────────
if (run_ORA) {

  source("ORA.R")

  dbs <- c("wikipathways", "reactome", "gobp")

  ORA_results <- lapply(setNames(dbs, dbs), function(db) {
    message("\n", strrep("=", 60))
    tryCatch(
      run_ora(
        dataset_path = data_path,
        pathway_db   = db,
        output_dir   = "results"
      ),
      error = function(e) {
        message("ERROR in ", db, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  message("\n", strrep("=", 60))
  message("ORA SUMMARY")
  message(strrep("=", 60))

  ora_summary <- data.frame(
    database          = dbs,
    enriched_pathways = sapply(ORA_results, function(r) {
      if (is.null(r)) NA else nrow(r$results)
    }),
    output_dir        = sapply(ORA_results, function(r) {
      if (is.null(r)) NA else r$dirs$plots
    }),
    stringsAsFactors = FALSE
  )

  print(ora_summary, row.names = FALSE)

} else {
  message("ORA skipped (run_ORA = FALSE)")
}


# ── 2. GSEA ──────────────────────────────────────────────────────────────────
if (run_GSEA) {

  source("GSEA.R")

  dbs <- c("wikipathways", "reactome", "gobp")

  GSEA_results <- lapply(setNames(dbs, dbs), function(db) {
    message("\n", strrep("=", 60))
    tryCatch(
      run_gsea(
        dataset_path = data_path,
        pathway_db   = db,
        output_dir   = "results"
      ),
      error = function(e) {
        message("ERROR in ", db, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  message("\n", strrep("=", 60))
  message("GSEA SUMMARY")
  message(strrep("=", 60))

  gsea_summary <- data.frame(
    database          = dbs,
    enriched_pathways = sapply(GSEA_results, function(r) {
      if (is.null(r)) NA else nrow(r$results)
    }),
    up_regulated      = sapply(GSEA_results, function(r) {
      if (is.null(r)) NA else sum(r$results$NES > 0)
    }),
    down_regulated    = sapply(GSEA_results, function(r) {
      if (is.null(r)) NA else sum(r$results$NES < 0)
    }),
    output_dir        = sapply(GSEA_results, function(r) {
      if (is.null(r)) NA else r$dirs$plots
    }),
    stringsAsFactors = FALSE
  )

  print(gsea_summary, row.names = FALSE)

} else {
  message("GSEA skipped (run_GSEA = FALSE)")
}

# ── 3. SPIA ──────────────────────────────────────────────────────────────────
if (run_SPIA) {

  source("SPIA.R")

  SPIA_result <- tryCatch(
    run_spia(
      dataset_path = data_path,
      output_dir   = "results"
    ),
    error = function(e) {
      message("ERROR in SPIA: ", conditionMessage(e))
      NULL
    }
  )

  message("\n", strrep("=", 60))
  message("SPIA SUMMARY")
  message(strrep("=", 60))

  if (!is.null(SPIA_result)) {
    sig <- SPIA_result$results_clean
    sig <- sig[!is.na(sig$pG) & sig$pG < 0.05, ]
    spia_summary <- data.frame(
      total_pathways    = nrow(SPIA_result$results_clean),
      significant_pG    = nrow(sig),
      activated         = sum(sig$tA > 0),
      inhibited         = sum(sig$tA < 0),
      output_dir        = SPIA_result$dirs$plots,
      stringsAsFactors  = FALSE
    )
    print(spia_summary, row.names = FALSE)
  }

} else {
  message("SPIA skipped (run_SPIA = FALSE)")
}
