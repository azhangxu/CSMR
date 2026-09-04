# scTenifoldKnk for one of the 12 study tasks. Run from the repository root:
# Rscript --vanilla scripts/virtual_knockout.R S8b
# In RStudio, edit TASK_ID below, then source this file.
# S8j reads its prepared count matrix and uses the supplied RPS23 run settings.
TASK_ID <- "S8b"
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) TASK_ID <- args[[1]]

library(Seurat)
library(Matrix)
library(scTenifoldKnk)

# These two filtering branches preserve differences in the supplied source scripts.
prepare_knockout_counts <- function(counts, hvg, task) {
  target <- task$target_gene
  required <- strsplit(task$must_include, ";", fixed = TRUE)[[1]]
  protected <- task$filter_mode == "protected"
  if (protected) required <- unique(c(target, required))
  required <- intersect(required, rownames(counts))
  high_expression <- rownames(counts)[
    Matrix::rowMeans(counts) > 0.1 &
      Matrix::rowSums(counts > 0) / ncol(counts) * 100 > 5
  ]
  genes <- unique(c(intersect(hvg, rownames(counts)), high_expression, required))
  x <- counts[genes, , drop = FALSE]
  expressed_in <- Matrix::rowSums(x > 0)
  if (protected) {
    threshold <- max(1L, min(25L, ceiling(ncol(x) * 0.05)))
    required_expressed <- intersect(required, names(expressed_in)[expressed_in > 0])
    keep <- expressed_in >= threshold | rownames(x) %in% required_expressed
    x <- x[keep, , drop = FALSE]
    if (nrow(x) > 5000) {
      required_top <- intersect(required_expressed, rownames(x))
      stopifnot(length(required_top) <= 5000)
      other <- setdiff(rownames(x), required_top)
      means <- Matrix::rowMeans(x[other, , drop = FALSE])
      selected <- names(sort(means, decreasing = TRUE))[
        seq_len(min(5000L - length(required_top), length(other)))
      ]
      x <- x[unique(c(required_top, selected)), , drop = FALSE]
    }
    stopifnot(target %in% rownames(x), sum(x[target, ] > 0) > 0)
  } else {
    # The legacy HLA-DQA2 source uses floor and does not rescue all marker genes.
    threshold <- min(25, floor(ncol(x) * 0.05))
    x <- x[expressed_in >= threshold, , drop = FALSE]
    if (nrow(x) > 5000) {
      selected <- names(sort(Matrix::rowMeans(x), decreasing = TRUE))[1:5000]
      selected <- unique(c(selected, target, required))
      x <- x[intersect(selected, rownames(x)), , drop = FALSE]
    }
    if (target %in% rownames(counts) && !target %in% rownames(x)) {
      x <- rbind(x, counts[target, , drop = FALSE])
    }
  }
  stopifnot(target %in% rownames(x), nrow(x) >= 3, ncol(x) >= 2)
  x
}

tasks <- read.csv("config/knockout_tasks.csv", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = "")
task <- tasks[tasks$task_id == TASK_ID, , drop = FALSE]
stopifnot(nrow(task) == 1, file.exists(task$input_file))
OUTDIR <- file.path("results", "virtual_knockout", TASK_ID)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (task$filter_mode == "precomputed") {
  # RPS23 starts from the selected count matrix used by the supplied run script.
  counts <- readRDS(task$input_file)
  stopifnot(is.matrix(counts) || inherits(counts, "Matrix"),
            task$target_gene %in% rownames(counts), nrow(counts) >= 3, ncol(counts) >= 2)
} else {
  obj <- readRDS(task$input_file)
  stopifnot(inherits(obj, "Seurat"))
  meta <- obj[[]]
  stopifnot(all(c(task$subtype_column, "group") %in% names(meta)))
  keep <- !is.na(meta[[task$subtype_column]]) & !is.na(meta$group) &
    meta[[task$subtype_column]] == task$subtype & meta$group == task$tumor_label
  stopifnot(sum(keep) >= 2)
  obj <- subset(obj, cells = rownames(meta)[keep])
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  stopifnot(task$target_gene %in% rownames(counts))
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  obj <- FindVariableFeatures(obj, assay = "RNA", nfeatures = task$hvg_n, verbose = FALSE)
  counts <- prepare_knockout_counts(counts, VariableFeatures(obj), task)
}
saveRDS(counts, file.path(OUTDIR, "count_subset.rds"))

# Use each recorded seed; HLA-DQA2 has no seed recorded in its source.
if (!is.na(task$seed)) set.seed(task$seed)
network_cells <- if (task$network_cells_mode == "floor_80_percent") {
  min(500, floor(ncol(counts) * 0.8))
} else {
  min(500L, ncol(counts))
}
write.csv(task, file.path(OUTDIR, "task_settings.csv"), row.names = FALSE)
write.csv(data.frame(genes = nrow(counts), cells = ncol(counts),
                     cells_per_network = network_cells),
          file.path(OUTDIR, "input_dimensions.csv"), row.names = FALSE)
result <- scTenifoldKnk(
  countMatrix = counts, qc = FALSE, gKO = task$target_gene,
  nc_nNet = 5, nc_nCells = network_cells, nc_nComp = 3,
  nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9,
  td_K = 3, td_maxIter = 1000, td_maxError = 1e-5, td_nDecimal = 3,
  ma_nDim = 2, nCores = 4
)
saveRDS(result, file.path(OUTDIR, "scTenifoldKnk_result.rds"))
deg <- as.data.frame(result$diffRegulation)
if (!"gene" %in% names(deg)) deg$gene <- rownames(deg)
deg <- deg[, c("gene", setdiff(names(deg), "gene")), drop = FALSE]
sig <- deg[!is.na(deg$p.adj) & deg$p.adj < 0.05, , drop = FALSE]
sig <- sig[order(abs(sig$FC), decreasing = TRUE), , drop = FALSE]
write.csv(deg, file.path(OUTDIR, "KO_all_genes.csv"), row.names = FALSE)
write.csv(sig, file.path(OUTDIR, "KO_significant_genes.csv"), row.names = FALSE)
# FC describes perturbation magnitude; it is not signed expression log-fold change.
write.csv(data.frame(task_id = TASK_ID, significant_genes = nrow(sig)),
          file.path(OUTDIR, "KO_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
message("Saved task ", TASK_ID, " in ", OUTDIR)
