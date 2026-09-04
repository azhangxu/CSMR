# Hierarchical annotation example: GSE166555, four compartments -> immune -> CD8.
# Run in a fresh R session from the repository root. See docs/SINGLE_CELL.md.
library(Seurat)
library(harmony)
library(data.table)

INPUT_DIR <- file.path("data", "scrna", "GSE166555")
OUTDIR <- file.path("results", "scrna", "GSE166555")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
sample_files <- list.files(INPUT_DIR, pattern = "\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(sample_files) >= 2)
objects <- lapply(sample_files, function(filename) {
  tab <- fread(filename, data.table = FALSE)
  rownames(tab) <- tab[, 1]
  counts <- tab[, -1, drop = FALSE]
  counts[is.na(counts)] <- 0
  sample_name <- gsub("\\.tsv\\.gz$", "", basename(filename))
  obj <- CreateSeuratObject(counts, project = sample_name, min.cells = 5, min.features = 300)
  obj$group <- if (grepl("_p[0-9]+n", filename, ignore.case = TRUE)) "Normal" else
    if (grepl("_p[0-9]+t", filename, ignore.case = TRUE)) "Cancer" else "Unknown"
  obj$sample <- sample_name
  obj
})
names(objects) <- gsub("\\.tsv\\.gz$", "", basename(sample_files))
sce <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce[["percent.rp"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HB[^(P)]")
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)

set.seed(10086)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(sce), npcs = 50)
sce <- RunHarmony(sce, "orig.ident")
sce <- FindNeighbors(sce, dims = 1:15, reduction = "harmony")
sce <- FindClusters(sce, resolution = seq(0.1, 2, by = 0.1), verbose = FALSE)
sce <- SetIdent(sce, value = "RNA_snn_res.0.6")
sce <- RunUMAP(sce, dims = 1:12, reduction = "harmony",
               umap.method = "uwot", metric = "cosine")
annotations <- read.csv("config/scrna_cluster_annotations.csv", check.names = FALSE)
rows <- subset(annotations, dataset == "GSE166555" & mapping == "celltype_annotation")
labels <- setNames(rows$label, rows$cluster)
sce$celltype <- unname(labels[as.character(Idents(sce))])
stopifnot(!anyNA(sce$celltype))
saveRDS(sce, file.path(OUTDIR, "sce.all_bindname.rds"))

# Immune includes lymphoid and myeloid cells; the original variable name was TNK.
Idents(sce) <- sce$celltype
immune <- JoinLayers(subset(sce, idents = "Immune"))
immune <- CreateSeuratObject(GetAssayData(immune, assay = "RNA", layer = "counts"),
                             meta.data = immune[[]])
set.seed(10086)
immune <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)
immune <- FindVariableFeatures(immune, selection.method = "vst", nfeatures = 2000)
immune <- ScaleData(immune, features = rownames(immune))
immune <- RunPCA(immune, features = VariableFeatures(immune), verbose = FALSE)
immune <- RunHarmony(immune, group.by.vars = "orig.ident")
immune <- FindNeighbors(immune, dims = 1:15, reduction = "harmony", verbose = FALSE)
immune <- FindClusters(immune, resolution = seq(0.1, 2, by = 0.1), verbose = FALSE)
immune <- SetIdent(immune, value = "RNA_snn_res.0.6")
immune <- RunUMAP(immune, dims = 1:15, reduction = "harmony",
                  umap.method = "uwot", metric = "cosine")
rows <- subset(annotations, dataset == "GSE166555" & mapping == "TNK_annotation")
labels <- setNames(rows$label, rows$cluster)
immune$TNK_celltype <- unname(labels[as.character(Idents(immune))])
stopifnot(!anyNA(immune$TNK_celltype))
saveRDS(immune, file.path(OUTDIR, "immune_annotated.rds"))

# Representative final lineage: CD8 T cells. Its graph uses 20 dimensions, not 15.
cd8 <- JoinLayers(subset(immune, subset = TNK_celltype == "CD8+T"))
set.seed(10086)
cd8 <- NormalizeData(cd8, normalization.method = "LogNormalize", scale.factor = 10000)
cd8 <- FindVariableFeatures(cd8, selection.method = "vst", nfeatures = 2000)
cd8 <- ScaleData(cd8, features = rownames(cd8))
cd8 <- RunPCA(cd8, features = VariableFeatures(cd8), verbose = FALSE)
cd8 <- RunHarmony(cd8, group.by.vars = "orig.ident")
cd8 <- FindNeighbors(cd8, dims = 1:20, reduction = "harmony", verbose = FALSE)
cd8 <- FindClusters(cd8, resolution = seq(0.1, 2, by = 0.1), verbose = FALSE)
cd8 <- SetIdent(cd8, value = "RNA_snn_res.0.6")
cd8 <- RunUMAP(cd8, dims = 1:20, reduction = "harmony",
               umap.method = "uwot", metric = "cosine")
rows <- subset(annotations, dataset == "GSE166555" & mapping == "cd8_annotation")
labels <- setNames(rows$label, rows$cluster)
cd8$cd8_subtype <- unname(labels[as.character(Idents(cd8))])
stopifnot(!anyNA(cd8$cd8_subtype))
saveRDS(cd8, file.path(OUTDIR, "CD8T_annotated.rds"))

# CD8_NC was not assigned in this source map; do not create an artificial comparison.
thresholds <- data.frame(name = c("strict", "moderate", "relaxed", "very_relaxed"),
                         min_pct = c(0.25, 0.1, 0.05, 0.01), logfc = c(0.5, 0.25, 0.15, 0.1))
group_counts <- list()
for (subtype in c("CD8_ET", "CD8_S100B")) {
  cells <- colnames(cd8)[cd8$cd8_subtype == subtype]
  obj <- JoinLayers(subset(cd8, cells = cells))
  group_n <- table(factor(obj$group, levels = c("Cancer", "Normal")))
  group_counts[[subtype]] <- data.frame(subtype, group = names(group_n), cells = as.integer(group_n))
  if (any(group_n < 3)) {
    message(subtype, " differential expression skipped: fewer than three cells in a group.")
    next
  }
  Idents(obj) <- obj$group
  for (i in seq_len(nrow(thresholds))) {
    deg <- FindMarkers(obj, ident.1 = "Cancer", ident.2 = "Normal", test.use = "wilcox",
                        min.pct = thresholds$min_pct[i], logfc.threshold = thresholds$logfc[i])
    deg$gene <- rownames(deg)
    deg$subtype <- subtype
    write.csv(deg, file.path(OUTDIR, paste0("DEG_", subtype, "_", thresholds$name[i], ".csv")), row.names = FALSE)
  }
}
write.csv(do.call(rbind, group_counts), file.path(OUTDIR, "CD8_group_counts.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
