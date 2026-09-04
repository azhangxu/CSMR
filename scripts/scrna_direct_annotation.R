# Direct annotation example: GSE117570 (P2 tumour/normal), then B-cell subtyping.
# Run in a fresh R session from the repository root. See docs/SINGLE_CELL.md.
library(Seurat)
library(harmony)
library(data.table)

INPUT_DIR <- file.path("data", "scrna", "GSE117570")
OUTDIR <- file.path("results", "scrna", "GSE117570")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
sample_files <- c("GSM3304009_P2_Tumor_processed_data.txt.gz",
                  "GSM3304010_P2_Normal_processed_data.txt.gz")
stopifnot(all(file.exists(file.path(INPUT_DIR, sample_files))))

# The original analysis used this pair, not every sample in the GEO series.
objects <- lapply(sample_files, function(filename) {
  tab <- fread(file.path(INPUT_DIR, filename), data.table = FALSE, check.names = FALSE)
  genes <- as.character(tab[[1]])
  genes[is.na(genes) | genes == ""] <- "UnknownGene"
  counts <- as.matrix(tab[, -1, drop = FALSE])
  mode(counts) <- "numeric"
  rownames(counts) <- make.unique(genes)
  counts[is.na(counts)] <- 0
  sample_name <- gsub("_processed_data.txt.gz", "",
                      gsub("^GSM[0-9]+_", "", filename))
  obj <- CreateSeuratObject(counts, min.cells = 5, min.features = 300)
  obj$sample <- sample_name
  obj$orig.ident <- sample_name
  obj
})
names(objects) <- c("P2_Tumor", "P2_Normal")
sce <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
sce$group <- ifelse(sce$sample == "P2_Tumor", "tumor", "normal")
Idents(sce) <- "group"
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce[["percent.rp"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HB[^(P)]")
# No mitochondrial-percentage cutoff was present at this step in the source.
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

set.seed(10086)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(sce), npcs = 50)
sce <- RunHarmony(sce, "orig.ident")
sce <- FindNeighbors(sce, dims = 1:15, reduction = "harmony")
sce <- FindClusters(sce, resolution = seq(0.1, 2, by = 0.1), verbose = FALSE)
sce <- SetIdent(sce, value = "RNA_snn_res.1.2")
sce <- RunUMAP(sce, dims = 1:8, reduction = "harmony",
               umap.method = "uwot", metric = "cosine")

annotations <- read.csv("config/scrna_cluster_annotations.csv", check.names = FALSE)
rows <- subset(annotations, dataset == "GSE117570" & mapping == "celltype_annotation")
labels <- setNames(rows$label, rows$cluster)
sce$celltype <- unname(labels[as.character(Idents(sce))])
stopifnot(!anyNA(sce$celltype))
saveRDS(sce, file.path(OUTDIR, "sce.all_bindname.rds"))

# Representative downstream lineage: B cells. Other lineages have their own maps.
Idents(sce) <- sce$celltype
b <- JoinLayers(subset(sce, idents = "B cells"))
b <- CreateSeuratObject(counts = GetAssayData(b, assay = "RNA", layer = "counts"),
                        meta.data = b[[]])
set.seed(10086)
b <- NormalizeData(b, normalization.method = "LogNormalize", scale.factor = 10000)
b <- FindVariableFeatures(b, selection.method = "vst", nfeatures = 2000)
b <- ScaleData(b, features = rownames(b))
b <- RunPCA(b, features = VariableFeatures(b), verbose = FALSE)
b <- RunHarmony(b, group.by.vars = "orig.ident")
b <- FindNeighbors(b, dims = 1:25, reduction = "harmony", verbose = FALSE)
b <- FindClusters(b, resolution = seq(0.1, 2, by = 0.1), verbose = FALSE)
b <- SetIdent(b, value = "RNA_snn_res.1")
b <- RunUMAP(b, dims = 1:25, reduction = "harmony",
             umap.method = "uwot", metric = "cosine")
rows <- subset(annotations, dataset == "GSE117570" & mapping == "b_annotation")
labels <- setNames(rows$label, rows$cluster)
b$b_subtype <- unname(labels[as.character(Idents(b))])
stopifnot(!anyNA(b$b_subtype))
saveRDS(b, file.path(OUTDIR, "B_annotated.rds"))

# The source tested Plasma only. Preserve its four thresholds; the paper uses relaxed.
thresholds <- data.frame(name = c("strict", "moderate", "relaxed", "very_relaxed"),
                         min_pct = c(0.25, 0.1, 0.05, 0.01), logfc = c(0.5, 0.25, 0.15, 0.1))
plasma <- JoinLayers(subset(b, subset = b_subtype == "Plasma"))
group_n <- table(factor(plasma$group, levels = c("tumor", "normal")))
write.csv(as.data.frame(group_n), file.path(OUTDIR, "Plasma_group_counts.csv"), row.names = FALSE)
if (all(group_n >= 3)) {
  Idents(plasma) <- plasma$group
  for (i in seq_len(nrow(thresholds))) {
    deg <- FindMarkers(plasma, ident.1 = "tumor", ident.2 = "normal", test.use = "wilcox",
                        min.pct = thresholds$min_pct[i], logfc.threshold = thresholds$logfc[i])
    deg$gene <- rownames(deg)
    deg$subtype <- "Plasma"
    # p_val is nominal; Seurat's p_val_adj is Bonferroni adjusted, not BH.
    write.csv(deg, file.path(OUTDIR, paste0("DEG_Plasma_", thresholds$name[i], ".csv")), row.names = FALSE)
  }
} else {
  message("Plasma differential expression skipped: fewer than three cells in a group.")
}
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
