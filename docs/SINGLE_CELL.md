# Single-cell workflows

## Representative analyses

Both scripts read real GEO input files and run QC, normalization, variable-gene selection, PCA, Harmony, clustering, manual source-derived annotation, a representative subtype analysis and tumour-versus-reference differential expression. They are executable examples, not 11 independent copies of the same workflow.

| Step | Direct example: GSE117570 | Hierarchical example: GSE166555 |
|---|---|---|
| Samples | P2 tumour and P2 normal only | Supplied `_p<number>t` and `_p<number>n` matrices |
| Object creation | min.cells 5; min.features 300 | min.cells 5; min.features 300 |
| Explicit subsequent QC | nFeature_RNA >200 and <6000 | nFeature_RNA >200 and <9000; percent.mt <25 |
| Initial graph dimensions | 1:15 | 1:15 |
| Initial selected resolution | 1.2 | 0.6 |
| Initial UMAP dimensions | 1:8 | 1:12 |
| Next stage | B cells; dimensions 1:25, resolution 1 | Immune cells; dimensions 1:15, resolution 0.6 |
| Final example lineage | Plasma | CD8; dimensions 1:20, resolution 0.6 |
| Comparison | tumor versus normal | Cancer versus Normal |

Creation and subsequent QC both apply; the >200 cutoff does not undo the earlier min.features=300 filter. Initial ScaleData retains the source default; downstream lineage stages scale all genes as in their source scripts. The full 0.1–2.0 resolution grid and recorded seeds are preserved. Graph and UMAP dimensions are different parameters.

Cluster labels are read from `config/scrna_cluster_annotations.csv`. Unknown labels remain unknown. GSE117570's B-cell source labels clusters 0–3 as Plasma and tests Plasma only. GSE166555's CD8 source assigns CD8_ET, CD8_S100B and unknown; it does not assign CD8_NC.

The examples save annotated Seurat objects and differential-expression tables. UMAP coordinates are retained within the objects; no figures are rendered or exported. Cluster labels come from the study's manual annotations. When adapting the workflows, assess cluster identity using marker evidence before reusing a numeric cluster map.

## Differential expression

Both examples export the source's four threshold sets: strict (min.pct 0.25, logfc 0.5), moderate (0.1, 0.25), relaxed (0.05, 0.15), and very_relaxed (0.01, 0.1). The supplied methods use relaxed. Comparisons require at least three cells per group.

FindMarkers uses the Wilcoxon test. Positive avg_log2FC refers to tumour/Cancer relative to the reference. `p_val` is nominal and `p_val_adj` is Seurat's Bonferroni-adjusted value; neither is replaced with BH-adjusted values. CSV exports replace the repeated Excel-workbook export blocks. No MR directional-concordance calculation is added here because its join to the MR evidence table is a separate step.

## Other datasets and special processing

- `config/scrna_datasets.csv`: all 11 GEO identifiers, 10 analysis units, initial parameters and dataset notes.
- `config/scrna_cluster_annotations.csv`: literal cluster maps and explicit cluster-selection labels extracted from the supplied analysis scripts, including downstream lineages. `source_file`, `mapping` and `mapping_line` identify the context. Maps with different contexts must not be concatenated or treated as a universal dictionary.
- `config/scrna_stage_settings.csv`: source-line references for selected resolutions, dimensions, subcluster settings and explicit QC. This is a reference inventory, not a complete machine-executable configuration.
GSE181294 directly annotates the initial cell populations and retains Healthy/Tumor while excluding AdjNormal. GSE114725 and GSE156632 need additional T/NK cleaning; simply changing the initial cluster map is insufficient.

GSE131907 reads `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` with its supplied `GSE131907_Lung_Cancer_cell_annotation.txt.gz`. The initial object uses min.cells=3 and selects Sample_Origin tLung/nLung. Source Cell_subtype labels select the lineages before reclustering; the two examples' initial de novo annotation is not a replacement for that input annotation.

The endometrial source uses tumour samples GSM5276933–GSM5276937 from GSE173682 and three reference samples from GSE179640, preprocesses the studies, then merges their counts and reconstructs study/sample/group metadata. It initially labels stromal, epithelial/cancer and immune cells, then subdivides immune cells. It should not be described as exactly the same four-compartment implementation as GSE166555.

The reference tables preserve recorded differences but do not automate every upstream cleaning operation or every branch in all datasets. The two representative scripts are the supported runnable single-cell examples in this release.
