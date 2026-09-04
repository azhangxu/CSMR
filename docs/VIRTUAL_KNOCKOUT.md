# Virtual knockout

Run one task from the repository root, for example:

```sh
Rscript --vanilla scripts/virtual_knockout.R S8b
```

In RStudio, set `TASK_ID` near the top of the same script and source it. There is no separate function file to manage.

| Task | Dataset | Target | Cell subtype |
|---|---|---|---|
| S8a | GSE181294 | HLA-DQA2 | Mono_C |
| S8b | GSE181294 | CASP8 | CD4_ET |
| S8c | GSE181294 | CASP8 | CD8_NC |
| S8d | GSE181294 | CASP8 | CD8_S100B |
| S8e | GSE181294 | CASP8 | NK |
| S8f | GSE181294 | NUCKS1 | CD8_NC |
| S8g | GSE181294 | PIP4K2A | NK |
| S8h | GSE181294 | RMDN3 | NK |
| S8i | GSE181294 | GCHFR | CD8_ET |
| S8j | GSE114725 | RPS23 | CD8_S100B |
| S8k | GSE114725 | APOBEC3A | Mono_C |
| S8l | GSE114725 | GDI2 | CD8_NC |

## Inputs

For tasks other than RPS23 (S8j), supply the actual annotated Seurat objects identified in `config/knockout_tasks.csv`. They must contain the RNA raw-count layer, the named subtype metadata column and `group`. Tumour labels are case-sensitive: Tumor for GSE181294 and tumor for GSE114725. The two representative single-cell examples do not generate these prostate/breast input objects.

RPS23 (S8j) instead reads `data/scrna/knockout_inputs/GSE114725/CD8_S100B_Tumor_count_subset_for_scTenifoldKnk.rds`, the saved gene-by-cell count matrix for tumour CD8_S100B cells. It uses this matrix unchanged; it does not repeat cell selection, HVG selection or gene filtering. Supply the matrix used by the intended analysis run.

For the 11 tasks that select genes, the `must_include` column contains the ordered, semicolon-separated gene list extracted from each source. It is the requested list; absent and zero-expression genes can be removed by the source-specific filtering rules. Do not replace these lists with one generic marker set.

## Preserved settings

All tasks use five networks, nc_nComp=3, nc_q=0.9, td_K=3, ma_nDim=2, qc=FALSE, four cores, and the source tensor-decomposition settings. Significance is `diffRegulation$p.adj < 0.05`.

Ten scripts use the `protected` branch: ceiling for the 5% low-expression cutoff (bounded by 1 and 25 cells), rescue of expressed prespecified genes, a 5,000-gene cap reserving their slots, up to 500 cells/network, and seed 123.

HLA-DQA2/Mono_C uses the `legacy` branch: floor for the 5% cutoff, no rescue of all marker genes after low-expression filtering, a top-5,000 selection with a subsequent union of still-present prespecified genes, target-gene reinsertion if needed, and min(500, floor(0.8*n_cells)) cells/network. Its source did not record a seed, so the task table leaves it blank.

RPS23 preserves seed 123 and min(500, floor(0.8*n_cells)) cells/network from the supplied run script. HVG count and `must_include` are blank for this task because its input is already selected. Upstream generation of that saved matrix is outside this entry point.

## Outputs

`results/virtual_knockout/<task>/` contains the selected count matrix, task settings, dimensions, scTenifoldKnk result object, all/significant gene CSVs, a summary and sessionInfo. FC is a perturbation-magnitude statistic, not a signed expression fold change.

