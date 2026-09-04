# CSMR study analysis code

Analysis code for cell-type specificity, disease specificity, MDTPS, single-cell analysis and virtual knockout, accompanying the existing [CSMR](https://github.com/azhangxu/CSMR) upstream scripts.

| Analysis | Entry point | Main output |
|---|---|---|
| Cell-type specificity | `scripts/analysis_final.R` | Excel workbook: gene-level Fisher results and gene specificity list |
| Disease specificity | `scripts/analysis_coloc_specificity.R` | Eight numerical result tables |
| MDTPS | `scripts/run_mdtps.R` | Excel workbook of target-prioritization scores |
| Direct single-cell annotation | `scripts/scrna_direct_annotation.R` | GSE117570 annotated objects and Plasma differential-expression tables |
| Hierarchical single-cell annotation | `scripts/scrna_hierarchical_annotation.R` | GSE166555 annotated objects and CD8 differential-expression tables |
| Virtual knockout | `scripts/virtual_knockout.R` | scTenifoldKnk objects and differential-regulation tables for a selected task |

## Run

Install the packages in [docs/SOFTWARE.md](docs/SOFTWARE.md), prepare the inputs described in [data/README.md](data/README.md), and run from the repository root. `research-code.Rproj` opens this directory in RStudio. Use a fresh R session for each analysis.

```sh
Rscript --vanilla scripts/analysis_final.R
Rscript --vanilla scripts/analysis_coloc_specificity.R
Rscript --vanilla scripts/run_mdtps.R
Rscript --vanilla scripts/scrna_direct_annotation.R
Rscript --vanilla scripts/scrna_hierarchical_annotation.R
Rscript --vanilla scripts/virtual_knockout.R S8b
```

Outputs are saved under `results/` and replaced when the same module is rerun. Study data and generated results are supplied separately. The code bundle contains analysis and table exports; figure rendering and Word report generation are omitted.

## Analysis scope

- [Specificity and MDTPS](docs/SPECIFICITY_AND_MDTPS.md): the cell-type entry point generates the two worksheets in `CellType_Specificity_Results.xlsx`. MDTPS uses prepared evidence, including GTEx and eQTLGen, through `R/calculate_mdtps.R`. A base-R example is available in `examples/mdtps_example.R`.
- [Single-cell workflows](docs/SINGLE_CELL.md): two runnable examples illustrate direct annotation (GSE117570) and hierarchical annotation (GSE166555). The 11 datasets for 10 cancers, cluster labels and stage settings are recorded in `config/`. Dataset-specific sample selection, cleaning and reclustering are described in the single-cell documentation; the examples do not automate all lineages in all datasets.
- [Virtual knockout](docs/VIRTUAL_KNOCKOUT.md): `config/knockout_tasks.csv` records the 12 gene–subtype–dataset combinations and their settings. Eleven tasks use annotated Seurat objects; RPS23 uses its prepared count matrix and seed 123.

The existing root-level CSMR scripts (`csMR.R`, `merge_cell.R`, `mr.R`, `coloc.R`, `find_diff.R`) provide the upstream analysis code. These downstream modules use prepared inputs and relative paths.
