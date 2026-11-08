# CSMR study analysis code

Analysis code for cell-type-specific eQTL instrument preparation, two-sample Mendelian randomization (MR), colocalization, cell-type and disease specificity, MDTPS, single-cell analysis and virtual knockout.

| Analysis | Script |
|---|---|
| OneK1K eQTL instrument preparation: significance filtering, exposure-data formatting and LD clumping | `csMR.R` |
| Merge cell-type eQTL exposure tables, add cell-type sample sizes and exclude instruments with F < 10 | `merge_cell.R` |
| Two-sample MR for 12 cancers (Wald ratio for one SNP, inverse-variance weighted for two or more SNPs), with Benjamini–Hochberg FDR correction within each cancer, Steiger directionality, Cochran's Q and MR-Egger intercept tests | `mr.R` |
| Bayesian colocalization of eQTL and GWAS signals using coloc.abf, for FDR-significant MR associations | `coloc.R` |
| Identify MR gene–cell–cancer combinations absent from the colocalization summary | `find_diff.R` |
| Cell-type specificity | `scripts/analysis_final.R` |
| Disease specificity | `scripts/analysis_coloc_specificity.R` |
| MDTPS | `scripts/run_mdtps.R` |
| Direct single-cell annotation | `scripts/scrna_direct_annotation.R` |
| Hierarchical single-cell annotation | `scripts/scrna_hierarchical_annotation.R` |
| Virtual knockout | `scripts/virtual_knockout.R` |

## Run MR and colocalization

Run the root-level scripts from the repository root in this order: `csMR.R`, `merge_cell.R`, `mr.R` and `coloc.R`. `find_diff.R` is an optional check. The scripts read:

- `data/onek1k/`: the 14 OneK1K tables (`<cell>_eqtl_table.tsv.gz`) and `eqtl_sample_size_table.xlsx` (columns `Prefix`, `Sample_size`);
- `data/outcome/`: one GWAS file per cancer in TwoSampleMR outcome format (`SNP`, `beta.outcome`, `se.outcome`, `pval.outcome`, alleles, `eaf.outcome`, `chr.outcome`, `pos.outcome`, `samplesize.outcome`), named `<Cancer_name>_full_outcome.rds` with spaces replaced by underscores, for example `Ovarian_cancer_full_outcome.rds`;
- `config/mr_outcomes.csv`: the 12 cancers, GWAS sources and case and control numbers.

LD clumping uses the OpenGWAS API. Save the token as `OPENGWAS_JWT=<token>` in `~/.Renviron`; do not write it in the scripts. Outputs are written to `results/instruments/`, `results/mr/` and `results/coloc/`.

## Run downstream modules

For the downstream modules, install the packages in [docs/SOFTWARE.md](docs/SOFTWARE.md), prepare the inputs described in [data/README.md](data/README.md), and run from the repository root. `research-code.Rproj` opens this directory in RStudio. Use a fresh R session for each analysis.

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
