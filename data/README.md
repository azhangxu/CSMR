# Input data specification

No real study workbooks are bundled in this code-only edition. Do not substitute the synthetic MDTPS example for the study inputs. The descriptions below state what the supplied scripts consume, rather than documenting the upstream methods that produced the files.

## Cell-type specificity

Both workbooks are read from their first worksheet.

| File under `cell_type/` | Required columns | Role |
|---|---|---|
| `merged_eqtl.xlsx` | `source_file`, `gene.exposure` | Count distinct eGenes within each cell type |
| `sig_coloc.xlsx` | `gene`, `cell` | Previously selected colocalization signals |

Recognized cell labels are `cd4nc`, `cd4et`, `cd4sox4`, `cd8nc`, `cd8et`, `cd8s100b`, `bin`, `bmem`, `plasma`, `nk`, `nkr`, `monoc`, `mononc`, and `dc`. Use matching identifiers in both input tables. Coloc genes must belong to the corresponding cell type's eGene background. These inputs generate the two result worksheets.

## Disease specificity

`disease/sig_coloc.xlsx` uses the first worksheet and requires `gene`, `cell`, and `cancer`.

`disease/supplement_total.xlsx` must contain a worksheet named `ST2.outcomes`. The first row is a title; the second row provides column headers. The script reads with `skip = 1` and assigns the following ten column names **by position**:

```text
Cancer, Link, N_Case, N_Control, N_Total,
ST3_Cancer_Name, N_FDR_sig, N_coloc, N_FDR_coloc, N_tested
```

`N_Case` must be numeric and positive for each analyzed cancer. `ST3_Cancer_Name` must match `cancer` in the coloc table, with one unambiguous case count for each cancer. The cancer labels present in the coloc input define the analysis universe; cancers absent from that file are not automatically added from the metadata workbook.

The disease module does not read its inputs from the cell-type folder. This allows the exact provenance of each module's inputs to be retained.

## MDTPS

`mdtps/gene_score_input.xlsx` must have these twelve columns on the worksheet selected by `INPUT_SHEET`:

| Columns | Meaning |
|---|---|
| `gene` | Gene identifier |
| `mr_sig`, `coloc_sig`, `gtex_sig`, `eqtlgen_sig`, `scrna_sig`, `ppi_sig` | Evidence support/significance flags |
| `mr_beta`, `gtex_beta`, `eqtlgen_beta`, `scrna_log2fc` | Signed effect values |
| `drug_stage` | Drug development stage |

Extra columns are retained. Include `cell` and `cancer` when a row represents a gene–cell–cancer combination. Define the intended row unit explicitly in the study documentation.

- Use `1` for supported/significant evidence, `0` for nonsignificant evidence, and `A` for an analysis not performed. The function also accepts several equivalent text labels, including Chinese labels documented in its source.
- Missing significance flags, blanks, and `A` are treated as zero support. Remaining evidence weights are not renormalized.
- Missing effect values, blanks, and `A` become `NA`. Supply signed betas or log fold changes, not untransformed odds ratios.
- Canonical drug-stage labels are `approved`, `clinical`, `preclinical`, and `none`. Unrecognized labels receive zero in the original implementation; check the input vocabulary before use.
- The function does not harmonize effect directions or determine significance thresholds. Those preparation steps must be performed and documented upstream.
- Duplicate identifiers are not automatically collapsed. Every supplied row is scored, and ranking is global across the supplied rows.

## Single-cell examples

Download the relevant source matrices from GEO and extract them before running the scripts:

```text
data/scrna/GSE117570/
  GSM3304009_P2_Tumor_processed_data.txt.gz
  GSM3304010_P2_Normal_processed_data.txt.gz
data/scrna/GSE166555/
  <original sample filename ending in .tsv.gz>
```

For GSE166555, each file is a wide gene-by-cell table with gene identifiers in its first column. Preserve its original filename containing the `_p<number>n` or `_p<number>t` group indicator. Store only the intended sample matrices in this input directory. GEO links for all datasets are in `config/scrna_datasets.csv`.

## Virtual knockout

Eleven tasks consume the actual final lineage-annotated Seurat objects:

```text
data/scrna/annotated/GSE181294/
  Mono_annotated.rds
  CD4T_annotated.rds
  CD8T_annotated.rds
  NK_annotated.rds
data/scrna/annotated/GSE114725/
  Mono_annotated.rds
  CD8T_annotated.rds
```

RPS23 (S8j) uses a prepared gene-by-cell count matrix instead of an annotated Seurat object:

```text
data/scrna/knockout_inputs/GSE114725/
  CD8_S100B_Tumor_count_subset_for_scTenifoldKnk.rds
```

This matrix must contain the selected genes and tumour CD8_S100B cells used for that run, including RPS23. It is passed unchanged to network reconstruction.

Use the exact filenames listed in `config/knockout_tasks.csv` if they differ from this outline. Required metadata and tumour-group labels are documented in `docs/VIRTUAL_KNOCKOUT.md`. These objects are intermediate scientific inputs and are not included in the code bundle. Raw matrices, RDS files and the original result workbook are excluded from Git tracking.
