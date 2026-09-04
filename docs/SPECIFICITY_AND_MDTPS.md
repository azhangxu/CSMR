# Specificity analyses and MDTPS

Run scripts from the repository root in separate R sessions. Required worksheets and columns are listed in [data/README.md](../data/README.md); package installation is described in [SOFTWARE.md](SOFTWARE.md).

## Cell-type specificity

```sh
Rscript --vanilla scripts/analysis_final.R
```

Inputs are `data/cell_type/merged_eqtl.xlsx` and `data/cell_type/sig_coloc.xlsx`. The coloc table must already contain the selected signals (PP.H4 > 0.7 in the study).

The script writes `results/cell_type/CellType_Specificity_Results.xlsx` with two worksheets:

| Worksheet | Contents |
|---|---|
| `Table1_GeneLevel_Fisher` | eGene and coloc-gene counts, proportions, Fisher odds ratios, 95% confidence intervals, P values and BH-adjusted P values for each cell type |
| `Table2_GeneList_by_Specificity` | Each coloc gene, the number of cell types in which it occurs, and their names |

The gene-level tests retain the study script's calculation: compare each cell type with the pooled remaining cell types, using a one-sided Fisher test (`greater` above the global coloc rate, `less` otherwise). Confidence intervals use the two-sided test. BH correction is applied across the 14 cell types to numeric P values before display rounding. Counts represent gene–cell-type pairs; a gene found in multiple cell types contributes to each corresponding cell type.

This entry point covers these two tables. It also saves `sessionInfo.txt`.

## Disease specificity

```sh
Rscript --vanilla scripts/analysis_coloc_specificity.R
```

Inputs are `data/disease/sig_coloc.xlsx` and `data/disease/supplement_total.xlsx`. The script exports eight CSV tables to `results/disease/`, covering sample-size associations, enrichment/risk ratios, gene sharing, Jaccard similarity and clustering. Filenames containing `Fig` identify the corresponding figure's source data; the files are numerical result tables.

## MDTPS

```sh
Rscript --vanilla scripts/run_mdtps.R
```

`R/calculate_mdtps.R` contains the scoring function; `scripts/run_mdtps.R` reads the prepared evidence table and writes the results. The input is `data/mdtps/gene_score_input.xlsx`; `INPUT_SHEET` selects the worksheet. Outputs are `results/mdtps/MDTPS_results.xlsx` (full and main results) and `sessionInfo.txt`.

```text
MDTPS = 75 B + 10 B (2 C - 1) + 15 B D
```

Evidence weights for `B` are MR 0.30, COLOC 0.20, GTEx 0.10, eQTLGen 0.10, scRNA 0.20 and PPI 0.10. Direction weights are MR 0.40 and GTEx/eQTLGen/scRNA 0.20 each, using significant evidence with nonzero, nonmissing effects. With fewer than two directional sources, `C = 0.5`; otherwise it is the absolute weighted mean direction. Drug-stage scores for `D` are approved 1.0, clinical 0.6, preclinical 0.3 and none 0.0. Scores are bounded to 0–100.

The function consumes prepared evidence flags and effects, including GTEx and eQTLGen; it does not generate those upstream analyses. Missing or unperformed evidence contributes no support, without reweighting the remaining evidence. All input rows are ranked together.

For a five-gene synthetic example requiring only base R:

```sh
Rscript --vanilla examples/mdtps_example.R
```
