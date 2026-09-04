# Cell-type specificity: the two tables in CellType_Specificity_Results.xlsx.
# Run from the repository root with prepared eQTL and colocalization inputs.
library(readxl)
library(dplyr)
library(openxlsx)

EQTL <- file.path("data", "cell_type", "merged_eqtl.xlsx")
COLOC <- file.path("data", "cell_type", "sig_coloc.xlsx")
OUTDIR <- file.path("results", "cell_type")
stopifnot(file.exists(EQTL), file.exists(COLOC))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

LINEAGE_ORDER <- c("CD4+ T", "CD8+ T", "B cell", "NK", "Myeloid")
LINEAGES <- list("CD4+ T"=c("cd4nc","cd4et","cd4sox4"),
                 "CD8+ T"=c("cd8nc","cd8et","cd8s100b"),
                 "B cell"=c("bin","bmem","plasma"),
                 "NK"=c("nk","nkr"),
                 "Myeloid"=c("monoc","mononc","dc"))
CELL_TO_LIN <- setNames(rep(names(LINEAGES), lengths(LINEAGES)), unlist(LINEAGES))

# Count distinct genes within each cell type.
eqtl <- read_excel(EQTL)
df <- eqtl %>%
  group_by(cell_type = source_file) %>%
  summarise(n_eGenes = n_distinct(gene.exposure), .groups="drop")
N_EGENE_TOTAL <- sum(df$n_eGenes)

coloc <- read_excel(COLOC)
coloc_gene <- coloc %>%
  group_by(cell_type = cell) %>%
  summarise(n_coloc_genes = n_distinct(gene), .groups="drop")
N_COLOC_GENE_TOTAL <- sum(coloc_gene$n_coloc_genes)

df <- df %>%
  left_join(coloc_gene, by="cell_type") %>%
  mutate(n_coloc_genes = coalesce(n_coloc_genes, 0L),
         gene_rate = n_coloc_genes / n_eGenes,
         lineage = factor(CELL_TO_LIN[cell_type], levels=LINEAGE_ORDER))

GLOBAL_RATE <- N_COLOC_GENE_TOTAL / N_EGENE_TOTAL

# Table 1: original gene-level Fisher tests; BH uses numeric P values.
fisher_res <- lapply(seq_len(nrow(df)), function(i) {
  a <- df$n_coloc_genes[i]; b <- df$n_eGenes[i] - a
  c <- N_COLOC_GENE_TOTAL - a; d <- (N_EGENE_TOTAL - df$n_eGenes[i]) - c
  mat <- matrix(c(a,b,c,d), nrow=2)
  alt <- ifelse(a/df$n_eGenes[i] > GLOBAL_RATE, "greater", "less")
  ft_p  <- fisher.test(mat, alternative=alt)
  ft_ci <- fisher.test(mat, conf.int=TRUE, conf.level=0.95)
  list(OR=as.numeric(ft_p$estimate), CI_l=ft_ci$conf.int[1],
       CI_h=ft_ci$conf.int[2], P=ft_p$p.value)
})

df$Fisher_OR  <- sapply(fisher_res, `[[`, "OR")
df$Fisher_CIl <- sapply(fisher_res, `[[`, "CI_l")
df$Fisher_CIh <- sapply(fisher_res, `[[`, "CI_h")
df$Fisher_P   <- sapply(fisher_res, `[[`, "P")
df$Fisher_FDR <- p.adjust(df$Fisher_P, method="BH")

tbl_fisher <- df[,c("cell_type","lineage","n_eGenes","n_coloc_genes","gene_rate",
                    "Fisher_OR","Fisher_CIl","Fisher_CIh","Fisher_P","Fisher_FDR")]
tbl_fisher$gene_rate  <- round(tbl_fisher$gene_rate, 4)
tbl_fisher$Fisher_OR  <- round(tbl_fisher$Fisher_OR, 2)
tbl_fisher$Fisher_P   <- signif(tbl_fisher$Fisher_P, 3)
tbl_fisher$Fisher_FDR <- signif(tbl_fisher$Fisher_FDR, 3)
tbl_fisher$Fisher_CIl[is.infinite(tbl_fisher$Fisher_CIl)] <- NA
tbl_fisher$Fisher_CIh[is.infinite(tbl_fisher$Fisher_CIh)] <- NA

# Table 2: number and names of cell types for each coloc gene.
gene_cells_list <- coloc %>%
  distinct(gene, cell) %>%
  group_by(gene) %>%
  summarise(
    N_cell_types = n(),
    Cell_types = paste(sort(cell), collapse = ", "),
    .groups = "drop") %>%
  arrange(N_cell_types, gene)

write.xlsx(
  list(Table1_GeneLevel_Fisher = tbl_fisher,
       Table2_GeneList_by_Specificity = gene_cells_list),
  file.path(OUTDIR, "CellType_Specificity_Results.xlsx"),
  rowNames = FALSE, overwrite = TRUE
)
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
message("Saved the two result tables in ", OUTDIR)
