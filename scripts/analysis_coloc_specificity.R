# ============================================================
# 疾病特异性分析 — 基于 Colocalization 结果 (PP.H4 > 0.7)
# 输入: sig_coloc.xlsx (coloc pairs) + supplement_total.xlsx (ncase)
# 输出: 8 个 CSV + 控制台报告
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)

# ---- 0. Paths: run from the repository root ----
INPUT_COLOC <- file.path("data", "disease", "sig_coloc.xlsx")
INPUT_ST2   <- file.path("data", "disease", "supplement_total.xlsx")
OUTDIR <- file.path("results", "disease")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(INPUT_COLOC), file.exists(INPUT_ST2))

message("INPUT_COLOC: ", INPUT_COLOC)
message("INPUT_ST2:   ", INPUT_ST2)

# ---- 1. 读取数据 ----
coloc <- read_excel(INPUT_COLOC, sheet = 1) %>%
  select(gene, cell, cancer) %>%
  mutate(cancer = trimws(cancer))            # 去掉前后空格

# ST2.outcomes: Row 1 = merged title, Row 2 = headers, Row 3+ = data
st2 <- read_excel(INPUT_ST2, sheet = "ST2.outcomes", skip = 1)
colnames(st2) <- c("Cancer","Link","N_Case","N_Control","N_Total",
                    "ST3_Cancer_Name","N_FDR_sig","N_coloc","N_FDR_coloc","N_tested")
st2 <- st2 %>% filter(!is.na(N_Case))  # drop header & empty rows

ncase_tbl <- st2 %>% select(cancer = ST3_Cancer_Name, ncase = N_Case) %>% distinct()

message(sprintf("Coloc: %d pairs | %d genes | %d cells | %d cancers",
                nrow(coloc), n_distinct(coloc$gene),
                n_distinct(coloc$cell), n_distinct(coloc$cancer)))

# ---- 2. 癌症简称 ----
cancer_short_map <- c(
  "Breast_Cancer"             = "Breast",
  "Clear_Cell_Renal_Cell_Carcinoma" = "ClearCell RCC",
  "Colorectal_Cancer"         = "Colorectal",
  "Endometrial cancer"        = "Endometrial",
  "Lung_Adenocarcinoma"       = "Lung Adeno",
  "Lung_Squamous_Cell_Carcinoma" = "Lung SCC",
  "Ovarian_cancer"            = "Ovarian_cancer",
  "Prostate_Cancer"           = "Prostate",
  "Small_Cell_Lung_Cancer"    = "SCLC",
  "Thyroid_cancer"            = "Thyroid",
  "lung_cancer"               = "Lung(overall)",
  "melanoma"                  = "Melanoma"
)

# ============================================================
# Step 1: Spearman — 诊断样本量偏倚
# ============================================================
message("\n=== Step 1: Spearman (ncase vs coloc genes) ===")

total_sig <- coloc %>%
  group_by(cancer) %>%
  summarise(total_sig_genes = n_distinct(gene), .groups = "drop") %>%
  left_join(ncase_tbl, by = "cancer")

s <- cor.test(log10(total_sig$ncase), log10(total_sig$total_sig_genes), method = "spearman")
message(sprintf("Spearman rho = %.3f, p = %.4f", s$estimate, s$p.value))

write.csv(
  total_sig %>% mutate(cancer_short = recode(cancer, !!!cancer_short_map)),
  file.path(OUTDIR, "Fig2_ncase_vs_genes.csv"), row.names = FALSE
)

# ============================================================
# Step 2: Poisson 回归 — 校正样本量, 识别富集/匮乏
# ============================================================
message("\n=== Step 2: Offset-corrected Poisson ===")

cell_cancer_counts <- coloc %>%
  group_by(cell, cancer) %>%
  summarise(sig = n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell, cancer, fill = list(sig = 0L)) %>%
  left_join(ncase_tbl, by = "cancer")

cell_types <- sort(unique(coloc$cell))

poisson_results <- lapply(cell_types, function(ct) {
  sub <- cell_cancer_counts %>% filter(cell == ct)
  fit <- glm(sig ~ 1, offset = log(sub$ncase), family = poisson(), data = sub)
  exp_val <- fitted(fit)
  rr <- sub$sig / exp_val; rr[is.nan(rr)] <- 0

  data.frame(
    cell_type = ct, cancer = sub$cancer,
    sig_genes = sub$sig, expected = exp_val,
    risk_ratio = rr, ncase = sub$ncase,
    stringsAsFactors = FALSE
  )
})

all_res <- bind_rows(poisson_results)
all_res$p_raw <- mapply(function(o, e) {
  if (is.na(e) || e == 0) return(1)
  if (o > e) ppois(o - 0.5, e, lower.tail = FALSE)
  else       ppois(o + 0.5, e, lower.tail = TRUE)
}, all_res$sig_genes, all_res$expected)
all_res$p_fdr <- p.adjust(all_res$p_raw, method = "BH")
all_res$cancer_short <- recode(all_res$cancer, !!!cancer_short_map)

fig1 <- all_res %>% select(cell_type, cancer, cancer_short, sig_genes,
                            expected, risk_ratio, p_raw, p_fdr, ncase)
write.csv(fig1, file.path(OUTDIR, "Fig1_cell_cancer_enrichment.csv"), row.names = FALSE)
write.csv(all_res %>% select(risk_ratio), file.path(OUTDIR, "Fig2_rr_distribution.csv"), row.names = FALSE)

enr <- filter(fig1, p_fdr < 0.05, risk_ratio > 2)
dep <- filter(fig1, p_fdr < 0.05, risk_ratio < 0.3)
message(sprintf("富集 (FDR<0.05, RR>2): %d", nrow(enr)))
message(sprintf("匮乏 (FDR<0.05, RR<0.3): %d", nrow(dep)))
if (nrow(enr) > 0) print(select(enr, cell_type, cancer_short, sig_genes, risk_ratio, p_fdr))

# ============================================================
# Step 3: Jaccard + 层次聚类 — 癌症间基因重叠模式
# ============================================================
message("\n=== Step 3: Jaccard + Clustering ===")

gene_cancer <- coloc %>% select(gene, cancer) %>% distinct()

# 基因跨癌症分布
gp <- gene_cancer %>% group_by(gene) %>%
  summarise(n_cancers = n_distinct(cancer),
            cancers = paste(sort(unique(cancer)), collapse = "; "), .groups = "drop")

pd <- gp %>% count(n_cancers, name = "n_genes") %>%
  arrange(n_cancers) %>% mutate(pct = round(n_genes / sum(n_genes) * 100, 1))

write.csv(pd, file.path(OUTDIR, "Fig4_gene_pancancer.csv"), row.names = FALSE)
write.csv(
  gp %>% mutate(specificity = case_when(
    n_cancers == 1 ~ "Cancer-specific", n_cancers %in% 2:3 ~ "Low pan-cancer",
    n_cancers %in% 4:5 ~ "Moderate", n_cancers >= 6 ~ "High pan-cancer"
  )) %>% arrange(desc(n_cancers), gene) %>% select(gene, n_cancers, cancers, specificity),
  file.path(OUTDIR, "Fig4_gene_detail.csv"), row.names = FALSE
)

# Jaccard 矩阵
sets <- split(gene_cancer$gene, gene_cancer$cancer)
cn   <- names(sets)
jmat <- matrix(NA, length(cn), length(cn), dimnames = list(cn, cn))
for (i in seq_along(cn)) for (j in seq_along(cn))
  jmat[i, j] <- length(intersect(sets[[i]], sets[[j]])) /
                length(union(sets[[i]], sets[[j]]))

cns <- recode(colnames(jmat), !!!cancer_short_map)
jdf <- cbind(Cancer = cns, as.data.frame(jmat))
colnames(jdf)[-1] <- cns
write.csv(jdf, file.path(OUTDIR, "Fig3_cancer_jaccard.csv"), row.names = FALSE)

# 层次聚类
colnames(jmat) <- rownames(jmat) <- cns
diag(jmat) <- NA
hc <- hclust(as.dist(1 - jmat), method = "ward.D2")
k2 <- cutree(hc, k = 2)

intra_j <- sapply(1:2, function(k) mean(jmat[k2==k, k2==k][upper.tri(jmat[k2==k, k2==k])]))
hi <- which.max(intra_j); lo <- which.min(intra_j)

cl_df <- data.frame(Cancer = cns, Cluster = k2,
  Group = ifelse(k2 == hi, "Cluster A (免疫共享型)", "Cluster B (独立型)"),
  Genes = sapply(cn, function(x) length(sets[[x]])),
  stringsAsFactors = FALSE
) %>% arrange(Group, desc(Genes))
write.csv(cl_df, file.path(OUTDIR, "Fig3_cancer_clusters.csv"), row.names = FALSE)

c1 <- cns[k2 == hi]; c2 <- cns[k2 == lo]
ic1 <- jmat[c1, c1][upper.tri(jmat[c1, c1])]
ic2 <- jmat[c2, c2][upper.tri(jmat[c2, c2])]
ix  <- as.vector(jmat[c1, c2])

cs_df <- data.frame(
  Comparison = c(paste0("Cluster A (n=", length(c1), ")"), paste0("Cluster B (n=", length(c2), ")"), "Inter-cluster"),
  Mean_J = c(mean(ic1), mean(ic2), mean(ix)),
  SD_J   = c(sd(ic1), sd(ic2), sd(ix)),
  stringsAsFactors = FALSE
)
write.csv(cs_df, file.path(OUTDIR, "Fig3_cluster_stats.csv"), row.names = FALSE)

message(sprintf("Cluster A: %s (intra-J=%.3f)", paste(c1, collapse=", "), mean(ic1)))
message(sprintf("Cluster B: %s (intra-J=%.3f)", paste(c2, collapse=", "), mean(ic2)))
message(sprintf("Inter-cluster J = %.3f", mean(ix)))

# ============================================================
# 完成
# ============================================================
cat("\n========== 完成 ==========\n")
for (f in sort(list.files(OUTDIR, pattern = "\\.csv$", full.names = TRUE)))
  cat(sprintf("  %s (%d rows)\n", basename(f), nrow(read.csv(f))))
