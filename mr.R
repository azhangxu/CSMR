library(TwoSampleMR)
library(dplyr)

exposure_file <- "results/instruments/merged_eqtl.xlsx"   # merge_cell.R 的输出
outcome_dir   <- "data/outcome"                          # <Cancer_name>_full_outcome.rds
out_dir       <- "results/mr"
outcomes      <- read.csv("config/mr_outcomes.csv")      # 12 个癌种及病例/对照数
dir.create(file.path(out_dir, "sig"), recursive = TRUE, showWarnings = FALSE)

exp <- openxlsx2::read_xlsx(exposure_file)

for (cancer in outcomes$cancer) {
  cancer_file <- gsub(" ", "_", cancer)
  full_outcome <- readRDS(file.path(outcome_dir, paste0(cancer_file, "_full_outcome.rds")))

  out <- TwoSampleMR::format_data(
    dat = full_outcome,
    type = "outcome",
    snps = exp$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )

  dat <- harmonise_data(exp, out, action = 2)
  res <- generate_odds_ratios(mr(dat, method_list = c("mr_wald_ratio", "mr_ivw")))
  res$fdr <- p.adjust(res$pval, method = "BH")   # 每个癌种内部做 Benjamini-Hochberg 校正

  # 敏感性分析：Steiger 方向检验；IVW Cochran's Q（≥2 个 SNP）；MR-Egger Q 和截距（≥3 个 SNP）
  steiger <- directionality_test(dat) |>
    select(id.exposure, id.outcome, correct_causal_direction, steiger_pval)
  het <- mr_heterogeneity(dat)
  het_ivw <- het |>
    filter(method == "Inverse variance weighted") |>
    select(id.exposure, id.outcome, het_Q_IVW = Q, het_Q_pval_IVW = Q_pval)
  het_egger <- het |>
    filter(method == "MR Egger") |>
    select(id.exposure, id.outcome, het_Q_MR_Egger = Q, het_Q_pval_MR_Egger = Q_pval)
  pleio <- mr_pleiotropy_test(dat) |>
    select(id.exposure, id.outcome, pleio_egger_intercept = egger_intercept,
           pleio_se = se, pleio_pval = pval)

  res <- res |>
    left_join(steiger, by = c("id.exposure", "id.outcome")) |>
    left_join(het_ivw, by = c("id.exposure", "id.outcome")) |>
    left_join(het_egger, by = c("id.exposure", "id.outcome")) |>
    left_join(pleio, by = c("id.exposure", "id.outcome"))
  res$cancer <- cancer

  writexl::write_xlsx(res, file.path(out_dir, paste0("res_mi_", cancer_file, "_full_outcome.xlsx")))
  writexl::write_xlsx(subset(res, fdr < 0.05),
                      file.path(out_dir, "sig", paste0("res_mi_sig_", cancer_file, "_full_outcome.xlsx")))
}
