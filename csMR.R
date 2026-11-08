library(data.table)
library(TwoSampleMR)
library(tidyverse)

# clump_data() 调用 OpenGWAS API。令牌写在 ~/.Renviron（OPENGWAS_JWT=<token>），不要写进代码
if (Sys.getenv("OPENGWAS_JWT") == "") stop("OPENGWAS_JWT is not set; add it to ~/.Renviron")

onek1k_dir     <- "data/onek1k"          # OneK1K 的 <cell>_eqtl_table.tsv.gz
instrument_dir <- "results/instruments"  # 每个细胞类型的工具变量
dir.create(instrument_dir, recursive = TRUE, showWarnings = FALSE)

# 可选：明确使用 TwoSampleMR 的 format_data
# install.packages("TwoSampleMR")
# install.packages(c("dplyr", "purrr"))
library(dplyr)
library(purrr)

#—— 可配置参数-----------------
fmt_args <- list(
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "RSID",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "A2_FREQ_ONEK1K",
  effect_allele_col = "A2",
  other_allele_col = "A1",
  pval_col = " P_VALUE",    # 注意：原列名前有空格
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  gene_col = "gene",
  id_col = "GENE",
  min_pval = 1e-200,
  z_col = "z",
  info_col = "info",
  chr_col = "chr",
  pos_col = "pos",
  log_pval = FALSE
)

#—— 工具函数：按基因分组后逐组调用 format_data -------
format_eqtl_by_gene <- function(df, fmt_args, prefer_best_per_snp = FALSE) {
  # 1) 处理列名里可能的前后空格，避免 " P_VALUE" 这类问题
  names(df) <- trimws(names(df))
  
  # 2) 兼容 gene 与 GENE：如果缺失其一但存在另一列，则自动补齐一份
  if (!("gene" %in% names(df)) && "GENE" %in% names(df)) {
    df <- df %>% mutate(gene = .data$GENE)
  }
  if (!("GENE" %in% names(df)) && "gene" %in% names(df)) {
    df <- df %>% mutate(GENE = .data$gene)
  }
  # 5) 按基因拆分
  pieces <- split(df, df[[fmt_args$id_col]], drop = TRUE)
  
  # 6) 对每个基因单独跑 format_data，再绑定
  res_list <- map(pieces, function(gdf) {
    TwoSampleMR::format_data(
      dat = gdf,
      type = fmt_args$type,
      snps = fmt_args$snps,
      header = fmt_args$header,
      phenotype_col = fmt_args$phenotype_col,
      snp_col = fmt_args$snp_col,
      beta_col = fmt_args$beta_col,
      se_col = fmt_args$se_col,
      eaf_col = fmt_args$eaf_col,
      effect_allele_col = fmt_args$effect_allele_col,
      other_allele_col = fmt_args$other_allele_col,
      pval_col = trimws(fmt_args$pval_col),  # 去掉空格后的列名
      units_col = fmt_args$units_col,
      ncase_col = fmt_args$ncase_col,
      ncontrol_col = fmt_args$ncontrol_col,
      gene_col = fmt_args$gene_col,
      id_col = fmt_args$id_col,
      min_pval = fmt_args$min_pval,
      z_col = fmt_args$z_col,
      info_col = fmt_args$info_col,
      chr_col = fmt_args$chr_col,
      pos_col = fmt_args$pos_col,
      log_pval = fmt_args$log_pval
    )
  })
  
  bind_rows(res_list)
}


files <- list.files(onek1k_dir, pattern = "eqtl_table\\.tsv\\.gz$", full.names = TRUE)
prefix <- str_extract(basename(files), "^[^_]+")





######
for (i in seq_along(files)) {
  f<-fread(files[i])
  f<-subset(f,P_VALUE<5e-08)
  f<- f |>
    separate(
      col = SNPID,
      into = c("chr", "pos", "allele"),
      sep = "[:_]"
    )
  p <- f$A2_FREQ_ONEK1K
  
  # 计算真正的 Z（绝对值，不带方向）
  Z <- abs(qnorm(f$P_VALUE / 2))
  
  # 计算 beta
  f$beta <- f$SPEARMANS_RHO / sqrt(2 * p * (1 - p))
  
  # 计算 se（取绝对值，防止负数）
  f$se <- abs(f$beta / Z)
  #—— 调用示例 ----------------------------------------
  # 如果希望每个基因内部，同一 SNP 仅保留最显著一条，设 TRUE
  exp <- format_eqtl_by_gene(f, fmt_args, prefer_best_per_snp = TRUE)
  exp1<-TwoSampleMR::clump_data(exp,clump_kb = 10000,clump_r2 = 0.001)
  writexl::write_xlsx(x = exp1,path = file.path(instrument_dir, paste0(prefix[i],".xlsx")))
}








