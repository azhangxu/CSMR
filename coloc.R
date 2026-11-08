library(coloc)
library(fs)
library(dplyr)
library(writexl)
library(tidyr)
library(purrr)
library(TwoSampleMR)
library(openxlsx2)
library(data.table)
###########################################################################################

onek1k_dir  <- "data/onek1k"      # <cell>_eqtl_table.tsv.gz 和 eqtl_sample_size_table.xlsx
outcome_dir <- "data/outcome"     # <Cancer_name>_full_outcome.rds
mr_sig_dir  <- "results/mr/sig"   # mr.R 输出的 FDR < 0.05 结果
out_dir     <- "results/coloc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outcomes   <- read.csv("config/mr_outcomes.csv")
samplesize <- read_xlsx(file.path(onek1k_dir, "eqtl_sample_size_table.xlsx"))

maf_from_eaf <- function(eaf) {
  eaf <- as.numeric(eaf)
  bad <- is.na(eaf) | eaf < 0 | eaf > 1
  maf <- pmin(eaf, 1 - eaf)
  maf[bad] <- NA_real_
  maf
}

# 建议先准备一个容器
all_summary <- list()
k <- 1
for (i in seq_len(nrow(outcomes))) {
  cancer      <- outcomes$cancer[i]
  cancer_file <- gsub(" ", "_", cancer)

  # 读sig：失败就跳过本 i
  sig <- tryCatch({
    openxlsx2::read_xlsx(file.path(mr_sig_dir, paste0("res_mi_sig_", cancer_file, "_full_outcome.xlsx")))
  }, error = function(e) {
    message(sprintf("[skip i=%d] read_xlsx error: %s", i, e$message))
    return(NULL)
  })
  if (is.null(sig)) next
  
  sig <- tryCatch({
    tidyr::separate(sig, id.exposure, into = c("gene", "cell"), sep = "_")
  }, error = function(e) {
    message(sprintf("[skip i=%d] separate error: %s", i, e$message))
    return(NULL)
  })
  if (is.null(sig)) next
  
  # 读RDS：失败就跳过本 i
  f <- tryCatch({
    readRDS(file.path(outcome_dir, paste0(cancer_file, "_full_outcome.rds")))
  }, error = function(e) {
    message(sprintf("[skip i=%d] readRDS error: %s", i, e$message))
    return(NULL)
  })
  if (is.null(f)) next
  
  # 每个癌种自己的病例/对照数（config/mr_outcomes.csv）
  f$ncase.outcome    <- outcomes$n_case[i]
  f$ncontrol.outcome <- outcomes$n_control[i]
  
  only <- unique(sig$cell)
  for (j in seq_along(only)) {
    
    # 读 eQTL：失败就跳过本 j
    cell <- tryCatch({
      data.table::fread(file.path(onek1k_dir, paste0(only[j], "_eqtl_table.tsv.gz")))
    }, error = function(e) {
      message(sprintf("[skip i=%d j=%d] fread error: %s", i, j, e$message))
      return(NULL)
    })
    if (is.null(cell) || nrow(cell) == 0) {
      message(sprintf("[skip i=%d j=%d] cell table empty.", i, j))
      next
    }
    
    sig_gene  <- subset(sig, cell == only[j])
    only_gene <- unique(sig_gene$gene)
    
    # 取 N（没有就跳过）
    N <- tryCatch({
      subset(samplesize, Prefix == only[j])$Sample_size[1]
    }, error = function(e) NA_real_)
    if (is.na(N)) {
      message(sprintf("[skip i=%d j=%d] no sample size for %s", i, j, only[j]))
      next
    }
    
    for (m in seq_along(only_gene)) {
      
      ok <- tryCatch({
        # ===== 核心计算开始 =====
        sig_cell <- dplyr::filter(cell, GENE == only_gene[m])
        if (nrow(sig_cell) == 0) stop("sig_cell empty")
        
        # 计算列
        sig_cell$maf          <- maf_from_eaf(sig_cell$A2_FREQ_ONEK1K)
        sig_cell[["z-score"]] <- qnorm(sig_cell$P_VALUE/2, lower.tail = FALSE)
        sig_cell[["BETA"]]    <- sig_cell$SPEARMANS_RHO /
          sqrt(2 * sig_cell$A2_FREQ_ONEK1K * (1 - sig_cell$A2_FREQ_ONEK1K))
        sig_cell[["se"]]      <- abs(sig_cell$BETA / sig_cell$`z-score`)
        sig_cell$verbeta      <- (sig_cell$se)^2
        sig_cell$N            <- N
        
        # 共定位区域：该基因在 OneK1K cis-eQTL 表中检验过的全部 SNP，与同一染色体上的 GWAS SNP 取交集
        gwas <- dplyr::filter(f, chr.outcome == sig_cell$CHR[1])
        if (nrow(gwas) == 0) stop("gwas window empty")
        
        gwas[["z-score"]] <- qnorm(gwas$pval.outcome/2, lower.tail = FALSE)
        gwas[["s"]]       <- gwas$ncase.outcome / gwas$samplesize.outcome
        gwas$verbeta      <- (gwas$se.outcome)^2
        
        # 匹配并去重
        sameSNP <- intersect(sig_cell$RSID, gwas$SNP)
        if (length(sameSNP) == 0) stop("no overlapping SNPs")
        
        sig_cell <- sig_cell |>
          dplyr::filter(RSID %in% sameSNP, !is.na(P_VALUE)) |>
          dplyr::arrange(RSID, P_VALUE) |>
          dplyr::distinct(RSID, .keep_all = TRUE)
        
        gwas <- gwas |>
          dplyr::filter(SNP %in% sameSNP) |>
          dplyr::arrange(SNP) |>
          tidyr::drop_na()
        
        common_snp <- intersect(sig_cell$RSID, gwas$SNP)
        if (length(common_snp) < 2) stop("too few common SNPs after re-alignment")
        
        # 按同一顺序对齐两张表
        # 用 match 保证顺序完全一致
        gwas     <- gwas[match(common_snp, gwas$SNP), ]
        sig_cell <- sig_cell[match(common_snp, sig_cell$RSID), ]
        
        # 再做一次长度检查，防止意外
        if (nrow(gwas) != nrow(sig_cell)) stop("row mismatch after SNP alignment")
        
        if (nrow(sig_cell) < 2 || nrow(gwas) < 2) stop("too few SNPs after QC")
        
        # coloc
        dataset1 <- list(
          snp     = gwas$SNP,
          pvalues = gwas$pval.outcome,
          type    = "cc",
          beta    = gwas$beta.outcome,
          verbeta = gwas$verbeta,
          s       = gwas$s[1],
          N       = gwas$samplesize.outcome[1]
        )
        dataset2 <- list(
          snp     = sig_cell$RSID,
          pvalues = sig_cell$P_VALUE,
          type    = "quant",
          N       = sig_cell$N[1],
          verbeta = sig_cell$verbeta,
          beta    = sig_cell$BETA
        )
        
        coloc_result <- coloc.abf(
          dataset1, dataset2,
          MAF = sig_cell$maf, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
        )
        
        Summary <- as.data.frame(t(coloc_result$summary))
        ss <- coloc_result$results |>
          dplyr::arrange(dplyr::desc(SNP.PP.H4))
        
        Summary$TOP_snp   <- ss$snp[1]
        Summary$SNP.PP.H4 <- ss$SNP.PP.H4[1]
        Summary$cell      <- only[j]
        Summary$gene      <- only_gene[m]
        Summary$cancer    <- cancer
        all_summary[[k]] <- Summary
        k <- k + 1
        TRUE
        # ===== 核心计算结束 =====
      }, error = function(e) {
        message(sprintf("[skip i=%d j=%d m=%d] %s", i, j, m, e$message))
        return(FALSE)
      })
      
      if (!ok) next
    } # m
  }   # j
}     # i

# 合并保存
final_summary <- dplyr::bind_rows(all_summary)
writexl::write_xlsx(final_summary, file.path(out_dir, "coloc_all_summary.xlsx"))


