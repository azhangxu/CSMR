library(coloc)
library(fs)
library(dplyr)
library(writexl)
library(tidyr)
library(purrr)
library(TwoSampleMR)
library(openxlsx2)
library(data.table)
setwd("/home/y410993/zx/OA/results/")
###########################################################################################

samplesize<-read_xlsx("/home/y410993/zx/onek1k/eqtl_sample_size_table.xlsx")

path <- "/home/y410993/zx/OA/"
rds_files <- list.files(path, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
file_names <- basename(rds_files)
rds_files <- rds_files[order(file_names)]
tumor_names <- sub("_full_outcome\\.rds$", "", file_names)

path1 <- "/home/y410993/zx/OA/results/sig/"
rds_files1 <- list.files(path1, pattern = "\\.xlsx$", full.names = TRUE, recursive = TRUE)
file_names1 <- basename(rds_files1)
rds_files1 <- rds_files1[order(file_names1)]

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
i=1
j=1
m=1
for (i in 1:9) {
  # 读sig：失败就跳过本 i
  sig <- tryCatch({
    openxlsx2::read_xlsx(rds_files1[i])
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
    readRDS(rds_files[i])
  }, error = function(e) {
    message(sprintf("[skip i=%d] readRDS error: %s", i, e$message))
    return(NULL)
  })
  if (is.null(f)) next
  
  f$ncase.outcome =	25509
  f$ncontrol.outcome = 40941
  
  only <- unique(sig$cell)
  for (j in seq_along(only)) {
    
    # 读 eQTL：失败就跳过本 j
    cell <- tryCatch({
      data.table::fread(paste0("/home/y410993/zx/onek1k/",
                               only[j], "_eqtl_table.tsv.gz"))
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
        sig_cell <- dplyr::filter(cell, GENE == sig_gene$gene[m])
        if (nrow(sig_cell) == 0) stop("sig_cell empty")
        
        # 计算列
        sig_cell$maf          <- maf_from_eaf(sig_cell$A2_FREQ_ONEK1K)
        sig_cell[["z-score"]] <- qnorm(sig_cell$P_VALUE/2, lower.tail = FALSE)
        sig_cell[["BETA"]]    <- sig_cell$SPEARMANS_RHO /
          sqrt(2 * sig_cell$A2_FREQ_ONEK1K * (1 - sig_cell$A2_FREQ_ONEK1K))
        sig_cell[["se"]]      <- abs(sig_cell$BETA / sig_cell$`z-score`)
        sig_cell$verbeta      <- (sig_cell$se)^2
        sig_cell$N            <- N
        
        # 取最小P的POS并圈定GWAS窗口
        POS <- sig_cell$POS[which.min(sig_cell$P_VALUE)]
        gwas <- dplyr::filter(
          f,
          chr.outcome == sig_cell$CHR[1]#,
          #pos.outcome >= (POS - 1e6),
          #pos.outcome <= (POS + 1e6)
        )
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
        Summary$cancer    <- tumor_names[i]
        #writexl::write_xlsx(Summary,paste0("/home/y410993/zx/onek1k/coloc_new/",sig_gene$gene[m],"_",only[j],"_",tumor_names[i],".xlsx")) 
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
writexl::write_xlsx(final_summary, "coloc_all_summary.xlsx")


