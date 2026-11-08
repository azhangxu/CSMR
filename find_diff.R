# 安装一次即可
# install.packages(c("readxl", "dplyr", "purrr", "stringr", "tidyr"))

library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

# 1. 设置这些 excel 所在的文件夹路径
data_dir <- "results/mr/sig/"   # mr.R 输出的 FDR < 0.05 结果
coloc_all_summary <- read_excel("results/coloc/coloc_all_summary.xlsx")   # coloc.R 的输出


file_paths <- list.files(
  path = data_dir,
  pattern = ".*_full_outcome\\.xlsx$",
  full.names = TRUE
)

read_one_file <- function(path) {
  read_excel(path, col_types = "text") |>   # ★ 所有列按字符读入
    mutate(
      file_name  = basename(path),
      cancer_type = file_name |>
        str_remove("^res_mi_sig_") |>
        str_remove("_full_outcome\\.xlsx$") |>
        str_replace_all("_", " ")   # 与 coloc 结果中的癌种名一致
    )
}

all_data <- map_dfr(file_paths, read_one_file)

all_data<- all_data  |>
  mutate(gene_cell_cancer = paste0(id.exposure,"_",cancer_type ))

coloc_all_summary <- coloc_all_summary |>
  unite(col = merged,gene , cell, cancer, sep = "_", remove = FALSE)

A_not_in_B <- setdiff(all_data$gene_cell_cancer, coloc_all_summary$merged)
