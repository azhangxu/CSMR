#------------------------------------------------------------
# 加载依赖包
#------------------------------------------------------------
# install.packages("readxl")
# install.packages("dplyr")
# install.packages("purrr")
library(readxl)
library(dplyr)
library(purrr)

#------------------------------------------------------------
# 设置路径（修改为你的文件所在目录）
#------------------------------------------------------------
input_dir <- "results/instruments"   # csMR.R 的输出
output_file <- file.path(input_dir, "merged_eqtl.xlsx")
samplesize_file <- "data/onek1k/eqtl_sample_size_table.xlsx"   # 列：Prefix, Sample_size

#------------------------------------------------------------
# 获取所有 .xlsx 文件
#------------------------------------------------------------
files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
files <- setdiff(files, output_file)   # 重跑时不读入上一次的合并结果
files
# 例如：
# [1] ".../bin.xlsx" ".../bmem.xlsx" ".../cd4et.xlsx" ...

#------------------------------------------------------------
# 定义读取函数：每个文件读入并添加文件名列
#------------------------------------------------------------
read_and_tag <- function(file) {
  message("📄 读取文件：", basename(file))
  dat <- readxl::read_xlsx(file)
  dat <- dat %>%
    mutate(source_file = tools::file_path_sans_ext(basename(file)))  # 去掉扩展名作为标签
  return(dat)
}

#------------------------------------------------------------
# 批量读取并合并
#------------------------------------------------------------
merged_data <- map_dfr(files, read_and_tag)

# 查看合并结果
dplyr::glimpse(merged_data)

#------------------------------------------------------------
# 可选：保存为 xlsx 或 csv
#------------------------------------------------------------
# 保存为 Excel
# install.packages("writexl")

merged_data <- merged_data %>%
  mutate(id.exposure = paste(id.exposure, source_file, sep = "_"))


# 每个细胞类型的 eQTL 样本量（Steiger 检验需要）
samplesize <- readxl::read_xlsx(samplesize_file)
merged_data$samplesize.exposure <- samplesize$Sample_size[match(merged_data$source_file, samplesize$Prefix)]

merged_data <- merged_data %>%
  mutate(
    # 当 se.exposure 为 NA 时，用 beta 和 pval 计算
    se.exposure = ifelse(
      is.na(se.exposure),
      abs(beta.exposure / qnorm(1 - pval.exposure / 2)),
      se.exposure
    )
  )

# 排除弱工具变量（F < 10）
merged_data$F_statistics<-(merged_data$beta.exposure)^2/(merged_data$se.exposure)^2
merged_data <- subset(merged_data, F_statistics >= 10)


# 检查结果
head(merged_data$id.exposure)
writexl::write_xlsx(merged_data, output_file)

# 或保存为压缩 CSV
# write.csv(merged_data, file.path(input_dir, "merged_eqtl.csv"), row.names = FALSE)
