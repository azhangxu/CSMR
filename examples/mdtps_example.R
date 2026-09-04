# Run from the repository root. These five rows are synthetic examples.
source("R/calculate_mdtps.R", encoding = "UTF-8")

# ============================================================
# 示例数据
# ============================================================

example_dat <- data.frame(

  gene = c(
    "GENE_A",
    "GENE_B",
    "GENE_C",
    "GENE_D",
    "GENE_E"
  ),

  mr_sig = c(
    1,
    1,
    0,
    "A",
    1
  ),

  mr_beta = c(
    0.25,
    -0.18,
    0.05,
    "A",
    0.20
  ),

  coloc_sig = c(
    1,
    1,
    0,
    "A",
    1
  ),

  gtex_sig = c(
    1,
    1,
    0,
    1,
    "A"
  ),

  gtex_beta = c(
    0.20,
    -0.15,
    0.03,
    0.12,
    "A"
  ),

  eqtlgen_sig = c(
    1,
    1,
    "A",
    0,
    1
  ),

  eqtlgen_beta = c(
    0.15,
    -0.12,
    "A",
    0.02,
    -0.10
  ),

  scrna_sig = c(
    1,
    1,
    1,
    1,
    1
  ),

  scrna_log2fc = c(
    0.60,
    0.45,
    -0.35,
    -0.70,
    0.40
  ),

  ppi_sig = c(
    1,
    1,
    0,
    "A",
    1
  ),

  drug_stage = c(
    "approved",
    "clinical",
    "none",
    "A",
    "preclinical"
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 运行示例
# ============================================================

result <- calculate_mdtps(
  dat = example_dat
)

print(result)


# ============================================================
# 查看主要结果
# ============================================================

result_main <- result[
  ,
  c(
    "gene",
    "mr_sig",
    "coloc_sig",
    "gtex_sig",
    "eqtlgen_sig",
    "scrna_sig",
    "ppi_sig",
    "n_supported_evidence",
    "B",
    "dir_mr",
    "dir_gtex",
    "dir_eqtlgen",
    "dir_scrna",
    "m_direction",
    "direction_weight_sum",
    "R",
    "C",
    "overall_direction",
    "drug_stage",
    "D",
    "BaseScore",
    "DirectionScore",
    "DrugBonus",
    "MDTPS",
    "Rank"
  )
]

print(result_main)


# ============================================================
# 从Excel读取并保存结果
#
# 如未安装：
# install.packages("readxl")
# install.packages("openxlsx")
# ============================================================

# library(readxl)
# library(openxlsx)

# 读取Excel数据
# dat <- read_excel("gene_score_input.xlsx")

# 计算评分
# result <- calculate_mdtps(dat = dat)

# 保存完整结果
# write.xlsx(
#   result,
#   file = "MDTPS完整评分结果.xlsx",
#   rowNames = FALSE
# )

# 保存主要结果
# result_main <- result[
#   ,
#   c(
#     "gene",
#     "n_supported_evidence",
#     "B",
#     "m_direction",
#     "R",
#     "C",
#     "overall_direction",
#     "D",
#     "BaseScore",
#     "DirectionScore",
#     "DrugBonus",
#     "MDTPS",
#     "Rank"
#   )
# ]

# write.xlsx(
#   result_main,
#   file = "MDTPS主要评分结果.xlsx",
#   rowNames = FALSE
# )
