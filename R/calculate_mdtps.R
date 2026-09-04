# ============================================================
# 多组学方向性靶点优先级评分
# Multi-omics Directional Target Prioritization Score
# 简称：MDTPS
#
# 最终评分公式：
#
# MDTPS_g =
# 75 × B_g
# + 10 × B_g × (2 × C_g - 1)
# + 15 × B_g × D_g
#
# 其中：
# B_g：基础多组学证据指数，范围为0~1
# C_g：方向一致性指数，范围为0~1
# D_g：药物转化指数，范围为0~1
#
# 显著性变量填写规则：
# 1 = 分析结果显著或提供支持
# 0 = 分析结果不显著
# A = 未进行该项分析，评分时按0处理
#
# 最终评分范围：0~100
# ============================================================


calculate_mdtps <- function(
    dat,

    # --------------------------------------------------------
    # 基础多组学证据权重
    #
    # MR：0.30
    # COLOC：0.20
    # GTEx：0.10
    # eQTLGen：0.10
    # 单细胞差异分析：0.20
    # PPI：0.10
    #
    # 权重之和必须等于1
    # --------------------------------------------------------
    evidence_weights = c(
      MR = 0.30,
      COLOC = 0.20,
      GTEx = 0.10,
      eQTLGen = 0.10,
      scRNA = 0.20,
      PPI = 0.10
    ),

    # --------------------------------------------------------
    # 方向一致性计算权重
    #
    # 只有以下四类结果参与方向计算：
    # MR、GTEx、eQTLGen、单细胞差异分析
    #
    # MR方向权重较高，为0.40
    # 其他三项各为0.20
    #
    # 权重之和必须等于1
    # --------------------------------------------------------
    direction_weights = c(
      MR = 0.40,
      GTEx = 0.20,
      eQTLGen = 0.20,
      scRNA = 0.20
    ),

    # 当综合方向绝对值小于0.20时，
    # 将总体方向标记为“方向冲突或不明确”
    direction_uncertain_cutoff = 0.20,

    # 最终结果保留的小数位数
    digits = 4
) {

  # ==========================================================
  # 第一部分：检查输入数据
  # ==========================================================

  required_cols <- c(
    "gene",
    "mr_sig",
    "mr_beta",
    "coloc_sig",
    "gtex_sig",
    "gtex_beta",
    "eqtlgen_sig",
    "eqtlgen_beta",
    "scrna_sig",
    "scrna_log2fc",
    "ppi_sig",
    "drug_stage"
  )

  missing_cols <- setdiff(
    required_cols,
    names(dat)
  )

  if (length(missing_cols) > 0) {
    stop(
      "输入数据中缺少以下变量：",
      paste(missing_cols, collapse = ", ")
    )
  }

  out <- as.data.frame(
    dat,
    stringsAsFactors = FALSE
  )


  # ==========================================================
  # 第二部分：检查权重
  # ==========================================================

  required_evidence_names <- c(
    "MR",
    "COLOC",
    "GTEx",
    "eQTLGen",
    "scRNA",
    "PPI"
  )

  if (!all(
    required_evidence_names %in% names(evidence_weights)
  )) {
    stop(
      "evidence_weights必须包含以下名称：",
      paste(required_evidence_names, collapse = ", ")
    )
  }

  evidence_weights <- evidence_weights[
    required_evidence_names
  ]

  if (abs(sum(evidence_weights) - 1) > 1e-8) {
    stop(
      "基础证据权重 evidence_weights 之和必须等于1。"
    )
  }

  required_direction_names <- c(
    "MR",
    "GTEx",
    "eQTLGen",
    "scRNA"
  )

  if (!all(
    required_direction_names %in% names(direction_weights)
  )) {
    stop(
      "direction_weights必须包含以下名称：",
      paste(required_direction_names, collapse = ", ")
    )
  }

  direction_weights <- direction_weights[
    required_direction_names
  ]

  if (abs(sum(direction_weights) - 1) > 1e-8) {
    stop(
      "方向证据权重 direction_weights 之和必须等于1。"
    )
  }


  # ==========================================================
  # 第三部分：定义显著性变量转换函数
  #
  # 允许输入：
  # 1、显著、Yes、Y、TRUE
  # 0、不显著、No、N、FALSE
  # A、未做、没有做、NOT DONE、NA、空白
  #
  # 其中A、未做和空白均按0处理
  # ==========================================================

  normalize_sig <- function(x, variable_name) {

    x <- trimws(
      as.character(x)
    )

    missing_position <- is.na(x)
    x <- toupper(x)
    x[missing_position] <- "0"

    x[x %in% c(
      "1",
      "YES",
      "Y",
      "TRUE",
      "显著",
      "SIGNIFICANT"
    )] <- "1"

    x[x %in% c(
      "0",
      "NO",
      "N",
      "FALSE",
      "不显著",
      "NONSIGNIFICANT",
      "NON-SIGNIFICANT",
      "A",
      "未做",
      "没有做",
      "未分析",
      "NOT DONE",
      "NA",
      "N/A",
      ""
    )] <- "0"

    invalid_values <- unique(
      x[!x %in% c("0", "1")]
    )

    if (length(invalid_values) > 0) {
      stop(
        variable_name,
        " 中存在无法识别的值：",
        paste(invalid_values, collapse = ", "),
        "。显著请填写1，不显著填写0，未分析填写A。"
      )
    }

    return(
      as.numeric(x)
    )
  }


  # ==========================================================
  # 第四部分：定义效应值转换函数
  #
  # mr_beta、gtex_beta、eqtlgen_beta和scrna_log2fc
  # 需要为数值型。
  #
  # 如果填写A、NA、空白或未做，则转换为真正的NA。
  # ==========================================================

  normalize_effect <- function(x, variable_name) {

    x_character <- trimws(
      as.character(x)
    )

    x_character[
      toupper(x_character) %in% c(
        "A",
        "NA",
        "N/A",
        "未做",
        "没有做",
        "未分析",
        "NOT DONE",
        ""
      )
    ] <- NA_character_

    x_numeric <- suppressWarnings(
      as.numeric(x_character)
    )

    invalid_position <- (
      !is.na(x_character) &
        is.na(x_numeric)
    )

    if (any(invalid_position)) {
      stop(
        variable_name,
        " 中存在无法转换为数值的内容：",
        paste(
          unique(x_character[invalid_position]),
          collapse = ", "
        )
      )
    }

    return(x_numeric)
  }


  # ==========================================================
  # 第五部分：标准化所有显著性变量
  # A会在这里统一转换为0
  # ==========================================================

  sig_cols <- c(
    "mr_sig",
    "coloc_sig",
    "gtex_sig",
    "eqtlgen_sig",
    "scrna_sig",
    "ppi_sig"
  )

  for (col in sig_cols) {
    out[[col]] <- normalize_sig(
      x = out[[col]],
      variable_name = col
    )
  }


  # ==========================================================
  # 第六部分：标准化所有效应值变量
  # A、NA和空白转换为NA
  # ==========================================================

  effect_cols <- c(
    "mr_beta",
    "gtex_beta",
    "eqtlgen_beta",
    "scrna_log2fc"
  )

  for (col in effect_cols) {
    out[[col]] <- normalize_effect(
      x = out[[col]],
      variable_name = col
    )
  }


  # ==========================================================
  # 第七部分：计算基础多组学证据指数B
  #
  # B_g =
  # 0.30 × I_MR
  # + 0.20 × I_COLOC
  # + 0.10 × I_GTEx
  # + 0.10 × I_eQTLGen
  # + 0.20 × I_scRNA
  # + 0.10 × I_PPI
  #
  # B的理论范围为0~1
  # ==========================================================

  out$B <- (
    evidence_weights["MR"] *
      out$mr_sig +

      evidence_weights["COLOC"] *
      out$coloc_sig +

      evidence_weights["GTEx"] *
      out$gtex_sig +

      evidence_weights["eQTLGen"] *
      out$eqtlgen_sig +

      evidence_weights["scRNA"] *
      out$scrna_sig +

      evidence_weights["PPI"] *
      out$ppi_sig
  )

  out$n_supported_evidence <- (
    out$mr_sig +
      out$coloc_sig +
      out$gtex_sig +
      out$eqtlgen_sig +
      out$scrna_sig +
      out$ppi_sig
  )


  # ==========================================================
  # 第八部分：计算各分析的效应方向
  #
  # sign(x)：
  # x > 0，返回+1
  # x < 0，返回-1
  # x = 0，返回0
  #
  # 如果效应值缺失，则方向设置为0
  # ==========================================================

  out$dir_mr <- ifelse(
    !is.na(out$mr_beta) &
      out$mr_beta != 0,
    sign(out$mr_beta),
    0
  )

  out$dir_gtex <- ifelse(
    !is.na(out$gtex_beta) &
      out$gtex_beta != 0,
    sign(out$gtex_beta),
    0
  )

  out$dir_eqtlgen <- ifelse(
    !is.na(out$eqtlgen_beta) &
      out$eqtlgen_beta != 0,
    sign(out$eqtlgen_beta),
    0
  )

  out$dir_scrna <- ifelse(
    !is.na(out$scrna_log2fc) &
      out$scrna_log2fc != 0,
    sign(out$scrna_log2fc),
    0
  )


  # ==========================================================
  # 第九部分：定义参与方向计算的有效指标H
  #
  # H = 1需要同时满足：
  # 1. 该分析结果显著，即sig=1
  # 2. 该分析具有明确的正负方向
  #
  # H = 0表示：
  # 不显著、未分析、效应值缺失或效应值为0
  # ==========================================================

  out$H_mr <- as.numeric(
    out$mr_sig == 1 &
      out$dir_mr != 0
  )

  out$H_gtex <- as.numeric(
    out$gtex_sig == 1 &
      out$dir_gtex != 0
  )

  out$H_eqtlgen <- as.numeric(
    out$eqtlgen_sig == 1 &
      out$dir_eqtlgen != 0
  )

  out$H_scrna <- as.numeric(
    out$scrna_sig == 1 &
      out$dir_scrna != 0
  )


  # ==========================================================
  # 第十部分：计算方向权重总和W
  #
  # W_g =
  # 0.40 × H_MR
  # + 0.20 × H_GTEx
  # + 0.20 × H_eQTLGen
  # + 0.20 × H_scRNA
  # ==========================================================

  out$direction_weight_sum <- (
    direction_weights["MR"] *
      out$H_mr +

      direction_weights["GTEx"] *
      out$H_gtex +

      direction_weights["eQTLGen"] *
      out$H_eqtlgen +

      direction_weights["scRNA"] *
      out$H_scrna
  )


  # ==========================================================
  # 第十一部分：计算加权方向分子N
  #
  # N_g =
  # 0.40 × H_MR × s_MR
  # + 0.20 × H_GTEx × s_GTEx
  # + 0.20 × H_eQTLGen × s_eQTLGen
  # + 0.20 × H_scRNA × s_scRNA
  # ==========================================================

  out$direction_numerator <- (
    direction_weights["MR"] *
      out$H_mr *
      out$dir_mr +

      direction_weights["GTEx"] *
      out$H_gtex *
      out$dir_gtex +

      direction_weights["eQTLGen"] *
      out$H_eqtlgen *
      out$dir_eqtlgen +

      direction_weights["scRNA"] *
      out$H_scrna *
      out$dir_scrna
  )


  # ==========================================================
  # 第十二部分：计算加权综合方向R
  #
  # R_g = N_g / W_g
  #
  # R范围为-1~1：
  # 接近+1：显著结果主要为正方向
  # 接近-1：显著结果主要为负方向
  # 接近0：不同分析方向相互冲突
  #
  # 如果没有任何方向证据，即W=0，则R=0
  # ==========================================================

  out$R <- ifelse(
    out$direction_weight_sum > 0,
    out$direction_numerator /
      out$direction_weight_sum,
    0
  )


  # ==========================================================
  # 第十三部分：计算有效方向证据数量m
  #
  # m_g =
  # H_MR + H_GTEx + H_eQTLGen + H_scRNA
  # ==========================================================

  out$m_direction <- (
    out$H_mr +
      out$H_gtex +
      out$H_eqtlgen +
      out$H_scrna
  )


  # ==========================================================
  # 第十四部分：计算方向一致性指数C
  #
  # 当m < 2时：
  # C = 0.5
  # 表示无法评价多方法一致性，不加不扣
  #
  # 当m >= 2时：
  # C = |R|
  #
  # C = 1：完全一致
  # C = 0.5：不加不扣
  # C < 0.5：存在明显方向冲突
  # C = 0：正负方向完全抵消
  # ==========================================================

  out$C <- ifelse(
    out$m_direction < 2,
    0.5,
    abs(out$R)
  )

  out$C <- pmin(
    pmax(out$C, 0),
    1
  )


  # ==========================================================
  # 第十五部分：判断总体方向
  # ==========================================================

  out$overall_direction <- ifelse(
    out$m_direction == 0,

    "无方向证据",

    ifelse(
      out$m_direction == 1,

      "方向证据不足",

      ifelse(
        abs(out$R) < direction_uncertain_cutoff,

        "方向冲突或不明确",

        ifelse(
          out$R > 0,
          "正方向",
          "负方向"
        )
      )
    )
  )


  # ==========================================================
  # 第十六部分：定义药物研发阶段转换函数
  #
  # approved / marketed / 上市：
  # D = 1.0
  #
  # clinical / phase I-II-III / 临床阶段：
  # D = 0.6
  #
  # preclinical / 临床前：
  # D = 0.3
  #
  # none / A / 无药物：
  # D = 0
  #
  # 如果一个基因对应多个药物，
  # 建议填写最高研发阶段
  # ==========================================================

  normalize_drug_stage <- function(x) {

    x <- tolower(
      trimws(
        as.character(x)
      )
    )

    x[is.na(x)] <- "none"

    score <- rep(
      0,
      length(x)
    )

    preclinical_pattern <- paste(
      c(
        "preclinical",
        "pre-clinical",
        "临床前"
      ),
      collapse = "|"
    )

    is_preclinical <- grepl(
      preclinical_pattern,
      x
    )

    score[is_preclinical] <- 0.3

    clinical_pattern <- paste(
      c(
        "clinical",
        "phase i",
        "phase ii",
        "phase iii",
        "phase 1",
        "phase 2",
        "phase 3",
        "临床阶段",
        "临床试验",
        "临床i期",
        "临床ii期",
        "临床iii期",
        "临床1期",
        "临床2期",
        "临床3期"
      ),
      collapse = "|"
    )

    is_clinical <- grepl(
      clinical_pattern,
      x
    )

    score[
      is_clinical &
        !is_preclinical
    ] <- 0.6

    approved_pattern <- paste(
      c(
        "approved",
        "marketed",
        "上市",
        "已上市",
        "批准",
        "已批准"
      ),
      collapse = "|"
    )

    is_approved <- grepl(
      approved_pattern,
      x
    )

    score[is_approved] <- 1.0

    return(score)
  }


  out$D <- normalize_drug_stage(
    out$drug_stage
  )


  # ==========================================================
  # 第十七部分：计算基础证据分
  #
  # BaseScore_g = 75 × B_g
  # ==========================================================

  out$BaseScore <- (
    75 *
      out$B
  )


  # ==========================================================
  # 第十八部分：计算方向一致性调整分
  #
  # DirectionScore_g =
  # 10 × B_g × (2 × C_g - 1)
  #
  # C = 1：最大加分
  # C = 0.5：不加不扣
  # C < 0.5：扣分
  # C = 0：最大扣分
  # ==========================================================

  out$DirectionScore <- (
    10 *
      out$B *
      (2 * out$C - 1)
  )


  # ==========================================================
  # 第十九部分：计算药物转化加分
  #
  # DrugBonus_g =
  # 15 × B_g × D_g
  #
  # 乘以B的目的是：
  # 药物信息只能增强已有多组学证据的候选基因
  # ==========================================================

  out$DrugBonus <- (
    15 *
      out$B *
      out$D
  )


  # ==========================================================
  # 第二十部分：计算最终MDTPS评分
  #
  # MDTPS_g =
  # BaseScore_g
  # + DirectionScore_g
  # + DrugBonus_g
  #
  # 即：
  # MDTPS_g =
  # 75B_g
  # + 10B_g(2C_g-1)
  # + 15B_gD_g
  # ==========================================================

  out$MDTPS <- (
    out$BaseScore +
      out$DirectionScore +
      out$DrugBonus
  )

  out$MDTPS <- pmin(
    pmax(out$MDTPS, 0),
    100
  )


  # ==========================================================
  # 第二十一部分：根据MDTPS进行排名
  # ==========================================================

  out$Rank <- rank(
    -out$MDTPS,
    ties.method = "min",
    na.last = "keep"
  )


  # ==========================================================
  # 第二十二部分：对数值结果进行四舍五入
  # ==========================================================

  round_cols <- c(
    "B",
    "direction_weight_sum",
    "direction_numerator",
    "R",
    "C",
    "D",
    "BaseScore",
    "DirectionScore",
    "DrugBonus",
    "MDTPS"
  )

  for (col in round_cols) {
    out[[col]] <- round(
      out[[col]],
      digits = digits
    )
  }


  # ==========================================================
  # 第二十三部分：按照最终排名排序并返回结果
  # ==========================================================

  out <- out[
    order(
      out$Rank,
      out$gene,
      na.last = TRUE
    ),
  ]

  rownames(out) <- NULL

  return(out)
}


