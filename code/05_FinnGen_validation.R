###############################################################################
# 补充分析 5: FinnGen外部验证
# 目的：回应编辑质疑——"缺乏独立外部验证"
# 方法：使用FinnGen R12/R13 GWAS数据重复正向MR分析
# 核心优势：FinnGen独立于UK Biobank，提供真正的外部验证
# 数据来源：FinnGen (finngen.fi) → 通过IEU OpenGWAS或本地下载
# 作者：MedResearch Pro 辅助生成
# 日期：2026-06-18
###############################################################################

# ==================== 0. 环境设置 ====================
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(ggplot2)

# ==================== 1. FinnGen表型ID查询 ====================

# FinnGen在IEU OpenGWAS中的格式: finngen_R12_{ICD10_code}
# 需要找到对应的焦虑/抑郁和颈肩痛表型

# 查询FinnGen中可用表型的方法:
# 1) 使用 ieugwasr::gwas_search() 搜索
# 2) 访问 finngen.fi 直接查询
# 3) 使用已知ICD-10/FinnGen编码

# 已知FinnGen编码 (R12版本):
# - 抑郁: F32_DEPRESS (Depressive episode)
# - 焦虑: F5_ANXIETY (Anxiety disorders) 或 F41_OTHER_ANXIETY
# - 颈肩痛: M53_NECKPAIN (Cervicalgia) 或 M79_OTHERPAIN

# 使用API搜索
cat("搜索FinnGen中可用的焦虑/抑郁和颈肩痛表型...\n")

# 方法: 使用 ieugwasr 包搜索
tryCatch({
  search_results <- ieugwasr::gwas_search("anxiety")
  cat("可用焦虑相关表型数量:", nrow(search_results), "\n")
  print(search_results[1:min(5, nrow(search_results)), ])
}, error = function(e) {
  cat("API搜索失败:", e$message, "\n")
  cat("请手动访问 finngen.fi 查询可用表型ID\n")
})

# ==================== 2. 提取FinnGen暴露数据 ====================

# 假设已确认FinnGen表型ID
# 以下为示例ID，需根据实际查询结果替换

# FinnGen 焦虑 (示例)
finngen_anx_id <- "finngen_R12_F5_ANXIETY"
finngen_dep_id <- "finngen_R12_F32_DEPRESS"

# FinnGen 颈肩痛 (需在FinnGen中确认对应表型)
finngen_neckpain_id <- "finngen_R12_M53_NECKPAIN"  # 需确认

# --- 提取FinnGen焦虑数据 ---
exposure_finngen_anx <- tryCatch({
  extract_instruments(finngen_anx_id, p1 = 5e-8, p2 = 5e-8)
}, error = function(e) {
  cat("FinnGen焦虑表型提取失败:", e$message, "\n")
  return(NULL)
})

# --- 提取FinnGen抑郁数据 ---
exposure_finngen_dep <- tryCatch({
  extract_instruments(finngen_dep_id, p1 = 5e-8, p2 = 5e-8)
}, error = function(e) {
  cat("FinnGen抑郁表型提取失败:", e$message, "\n")
  return(NULL)
})

# --- 提取FinnGen结局数据 (颈肩痛) ---
outcome_finngen_neck <- tryCatch({
  extract_outcome_data(
    snps = c(exposure_finngen_anx$SNP, exposure_finngen_dep$SNP),
    outcomes = finngen_neckpain_id
  )
}, error = function(e) {
  cat("FinnGen颈肩痛表型提取失败:", e$message, "\n")
  cat("可考虑使用UKB颈肩痛作为结局，但暴露用FinnGen\n")
  return(NULL)
})

# ==================== 3. 替代方案: FinnGen暴露 + UKB结局 ====================

# 如果FinnGen中颈肩痛表型不可用，用FinnGen焦虑/抑郁作为暴露，
# UKB CNSP (ukb-b-16118) 作为结局
# 这仍提供独立验证，因为暴露数据来自不同人群

if(is.null(outcome_finngen_neck)) {
  cat("\n采用替代方案: FinnGen暴露 + UKB结局\n")
  cat("这提供了暴露层面的独立验证\n")
  
  outcome_ukb <- extract_outcome_data(
    snps = c(exposure_finngen_anx$SNP, exposure_finngen_dep$SNP),
    outcomes = "ukb-b-16118"
  )
}

# ==================== 4. LD Clumping (FinnGen暴露) ====================

if(!is.null(exposure_finngen_anx)) {
  exposure_finngen_anx_clumped <- clump_data(
    exposure_finngen_anx,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR"
  )
  cat("FinnGen焦虑工具变量数量:", nrow(exposure_finngen_anx_clumped), "\n")
}

if(!is.null(exposure_finngen_dep)) {
  exposure_finngen_dep_clumped <- clump_data(
    exposure_finngen_dep,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR"
  )
  cat("FinnGen抑郁工具变量数量:", nrow(exposure_finngen_dep_clumped), "\n")
}

# ==================== 5. Harmonize ====================

if(!is.null(exposure_finngen_anx) && !is.null(outcome_ukb)) {
  dat_finngen_anx <- harmonise_data(
    exposure_finngen_anx_clumped,
    outcome_ukb[outcome_ukb$SNP %in% exposure_finngen_anx_clumped$SNP, ]
  )
}

if(!is.null(exposure_finngen_dep) && !is.null(outcome_ukb)) {
  dat_finngen_dep <- harmonise_data(
    exposure_finngen_dep_clumped,
    outcome_ukb[outcome_ukb$SNP %in% exposure_finngen_dep_clumped$SNP, ]
  )
}

# ==================== 6. F-statistic检查 ====================

if(exists("dat_finngen_anx")) {
  dat_finngen_anx$F_stat <- (dat_finngen_anx$beta.exposure / dat_finngen_anx$se.exposure)^2
  cat("FinnGen焦虑: Mean F =", mean(dat_finngen_anx$F_stat), 
      "Min F =", min(dat_finngen_anx$F_stat), "\n")
}

if(exists("dat_finngen_dep")) {
  dat_finngen_dep$F_stat <- (dat_finngen_dep$beta.exposure / dat_finngen_dep$se.exposure)^2
  cat("FinnGen抑郁: Mean F =", mean(dat_finngen_dep$F_stat), 
      "Min F =", min(dat_finngen_dep$F_stat), "\n")
}

# ==================== 7. MR分析 (FinnGen验证) ====================

# --- FinnGen焦虑 → CNSP ---
if(exists("dat_finngen_anx")) {
  res_finngen_anx <- mr(
    dat_finngen_anx,
    method_list = c("mr_ivw_fe", "mr_weighted_median", 
                    "mr_egger_regression_bootstrap")
  )
  
  het_finngen_anx <- mr_heterogeneity(dat_finngen_anx)
  pleio_finngen_anx <- mr_pleiotropy_test(dat_finngen_anx)
  
  cat("\nFinnGen焦虑 → CNSP (UKB) MR结果:\n")
  anx_ivw <- res_finngen_anx[res_finngen_anx$method == "Inverse variance weighted (fixed effects)", ]
  cat("  IVW OR:", round(exp(anx_ivw$estimate), 3), 
      "P:", anx_ivw$pval, "\n")
}

# --- FinnGen抑郁 → CNSP ---
if(exists("dat_finngen_dep")) {
  res_finngen_dep <- mr(
    dat_finngen_dep,
    method_list = c("mr_ivw_fe", "mr_weighted_median", 
                    "mr_egger_regression_bootstrap")
  )
  
  het_finngen_dep <- mr_heterogeneity(dat_finngen_dep)
  pleio_finngen_dep <- mr_pleiotropy_test(dat_finngen_dep)
  
  cat("\nFinnGen抑郁 → CNSP (UKB) MR结果:\n")
  dep_ivw <- res_finngen_dep[res_finngen_dep$method == "Inverse variance weighted (fixed effects)", ]
  cat("  IVW OR:", round(exp(dep_ivw$estimate), 3), 
      "P:", dep_ivw$pval, "\n")
}

# ==================== 8. MR-PRESSO (FinnGen验证) ====================

if(exists("dat_finngen_anx")) {
  presso_finngen_anx <- mr_presso(
    BetaOutcome = "beta.outcome",
    BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome",
    SdExposure = "se.exposure",
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = dat_finngen_anx
  )
}

if(exists("dat_finngen_dep")) {
  presso_finngen_dep <- mr_presso(
    BetaOutcome = "beta.outcome",
    BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome",
    SdExposure = "se.exposure",
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = dat_finngen_dep
  )
}

# ==================== 9. 与原UKB结果对比 ====================

# 构建对比表
validation_comparison <- data.frame(
  Source = c("UKB (original)", "FinnGen anxiety → UKB CNSP", "FinnGen depression → UKB CNSP"),
  Exposure = c("GP-ANXD (ukb-b-6991)", "FinnGen Anxiety", "FinnGen Depression"),
  NSNP = c(192, 
           ifelse(exists("dat_finngen_anx"), nrow(dat_finngen_anx), NA),
           ifelse(exists("dat_finngen_dep"), nrow(dat_finngen_dep), NA)),
  IVW_OR = c(1.18, 
             ifelse(exists("res_finngen_anx"), round(exp(anx_ivw$estimate), 3), NA),
             ifelse(exists("res_finngen_dep"), round(exp(dep_ivw$estimate), 3), NA)),
  IVW_P = c(5.37e-7,
            ifelse(exists("res_finngen_anx"), anx_ivw$pval, NA),
            ifelse(exists("res_finngen_dep"), dep_ivw$pval, NA))
)

cat("\n========================================\n")
cat("UKB vs FinnGen验证结果对比\n")
cat("========================================\n")
print(validation_comparison)

# 判断: 结果是否可重现
if(exists("res_finngen_anx") && exists("res_finngen_dep")) {
  if((exp(anx_ivw$estimate) > 1 && anx_ivw$pval < 0.05) ||
     (exp(dep_ivw$estimate) > 1 && dep_ivw$pval < 0.05)) {
    cat("\n✅ FinnGen验证结果方向一致 → 结论具有外部位复制性\n")
  } else {
    cat("\n⚠️ FinnGen验证结果不一致 → 需深入调查异质性来源\n")
  }
}

# ==================== 10. 森林图: UKB vs FinnGen对比 ====================

# 如果FinnGen结果可用，生成对比森林图
if(exists("res_finngen_anx") || exists("res_finngen_dep")) {
  
  forest_dat <- data.frame(
    Dataset = character(),
    Phenotype = character(),
    OR = numeric(),
    lower = numeric(),
    upper = numeric()
  )
  
  # UKB original
  forest_dat <- rbind(forest_dat, data.frame(
    Dataset = "UK Biobank",
    Phenotype = "GP-ANXD",
    OR = 1.18,
    lower = 1.11,
    upper = 1.27
  ))
  
  # FinnGen anxiety
  if(exists("res_finngen_anx")) {
    forest_dat <- rbind(forest_dat, data.frame(
      Dataset = "FinnGen",
      Phenotype = "Anxiety",
      OR = exp(anx_ivw$estimate),
      lower = exp(anx_ivw$estimate - 1.96 * anx_ivw$se),
      upper = exp(anx_ivw$estimate + 1.96 * anx_ivw$se)
    ))
  }
  
  # FinnGen depression
  if(exists("res_finngen_dep")) {
    forest_dat <- rbind(forest_dat, data.frame(
      Dataset = "FinnGen",
      Phenotype = "Depression",
      OR = exp(dep_ivw$estimate),
      lower = exp(dep_ivw$estimate - 1.96 * dep_ivw$se),
      upper = exp(dep_ivw$estimate + 1.96 * dep_ivw$se)
    ))
  }
  
  ggplot(forest_dat, aes(x = OR, y = interaction(Phenotype, Dataset), 
                         color = Dataset)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), 
                   height = 0.2, position = position_dodge(0.5)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_x_log10() +
    scale_color_manual(values = c("#2C3E50", "#E74C3C")) +
    labs(
      title = "External Validation: UK Biobank vs FinnGen",
      subtitle = "Causal Effect of Mental Health on Chronic Neck/Shoulder Pain",
      x = "Odds Ratio (log scale)",
      y = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# ==================== 11. 保存结果 ====================

if(exists("res_finngen_anx")) {
  write.csv(res_finngen_anx, "FinnGen_anxiety_to_CNSP_MR_results.csv", row.names = FALSE)
}
if(exists("res_finngen_dep")) {
  write.csv(res_finngen_dep, "FinnGen_depression_to_CNSP_MR_results.csv", row.names = FALSE)
}
write.csv(validation_comparison, "validation_comparison_UKB_vs_FinnGen.csv", row.names = FALSE)

cat("\nFinnGen外部验证分析完成。\n")
cat("========================================\n")
cat("外部验证的意义:\n")
cat("========================================\n")
cat("编辑要求'independent replication'，FinnGen是最可行的选择:\n")
cat("  1. FinnGen独立于UKB，样本无重叠\n")
cat("  2. 与UKB同为欧洲人群(可比性高)\n")
cat("  3. 表型定义可能略有不同，但方向一致即支持结论稳定性\n")
cat("  4. 如果FinnGen颈肩痛表型也可用，则可做完全独立的MR\n")
cat("\n如果FinnGen数据不可用，可考虑:\n")
cat("  - 使用Neuroticism (UK Biobank field 20127) 作为替代暴露重复分析\n")
cat("  - 使用更大样本量的焦虑/抑郁GWAS (如23andMe数据) 作为验证\n")
cat("  - 在UKB内部做10-fold交叉验证 (虽非真正外部验证，但可展示稳定性)\n")
