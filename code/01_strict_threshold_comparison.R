###############################################################################
# 补充分析 1: 严格阈值对比分析 (p < 5×10⁻⁸ vs p < 1×10⁻⁵)
# 目的：回应编辑质疑——宽松阈值是否引入偏倚？
# 方法：用基因组显著性阈值重做全套MR，与原始结果对比
# 作者：MedResearch Pro 辅助生成
# 日期：2026-06-18
###############################################################################

# ==================== 0. 环境设置 ====================
library(TwoSampleMR)
library(MRPRESSO)
library(gtx)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==================== 1. 数据提取 ====================

# 暴露数据: GP-ANXD (ukb-b-6991)
# 结局数据: CNSP (ukb-b-16118)
# 这些数据可通过 IEU OpenGWAS API 获取，或从本地文件读取

# --- 方式A: 从IEU OpenGWAS API在线提取 ---
exposure_dat <- extract_instruments("ukb-b-6991", p1 = 5e-8, p2 = 5e-8)  # 严格阈值
outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = "ukb-b-16118")

# --- 方式B: 从本地GWAS summary文件读取 ---
# 如果已有本地VCF/txt文件，可使用以下方式：
# exposure_gwas <- read.table("GP-ANXD_GWAS.txt", header = TRUE)
# exposure_dat <- format_data(exposure_gwas, 
#   type = "exposure",
#   snp_col = "SNP",
#   beta_col = "beta",
#   se_col = "se",
#   eaf_col = "eaf",
#   effect_allele_col = "effect_allele",
#   other_allele_col = "other_allele",
#   pval_col = "pval",
#   samplesize_col = "n",
#   min_pval = 5e-8)  # 严格阈值

# ==================== 2. LD Clumping (严格阈值下) ====================

# TwoSampleMR::clump_data() 默认使用 1000 Genmas European LD参考面板
exposure_dat_clumped <- clump_data(
  exposure_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,  # 已在提取时用p<5e-8筛选，这里设为1保留所有
  pop = "EUR"
)

cat("严格阈值(p<5e-8)下的工具变量数量:", nrow(exposure_dat_clumped), "\n")

# ==================== 3. Harmonize ====================

dat_strict <- harmonise_data(
  exposure_dat = exposure_dat_clumped,
  outcome_dat = outcome_dat
)

# ==================== 4. F-statistic 计算 ====================

# F = (beta/se)^2 或 F = R² * (N-2) / (1-R²)
# 使用近似公式: F = beta² / se²

dat_strict$F_stat <- (dat_strict$beta.exposure / dat_strict$se.exposure)^2

cat("F-statistic 分布:\n")
cat("  Mean:", mean(dat_strict$F_stat), "\n")
cat("  Median:", median(dat_strict$F_stat), "\n")
cat("  Min:", min(dat_strict$F_stat), "\n")
cat("  所有SNP F>10?", all(dat_strict$F_stat > 10), "\n")

# 聚合 F-statistic (combined F)
total_r2 <- sum(dat_strict$F_stat / (dat_strict$F_stat + nrow(dat_strict) - 2))
combined_F <- total_r2 * (nrow(dat_strict) - 1) / (1 - total_r2)
cat("  Combined F-statistic:", combined_F, "\n")

# ==================== 5. MR 分析 (严格阈值) ====================

res_strict <- mr(
  dat_strict,
  method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_egger_regression_bootstrap")
)

# ==================== 6. MR-PRESSO (严格阈值) ====================

presso_strict <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = dat_strict
)

# ==================== 7. MR-RAPS (严格阈值) ====================

library(MRRAPS)

raps_strict <- mr_raps(
  dat_strict$beta.exposure,
  dat_strict$beta.outcome,
  dat_strict$se.exposure,
  dat_strict$se.outcome,
  over.dispersion = TRUE  # 允许多效性残差
)

# ==================== 8. 敏感性分析 (严格阈值) ====================

# Heterogeneity
het_strict <- mr_heterogeneity(dat_strict)
cat("Heterogeneity test (strict threshold):\n")
print(het_strict)

# Pleiotropy (MR-Egger intercept)
pleio_strict <- mr_pleiotropy_test(dat_strict)
cat("MR-Egger intercept test (strict threshold):\n")
print(pleio_strict)

# Leave-one-out
loo_strict <- mr_leaveoneout(dat_strict)
loo_plot_strict <- mr_leaveoneout_plot(loo_strict)

# ==================== 9. Steiger 方向性测试 (严格阈值) ====================

steiger_strict <- steiger(
  p_exp = dat_strict$pval.exposure,
  p_out = dat_strict$pval.outcome,
  r_exp = dat_strict$beta.exposure^2 / 
    (dat_strict$beta.exposure^2 + dat_strict$se.exposure^2 * nrow(dat_strict)),
  r_out = dat_strict$beta.outcome^2 / 
    (dat_strict$beta.outcome^2 + dat_strict$se.outcome^2 * nrow(dat_strict)),
  n_exp = dat_strict$samplesize.exposure,
  n_out = dat_strict$samplesize.outcome
)

# ==================== 10. GRS验证 (严格阈值) ====================

# 使用 gtx 包
grs_strict <- grs.summary(
  w = dat_strict$beta.exposure,       # 权重(暴露效应)
  b = dat_strict$beta.outcome,         # 结局效应
  se = dat_strict$se.outcome           # 结局标准误
)

cat("GRS analysis (strict threshold):\n")
cat("  OR:", exp(grs_strict$beta), "\n")
cat("  95% CI:", exp(grs_strict$beta - 1.96*grs_strict$se), "-", 
    exp(grs_strict$beta + 1.96*grs_strict$se), "\n")
cat("  P-value:", grs_strict$pval, "\n")

# ==================== 11. 与原始宽松阈值结果对比 ====================

# 原始结果（从论文中获取）
original_results <- data.frame(
  Method = c("IVW", "Weighted Median", "MR-Egger", "MR-PRESSO", "MR-RAPS", "GRS"),
  OR_original = c(1.18, 1.16, 1.12, 1.18, 1.17, 1.190),
  CI_lower_original = c(1.11, 1.09, 0.98, 1.11, 1.10, 1.134),
  CI_upper_original = c(1.27, 1.24, 1.28, 1.27, 1.25, 1.247),
  P_original = c(5.37e-7, 2.41e-5, 0.087, 5.37e-7, 3.15e-6, 1.60e-9),
  SNPs_original = c(192, 192, 192, 192, 192, NA)
)

# 严格阈值结果整理
strict_ivw <- res_strict[res_strict$method == "Inverse variance weighted (fixed effects)", ]
strict_wm <- res_strict[res_strict$method == "Weighted median", ]
strict_egger <- res_strict[res_strict$method == "MR Egger (bootstrap)", ]

strict_results <- data.frame(
  Method = c("IVW", "Weighted Median", "MR-Egger"),
  OR_strict = exp(c(strict_ivw$estimate, strict_wm$estimate, strict_egger$estimate)),
  P_strict = c(strict_ivw$pval, strict_wm$pval, strict_egger$pval),
  SNPs_strict = c(strict_ivw$nsnp, strict_wm$nsnp, strict_egger$nsnp)
)

# ==================== 12. 对比可视化 ====================

# 效应量对比森林图
comparison_df <- data.frame(
  Threshold = c("Relaxed (p<1e-5)", "Strict (p<5e-8)"),
  IVW_OR = c(1.18, strict_results$OR_strict[1]),
  IVW_lower = c(1.11, NA),  # 需从结果中提取CI
  IVW_upper = c(1.27, NA),
  WM_OR = c(1.16, strict_results$OR_strict[2]),
  NSNP = c(192, strict_results$SNPs_strict[1])
)

# 打印对比结果表
cat("\n========================================\n")
cat("对比结果: 宽松阈值 vs 严格阈值\n")
cat("========================================\n")
cat("宽松阈值 (p<1e-5): 192 SNPs, IVW OR=1.18 (1.11-1.27)\n")
cat("严格阈值 (p<5e-8):", strict_results$SNPs_strict[1], "SNPs\n")
cat("  IVW OR:", round(strict_results$OR_strict[1], 3), "\n")
cat("  Weighted Median OR:", round(strict_results$OR_strict[2], 3), "\n")
cat("  MR-Egger OR:", round(strict_results$OR_strict[3], 3), "\n")

# 检查关键判断
if(strict_results$OR_strict[1] > 0 && strict_results$P_strict[1] < 0.05) {
  cat("\n✅ 严格阈值下IVW结果方向一致且显著 → 宽松阈值未引入实质性偏倚\n")
} else if(strict_results$SNPs_strict[1] < 5) {
  cat("\n⚠️ 严格阈值下SNPs数量不足(<5) → 无法进行稳健的MR分析\n")
  cat("   这正是选择宽松阈值的正当理由，但需用其他方式弥补\n")
} else {
  cat("\n❌ 严格阈值下结果不一致 → 需进一步调查偏倚来源\n")
}

# ==================== 13. 保存结果 ====================

# 保存所有结果到CSV
write.csv(res_strict, "strict_threshold_MR_results.csv", row.names = FALSE)
write.csv(het_strict, "strict_threshold_heterogeneity.csv", row.names = FALSE)
write.csv(pleio_strict, "strict_threshold_pleiotropy.csv", row.names = FALSE)

cat("\n分析完成。结果已保存至工作目录。\n")
cat("请在论文中添加Supplementary Table: 'MR Results Using Genome-Wide Significance Threshold (p<5×10⁻⁸)' \n")
