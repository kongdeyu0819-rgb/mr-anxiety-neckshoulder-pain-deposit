###############################################################################
# 补充分析 4: 逆向MR统计功效计算
# 目的：回应编辑质疑——"反向MR零结果可能是统计功效不足而非真正零效应"
# 方法：使用Brion et al. 2013的MR功效公式计算最小可检测效应量
# 核心判断：在当前工具变量数量/R²水平下，逆向MR能否检测到合理大小的效应？
# 作者：MedResearch Pro 辅助生成
# 日期：2026-06-18
###############################################################################

# ==================== 0. 环境设置 ====================
library(TwoSampleMR)
library(pwr)  # 统计功效计算
library(ggplot2)
library(dplyr)

# ==================== 1. MR功效计算原理 ====================

# Brion et al. (2013) "Power and sample size approximations for 
# Mendelian randomization studies using one genetic instrument"
# 
# 核心公式:
# 对于单SNP的MR功效计算:
#   Power = Φ( Z_α/2 + √(n * β² * R² / (1-R²)) ) - Φ( Z_α/2 - √(n * β² * R² / (1-R²)) )
# 其中:
#   n = 样本量 (结局GWAS)
#   β = 因果效应 (log OR)
#   R² = SNPs解释的暴露方差比例
#   α = I类错误率 (通常0.05)
#   Φ = 标准正态分布累积函数
#
# 对于多SNP的IVW MR:
#   近似功效 = Power(n_effective, β_effective, R²_effective)
#   n_effective ≈ N * β² / se²

# ==================== 2. 从实际MR结果中提取参数 ====================

# --- 逆向MR (CNSP → GP-ANXD) 的实际参数 ---
# 从你的论文结果中提取:
# - 工具变量数量: 约10-15个(从CNSP GWAS中提取)
# - 每个SNP的F值: 需从实际分析中计算
# - 总R²: 工具变量解释的CNSP方差比例
# - 样本量: CNSP GWAS (72,887 cases + 32,509 controls)

# 由于你未提供SNP水平的详细数据，以下用合理的近似值
# 从论文Table 2中提到的逆向MR结果可知:
# - 逆向MR不显著
# - 异质性显著 (Q p=0.014)
# - MR-Egger intercept不显著
# - 说明工具变量数量有限

# 读取正向MR结果中的SNP数据(如果保存了)
# snp_data <- read.csv("forward_MR_SNPs_data.csv")

# 如果没有保存，基于论文摘要估算:
N_cnsp_exposure <- 72887 + 32509  # CNSP GWAS总样本量
N_gpanxd_outcome <- 158565 + 300995  # GP-ANXD GWAS总样本量(作为逆向结局)

# 逆向MR参数估算
n_snps_reverse <- 10  # 假设CNSP GWAS在p<1e-5阈值下有10个独立SNPs
# 注意: 这需要从实际CNSP GWAS中提取后确认

# 每个SNP的R²估算(基于单个SNP的variance explained)
# 对于二分结局的GWAS, R² ≈ 2 * MAF * (1-MAF) * beta² / sqrt(N)
# 近似: 每个SNP的R² ≈ 0.001 - 0.01 (常见范围)

# 总R² (所有SNPs合计)
R2_total_reverse <- 0.005  # 保守估计: 所有SNPs共解释0.5%的方差
# 这个值需要从实际CNSP GWAS计算中确认

# ==================== 3. 精确功效计算函数 (Brion et al. 2013) ====================

# 单SNP MR功效
mr_power_single_snp <- function(n_outcome, beta, R2, alpha = 0.05) {
  # n_outcome: 结局GWAS样本量
  # beta: 目标因果效应 (log OR)
  # R2: 该SNP解释的暴露方差比例
  # alpha: I类错误率
  
  z_alpha <- qnorm(1 - alpha/2)
  
  # 非中心参数
  ncp <- sqrt(n_outcome * beta^2 * R2 / (1 - R2))
  
  # 功效
  power <- pnorm(ncp - z_alpha) + pnorm(-ncp - z_alpha, lower.tail = TRUE)
  
  return(power)
}

# 多SNP IVW MR功效 (近似)
mr_power_ivw <- function(n_outcome, beta, R2_total, n_snps, alpha = 0.05) {
  # R2_total: 所有SNPs合计解释的暴露方差
  # n_snps: 工具变量数量
  # 注: 多SNP的近似功效基于聚合样本量
  
  # 有效样本量
  n_eff <- n_outcome * R2_total
  
  z_alpha <- qnorm(1 - alpha/2)
  ncp <- sqrt(n_eff * beta^2)
  
  power <- pnorm(ncp - z_alpha) + pnorm(-ncp - z_alpha, lower.tail = TRUE)
  
  return(power)
}

# ==================== 4. 计算逆向MR的最小可检测OR ====================

# 给定功效阈值(如80%)，反推最小可检测效应量
detectable_or <- function(n_outcome, R2_total, n_snps, power_target = 0.8, alpha = 0.05) {
  # 使用二分法搜索
  log_or_low <- 0
  log_or_high <- 2  # 对应OR ~7.4
  log_or_mid <- (log_or_low + log_or_high) / 2
  
  while((log_or_high - log_or_low) > 0.001) {
    pwr <- mr_power_ivw(n_outcome, log_or_mid, R2_total, n_snps, alpha)
    if(pwr < power_target) {
      log_or_low <- log_or_mid
    } else {
      log_or_high <- log_or_mid
    }
    log_or_mid <- (log_or_low + log_or_high) / 2
  }
  
  return(exp(log_or_mid))
}

# ==================== 5. 逆向MR关键参数计算 ====================

# 从实际CNSP GWAS数据中提取工具变量(需替换为实际数据)
# 以下为近似计算——实际运行时应替换为你实际提取的SNPs

# 假设CNSP GWAS在p<1e-5阈值下有n个SNPs
# 并计算每个SNP的R²

# 如果需要从实际GWAS数据中提取，使用:
# cnsp_exposure <- extract_instruments("ukb-b-16118", p1 = 1e-5, p2 = 1e-5)
# cnsp_clumped <- clump_data(cnsp_exposure)
# cat("CNSP工具变量数量:", nrow(cnsp_clumped), "\n")
# 
# 计算每个SNP的R²
# cnsp_clumped$R2 <- 2 * cnsp_clumped$beta.exposure^2 * cnsp_clumped$eaf * (1-cnsp_clumped$eaf)

# 临时使用近似值:
n_snps_reverse_actual <- 12  # 需从实际CNSP GWAS确认
R2_total_reverse_actual <- 0.008  # 所有SNPs合计R² (需从实际计算确认)
N_outcome_reverse <- 158565 + 300995  # GP-ANXD GWAS总样本量

# 计算最小可检测OR (80%功效)
min_detectable_or <- detectable_or(
  n_outcome = N_outcome_reverse,
  R2_total = R2_total_reverse_actual,
  n_snps = n_snps_reverse_actual,
  power_target = 0.8
)

cat("========================================\n")
cat("逆向MR统计功效计算结果\n")
cat("========================================\n")
cat("工具变量数量:", n_snps_reverse_actual, "\n")
cat("总R² (所有SNPs合计):", R2_total_reverse_actual, "\n")
cat("结局GWAS样本量:", N_outcome_reverse, "\n")
cat("\n最小可检测OR (80%功效):", round(min_detectable_or, 3), "\n")
cat("最小可检测OR (50%功效):", 
    round(detectable_or(N_outcome_reverse, R2_total_reverse_actual, 
                        n_snps_reverse_actual, power_target = 0.5), 3), "\n")

# ==================== 6. 功效曲线可视化 ====================

# 构建不同OR值下的功效曲线
or_range <- seq(1.0, 2.0, by = 0.01)
power_curve <- sapply(log(or_range), function(log_or) {
  mr_power_ivw(N_outcome_reverse, log_or, R2_total_reverse_actual, 
               n_snps_reverse_actual)
})

power_df <- data.frame(
  OR = or_range,
  Power = power_curve
)

ggplot(power_df, aes(x = OR, y = Power)) +
  geom_line(color = "#E74C3C", linewidth = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = min_detectable_or, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = min_detectable_or + 0.05, y = 0.85, 
           label = paste("Minimum detectable OR =", round(min_detectable_or, 2)),
           hjust = 0) +
  labs(
    title = "Statistical Power of Reverse MR Analysis (CNSP → GP-ANXD)",
    subtitle = paste("Based on", n_snps_reverse_actual, "instruments, R² =", R2_total_reverse_actual),
    x = "Odds Ratio (causal effect size)",
    y = "Statistical Power"
  ) +
  ylim(0, 1) +
  theme_minimal()

# ==================== 7. 不同工具变量数量下的功效比较 ====================

# 展示: 如果CNSP GWAS样本量更大(如150K cases)，功效会如何？
# 这是为未来研究提供建议的量化依据

sample_range <- c(50000, 100000, 200000, 500000)
power_sample_compare <- data.frame()

for(N in sample_range) {
  for(or_val in seq(1.05, 1.5, by = 0.05)) {
    pwr_val <- mr_power_ivw(N, log(or_val), R2_total_reverse_actual, n_snps_reverse_actual)
    power_sample_compare <- rbind(power_sample_compare, data.frame(
      SampleSize = N,
      OR = or_val,
      Power = pwr_val
    ))
  }
}

ggplot(power_sample_compare, aes(x = OR, y = Power, color = factor(SampleSize))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_viridis_d(name = "Outcome GWAS\nSample Size") +
  labs(
    title = "Reverse MR Power: Impact of GWAS Sample Size",
    x = "Odds Ratio",
    y = "Statistical Power"
  ) +
  theme_minimal()

# ==================== 8. 敏感度分析结果解读 ====================

cat("\n========================================\n")
cat("功效计算结果解读\n")
cat("========================================\n")

if(min_detectable_or > 1.2) {
  cat("\n⚠️ 关键发现: 逆向MR的最小可检测OR >", round(min_detectable_or, 2), "\n")
  cat("   这意味着: 如果真实逆向效应OR <", round(min_detectable_or, 2), "\n")
  cat("   我们当前的分析无法检测到 (统计功效 < 80%)\n")
  cat("   因此，逆向MR的零结果很可能是统计功效不足导致的，\n")
  cat("   而非真正的零效应。\n")
  cat("\n   建议在论文中修改结论表述:\n")
  cat("   'The non-significant reverse MR finding should be interpreted \n")
  cat("    as inconclusive rather than evidence against reverse causation, \n")
  cat("    given that the current analysis had limited power to detect \n")
  cat("    inverse effects smaller than OR = ", round(min_detectable_or, 2), ".'\n", sep = "")
  
} else if(min_detectable_or > 1.1) {
  cat("\n✅ 逆向MR的功效尚可: 最小可检测OR =", round(min_detectable_or, 2), "\n")
  cat("   这意味着: 如果真实效应OR >", round(min_detectable_or, 2), "，我们的分析有80%功效检测\n")
  cat("   但较小效应(OR <", round(min_detectable_or, 2), ")仍可能未被检测到\n")
  cat("   结论表述需谨慎，不宜断言'完全无逆向效应'\n")
  
} else {
  cat("\n✅ 逆向MR功效充足: 最小可检测OR =", round(min_detectable_or, 2), "\n")
  cat("   这意味着: 即使较小的逆向效应也能被检测到\n")
  cat("   在此情况下，零结果更可能是真实零效应\n")
  cat("   但仍需讨论异质性(Cochran's Q p=0.014)的影响\n")
}

# ==================== 9. 定量讨论: 在论文中如何表述 ====================

cat("\n========================================\n")
cat("论文修改建议: 量化逆向MR功效不足\n")
cat("========================================\n")

cat("\n在Discussion的逆向分析段落中，添加以下定量讨论:\n")
cat("\n---\n")
cat("To quantify the power of the reverse MR analysis, we estimated the \n")
cat("minimum detectable causal effect using the framework of Brion et al. \n")
cat("(2013). Based on the", n_snps_reverse_actual, "instruments with a \n")
cat("combined variance explained (R²) of", R2_total_reverse_actual, "in the \n")
cat("exposure (CNSP), the reverse MR analysis had 80% power to detect an \n")
cat("OR of", round(min_detectable_or, 2), "or greater. The non-significant \n")
cat("reverse finding therefore does not definitively exclude small-to-\n")
cat("moderate reverse effects (OR <", round(min_detectable_or, 2), "). This \n")
cat("limitation should be considered when interpreting the apparent \n")
cat("unidirectional pattern. Future reverse MR studies using larger CNSP \n")
cat("GWAS datasets (e.g., with >150,000 cases) are needed to achieve \n")
cat("sufficient power to conclusively evaluate the reverse causal \n")
cat("hypothesis.\n")
cat("---\n")

# ==================== 10. 保存结果 ====================

power_results <- data.frame(
  Analysis = "Reverse MR (CNSP → GP-ANXD)",
  N_SNPs = n_snps_reverse_actual,
  R2_total = R2_total_reverse_actual,
  Outcome_sample_size = N_outcome_reverse,
  Min_detectable_OR_50power = detectable_or(N_outcome_reverse, R2_total_reverse_actual, n_snps_reverse_actual, 0.5),
  Min_detectable_OR_80power = min_detectable_or,
  Interpretation = ifelse(min_detectable_or > 1.2, "Power limited, null finding inconclusive", "Power adequate")
)

write.csv(power_results, "reverse_MR_power_analysis.csv", row.names = FALSE)
write.csv(power_df, "reverse_MR_power_curve_data.csv", row.names = FALSE)

cat("\n功效分析完成。\n")
cat("关键输出文件:\n")
cat("  - reverse_MR_power_analysis.csv: 数字摘要\n")
cat("  - reverse_MR_power_curve_data.csv: 功效曲线数据\n")
cat("  - 功效曲线图 (由ggplot生成)\n")
cat("\n请将功效计算结果写入论文的Supplementary Material中。\n")
