###############################################################################
# 补充分析 2: 表型对比分析 (PGC-MDD vs GP-ANXD)
# 目的：回应编辑质疑——GP咨询表型是否适用于临床诊断焦虑/抑郁？
# 方法：用PGC-MDD临床诊断抑郁症GWAS做正向MR→CNSP，与GP-ANXD结果对比
# 数据来源：
#   - PGC-MDD: Howard et al. 2019, Nat Neurosci (807,553人)
#   - PGC-ANX: Otowa et al. 2022 (需确认公开可用性)
#   - FinnGen anxiety/depression (替代方案)
# 作者：MedResearch Pro 辅助生成
# 日期：2026-06-18
###############################################################################

# ==================== 0. 环境设置 ====================
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(ggplot2)

# ==================== 1. 获取PGC-MDD GWAS数据 ====================

# --- 方式A: 从IEU OpenGWAS获取PGC-MDD ---
# PGC-MDD在OpenGWAS中的ID (需确认最新版本):
# ieu-a-797 (PGC-MDD 2019, 23andMe补充版)
# ieu-a-300 (PGC-MDD 2018)
# pgc-mdd-2019-utrecht (需搜索确认)

# 尝试多个可能的ID
pgc_ids <- c("ieu-a-797", "ieu-a-300", "pgc-mdd-2019-utrecht")

for(id in pgc_ids) {
  tryCatch({
    exposure_pgc <- extract_instruments(id, p1 = 5e-8, p2 = 5e-8)
    if(nrow(exposure_pgc) > 5) {
      cat("找到PGC-MDD数据:", id, "SNPs数量:", nrow(exposure_pgc), "\n")
      break
    }
  }, error = function(e) {
    cat("ID", id, "不可用:", e$message, "\n")
  })
}

# --- 方式B: 从本地文件读取PGC-MDD ---
# 如果从PGC官网下载了GWAS summary statistics:
# pgc_url <- "https://pgc.unc.edu/"
# 下载后用 format_data() 处理
#
# pgc_gwas <- read.table("PGC_MDD_GWAS.txt", header = TRUE)
# exposure_pgc <- format_data(pgc_gwas,
#   type = "exposure",
#   snp_col = "SNP",
#   beta_col = "beta",
#   se_col = "se",
#   eaf_col = "eaf",
#   effect_allele_col = "effect_allele",
#   other_allele_col = "other_allele",
#   pval_col = "pval",
#   min_pval = 5e-8)

# ==================== 2. 获取PGC-ANX (焦虑症) GWAS数据 ====================

# PGC-ANX: Otowa et al. 2022
# IEU OpenGWAS可能的ID: ieu-a-302, pgc-anx-2022
# 如果PGC-ANX不可用，可使用FinnGen焦虑表型作为替代

# --- FinnGen替代 ---
# FinnGen R12: 焦虑相关表型
# 需从FinnGen下载对应summary statistics
# finngen_url <- "https://finngen.fi/"

# 尝试IEU OpenGWAS中的焦虑相关数据
anx_ids <- c("ieu-a-302", "finngen_R12_F5_ANXIETY", "finngen_R12_F32_DEPRESS")

for(id in anx_ids) {
  tryCatch({
    exposure_anx <- extract_instruments(id, p1 = 5e-8, p2 = 5e-8)
    if(nrow(exposure_anx) > 5) {
      cat("找到焦虑/抑郁GWAS数据:", id, "SNPs数量:", nrow(exposure_anx), "\n")
      break
    }
  }, error = function(e) {
    cat("ID", id, "不可用:", e$message, "\n")
  })
}

# ==================== 3. LD Clumping ====================

exposure_pgc_clumped <- clump_data(
  exposure_pgc,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pop = "EUR"
)

exposure_anx_clumped <- clump_data(
  exposure_anx,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pop = "EUR"
)

cat("PGC-MDD工具变量数量:", nrow(exposure_pgc_clumped), "\n")
cat("PGC-ANX工具变量数量:", nrow(exposure_anx_clumped), "\n")

# ==================== 4. 提取结局数据 ====================

outcome_dat_pgc <- extract_outcome_data(
  snps = exposure_pgc_clumped$SNP, 
  outcomes = "ukb-b-16118"  # CNSP
)

outcome_dat_anx <- extract_outcome_data(
  snps = exposure_anx_clumped$SNP,
  outcomes = "ukb-b-16118"  # CNSP
)

# ==================== 5. Harmonize ====================

dat_pgc <- harmonise_data(exposure_pgc_clumped, outcome_dat_pgc)
dat_anx <- harmonise_data(exposure_anx_clumped, outcome_dat_anx)

# ==================== 6. F-statistic检查 ====================

dat_pgc$F_stat <- (dat_pgc$beta.exposure / dat_pgc$se.exposure)^2
dat_anx$F_stat <- (dat_anx$beta.exposure / dat_anx$se.exposure)^2

cat("PGC-MDD: Mean F=", mean(dat_pgc$F_stat), "Min F=", min(dat_pgc$F_stat), "\n")
cat("PGC-ANX: Mean F=", mean(dat_anx$F_stat), "Min F=", min(dat_anx$F_stat), "\n")

# ==================== 7. MR分析 ====================

# --- PGC-MDD → CNSP ---
res_pgc <- mr(
  dat_pgc,
  method_list = c("mr_ivw_fe", "mr_weighted_median", 
                  "mr_egger_regression_bootstrap")
)

# --- PGC-ANX → CNSP ---
res_anx <- mr(
  dat_anx,
  method_list = c("mr_ivw_fe", "mr_weighted_median",
                  "mr_egger_regression_bootstrap")
)

# ==================== 8. MR-PRESSO ====================

presso_pgc <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = dat_pgc
)

presso_anx <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = dat_anx
)

# ==================== 9. 敏感性分析 ====================

# --- PGC-MDD ---
het_pgc <- mr_heterogeneity(dat_pgc)
pleio_pgc <- mr_pleiotropy_test(dat_pgc)
loo_pgc <- mr_leaveoneout(dat_pgc)

# --- PGC-ANX ---
het_anx <- mr_heterogeneity(dat_anx)
pleio_anx <- mr_pleiotropy_test(dat_anx)
loo_anx <- mr_leaveoneout(dat_anx)

# ==================== 10. Steiger方向性测试 ====================

steiger_pgc <- steiger(
  p_exp = dat_pgc$pval.exposure,
  p_out = dat_pgc$pval.outcome,
  r_exp = dat_pgc$beta.exposure^2 / 
    (dat_pgc$beta.exposure^2 + dat_pgc$se.exposure^2 * nrow(dat_pgc)),
  r_out = dat_pgc$beta.outcome^2 / 
    (dat_pgc$beta.outcome^2 + dat_pgc$se.outcome^2 * nrow(dat_pgc)),
  n_exp = dat_pgc$samplesize.exposure,
  n_out = dat_pgc$samplesize.outcome
)

steiger_anx <- steiger(
  p_exp = dat_anx$pval.exposure,
  p_out = dat_anx$pval.outcome,
  r_exp = dat_anx$beta.exposure^2 / 
    (dat_anx$beta.exposure^2 + dat_anx$se.exposure^2 * nrow(dat_anx)),
  r_out = dat_anx$beta.outcome^2 / 
    (dat_anx$beta.outcome^2 + dat_anx$se.exposure^2 * nrow(dat_anx)),
  n_exp = dat_anx$samplesize.exposure,
  n_out = dat_anx$samplesize.outcome
)

# ==================== 11. 遗传相关性计算 ====================

# 计算GP-ANXD与PGC-MDD之间的遗传相关性
# 使用LDSC (LD Score Regression)

# 注意：LDSC需要特定的输入格式
# 如果数据格式兼容，可直接用 ieugwasr::ldsc()
# 否则需单独运行LDSC软件

# 简化版：使用IEU OpenGWAS的在线遗传相关性计算
# 需确保两个数据集的ID正确

tryCatch({
  rg_anxd_mdd <- ldsc(
    trait1 = "ukb-b-6991",  # GP-ANXD
    trait2 = "ieu-a-797"    # PGC-MDD
  )
  cat("GP-ANXD与PGC-MDD遗传相关性:", rg_anxd_mdd$rg, "\n")
  cat("P-value:", rg_anxd_mdd$p, "\n")
}, error = function(e) {
  cat("LDSC在线计算失败，需手动运行LDSC软件\n")
})

# ==================== 12. 三种表型结果对比汇总 ====================

# 整理三种暴露表型的IVW结果
comparison_table <- data.frame(
  Phenotype = c("GP-ANXD (ukb-b-6991)", 
                "PGC-MDD (clinical diagnosis)", 
                "PGC-ANX (clinical diagnosis)"),
  Source = c("UK Biobank GP consultation", 
             "PGC consortium", 
             "PGC consortium / FinnGen"),
  NSNP = c(192, 
           nrow(dat_pgc[dat_pgc$F_stat > 10, ]),
           nrow(dat_anx[dat_anx$F_stat > 10, ])),
  Threshold = c("p<1e-5", "p<5e-8", "p<5e-8"),
  # IVW结果需从res_pgc和res_anx中提取
  IVW_OR = c(1.18, NA, NA),
  IVW_P = c(5.37e-7, NA, NA)
)

# 从PGC结果中填充
pgc_ivw <- res_pgc[res_pgc$method == "Inverse variance weighted (fixed effects)", ]
anx_ivw <- res_anx[res_anx$method == "Inverse variance weighted (fixed effects)", ]

comparison_table$IVW_OR[2] <- exp(pgc_ivw$estimate)
comparison_table$IVW_P[2] <- pgc_ivw$pval
comparison_table$IVW_OR[3] <- exp(anx_ivw$estimate)
comparison_table$IVW_P[3] <- anx_ivw$pval

cat("\n========================================\n")
cat("三种暴露表型对CNSP的因果效应对比\n")
cat("========================================\n")
print(comparison_table)

# 关键判断逻辑
if(all(comparison_table$IVW_P < 0.05) && 
   all(comparison_table$IVW_OR > 1)) {
  cat("\n✅ 三种表型均显示正向因果效应 → GP-ANXD结果可推广到临床诊断表型\n")
} else {
  cat("\n⚠️ 表型间存在差异 → 需在Discussion中详细讨论差异来源\n")
}

# ==================== 13. 可视化：三种表型森林图对比 ====================

# 构建森林图数据
forest_data <- data.frame(
  Phenotype = comparison_table$Phenotype,
  OR = comparison_table$IVW_OR,
  # 需从MR结果中提取CI
  lower = NA,
  upper = NA
)

# 从结果中提取CI (需根据实际MR结果补充)
# forest_data$lower[2] <- exp(pgc_ivw$estimate - 1.96 * pgc_ivw$se)
# forest_data$upper[2] <- exp(pgc_ivw$estimate + 1.96 * pgc_ivw$se)

ggplot(forest_data, aes(x = OR, y = Phenotype)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_x_log10() +
  labs(
    title = "Phenotype Comparison: Causal Effect on Chronic Neck/Shoulder Pain",
    x = "Odds Ratio (log scale)",
    y = "Exposure Phenotype"
  ) +
  theme_minimal()

# ==================== 14. 保存结果 ====================

write.csv(res_pgc, "PGC_MDD_to_CNSP_MR_results.csv", row.names = FALSE)
write.csv(res_anx, "PGC_ANX_to_CNSP_MR_results.csv", row.names = FALSE)
write.csv(comparison_table, "phenotype_comparison_summary.csv", row.names = FALSE)

cat("\n表型对比分析完成。请在论文中添加: \n")
cat("'Supplementary Analysis: MR Analysis Using Clinically Diagnosed Depression/Anxiety Phenotypes' \n")
cat("并在Discussion中详细讨论GP-ANXD与临床诊断表型的一致性与差异。\n")
