###############################################################################
# 补充分析 6: 中介MR (Mediation MR) — 追踪因果通路
# 目的：将机制讨论从"推测"升级为"半验证"，回应编辑对机制讨论的质疑
# 方法：两步MR (two-step MR) 检验"焦虑/抑郁 → 中介变量 → 颈肩痛"通路
# 中介变量选择：
#   1) 炎症标记物 (CRP, IL-6 available via GWAS)
#   2) BMI (强遗传工具，中介效应明确)
#   3) 睡眠 (失眠 GWAS available)
#   4) 身体活动 (UKB physical activity GWAS)
# 作者：MedResearch Pro 辅助生成
# 日期：2026-06-18
###############################################################################

# ==================== 0. 环境设置 ====================
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(ggplot2)
library(igraph)  # 路径图可视化

# ==================== 1. 中介MR原理 ====================

# 两步MR (two-step MR) 框架:
# 
# 假设: Exposure (X) → 中介 (M) → Outcome (Y)
# 
# Step 1: X → M (MR1): 检验X对M的因果效应
# Step 2: M → Y (MR2): 检验M对Y的因果效应
# 
# 额外: X → Y (MR total): 总效应 (已有)
# 
# 中介比例估计:
#   Proportion mediated = (β_XM * β_MY) / β_XY
#   其中:
#     β_XM = X→M的MR估计
#     β_MY = M→Y的MR估计
#     β_XY = X→Y的总效应MR估计
#
# 注意: 中介MR比标准MR更复杂，需要:
#   - 所有三个MR (X→M, M→Y, X→Y) 方向正确
#   - 中介工具变量与X的独立性
#   - 无关于M的多效性

# ==================== 2. 中介变量GWAS数据来源 ====================

# 1) CRP (C-reactive protein)
#    IEU OpenGWAS: ukb-b-20175 (CRP in UKB)
#    Available: 较多欧洲人群GWAS

# 2) BMI
#    IEU OpenGWAS: ieu-a-835 (BMI GIANT + UKB)
#    强工具变量 (n SNPs > 100 at p<5e-8)

# 3) 失眠 (Insomnia)
#    IEU OpenGWAS: ukb-b-5657 (Insomnia in UKB)
#    PGC insomnia: ieu-a-1093

# 4) 身体活动 (Physical activity)
#    IEU OpenGWAS: ukb-b-13790 (Moderate/vigorous PA)
#    Broad PA: ieu-a-903

# 5) 疼痛敏感性 (Pain sensitivity)
#    可能可用 UKB pain ratings

# ==========================================
# 中介MR Step 1: GP-ANXD → 中介变量
# ==========================================

# --- 2a. CRP作为中介 ---
cat("========================================\n")
cat("中介分析 1: CRP (炎症标记物)\n")
cat("========================================\n")

# 提取GP-ANXD SNPs
exposure_anxd <- extract_instruments("ukb-b-6991", p1 = 1e-5, p2 = 1e-5)
exposure_anxd_clumped <- clump_data(exposure_anxd)

# 提取这些SNPs在CRP中的效应
outcome_crp <- extract_outcome_data(
  snps = exposure_anxd_clumped$SNP,
  outcomes = "ukb-b-20175"  # CRP
)

# Harmonize
dat_x_to_m_crp <- harmonise_data(exposure_anxd_clumped, outcome_crp)

# MR: ANXD → CRP
res_x_to_m_crp <- mr(
  dat_x_to_m_crp,
  method_list = c("mr_ivw_fe", "mr_weighted_median")
)

crp_beta <- res_x_to_m_crp[res_x_to_m_crp$method == "Inverse variance weighted (fixed effects)", ]$estimate
cat("GP-ANXD → CRP: beta =", crp_beta, "\n")

# ==========================================
# 中介MR Step 2: 中介变量 → CNSP
# ==========================================

# 提取CRP SNPs
exposure_crp <- extract_instruments("ukb-b-20175", p1 = 5e-8, p2 = 5e-8)
exposure_crp_clumped <- clump_data(exposure_crp)

# 提取这些SNPs在CNSP中的效应
outcome_cnsp_for_crp <- extract_outcome_data(
  snps = exposure_crp_clumped$SNP,
  outcomes = "ukb-b-16118"  # CNSP
)

# Harmonize
dat_m_to_y_crp <- harmonise_data(exposure_crp_clumped, outcome_cnsp_for_crp)

# MR: CRP → CNSP
res_m_to_y_crp <- mr(
  dat_m_to_y_crp,
  method_list = c("mr_ivw_fe", "mr_weighted_median")
)

crp_to_cnsp_beta <- res_m_to_y_crp[res_m_to_y_crp$method == "Inverse variance weighted (fixed effects)", ]$estimate
cat("CRP → CNSP: beta =", crp_to_cnsp_beta, "\n")

# ==========================================
# 中介比例计算 (CRP)
# ==========================================

# 总效应: ANXD → CNSP (已有)
beta_xy <- log(1.18)  # 约0.165

# 直接效应检验 (需MVMR，此处简化)
# 中介比例 = (beta_XM * beta_MY) / beta_XY
beta_xm <- crp_beta
beta_my <- crp_to_cnsp_beta

# 需要标准化 (因为scale不同)
# 使用标准化系数或用IVW OR近似计算

# 使用OR近似 (更直观)
or_xm <- exp(beta_xm)  # X→M的OR
or_my <- exp(beta_my)  # M→Y的OR
or_xy <- 1.18          # X→Y的OR

# 中介比例 (log scale)
prop_mediated_crp <- (beta_xm * beta_my) / beta_xy
cat("\nCRP中介比例 (approx):", round(prop_mediated_crp * 100, 1), "%\n")

# ==================== 3. BMI作为中介 ====================

cat("\n========================================\n")
cat("中介分析 2: BMI\n")
cat("========================================\n")

# ANXD → BMI
exposure_anxd_bmi <- extract_instruments("ukb-b-6991", p1 = 1e-5, p2 = 1e-5)
exposure_anxd_clumped_bmi <- clump_data(exposure_anxd_bmi)

outcome_bmi <- extract_outcome_data(
  snps = exposure_anxd_clumped_bmi$SNP,
  outcomes = "ieu-a-835"  # BMI
)

dat_x_to_m_bmi <- harmonise_data(exposure_anxd_clumped_bmi, outcome_bmi)
res_x_to_m_bmi <- mr(dat_x_to_m_bmi)
beta_xm_bmi <- res_x_to_m_bmi[res_x_to_m_bmi$method == "Inverse variance weighted (fixed effects)", ]$estimate

# BMI → CNSP
exposure_bmi <- extract_instruments("ieu-a-835", p1 = 5e-8, p2 = 5e-8)
exposure_bmi_clumped <- clump_data(exposure_bmi)

outcome_cnsp_for_bmi <- extract_outcome_data(
  snps = exposure_bmi_clumped$SNP,
  outcomes = "ukb-b-16118"
)

dat_m_to_y_bmi <- harmonise_data(exposure_bmi_clumped, outcome_cnsp_for_bmi)
res_m_to_y_bmi <- mr(dat_m_to_y_bmi)
beta_my_bmi <- res_m_to_y_bmi[res_m_to_y_bmi$method == "Inverse variance weighted (fixed effects)", ]$estimate

prop_mediated_bmi <- (beta_xm_bmi * beta_my_bmi) / beta_xy
cat("BMI中介比例 (approx):", round(prop_mediated_bmi * 100, 1), "%\n")

# ==================== 4. 失眠作为中介 ====================

cat("\n========================================\n")
cat("中介分析 3: 失眠 (Insomnia)\n")
cat("========================================\n")

# ANXD → 失眠
outcome_insomnia <- extract_outcome_data(
  snps = exposure_anxd_clumped$SNP,
  outcomes = "ukb-b-5657"  # Insomnia
)

dat_x_to_m_insomnia <- harmonise_data(exposure_anxd_clumped, outcome_insomnia)
res_x_to_m_insomnia <- mr(dat_x_to_m_insomnia)
beta_xm_insomnia <- res_x_to_m_insomnia[res_x_to_m_insomnia$method == "Inverse variance weighted (fixed effects)", ]$estimate

# 失眠 → CNSP
exposure_insomnia <- extract_instruments("ukb-b-5657", p1 = 5e-8, p2 = 5e-8)
exposure_insomnia_clumped <- clump_data(exposure_insomnia)

outcome_cnsp_for_insomnia <- extract_outcome_data(
  snps = exposure_insomnia_clumped$SNP,
  outcomes = "ukb-b-16118"
)

dat_m_to_y_insomnia <- harmonise_data(exposure_insomnia_clumped, outcome_cnsp_for_insomnia)
res_m_to_y_insomnia <- mr(dat_m_to_y_insomnia)
beta_my_insomnia <- res_m_to_y_insomnia[res_m_to_y_insomnia$method == "Inverse variance weighted (fixed effects)", ]$estimate

prop_mediated_insomnia <- (beta_xm_insomnia * beta_my_insomnia) / beta_xy
cat("失眠中介比例 (approx):", round(prop_mediated_insomnia * 100, 1), "%\n")

# ==================== 5. 中介结果汇总 ====================

mediation_summary <- data.frame(
  Mediator = c("CRP (inflammation)", "BMI", "Insomnia"),
  X_to_M_OR = c(exp(beta_xm), exp(beta_xm_bmi), exp(beta_xm_insomnia)),
  X_to_M_P = c(
    res_x_to_m_crp[res_x_to_m_crp$method == "Inverse variance weighted (fixed effects)", ]$pval,
    res_x_to_m_bmi[res_x_to_m_bmi$method == "Inverse variance weighted (fixed effects)", ]$pval,
    res_x_to_m_insomnia[res_x_to_m_insomnia$method == "Inverse variance weighted (fixed effects)", ]$pval
  ),
  M_to_Y_OR = c(exp(beta_my), exp(beta_my_bmi), exp(beta_my_insomnia)),
  M_to_Y_P = c(
    res_m_to_y_crp[res_m_to_y_crp$method == "Inverse variance weighted (fixed effects)", ]$pval,
    res_m_to_y_bmi[res_m_to_y_bmi$method == "Inverse variance weighted (fixed effects)", ]$pval,
    res_m_to_y_insomnia[res_m_to_y_insomnia$method == "Inverse variance weighted (fixed effects)", ]$pval
  ),
  Proportion_mediated = c(prop_mediated_crp, prop_mediated_bmi, prop_mediated_insomnia) * 100
)

cat("\n========================================\n")
cat("中介MR结果汇总\n")
cat("========================================\n")
print(mediation_summary)

# ==================== 6. 中介结果解读 ====================

cat("\n========================================\n")
cat("中介MR结果解读\n")
cat("========================================\n")

# 判断各中介通路是否显著
for(i in 1:nrow(mediation_summary)) {
  med <- mediation_summary$Mediator[i]
  p_xm <- mediation_summary$X_to_M_P[i]
  p_my <- mediation_summary$M_to_Y_P[i]
  prop <- mediation_summary$Proportion_mediated[i]
  
  cat("\n中介变量:", med, "\n")
  
  if(p_xm < 0.05 && p_my < 0.05) {
    cat("  ✅ X→M 和 M→Y 均显著 → 中介效应存在\n")
    cat("   中介比例约:", round(prop, 1), "%\n")
    cat("   这是机制讨论的实证基础，可写入Discussion\n")
  } else if(p_xm < 0.05 && p_my >= 0.05) {
    cat("  ⚠️ X→M显著，但M→Y不显著 → 中介效应不支持\n")
    cat("   但M→Y可能因统计功效不足不显著，需谨慎解读\n")
  } else if(p_xm >= 0.05 && p_my < 0.05) {
    cat("  ⚠️ X→M不显著，但M→Y显著 → 不支持中介\n")
  } else {
    cat("  ❌ X→M和M→Y均不显著 → 不支持中介\n")
  }
}

# ==================== 7. 限制性与敏感性讨论 ====================

cat("\n========================================\n")
cat("中介MR的局限性 (须在论文中说明)\n")
cat("========================================\n")
cat("1. 两步MR假设无水平多效性，即M的工具变量\n")
cat("   仅通过M影响Y，不直接影响Y\n")
cat("2. 中介MR比标准MR需要更强的工具变量\n")
cat("3. 比例中介估计需要MVMR才能得到直接效应\n")
cat("4. 更复杂的因果结构(如多中介、时序依赖)需要\n")
cat("   结构方程模型或多变量MR框架\n")
cat("5. 本分析仅提供初步证据，非机制验证\n")

# ==================== 8. 中介路径图可视化 ====================

# 创建因果路径图
if(any(mediation_summary$X_to_M_P < 0.05 & mediation_summary$M_to_Y_P < 0.05)) {
  cat("\n生成中介路径图...\n")
  
  # 使用igraph绘制路径图
  # 此处生成简化的文字描述，实际绘图可用igraph或ggplot
  cat("路径图应包含:\n")
  cat("  - X (GP-ANXD) → M (中介) → Y (CNSP)\n")
  cat("  - 标注各路径的OR和P值\n")
  cat("  - 用箭头粗细表示效应大小\n")
}

# ==================== 9. 保存结果 ====================

write.csv(mediation_summary, "mediation_MR_summary.csv", row.names = FALSE)

# 保存各中介MR详细结果
if(exists("res_x_to_m_crp")) {
  write.csv(res_x_to_m_crp, "mediation_ANXD_to_CRP_MR.csv", row.names = FALSE)
  write.csv(res_m_to_y_crp, "mediation_CRP_to_CNSP_MR.csv", row.names = FALSE)
}
if(exists("res_x_to_m_bmi")) {
  write.csv(res_x_to_m_bmi, "mediation_ANXD_to_BMI_MR.csv", row.names = FALSE)
  write.csv(res_m_to_y_bmi, "mediation_BMI_to_CNSP_MR.csv", row.names = FALSE)
}
if(exists("res_x_to_m_insomnia")) {
  write.csv(res_x_to_m_insomnia, "mediation_ANXD_to_Insomnia_MR.csv", row.names = FALSE)
  write.csv(res_m_to_y_insomnia, "mediation_Insomnia_to_CNSP_MR.csv", row.names = FALSE)
}

cat("\n中介MR分析完成。\n")
cat("========================================\n")
cat("中介MR在论文中的呈现方式:\n")
cat("========================================\n")
cat("1. 新增Results小节: 'Mediation Mendelian Randomization Analysis'\n")
cat("2. 展示各中介通路的MR结果(OR, 95% CI, P)\n")
cat("3. 用路径图(Figure)展示显著的中介通路\n")
cat("4. 在Discussion中，将机制讨论从'假设'升级为'MR支持的潜在通路'\n")
cat("5. 明确标注中介MR的局限性(见上文)\n")
cat("\n注意: 中介MR比标准MR更复杂，建议:\n")
cat("  - 咨询该领域专家或统计学家\n")
cat("  - 考虑使用更严谨的方法如MVMR-based mediation\n")
cat("  - 如果结果不支持中介，可改为'Exploratory Mediation Analysis'\n")
