cat("=== 安装 MR 分析核心包 (R 4.6.0) ===\n")

# 设置 CRAN 镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 安装 remotes 包（用于 GitHub 安装）
cat("1. 安装 remotes 包...\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  cat("  ✅ remotes 安装成功\n")
} else {
  cat("  ✅ remotes 已安装\n")
}

# 安装基础依赖
cat("\n2. 安装基础依赖包...\n")
basic_pkgs <- c("ggplot2", "dplyr", "data.table", "meta", "jsonlite", "httr",
                "curl", "tidyr", "rlang", "vctrs", "stringr", "purrr")
for (pkg in basic_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("   安装 %s...\n", pkg))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      cat(sprintf("   ✅ %s 安装成功\n", pkg))
    }, error = function(e) {
      cat(sprintf("   ❌ %s 安装失败: %s\n", pkg, e$message))
    })
  } else {
    cat(sprintf("   ✅ %s 已安装\n", pkg))
  }
}

# 安装 ieugwasr (OpenGWAS API 接口)
cat("\n3. 安装 ieugwasr...\n")
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  tryCatch({
    install.packages("ieugwasr", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"), quiet = TRUE)
    cat("  ✅ ieugwasr 安装成功\n")
  }, error = function(e) {
    cat(sprintf("  ❌ CRAN安装失败: %s\n", e$message))
    cat("  尝试从 GitHub 安装...\n")
    tryCatch({
      remotes::install_github("MRCIEU/ieugwasr", quiet = TRUE, upgrade = "never")
      cat("  ✅ ieugwasr 从 GitHub 安装成功\n")
    }, error = function(e2) {
      cat(sprintf("  ❌ GitHub 安装也失败: %s\n", e2$message))
    })
  })
}

# 安装 TwoSampleMR
cat("\n4. 安装 TwoSampleMR...\n")
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  tryCatch({
    install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"), quiet = TRUE)
    cat("  ✅ TwoSampleMR 安装成功\n")
  }, error = function(e) {
    cat(sprintf("  ❌ CRAN安装失败: %s\n", e$message))
    cat("  尝试从 GitHub 安装...\n")
    tryCatch({
      remotes::install_github("MRCIEU/TwoSampleMR", quiet = TRUE, upgrade = "never")
      cat("  ✅ TwoSampleMR 从 GitHub 安装成功\n")
    }, error = function(e2) {
      cat(sprintf("  ❌ GitHub 安装也失败: %s\n", e2$message))
    })
  })
}

# 安装 MendelianRandomization
cat("\n5. 安装 MendelianRandomization...\n")
if (!requireNamespace("MendelianRandomization", quietly = TRUE)) {
  tryCatch({
    install.packages("MendelianRandomization", repos = "https://cloud.r-project.org", quiet = TRUE)
    cat("  ✅ MendelianRandomization 安装成功\n")
  }, error = function(e) {
    cat(sprintf("  ❌ MendelianRandomization 安装失败: %s\n", e$message))
  })
}

# 安装 MRPRESSO
cat("\n6. 检查 MRPRESSO...\n")
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  tryCatch({
    remotes::install_github("rondolab/MR-PRESSO", quiet = TRUE, upgrade = "never")
    cat("  ✅ MRPRESSO 安装成功\n")
  }, error = function(e) {
    cat(sprintf("  ❌ MRPRESSO 安装失败: %s\n", e$message))
  })
} else {
  cat("  ✅ MRPRESSO 已安装\n")
}

# 安装 MVMR 包
cat("\n7. 安装 MVMR (多变量MR)...\n")
if (!requireNamespace("MVMR", quietly = TRUE)) {
  tryCatch({
    remotes::install_github("WSpiller/MVMR", quiet = TRUE, upgrade = "never")
    cat("  ✅ MVMR 安装成功\n")
  }, error = function(e) {
    cat(sprintf("  ❌ MVMR 安装失败: %s\n", e$message))
  })
}

# 最终检查
cat("\n=== 安装后包检查 ===\n")
final_pkgs <- c("TwoSampleMR", "ieugwasr", "MendelianRandomization",
                "MRPRESSO", "MVMR", "meta", "ggplot2", "dplyr",
                "data.table", "jsonlite", "httr")
for (pkg in final_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- as.character(packageVersion(pkg))
    cat(sprintf("  ✅ %s (版本 %s)\n", pkg, ver))
  } else {
    cat(sprintf("  ❌ %s - 未安装\n", pkg))
  }
}

cat("\n=== 完成 ===\n")
