# GitHub + Zenodo DOI 操作指引

## 当前进度（2026-09-15 自动执行完成）

- ✅ Token 验证：`ghp_...`（classic，scope=repo，账号 kongdeyu0819-rgb）可用；另一枚 fine-grained token 无效（401）
- ✅ 仓库已创建并推送：https://github.com/kongdeyu0819-rgb/mr-anxiety-neckshoulder-pain-deposit （36 个文件，含 code/、results/、README、LICENSE、.zenodo.json）
- ✅ Release v1.0.0 已发布：https://github.com/kongdeyu0819-rgb/mr-anxiety-neckshoulder-pain-deposit/releases/tag/v1.0.0
- ✅ 稿件 Declarations 已回填 GitHub 仓库 URL
- ⬜ 仅剩 Zenodo DOI：见下方"剩余唯一一步"
- 说明：`METRN_IVDD_manuscript_repo` 为 METRN 论文项目（对应已存在的 `METRN_IVDD_transcriptomic_analysis` 库），与本篇稿件无关，未在本库使用

## 剩余唯一一步：Zenodo DOI（二选一）

**方式 A（网页，约 2 分钟）**：
1. 打开 https://zenodo.org → 右上角 **Log in** → 选 **Log in with GitHub**（用 kongdeyu0819-rgb 账号授权）
2. 首次登录后访问 https://zenodo.org/account/settings/github/
3. 找到 `mr-anxiety-neckshoulder-pain-deposit`，把开关拨到 **On**
4. Zenodo 检测到已发布的 v1.0.0 Release，几分钟内自动生成存档记录
5. 在该页面记录处查看**版本 DOI**（形如 10.5281/zenodo.XXXXXXX），把它发给我，我来回填稿件并重建 docx

**方式 B（给我一个 Zenodo token，全程免操作）**：
1. 登录 https://zenodo.org 后访问 https://zenodo.org/account/settings/applications/tokens/new/
2. 创建 token 时勾选 **deposit:actions** 和 **deposit:write**
3. 把 token 发给我，我通过 API 直接创建并发布 deposit，立即回传 DOI

---


## 安全检查记录（已执行）

- `code/` 仅收录 8 个核心分析脚本；**test_\*、verify_\*、debug_\*、check_\*、cherry_\* 类文件因含硬编码 API token/密钥，已全部排除**。
- `results/` 仅含派生结果 CSV/JSON，无个体级数据。
- 原始 FinnGen gz 文件（约 800 MB/个）不入库，README 已写明下载途径。

## 需要您提供的输入汇总

1. GitHub 操作（第 1–3 步，需您的账号）。
2. DOI 与仓库 URL 回填（可交给我执行）。
3. ORCID（如有；无则保留占位或从 Declarations 删除）。
