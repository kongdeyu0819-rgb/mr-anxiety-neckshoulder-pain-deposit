# GitHub + Zenodo DOI 操作指引（需您本人操作的部分已标注）

本机探测结果：git 2.54 已安装，但 **gh CLI 未安装、无 GitHub 凭据、无全局 git 身份配置**。因此创建仓库、推送、授权 Zenodo 三步需要您本人操作（约 15 分钟）。其余文件已全部备好。

---

## 第 0 步：一次性准备（二选一）

**方式 A（推荐，纯网页，无需命令行）**：直接在 GitHub 网页上传文件夹。跳过第 1、2 步的 git 命令，用第 1 步的"网页上传替代"操作。

**方式 B（命令行）**：先配置 git 身份（在任意终端执行，替换为您的 GitHub 用户名和邮箱）：

```
git config --global user.name "你的GitHub用户名"
git config --global user.email "你的GitHub注册邮箱"
```

---

## 第 1 步：创建 GitHub 仓库并上传

1. 登录 https://github.com → 右上角 **+** → **New repository**。
2. Repository name 建议：`mr-anxiety-neckshoulder-pain`（可见性必须选 **Public**，否则 Zenodo 无法存档）。不要勾选 README/LICENSE 初始化（已备好）。
3. 点击 **Create repository**。

**网页上传替代（方式 A）**：在新建的空仓库页面点击 **uploading an existing file**，把本文件夹 `repository_deposit` 里的内容（code 文件夹、results 文件夹、README.md、LICENSE、.zenodo.json）逐个拖入页面。GitHub 网页上传不支持拖入整个文件夹树时，需分两次：先拖 code 和 results 文件夹，再拖三个文件。提交信息填 `Initial deposit: analysis code and results` → **Commit changes**。

**命令行替代（方式 B）**：

```
cd "D:\workburry数据\AGENT\医学科研专家\MR_NeckPain_Anxiety\PRI大修_2026-10\repository_deposit"
git init -b main
git add .
git commit -m "Initial deposit: analysis code and results"
git remote add origin https://github.com/<你的用户名>/mr-anxiety-neckshoulder-pain.git
git push -u origin main
```

推送时会弹出浏览器要求登录 GitHub 授权一次。

---

## 第 2 步：给仓库打版本号（Zenodo 只存档正式 Release）

网页上：仓库页 → **Releases** → **Create a new release** → **Choose a tag** 输入 `v1.0.0` → Release title 填 `Version 1.0.0 (manuscript submission)` → **Publish release**。

---

## 第 3 步：Zenodo 获取 DOI

1. 打开 https://zenodo.org → 用 **GitHub 账号登录**（Log in with GitHub）。
2. 首次登录后访问 https://zenodo.org/account/settings/github/ ：找到你的仓库 `mr-anxiety-neckshoulder-pain`，把开关拨到 **On**（开启存档）。
3. Zenodo 只存档 **打了 Release 的提交**：回到 GitHub 把第 2 步的 Release 发布后，Zenodo 会自动在几分钟内生成一条存档记录。
4. 在 https://zenodo.org/account/settings/github/ 对应记录处可看到两个 DOI：**概念 DOI**（代表所有版本）和**版本 DOI**（代表本版本）。论文中引用**版本 DOI**。

---

## 第 4 步：回填三处占位符（DOI 拿到后）

拿到 DOI（形如 `10.5281/zenodo.17182901`）后，需要回填以下三处，然后**重跑 build_docx.py 重新生成 docx**（或直接告知我，我来回填并重建）：

| 文件 | 位置 | 占位符 |
|---|---|---|
| Revised_Manuscript_v11.md | Declarations → Availability of data and materials | `[GitHub repository URL - to be added upon deposit]` 与 `[Zenodo DOI - to be added upon deposit]` |
| Cover_Letter_PRI_Revision.md | 正文提及数据可得性处 | 同上（如该文件无占位符则无需改） |
| .zenodo.json | creators[0].orcid | `0000-0000-0000-0000`（换成真实 ORCID，如有） |

---

## 安全检查记录（已执行）

- `code/` 仅收录 8 个核心分析脚本；**test_\*、verify_\*、debug_\*、check_\*、cherry_\* 类文件因含硬编码 API token/密钥，已全部排除**。
- `results/` 仅含派生结果 CSV/JSON，无个体级数据。
- 原始 FinnGen gz 文件（约 800 MB/个）不入库，README 已写明下载途径。

## 需要您提供的输入汇总

1. GitHub 操作（第 1–3 步，需您的账号）。
2. DOI 与仓库 URL 回填（可交给我执行）。
3. ORCID（如有；无则保留占位或从 Declarations 删除）。
