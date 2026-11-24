# GitHub + Zenodo 发布步骤（代码 DOI）

代码包已整理在 `06_投稿材料\Code_MGG\`，路径、姓名、邮箱均已清理；本机未安装 git，也未登录 GitHub，因此仓库与 DOI 需要你本人操作。

## 一、GitHub 建仓（约 5 分钟）

1. 新建仓库：账号 → New repository → 名称建议 `FECD-disulfidptosis-transcriptome` → **Public** → 不勾选 README/GitHub，直接 Create。
2. 上传：网页版 “uploading an existing file”，把 `Code_MGG` 里的 `scripts/`、`README.md`、`LICENSE`、`sessionInfo.txt`、`zenodo_metadata.json`、`.gitignore` 全部拖入；提交信息写 `v1.0.0 analysis code`。
3. 记录仓库地址，例如 `https://github.com/<账号>/FECD-disulfidptosis-transcriptome`。

## 二、Zenodo 取 DOI（约 10 分钟）

1. 用 GitHub 账号登录 <https://zenodo.org>（Settings → GitHub 授权）。
2. 在 GitHub 仓库页 → Settings → Integrations → Zenodo → 打开开关；回到 Zenodo → GitHub 面板 → 找到该仓库 → **Enable**。
3. 在 GitHub 仓库 **Releases** → Draft a new release → tag `v1.0.0` → Publish release。
4. Zenodo 自动归档并生成 DOI（形如 `10.5281/zenodo.xxxxxxx`）。把 DOI 与仓库地址发我，我写入正文 Code Availability 与 Cover Letter。

> 注意：Zenodo 归档的是 release 的快照，所以必须先发布 release；仅 push 代码不会触发归档。

## 三、需要回填到稿件的两处

- `MGG_投稿稿_全文.docx` → Statements and Declarations → Code availability。
- `CoverLetter_MGG.docx` → 末尾 Code availability 一句。

## 四、上传前自检

- 仓库内不得出现本机绝对路径、个人邮箱或与代码无关的姓名；发布前用文本检索再查一遍。
- 不要上传 `02_原始数据/` 与 `04_分析结果/` 原始大文件（`.gitignore` 已排除），只上传脚本与说明。
- `README.md` 中的数据下载地址保持 GEO 官方 accession，不再分发原始矩阵。
