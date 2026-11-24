## =====================================================================
## S00 创建项目输出目录
## 目的：让后续脚本可以在一个只含 02_原始数据 的新目录中直接运行
## 用法：Rscript scripts/S00_setup_directories.R
## =====================================================================

root <- Sys.getenv("FECD_PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]))), ".."), mustWork = FALSE)
  } else {
    root <- normalizePath(getwd(), mustWork = FALSE)
  }
}

dirs <- c(
  "02_原始数据/03_处理后数据",
  "04_分析结果/01_结果表格",
  "04_分析结果/02_图片/主图",
  "04_分析结果/02_图片/附图",
  "04_分析结果/03_中间文件",
  "04_分析结果/03_独立验证_GSE74123",
  "04_分析结果/04_分析报告",
  "05_文稿"
)
for (d in dirs) {
  dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
}

cat("Project root:", root, "\n")
cat("Output directories are ready.\n")
