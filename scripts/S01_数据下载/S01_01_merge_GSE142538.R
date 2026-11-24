## =====================================================================
## S01_01 合并 GSE142538 的 25 个 Cufflinks 文件为基因级计数矩阵
## 输入：02_原始数据/01_GEO原始下载/GSE142538/*.txt.gz
## 每个文件 4 列：tracking_id(XLOC) | gene_name | raw_frags | FPKM
## 处理：gene_name 若含逗号（多基因共用转录本）取第一个；同一基因的多转录本 raw_frags 求和
## 输出：02_原始数据/03_处理后数据/GSE142538_counts_gene.csv
##       02_原始数据/03_处理后数据/GSE142538_sampleInfo.csv
## =====================================================================

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("FECD_PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]))), "..", ".."), mustWork = FALSE)
  } else {
    root <- normalizePath(getwd(), mustWork = FALSE)
  }
}
in_dir  <- file.path(root, "02_原始数据/01_GEO原始下载/GSE142538")
out_dir <- file.path(root, "02_原始数据/03_处理后数据")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
cat("文件数：", length(files), "\n")

## 从文件名解析样本名与分组
parse_name <- function(f) {
  b <- basename(f)
  b <- sub("\\.txt\\.gz$", "", b)
  b <- sub("^GSM[0-9]+_", "", b)
  b
}
samples <- sapply(files, parse_name)
group <- ifelse(grepl("^Control", samples), "Control",
         ifelse(grepl("^Pre_S", samples),   "Pre_S",
         ifelse(grepl("^FECD_REP", samples), "FECD_REP",
         ifelse(grepl("^FECD_NR", samples),  "FECD_NR", "Other"))))
cat("\n分组分布：\n"); print(table(group))

## 逐文件聚合到基因级
gene_counts_list <- list()
for (i in seq_along(files)) {
  d <- fread(files[i], header = TRUE, data.table = FALSE, showProgress = FALSE)
  ## 取 gene_name 第一个（逗号前）
  g <- sub(",.*$", "", d$gene_name)
  ## 丢弃空基因名
  keep <- !is.na(g) & g != "" & g != "-"
  d <- d[keep, ]; g <- g[keep]
  agg <- tapply(d$raw_frags, g, sum, na.rm = TRUE)
  gene_counts_list[[samples[i]]] <- agg
  cat(sprintf("  [%2d/25] %-14s 基因数=%5d\n", i, samples[i], length(agg)))
}

## 合并为矩阵（取并集基因，缺失填 0）
all_genes <- unique(unlist(lapply(gene_counts_list, names)))
cat("\n并集基因数：", length(all_genes), "\n")

mat <- matrix(0L, nrow = length(all_genes), ncol = length(samples),
              dimnames = list(all_genes, samples))
for (s in samples) {
  v <- gene_counts_list[[s]]
  mat[names(v), s] <- as.integer(v)
}

cat("矩阵维度：", nrow(mat), "genes x", ncol(mat), "samples\n")
cat("总 counts：", format(sum(mat), big.mark = ","), "\n")

## 保存
write.csv(mat, file.path(out_dir, "GSE142538_counts_gene.csv"), quote = FALSE)
info <- data.frame(sample = samples, group = group,
                   gsm = sub("_.*$", "", basename(files)),
                   stringsAsFactors = FALSE)
write.csv(info, file.path(out_dir, "GSE142538_sampleInfo.csv"), row.names = FALSE)

cat("\n=== 完成 ===\n")
cat("GSE142538_counts_gene.csv\nGSE142538_sampleInfo.csv\n")
print(info, row.names = FALSE)
