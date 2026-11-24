## =====================================================================
## S04_02 全部 19 个双硫死亡基因投影到 WGCNA 六大模块
## 背景：S04_01 中 MAD top 25% 只保留了 5 个双硫死亡基因进入网络构建，
##       其余 14 个被筛掉。本脚本用已构建的模块 eigengene（T05c），
##       把全部 19 个基因投影回去，计算每个基因与 6 个模块的 kME。
## 输出：T08_disulfidptosisGene_allModule_kME.csv
##       Fig2D_disulfidptosisGene_kME_heatmap.png / .pdf
## =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

root <- Sys.getenv("FECD_PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]))), "..", ".."), mustWork = FALSE)
  } else {
    root <- normalizePath(getwd(), mustWork = FALSE)
  }
}
out_tab <- file.path(root, "04_分析结果/01_结果表格")
out_fig <- file.path(root, "04_分析结果/02_图片/主图")

set.seed(20261121)

## ---------- 1. 读模块 eigengene ----------
me <- as.matrix(read.csv(file.path(out_tab, "T05c_moduleEigengenes.csv"), row.names = 1, check.names = FALSE))
## 去掉 MEgrey（若存在）
me <- me[, !colnames(me) %in% "MEgrey", drop = FALSE]
cat("模块 eigengene：", nrow(me), "样本 x", ncol(me), "模块\n")
print(colnames(me))

## ---------- 2. 读表达矩阵并计算 logCPM ----------
counts <- as.matrix(read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_counts_gene.csv"),
                             row.names = 1, check.names = FALSE))
info <- read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_sampleInfo.csv"), stringsAsFactors = FALSE)

samples <- rownames(me)
stopifnot(all(samples %in% colnames(counts)))

cnt <- counts[, samples, drop = FALSE]
lib <- colSums(cnt); lib[lib <= 0] <- 1
logcpm <- log2(sweep(cnt, 2, lib, "/") * 1e6 + 1)

## ---------- 3. 双硫死亡基因 ----------
dsg <- read.csv(file.path(root, "02_原始数据/02_基因集/双硫死亡基因集/DisulfidptosisGenes_Liu2023_v01.csv"),
                stringsAsFactors = FALSE)
dsg_all <- unique(dsg$gene)
hit <- intersect(dsg_all, rownames(logcpm))
miss <- setdiff(dsg_all, rownames(logcpm))
cat("\n命中：", length(hit), "/", length(dsg_all), "\n")
if (length(miss)) cat("缺失：", paste(miss, collapse = ", "), "\n")

## ---------- 4. 计算 kME（基因与模块 eigengene 的 Pearson 相关） ----------
expr <- t(logcpm[hit, , drop = FALSE])       # 样本 x 基因
kme <- cor(expr, me, use = "p")              # 基因 x 模块
kme_df <- as.data.frame(kme)
kme_df$gene <- rownames(kme_df)

## 最佳模块
best <- data.frame(
  gene = rownames(kme),
  best_module = colnames(kme)[apply(abs(kme), 1, which.max)],
  best_kME = apply(kme, 1, function(x) x[which.max(abs(x))]),
  stringsAsFactors = FALSE
)
best$category <- dsg$category[match(best$gene, dsg$gene)]
best <- best[order(-abs(best$best_kME)), ]

cat("\n=== 19 个双硫死亡基因的最佳模块 ===\n")
print(best, row.names = FALSE)

cat("\n各模块归属统计：\n")
print(table(best$best_module))

## 方向判断：turquoise 与 FECD 正相关，blue 负相关
cat("\n方向解读（turquoise+/blue-）：\n")
dir_map <- c(MEturquoise = "up in FECD", MEblue = "down in FECD",
             MEmagenta = "up (weak)", MEgreenyellow = "down",
             MElightcyan = "up (weak)", MEblack = "neutral")
best$direction <- ifelse(best$best_kME > 0,
                         paste0(best$best_module, " (same direction)"),
                         paste0(best$best_module, " (opposite direction)"))
print(best[, c("gene", "best_module", "best_kME", "direction", "category")], row.names = FALSE)

## ---------- 5. 保存 ----------
write.csv(kme_df, file.path(out_tab, "T08_disulfidptosisGene_allModule_kME.csv"), row.names = FALSE)
write.csv(best, file.path(out_tab, "T08b_disulfidptosisGene_bestModule.csv"), row.names = FALSE)

## ---------- 6. 热图 ----------
long <- melt(as.matrix(kme_df[, setdiff(colnames(kme_df), "gene")]))
names(long) <- c("gene", "module", "kME")
long$gene <- factor(long$gene, levels = rev(best$gene))
long$module <- factor(long$module, levels = colnames(kme))

p <- ggplot(long, aes(x = module, y = gene, fill = kME)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(abs(kME) >= 0.5, sprintf("%.2f", kME), "")),
            size = 2.2, color = ifelse(abs(long$kME) >= 0.7, "white", "grey20")) +
  ## PDLIM1 的 |kME| 都低于 0.5，色块接近白色；单独标注其两个投影值，避免看起来是空白行
  geom_text(data = long[long$gene == "PDLIM1" & abs(long$kME) >= 0.15, ],
            aes(label = sprintf("%.2f", kME)), size = 1.9, color = "grey20") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "kME") +
  labs(x = NULL, y = NULL,
       title = "Disulfidptosis genes across WGCNA modules",
       subtitle = "GSE142538 | kME = correlation with module eigengene") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8))

ggsave(file.path(out_fig, "Fig2D_disulfidptosisGene_kME_heatmap.png"), p,
       width = 4.6, height = 6.2, dpi = 300)
ggsave(file.path(out_fig, "Fig2D_disulfidptosisGene_kME_heatmap.pdf"), p,
       width = 4.6, height = 6.2)

cat("\n=== 完成 ===\n")
cat("T08 / T08b / Fig2D\n")
