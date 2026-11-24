## =====================================================================
## S02_01 双硫死亡基因在 GSE142538 的检出率与表达量级核查
## 目的：接受 KC 项目的教训，正式分析前先判断"基因集在该组织里有没有料"
##       判定标准：24 个基因中，log2CPM > 4 的至少 8 个，否则该基因集不适用
## 输入：02_原始数据/03_处理后数据/GSE142538_counts_gene.csv + sampleInfo
##       02_原始数据/02_基因集/双硫死亡基因集/DisulfidptosisGenes_Liu2023_v01.csv
## 输出：T01_disulfidptosisGene_detection.csv（24 基因检出与表达）
##       T02_sample_group_summary.csv
##       FigS1A_detectionRate.png / FigS1B_expressionLevel.png
## 参数：log2(CPM+1)；TMM 不用于本核查（仅做量级判断）
## =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
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
out_fig <- file.path(root, "04_分析结果/02_图片/附图")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

set.seed(20261121)

## ---------- 1. 读数据 ----------
cat("=== 读入计数矩阵 ===\n")
counts <- as.matrix(read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_counts_gene.csv"),
                             row.names = 1, check.names = FALSE))
info <- read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_sampleInfo.csv"),
                 stringsAsFactors = FALSE)
cat("矩阵：", nrow(counts), "genes x", ncol(counts), "samples\n")
print(table(info$group))

## 样本列与 info 对齐
stopifnot(identical(colnames(counts), info$sample))

## ---------- 2. 双硫死亡基因 ----------
lact <- read.csv(file.path(root, "02_原始数据/02_基因集/双硫死亡基因集/DisulfidptosisGenes_Liu2023_v01.csv"),
                 stringsAsFactors = FALSE)
dsg <- unique(lact$gene)
cat("\n双硫死亡基因集：", length(dsg), "个\n")

hit <- dsg[dsg %in% rownames(counts)]
miss <- setdiff(dsg, rownames(counts))
cat("命中：", length(hit), "/", length(dsg), "\n")
if (length(miss) > 0) cat("未命中：", paste(miss, collapse = ", "), "\n")

## ---------- 3. CPM 归一化 ----------
lib <- colSums(counts)
cpm <- sweep(counts, 2, lib, "/") * 1e6
logcpm <- log2(cpm + 1)

sub <- logcpm[hit, , drop = FALSE]

## ---------- 4. 检出率与表达量级 ----------
det <- data.frame(
  gene = hit,
  category = lact$category[match(hit, lact$gene)],
  pct_detected = rowMeans(sub > 0) * 100,
  mean_logcpm = rowMeans(sub),
  median_logcpm = apply(sub, 1, median),
  max_logcpm = apply(sub, 1, max),
  stringsAsFactors = FALSE
)
det <- det[order(-det$mean_logcpm), ]

cat("\n=== 24 基因表达量级排序 ===\n")
print(det[, c("gene", "category", "pct_detected", "mean_logcpm", "max_logcpm")], row.names = FALSE)

## 分层
det$level <- ifelse(det$mean_logcpm >= 8, "High (>=8)",
             ifelse(det$mean_logcpm >= 6, "Moderate (6-8)",
             ifelse(det$mean_logcpm >= 4, "Low-mid (4-6)",
             ifelse(det$mean_logcpm >= 2, "Low (2-4)", "Very low (<2)"))))
cat("\n表达量级分布：\n"); print(table(det$level))

n_valid <- sum(det$mean_logcpm >= 4)
cat("\n>>> log2CPM >= 4 的基因数：", n_valid, "/", nrow(det), "\n")
cat(">>> 判定标准（>=8 个）：", ifelse(n_valid >= 8, "PASS — 可继续", "FAIL — 信号不足，需换基因集或换方向"), "\n")

## ---------- 5. 组间比较（FECD_REP vs Control） ----------
cat("\n=== FECD_REP vs Control 的组间差异（探索性） ===\n")
grp_lv <- c("Control", "FECD_REP")
sel <- info$group %in% grp_lv
g2 <- factor(info$group[sel], levels = grp_lv)

cmp <- data.frame()
for (g in hit) {
  v <- sub[g, sel]
  wt <- wilcox.test(v[g2 == "FECD_REP"], v[g2 == "Control"])
  cmp <- rbind(cmp, data.frame(
    gene = g, Control_mean = mean(v[g2 == "Control"]), FECD_mean = mean(v[g2 == "FECD_REP"]),
    diff = mean(v[g2 == "FECD_REP"]) - mean(v[g2 == "Control"]),
    p_value = wt$p.value, stringsAsFactors = FALSE))
}
cmp <- cmp[order(cmp$p_value), ]
print(cmp, row.names = FALSE)
cat("\np<0.05 的基因数：", sum(cmp$p_value < 0.05), "/", nrow(cmp), "\n")

## ---------- 6. Fig S1A 检出率 ----------
det$gene <- factor(det$gene, levels = rev(det$gene))
det$category <- factor(det$category)

p1 <- ggplot(det, aes(x = gene, y = pct_detected, fill = category)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Cystine transport"   = "#B2182B",
    "ER glycosylation"    = "#D6604D",
    "Cytoskeleton"        = "#4393C3",
    "Metabolism"          = "#2166AC",
    "Mitochondria"        = "#5AAE61")) +
  coord_flip() +
  labs(x = NULL, y = "Detection rate (% samples)",
       title = "Disulfidptosis gene detection in FECD endothelium",
       subtitle = "GSE142538 | 25 samples (9 Control / 6 Pre_S / 6 FECD_REP / 4 FECD_NR)") +
  theme_classic(base_size = 9) +
  theme(legend.title = element_blank(), legend.position = "right",
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8))

ggsave(file.path(out_fig, "FigS1A_detectionRate.png"), p1, width = 6.2, height = 5.0, dpi = 300)
ggsave(file.path(out_fig, "FigS1A_detectionRate.pdf"), p1, width = 6.2, height = 5.0)

## ---------- 7. Fig S1B 表达量级 ----------
p2 <- ggplot(det, aes(x = gene, y = mean_logcpm, fill = level)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "#B2182B", linewidth = 0.5) +
  annotate("text", x = 0.6, y = 4.3, label = "log2CPM = 4 (signal floor)",
           hjust = 0, size = 2.6, color = "#B2182B") +
  scale_fill_manual(values = c("High (>=8)" = "#08519C", "Moderate (6-8)" = "#4292C6",
                               "Low-mid (4-6)" = "#9ECAE1", "Low (2-4)" = "#FEE0D2",
                               "Very low (<2)" = "#FC9272")) +
  coord_flip() +
  labs(x = NULL, y = expression(log[2]*"(CPM+1) mean across 25 samples"),
       title = "Expression magnitude of disulfidptosis genes",
       subtitle = sprintf("GSE142538 | genes with log2CPM >= 4: %d/%d", n_valid, nrow(det))) +
  theme_classic(base_size = 9) +
  theme(legend.title = element_blank(), legend.position = "right",
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8))

ggsave(file.path(out_fig, "FigS1B_expressionLevel.png"), p2, width = 6.2, height = 5.0, dpi = 300)
ggsave(file.path(out_fig, "FigS1B_expressionLevel.pdf"), p2, width = 6.2, height = 5.0)

## ---------- 8. 保存 ----------
write.csv(det, file.path(out_tab, "T01_disulfidptosisGene_detection.csv"), row.names = FALSE)
write.csv(cmp, file.path(out_tab, "T02_geneFECDvsControl.csv"), row.names = FALSE)
write.csv(data.frame(sample = info$sample, group = info$group, lib_size = lib,
                     stringsAsFactors = FALSE),
          file.path(out_tab, "T02b_sample_libsize.csv"), row.names = FALSE)

cat("\n=== 完成 ===\n")
