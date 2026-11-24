## =====================================================================
## S06_05  外部验证失败的归因诊断（GSE112201）
## ---------------------------------------------------------------------
## 背景：
##   S06_04 显示 9 个候选基因在 GSE112201 中方向一致率仅 5/9（55.6%），
##   4 个统计核心只有 2 个一致，且无一个 padj<0.05。
##   验证集 logFC 全部落在 -0.56 ~ +0.19，而训练集是 -0.76 ~ +0.93，
##   组间差异整体塌缩。
##
## 必须回答的问题：这是"我们的基因不稳"，还是"这个数据集没有信号"？
##   两者的处理方式完全不同：
##     · 若数据集本身两组分不开 → 该数据集不适合作为验证集，方法学中说明，
##       候选基因的可靠性以训练集内部的多方法互证为准
##     · 若数据集两组分得很开、只是我们的 9 个基因不跟随 → 说明候选基因
##       确实不稳健，需要重新审视 S06 的筛选逻辑
##
## 三个诊断指标：
##   ① 全转录组 DEG 数量：组间差异到底有多大
##      参照系：训练集 GSE142538 是 3,964 个（padj<0.05 且 |logFC|>1）
##   ② PCA：4 例对照与 24 例 FECD 在全局转录组上是否分离
##      若对照点散落在 FECD 点中间，说明无全局分离
##   ③ 9 个候选基因的组内分布：箱线图直接看两组分布是否重叠
##
## 输入：GSE112201_GeneCounts.tsv.gz / T16_candidateGene_intersection.csv
## 输出：FigS3A_S3B_GSE112201_diagnostics.png/.pdf
##       控制台打印全部诊断数字
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

options(connectionObserver = NULL)

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
raw_dir <- file.path(root, "02_原始数据/01_GEO原始下载")

## ---------- 1. 重建表达矩阵 ----------
cat("\n========== 1. 重建 GSE112201 全矩阵 ==========\n")
raw <- fread(file.path(raw_dir, "GSE112201_GeneCounts.tsv.gz"))
sample_cols <- names(raw)[-(1:5)]
sample_cols <- gsub("[[:space:]]", "", sample_cols)
setnames(raw, old = names(raw)[-(1:5)], new = sample_cols)

ctrl_ids <- c("RNA115", "RNA136", "RNA84", "RNA89")
grp <- ifelse(sample_cols %in% ctrl_ids, "Control", "FECD")

allcnt <- as.matrix(raw[, ..sample_cols])
rownames(allcnt) <- raw$GeneID
if (any(duplicated(rownames(allcnt)))) {
  allcnt <- rowsum(allcnt, group = rownames(allcnt))
}

dge <- DGEList(counts = allcnt, group = grp)
keep <- filterByExpr(dge, group = grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 2)
cat("过滤后矩阵：", nrow(logcpm), "基因 × ", ncol(logcpm), "样本\n")

## ---------- 2. 诊断①：全转录组 DEG 数量 ----------
cat("\n========== 2. 诊断① 全转录组 DEG ==========\n")
design <- model.matrix(~ 0 + factor(grp, levels = c("Control", "FECD")))
colnames(design) <- c("Control", "FECD")
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
cont <- makeContrasts(FECD - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
tt <- topTable(fit2, number = Inf, sort.by = "P")

cat("padj<0.05           ：", sum(tt$adj.P.Val < 0.05), "\n")
cat("padj<0.05 且 |logFC|>1：",
    sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1), "\n")
cat("padj<0.01           ：", sum(tt$adj.P.Val < 0.01), "\n")
cat("padj<0.05 且 |logFC|>0.5：",
    sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5), "\n")
cat("\nTop 15 DEG：\n")
print(head(tt[, c("logFC", "AveExpr", "P.Value", "adj.P.Val")], 15))

## 对照：训练集 GSE142538 的 DEG 规模（读 S03 的结果，用于并列比较）
t03b <- file.path(out_tab, "T03b_DEG_significant.csv")
if (file.exists(t03b)) {
  n_train <- nrow(fread(t03b))
  cat("\n参照——训练集 GSE142538（9 vs 6）：", n_train,
      "个 DEG（padj<0.05 且 |logFC|>1）\n")
}

## ---------- 3. 诊断②：PCA ----------
cat("\n========== 3. 诊断② PCA ==========\n")
## 取变异最大的 5000 个基因做 PCA
v_by_gene <- apply(logcpm, 1, var)
top5000 <- logcpm[order(-v_by_gene), ][1:min(5000, nrow(logcpm)), ]
pca <- prcomp(t(top5000), scale. = TRUE)

pca_df <- data.frame(
  sample = colnames(logcpm),
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  group  = grp,
  stringsAsFactors = FALSE
)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
cat("PC1 解释方差：", var_exp[1], "%   PC2：", var_exp[2], "%\n")

## 组间在 PC1 上的差异有多大？（若很小，说明无全局分离）
pc1_ctrl <- pca_df$PC1[pca_df$group == "Control"]
pc1_fecd <- pca_df$PC1[pca_df$group == "FECD"]
cat("PC1 均值：Control =", round(mean(pc1_ctrl), 2),
    "  FECD =", round(mean(pc1_fecd), 2),
    "  差值 =", round(mean(pc1_fecd) - mean(pc1_ctrl), 2), "\n")
cat("PC1 标准差：Control =", round(sd(pc1_ctrl), 2),
    "  FECD =", round(sd(pc1_fecd), 2), "\n")
cat("对照样本在 PC1 上的位次（共 ", nrow(pca_df), " 例）：\n", sep = "")
ctrl_rank <- rank(pca_df$PC1)[pca_df$group == "Control"]
for (i in seq_along(ctrl_rank)) {
  cat("  ", pca_df$sample[pca_df$group == "Control"][i],
      " -> 第 ", as.integer(ctrl_rank[i]), " 位\n", sep = "")
}

## ---------- 4. 诊断③：候选基因箱线图 ----------
cat("\n========== 4. 诊断③ 候选基因分布 ==========\n")
t16  <- fread(file.path(out_tab, "T16_candidateGene_intersection.csv"))
cand <- t16[pass_all == TRUE, gene]
t18  <- fread(file.path(out_tab, "T18_stability_selection.csv"))
core <- t18[freq_enet >= 0.5, gene]

box_df <- data.frame(
  gene  = rep(rownames(logcpm)[rownames(logcpm) %in% cand], times = ncol(logcpm)),
  value = as.numeric(logcpm[rownames(logcpm) %in% cand, ]),
  group = rep(grp, each = sum(rownames(logcpm) %in% cand)),
  stringsAsFactors = FALSE
)
## 基因按训练集 logFC 排序，保证图面有序
ord <- t16[gene %in% cand][order(-logFC), gene]
box_df$gene <- factor(box_df$gene, levels = ord)
box_df$group <- factor(box_df$group, levels = c("Control", "FECD"))

## 图上标出该基因的训练集 logFC，便于一眼对照
lab_df <- t16[gene %in% cand][, .(gene, logFC)]
lab_df$gene <- factor(lab_df$gene, levels = ord)
lab_df$tag <- ifelse(lab_df$gene %in% core, "★", "")

p_box <- ggplot(box_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35,
               colour = "grey30") +
  geom_jitter(width = 0.14, size = 0.7, colour = "grey25", alpha = 0.6) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Control" = "#2C6FAC", "FECD" = "#C0392B"), name = NULL) +
  labs(title = "GSE112201: candidate gene distribution across groups",
       subtitle = "4 control vs 24 FECD; nine candidate genes ordered by discovery log2FC",
       x = NULL, y = "log2CPM") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text     = element_text(size = 7.5, face = "bold"),
        axis.text.x    = element_text(size = 7, angle = 0),
        axis.text.y    = element_text(size = 6.5),
        axis.title.y   = element_text(size = 7.5),
        legend.position = "bottom",
        legend.text    = element_text(size = 7),
        plot.title     = element_text(size = 9.5, face = "bold"),
        plot.subtitle  = element_text(size = 7.2, colour = "grey35"),
        plot.margin    = margin(8, 8, 6, 6))

## ---------- 5. 拼图：PCA + 箱线 ----
## 用 patchwork 若不可用则分开保存
cat("\n========== 5. 出图 ==========\n")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_text(aes(label = ifelse(group == "Control", sample, "")),
            size = 2.1, vjust = -0.9, colour = "#2C6FAC") +
  scale_colour_manual(values = c("Control" = "#2C6FAC", "FECD" = "#C0392B"), name = NULL) +
  labs(title = "PCA of GSE112201 (top 5000 variable genes)",
       subtitle = paste0("PC1 ", var_exp[1], "% | PC2 ", var_exp[2], "%"),
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        axis.title    = element_text(size = 8),
        axis.text     = element_text(size = 7.3, colour = "black"),
        axis.line     = element_line(linewidth = 0.35),
        legend.position = "bottom",
        legend.text   = element_text(size = 7),
        plot.margin   = margin(8, 10, 6, 6))

ggsave(file.path(out_fig, "FigS3A_GSE112201_PCA.png"), p_pca,
       width = 4.8, height = 4.2, dpi = 300)
ggsave(file.path(out_fig, "FigS3A_GSE112201_PCA.pdf"), p_pca,
       width = 4.8, height = 4.2)
cat("FigS3A (PCA) 已输出\n")

ggsave(file.path(out_fig, "FigS3B_GSE112201_candidateBox.png"), p_box,
       width = 6.6, height = 5.6, dpi = 300)
ggsave(file.path(out_fig, "FigS3B_GSE112201_candidateBox.pdf"), p_box,
       width = 6.6, height = 5.6)
cat("FigS3B (箱线图) 已输出\n")

cat("\n========== 完成 ==========\n")
