## =====================================================================
## S06_04  外部验证（GSE112201：24 FECD vs 4 Control）
## ---------------------------------------------------------------------
## 目的：检验 S06_01 三层证据筛出的 9 个候选基因，在一个完全独立的
##       队列中表达方向是否与训练集（GSE142538）一致。
##       这是整个 S06 里最能"顶住审稿"的一步——因为它是真正的外部数据。
##
## 样本分组依据（已逐条核对 GEO，2026-11-21）：
##   GSE112201 共 28 例。
##   GSM4932380  2011-376 RNA115 [control]   Characteristics: group: Control
##   GSM4932381  2011-409 RNA136 [control]   Characteristics: group: Control
##   GSM4932382  2011-325 RNA84  [control]   Characteristics: group: Control
##   GSM4932383  2011-305 RNA89  [control]   Characteristics: group: Control
##   其余 GSM3060909–GSM3060932 共 24 例，标题中无 [control] 标注，为 FECD。
##   注意：这 4 例对照是 2020-11-24 补交的，早前 24 例 FECD 于 2018 年提交。
##   原始矩阵列名与 GEO 标题一一对应（GSM3060920 的 "RNA 120" 在文件中
##   写作 "RNA120"，已去空格）。
##
## 处理口径：
##   · 原始 counts（HTSeq-count 产出）→ TMM 归一化 → log2CPM（prior.count=2）
##   · limma 线性模型，对比 FECD - Control
##   · 与训练集 log2FC 比较**符号**而非大小——两个队列平台、深度、样本量
##     都不同，效应量不可直接比较，方向一致性才是可验证的部分
##
## 预期：若 9 个基因中方向一致数 ≥7，说明候选基因在独立队列中可复现；
##       若仅 4–5 个一致，则需回头检查训练集信号是否稳定。
##
## 输入：GSE112201_GeneCounts.tsv.gz
##       T16_candidateGene_intersection.csv
## 输出：T19_external_validation.csv
##       FigS3C_external_validation.png/.pdf
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

## ---------- 1. 读入 GSE112201 ----------
cat("\n========== 1. 读入 GSE112201 ==========\n")
raw <- fread(file.path(raw_dir, "GSE112201_GeneCounts.tsv.gz"))
cat("原始矩阵：", nrow(raw), "行 × ", ncol(raw), "列\n")
cat("前 5 列：", paste(names(raw)[1:5], collapse = " | "), "\n")

## 第 6 列起是样本；列名去掉可能存在的空格（GEO 标题里 "RNA 120" → 文件内 "RNA120"）
sample_cols <- names(raw)[-(1:5)]
sample_cols <- gsub("[[:space:]]", "", sample_cols)
setnames(raw, old = names(raw)[-(1:5)], new = sample_cols)
cat("样本列数：", length(sample_cols), "\n")

## ---------- 2. 分组（依据 GEO 逐条核对结果） ----------
ctrl_ids <- c("RNA115", "RNA136", "RNA84", "RNA89")
grp <- ifelse(sample_cols %in% ctrl_ids, "Control", "FECD")
cat("\n分组：\n")
print(table(grp))
cat("对照组样本：", paste(sample_cols[sample_cols %in% ctrl_ids], collapse = ", "), "\n")

## 反向检查：确认 4 个对照全部命中，防止拼写不一致导致漏配
stopifnot(sum(grp == "Control") == 4)

## ---------- 3. 取候选基因 ----------
cat("\n========== 3. 提取候选基因 ==========\n")
t16  <- fread(file.path(out_tab, "T16_candidateGene_intersection.csv"))
cand <- t16[pass_all == TRUE, gene]
cat("候选基因：", length(cand), "个\n")

sub <- raw[GeneID %in% cand]
cat("在 GSE112201 中命中：", nrow(sub), " / ", length(cand), "\n")
miss <- setdiff(cand, sub$GeneID)
if (length(miss) > 0) cat("缺失基因：", paste(miss, collapse = ", "), "\n")

## 顺带看一眼 19 个可检测基因里还有多少能在本数据集找到
all19 <- t16$gene
cat("19 个可检测基因中，本数据集命中：", sum(all19 %in% raw$GeneID), "\n")

## ---------- 4. 归一化 + 差异表达 ----------
## 关键：TMM 的文库大小必须用**全部基因**的列和来算，不能用这 9 个候选基因
## 的列和。否则归一化基准被 9 个基因的绝对表达量绑架，log2CPM 会虚高
## （第一版就踩了这个坑：文库只有 0.05-0.27 M reads、log2CPM 冲到 19.9）
cat("\n========== 4. TMM 归一化与 limma ==========\n")

allcnt <- as.matrix(raw[, ..sample_cols])
rownames(allcnt) <- raw$GeneID
## 同一 GeneID 若出现多行（不同转录本合并残留），保留 counts 之和
if (any(duplicated(rownames(allcnt)))) {
  cat("GeneID 重复：", sum(duplicated(rownames(allcnt))), "行，合并求和\n")
  allcnt <- rowsum(allcnt, group = rownames(allcnt))
}
cat("全基因矩阵：", nrow(allcnt), "基因 × ", ncol(allcnt), "样本\n")

dge_all <- DGEList(counts = allcnt, group = grp)
## 低表达基因过滤（与 S03 的处理逻辑一致，避免噪声基因进入方差收缩）
keep <- filterByExpr(dge_all, group = grp)
cat("filterByExpr 保留：", sum(keep), " / ", nrow(allcnt),
    " (", round(100 * sum(keep) / nrow(allcnt), 1), "%)\n", sep = "")
dge_all <- dge_all[keep, , keep.lib.sizes = FALSE]

dge_all <- calcNormFactors(dge_all, method = "TMM")
cat("文库大小（百万 reads）：\n")
print(round(dge_all$samples$lib.size / 1e6, 2))
cat("TMM 归一化因子范围：", round(min(dge_all$samples$norm.factors), 3),
    "~", round(max(dge_all$samples$norm.factors), 3), "\n")

## 全矩阵 log2CPM（用于看候选基因的表达量级）
logcpm_all <- cpm(dge_all, log = TRUE, prior.count = 2)
logcpm <- logcpm_all[rownames(logcpm_all) %in% cand, , drop = FALSE]
cat("\n候选基因 log2CPM 范围：", round(min(logcpm), 2), "~", round(max(logcpm), 2), "\n")
cat("各基因均值（按组）：\n")
for (g in rownames(logcpm)) {
  cat(sprintf("  %-8s Control=%6.2f  FECD=%6.2f  diff=%+.2f\n",
              g,
              mean(logcpm[g, grp == "Control"]),
              mean(logcpm[g, grp == "FECD"]),
              mean(logcpm[g, grp == "FECD"]) - mean(logcpm[g, grp == "Control"])))
}

## limma：4 例对照 + 24 例 FECD。
## 注意 eBayes 的方差收缩必须基于全转录组的方差分布，不能只用 9 个基因，
## 否则 moderation 失去参照。故在全矩阵上拟合，最后取候选基因的结果。
design <- model.matrix(~ 0 + factor(grp, levels = c("Control", "FECD")))
colnames(design) <- c("Control", "FECD")
v <- voom(dge_all, design, plot = FALSE)
fit <- lmFit(v, design)
cont <- makeContrasts(FECD - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
tt_all <- topTable(fit2, number = Inf, sort.by = "none")
tt_all$gene <- rownames(tt_all)
tt <- tt_all[tt_all$gene %in% cand, ]
cat("\nlimma 结果（候选基因）：\n")
print(tt[, c("gene", "logFC", "P.Value", "adj.P.Val")], row.names = FALSE)

## ---------- 5. 与训练集比对方向 ----------
cat("\n========== 5. 方向一致性 ==========\n")
val <- data.table(gene      = tt$gene,
                  logFC_val = tt$logFC,
                  p_val     = tt$P.Value,
                  padj_val  = tt$adj.P.Val)

m <- merge(t16[, .(gene, logFC_train = logFC, padj_train = adj.P.Val,
                   best_module, best_kME)],
           val, by = "gene", all.x = TRUE)
m <- m[order(-abs(logFC_train))]
m[, consistent := sign(logFC_train) == sign(logFC_val)]

print(m[, .(gene, logFC_train, logFC_val, padj_train, padj_val,
            best_module, consistent)], row.names = FALSE)

n_total <- sum(!is.na(m$consistent))
n_cons  <- sum(m$consistent, na.rm = TRUE)
cat("\n方向一致：", n_cons, " / ", n_total, "  (",
    round(100 * n_cons / n_total, 1), "%)\n", sep = "")

## 4 个统计核心（稳定性选择 ≥0.5）的方向一致性单独统计
t18 <- fread(file.path(out_tab, "T18_stability_selection.csv"))
core4 <- t18[freq_enet >= 0.5, gene]
cat("统计核心（稳定性频率>=0.5）：", paste(core4, collapse = ", "), "\n")
core_m <- m[gene %in% core4]
cat("核心基因方向一致：", sum(core_m$consistent, na.rm = TRUE), " / ", nrow(core_m), "\n")

## 显著性：外部队列 padj<0.05 的基因
cat("\n外部队列 padj<0.05 的候选基因：\n")
sig_val <- m[!is.na(padj_val) & padj_val < 0.05]
if (nrow(sig_val) == 0) {
  cat("  （无）\n")
} else {
  print(sig_val[, .(gene, logFC_val, padj_val, consistent)], row.names = FALSE)
}

## ---------- 6. 保存 ----------
write.csv(m, file.path(out_tab, "T19_external_validation.csv"), row.names = FALSE)
cat("\nT19 已保存\n")

## ---------- 7. FigS3C：两队列 logFC 对照散点 ----------
cat("\n========== 7. 绘制 FigS3C ==========\n")

plot_df <- m[!is.na(logFC_val)]
plot_df$label_dir <- ifelse(plot_df$consistent, "same direction", "opposite")

p_sc <- ggplot(plot_df, aes(x = logFC_train, y = logFC_val)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  ## 一二三象限的浅色提示区（同向区）
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#C0392B", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#C0392B", alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = label_dir), size = 2.0) +
  geom_text(aes(label = gene), size = 2.2, vjust = -0.85, colour = "grey20") +
  scale_colour_manual(values = c("same direction" = "#C0392B",
                                 "opposite" = "#2C6FAC"), name = NULL) +
  labs(title = "Exploratory comparison of candidate genes (GSE112201)",
       subtitle = paste0("x: log2FC in GSE142538 (9 control vs 6 FECD) | ",
                         "y: log2FC in GSE112201 (4 control vs 24 FECD)"),
       x = "log2FC — discovery set (GSE142538)", y = "log2FC — exploratory set (GSE112201)") +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.0, colour = "grey35"),
        axis.title    = element_text(size = 8),
        axis.text     = element_text(size = 7.3, colour = "black"),
        axis.line     = element_line(linewidth = 0.35),
        legend.position = "bottom",
        legend.text   = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"),
        plot.margin   = margin(8, 10, 6, 6))

ggsave(file.path(out_fig, "FigS3C_external_validation.png"), p_sc,
       width = 5.4, height = 5.0, dpi = 600)
ggsave(file.path(out_fig, "FigS3C_external_validation.pdf"), p_sc,
       width = 5.4, height = 5.0)
cat("FigS3C 已输出\n")

cat("\n========== 完成 ==========\n")
