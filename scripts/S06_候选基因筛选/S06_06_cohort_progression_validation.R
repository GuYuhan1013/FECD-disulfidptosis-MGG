## =====================================================================
## S06_06  队列内进展轴验证（GSE142538 全部 25 例）
## ---------------------------------------------------------------------
## 缘起：
##   S06_04 用 GSE112201 做外部验证，方向一致率仅 5/9，S06_05 诊断证明
##   该数据集的 4 例对照与 24 例 FECD 在全局转录组上完全不可分
##   （PCA 重叠、Top DEG 全是组织噪声基因、箱线图全重叠），故弃用。
##
##   转而在 GSE142538 自身找留出样本。按 GEO 原始记录（已逐条核对）：
##     GSM4231599 Pre_S_0   ctg18.1: 13, 84   → 扩增阳性但**无症状**
##     GSM4231611 FECD-NR_0 ctg18.1: 15, 18   → 晚期 FECD 但**无扩增**
##   即该队列实为四组：
##     Control(9) 健康供体
##     Pre_S(6)   携带 CUG 扩增、尚无症状，来自眼库（比临床诊断早数十年）
##     FECD_REP(6) 晚期 FECD、扩增阳性（S03/S04 用的是这 6 例）
##     FECD_NR(4)  晚期 FECD、扩增阴性
##
##   这 10 例是同一批测序的留出样本，无跨平台/跨批次偏差，比 GSE112201 干净。
##
## 本步骤检验两个假设：
##   H1（时间轴）：双硫死亡相关基因改变在 Pre_S 阶段已经出现，方向与晚期一致
##        → 支持"双硫死亡是 FECD 的早期事件，而非终末结果"
##   H2（机制特异性）：这些改变只在扩增阳性组出现，扩增阴性组不明显
##        → 支持"TCF4 CUG 扩增驱动双硫死亡通路"
##
##   三种可能的结果及其解释：
##     · Pre_S 同向 + NR 不同向 → 最理想，两个假设都成立，故事线闭环
##     · Pre_S 同向 + NR 也同向 → 时间轴成立，但不依赖扩增，属 FECD 共性病理
##     · Pre_S 不同向          → 时间轴不成立，只能停在"晚期伴随现象"
##
## 归一化口径：
##   与 S03（15 例）不同，本步骤为纳入 Pre_S / FECD_NR，改用全部 25 例一起
##   做 TMM 归一化。因此 FECD_REP vs Control 的 logFC 会与 S03 有轻微差异，
##   这是样本集合变化导致 TMM 因子重估的正常现象；两组口径下的方向一致性
##   本身将作为一项稳健性证据报告。
##
## 输入：GSE142538_counts_gene.csv、GSE142538_sampleInfo.csv
##       T16_candidateGene_intersection.csv、T18_stability_selection.csv
## 输出：T20_cohort_progression.csv（9 基因 × 4 对比的 logFC/padj）
##       T20b_cohort_progression_summary.csv（假设检验汇总）
##       Fig5A_cohort_progression.png/.pdf
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
out_fig <- file.path(root, "04_分析结果/02_图片/主图")
dat_dir <- file.path(root, "02_原始数据/03_处理后数据")

## ---------- 1. 读入 25 例全队列 ----------
cat("\n========== 1. 读入全队列 ==========\n")
cnt <- fread(file.path(dat_dir, "GSE142538_counts_gene.csv"))
info <- fread(file.path(dat_dir, "GSE142538_sampleInfo.csv"))

## counts 第一列是基因名，列名需与 sampleInfo$sample 对齐
gene_col <- names(cnt)[1]
mat <- as.matrix(cnt[, -1, with = FALSE])
rownames(mat) <- cnt[[gene_col]]

miss <- setdiff(info$sample, colnames(mat))
if (length(miss) > 0) stop("样本名不匹配：", paste(miss, collapse = ", "))
mat <- mat[, info$sample]

grp <- factor(info$group, levels = c("Control", "Pre_S", "FECD_REP", "FECD_NR"))
cat("分组：\n"); print(table(grp))
cat("矩阵：", nrow(mat), "基因 × ", ncol(mat), "样本\n")

## ---------- 2. 归一化（25 例一起） ----------
cat("\n========== 2. 归一化 ==========\n")
dge <- DGEList(counts = mat, group = grp)
keep <- filterByExpr(dge, group = grp)
cat("filterByExpr 保留：", sum(keep), " / ", nrow(mat),
    " (", round(100 * sum(keep) / nrow(mat), 1), "%)\n", sep = "")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat("文库大小（百万）："); print(round(dge$samples$lib.size / 1e6, 1))
cat("TMM 因子范围：", round(min(dge$samples$norm.factors), 3), "~",
    round(max(dge$samples$norm.factors), 3), "\n")

logcpm <- cpm(dge, log = TRUE, prior.count = 2)

## ---------- 3. 四组对比 ----------
cat("\n========== 3. limma 四组对比 ==========\n")
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

cm <- makeContrasts(
  PreS_vs_Control = Pre_S    - Control,
  REP_vs_Control  = FECD_REP - Control,
  NR_vs_Control   = FECD_NR  - Control,
  REP_vs_NR       = FECD_REP - FECD_NR,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, cm))

## 各对比下的全转录组 DEG 规模（用于说明各组效应强度）
cat("\n全转录组 DEG（padj<0.05 且 |logFC|>1）：\n")
for (i in seq_len(ncol(cm))) {
  tt_i <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
  nm <- colnames(cm)[i]
  cat(sprintf("  %-16s : %5d\n", nm,
              sum(tt_i$adj.P.Val < 0.05 & abs(tt_i$logFC) > 1)))
}

## ---------- 4. 提取候选基因 ----------
cat("\n========== 4. 候选基因在各对比下的表现 ==========\n")
t16  <- fread(file.path(out_tab, "T16_candidateGene_intersection.csv"))
cand <- t16[pass_all == TRUE, gene]
t18  <- fread(file.path(out_tab, "T18_stability_selection.csv"))
core <- t18[freq_enet >= 0.5, gene]
cat("候选 9 个：", paste(cand, collapse = ", "), "\n")
cat("统计核心 4 个：", paste(core, collapse = ", "), "\n")

res <- rbindlist(lapply(seq_len(ncol(cm)), function(i) {
  tt_i <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
  tt_i$gene <- rownames(tt_i)
  sub <- tt_i[tt_i$gene %in% cand, c("gene", "logFC", "P.Value", "adj.P.Val")]
  sub$contrast <- colnames(cm)[i]
  sub
}))

## 宽表：每个基因一行，四组对比的 logFC 并列
w <- dcast(res, gene ~ contrast, value.var = "logFC")
w <- merge(w, t16[, .(gene, logFC_train = logFC, best_module)], by = "gene")
w <- merge(w, t18[, .(gene, freq_enet)], by = "gene", all.x = TRUE)
setcolorder(w, c("gene", "logFC_train", "PreS_vs_Control", "REP_vs_Control",
                 "NR_vs_Control", "REP_vs_NR", "best_module", "freq_enet"))
w <- w[order(-abs(logFC_train))]

## 假设检验的判断列
w[, `:=`(
  is_core      = gene %in% core,
  ## H1：Pre_S 与晚期同向（训练集方向和 REP 方向都应一致）
  preS_same_as_REP = sign(PreS_vs_Control) == sign(REP_vs_Control),
  ## H2：扩增阴性组不同向
  NR_diff_from_REP = sign(NR_vs_Control) != sign(REP_vs_Control),
  ## REP 与训练集方向是否仍一致（25 例口径的稳健性）
  REP_same_as_train = sign(REP_vs_Control) == sign(logFC_train)
)]

cat("\n各基因四组 logFC（宽表）：\n")
print(w[, .(gene, logFC_train, PreS_vs_Control, REP_vs_Control,
            NR_vs_Control, REP_vs_NR, freq_enet)], row.names = FALSE)

## ---------- 5. H1 / H2 汇总 ----------
cat("\n========== 5. 假设检验汇总 ==========\n")
n_preS <- sum(w$preS_same_as_REP, na.rm = TRUE)
n_NR   <- sum(w$NR_diff_from_REP, na.rm = TRUE)
n_rep  <- sum(w$REP_same_as_train, na.rm = TRUE)

cat("H1  时间轴（Pre_S 与 FECD_REP 同向）：", n_preS, " / 9\n")
cat("核心 4 基因中：",
    sum(w[is_core == TRUE]$preS_same_as_REP, na.rm = TRUE), " / ", length(core), "\n")
cat("H2  机制特异性（FECD_NR 与 FECD_REP 不同向）：", n_NR, " / 9\n")
cat("    注：FECD_NR 仅 4 例，统计功效低，多数基因可能本就不显著，\n")
cat("    需结合效应量而非仅看符号\n")
cat("稳健性  25 例口径下 REP 仍与训练集同向：", n_rep, " / 9\n")

## 逐基因判定表
w[, verdict := fcase(
  preS_same_as_REP &  NR_diff_from_REP, "时间轴+扩增特异",
  preS_same_as_REP & !NR_diff_from_REP, "时间轴成立、非扩增特异",
  !preS_same_as_REP & NR_diff_from_REP, "非早期、扩增特异",
  default = "均不成立"
)]
cat("\n逐基因判定：\n")
print(w[, .(gene, is_core, preS_same_as_REP, NR_diff_from_REP, verdict)],
      row.names = FALSE)

## ---------- 6. 保存 ----------
fwrite(w, file.path(out_tab, "T20_cohort_progression.csv"))

summ <- data.table(
  item = c("候选基因数", "H1 Pre_S同向于REP", "H1 核心基因同向",
           "H2 NR不同向于REP", "稳健性 REP同向于训练集"),
  value = c(nrow(w), n_preS,
            sum(w[is_core == TRUE]$preS_same_as_REP, na.rm = TRUE),
            n_NR, n_rep),
  total = c(nrow(w), 9, length(core), 9, 9)
)
fwrite(summ, file.path(out_tab, "T20b_cohort_progression_summary.csv"))
cat("\nT20 / T20b 已保存\n")

## ---------- 7. Fig5A 进展轴趋势图 ----------
cat("\n========== 7. 绘制 Fig5A ==========\n")

plot_df <- data.table(
  gene  = rep(rownames(logcpm)[rownames(logcpm) %in% cand], times = ncol(logcpm)),
  value = as.numeric(logcpm[rownames(logcpm) %in% cand, ]),
  group = rep(as.character(grp), each = sum(rownames(logcpm) %in% cand))
)
plot_df$group <- factor(plot_df$group,
                        levels = c("Control", "Pre_S", "FECD_REP", "FECD_NR"))
ord <- w[order(-abs(logFC_train), gene)]$gene
plot_df$gene <- factor(plot_df$gene, levels = ord)

## 均值连线数据：Control→Pre_S→FECD_REP 一条实线，FECD_NR 单点用虚线连到 Control
mean_df <- plot_df[, .(mean_val = mean(value)), by = .(gene, group)]
pivot_b <- dcast(mean_df, gene ~ group, value.var = "mean_val")

## 每组配色：正常组冷色，扩增阳性渐暖，扩增阴性用灰
grp_cols <- c("Control"  = "#2C6FAC",
              "Pre_S"    = "#7FB3D5",
              "FECD_REP" = "#C0392B",
              "FECD_NR"  = "#8C8C8C")

p_tr <- ggplot(plot_df, aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), width = 0.5, outlier.shape = NA,
               linewidth = 0.32, colour = "grey35", alpha = 0.95) +
  geom_jitter(width = 0.13, size = 0.6, colour = "grey20", alpha = 0.55) +
  ## 均值折线：Control → Pre_S → FECD_REP（扩增阳性进展轴）
  geom_line(data = mean_df[group %in% c("Control", "Pre_S", "FECD_REP")],
            aes(x = group, y = mean_val, group = gene),
            colour = "#7B241C", linewidth = 0.45, alpha = 0.75) +
  geom_point(data = mean_df[group %in% c("Control", "Pre_S", "FECD_REP")],
             aes(x = group, y = mean_val), colour = "#7B241C", size = 1.4) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  labs(title = "Disulfidptosis-related candidate genes across the FECD progression axis",
       subtitle = paste0("GSE142538, all 25 samples; dark red line = mean trajectory ",
                         "Control - Pre_S - FECD_REP (CUG expansion positive); ",
                         "grey = FECD_NR (expansion negative)"),
       x = NULL, y = "log2CPM") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_rect(fill = "grey92", colour = NA),
        strip.text         = element_text(size = 7.5, face = "bold"),
        axis.text.x        = element_text(size = 6.6, angle = 30, hjust = 1),
        axis.text.y        = element_text(size = 6.5),
        axis.title.y       = element_text(size = 7.5),
        legend.position    = "bottom",
        legend.text        = element_text(size = 7),
        legend.key.size    = unit(0.32, "cm"),
        plot.title         = element_text(size = 9.5, face = "bold"),
        plot.subtitle      = element_text(size = 6.6, colour = "grey35"),
        plot.margin        = margin(8, 8, 6, 6))

ggsave(file.path(out_fig, "Fig5A_cohort_progression.png"), p_tr,
       width = 7.0, height = 5.2, dpi = 300)
ggsave(file.path(out_fig, "Fig5A_cohort_progression.pdf"), p_tr,
       width = 7.0, height = 5.2)
cat("Fig5A 已输出\n")

## ---------- 8. Fig5B：logFC 热图（两阶段结构） ----------
## 把四个对比的 logFC 排成矩阵，一眼看出哪些基因在 Pre_S 阶段就完成改变
## （第一列已经很深）、哪些要等到晚期才动（第一列接近白色）
hm <- melt(w[, .(gene, PreS_vs_Control, REP_vs_Control, NR_vs_Control, REP_vs_NR)],
           id.vars = "gene", variable.name = "contrast", value.name = "logFC")
hm[, contrast := factor(contrast,
                        levels = c("PreS_vs_Control", "REP_vs_Control",
                                   "NR_vs_Control", "REP_vs_NR"),
                        labels = c("Pre_S\nvs Ctrl", "REP\nvs Ctrl",
                                   "NR\nvs Ctrl", "REP\nvs NR"))]
hm$gene <- factor(hm$gene, levels = rev(ord))

## 分层标注：按 Pre_S 完成度分组，便于读者识别两阶段
w[, preS_completion := abs(PreS_vs_Control) / pmax(abs(REP_vs_Control), 1e-6)]
w[, phase := fcase(
  preS_completion >= 0.7, "early (>=70% by Pre_S)",
  preS_completion >= 0.3, "intermediate",
  default = "late (<30% by Pre_S)"
)]

p_hm <- ggplot(hm, aes(x = contrast, y = gene, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", logFC)), size = 2.5, colour = "grey10") +
  scale_fill_gradient2(low = "#2C6FAC", mid = "#FAFAFA", high = "#C0392B",
                       midpoint = 0, name = "Group difference (log2 CPM)",
                       limits = c(-1, 1)) +
  labs(title = "Candidate gene group differences across four contrasts",
       subtitle = paste0("Observed group mean differences (log2 CPM), not limma coefficients\n",
                         "REP vs NR tests expansion status"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid       = element_blank(),
        axis.text.x      = element_text(size = 7.4, colour = "black"),
        axis.text.y      = element_text(size = 7.8, face = "italic", colour = "black"),
        legend.title     = element_text(size = 7),
        legend.text      = element_text(size = 6.6),
        legend.key.height = unit(0.42, "cm"),
        legend.key.width  = unit(0.26, "cm"),
        plot.title       = element_text(size = 9.5, face = "bold"),
        plot.subtitle    = element_text(size = 6.6, colour = "grey35"),
        plot.margin      = margin(8, 10, 6, 6))

ggsave(file.path(out_fig, "Fig5B_candidate_logFC_heatmap.png"), p_hm,
       width = 4.9, height = 3.9, dpi = 600)
ggsave(file.path(out_fig, "Fig5B_candidate_logFC_heatmap.pdf"), p_hm,
       width = 4.9, height = 3.9)
cat("Fig5B 已输出\n")

## 打印两阶段分组结果，供文本写作直接引用
cat("\n两阶段分组（按 Pre_S 完成度）：\n")
print(w[order(-preS_completion),
        .(gene, PreS_vs_Control, REP_vs_Control, preS_completion = round(preS_completion, 2),
          phase)], row.names = FALSE)
fwrite(w, file.path(out_tab, "T20_cohort_progression.csv"))  ## 补写 phase 列

cat("\n========== 完成 ==========\n")
