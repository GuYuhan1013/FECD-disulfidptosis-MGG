## =====================================================================
## S06_07  通路轴评分（Pathway axis scoring，GSE142538 全 25 例）
## ---------------------------------------------------------------------
## 目的：
##   S06_06 在单基因层面发现了"两阶段"结构（能量端早期完成、骨架端晚期完成）。
##   但单基因的时序可能是巧合——如果整体轴水平也呈现相反方向的时序，
##   才说明这是一个通路级别的现象，而不是几个基因的偶然。
##
## 设计：
##   按 Liu et al. 2023 的功能分类，把命中的 19 个基因归入两条轴：
##
##   Energy / anabolic axis（4 个，方向为下调）
##     GYS1    糖原合成酶 1 —— 糖原是 NADPH 的重要来源，其下调意味着还原力供应受损
##     OXSM    线粒体 FAS II —— 与 NADPH 消耗竞争
##     LRPPRC  线粒体 RNA 稳定性 —— 线粒体翻译与氧化磷酸化
##     RPN1    ER N-糖基化 —— 蛋白合成代谢负荷
##     （NUBPL 因在训练集中方向相反且不显著，不纳入主评分，仅作敏感性分析）
##
##   Cytoskeleton axis（12 个，方向为上调）
##     FLNA FLNB IQGAP1 NCKAP1 ACTB ACTN4 CAPZB CD2AP DSTN MYH9 MYH10 PDLIM1
##
##   评分方法：每个基因在全部 25 例内做 z-score 标准化（消除基因间的
##   表达量级差异），再对轴内基因取均值。这样得到的评分保留方向：
##   正值 = 相对上调，负值 = 相对下调。
##
##   为什么不用 ssGSEA：
##     19 个基因的基因集对排名型富集评分而言过小，评分对基因集内部
##     排序高度敏感、重复性差。z-score 均值法透明、可复核，且不依赖
##     基因集大小。若审稿人要求，可补充 GSVA 结果作为敏感性分析。
##
## 统计：
##   Kruskal-Wallis（四组整体差异）
##   Wilcoxon 两两比较（Control vs Pre_S / Control vs REP / Control vs NR / REP vs NR）
##   由于 n 小，p 值仅作参考，效应量（组均值差）为主要依据
##
## 输入：GSE142538_counts_gene.csv、GSE142538_sampleInfo.csv
## 输出：T21_axis_score_perSample.csv（25 行）
##       T21b_axis_score_summary.csv（组均值与检验）
##       Fig5C_pathway_axis_score.png/.pdf
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
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

## ---------- 1. 重建表达矩阵（与 S06_06 完全一致的口径） ----------
cat("\n========== 1. 重建矩阵 ==========\n")
cnt  <- fread(file.path(dat_dir, "GSE142538_counts_gene.csv"))
info <- fread(file.path(dat_dir, "GSE142538_sampleInfo.csv"))

gene_col <- names(cnt)[1]
mat <- as.matrix(cnt[, -1, with = FALSE])
rownames(mat) <- cnt[[gene_col]]
mat <- mat[, info$sample]

grp <- factor(info$group, levels = c("Control", "Pre_S", "FECD_REP", "FECD_NR"))

dge <- DGEList(counts = mat, group = grp)
keep <- filterByExpr(dge, group = grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 2)
cat("矩阵：", nrow(logcpm), " × ", ncol(logcpm), "\n")

## ---------- 2. 轴定义 ----------
cat("\n========== 2. 轴定义 ==========\n")
axis_energy <- c("GYS1", "OXSM", "LRPPRC", "RPN1")
axis_cyto   <- c("NCKAP1", "ACTB", "ACTN4", "CAPZB", "CD2AP",
                 "DSTN", "FLNA", "FLNB", "IQGAP1", "MYH9", "MYH10", "PDLIM1")

cat("能量轴（", length(axis_energy), " 个）：", paste(axis_energy, collapse = ", "), "\n")
cat("骨架轴（", length(axis_cyto), " 个）：", paste(axis_cyto, collapse = ", "), "\n")

stopifnot(all(c(axis_energy, axis_cyto) %in% rownames(logcpm)))

## ---------- 3. z-score 与轴评分 ----------
cat("\n========== 3. 计算轴评分 ==========\n")
## t(scale(t(x))) = 对每个基因（行）在 25 例内做 z 标准化
z <- t(scale(t(logcpm)))

score <- data.table(
  sample = colnames(logcpm),
  group  = as.character(grp),
  energy = colMeans(z[axis_energy, , drop = FALSE]),
  cyto   = colMeans(z[axis_cyto, , drop = FALSE])
)

## 单基因 z 方向核查：确认每个轴内基因的方向与预期一致，
## 若某个基因严重逆向抵消，需在文本中说明
cat("\n能量轴各基因的 z 均值（FECD_REP 组）：\n")
for (g in axis_energy) {
  cat(sprintf("  %-8s %+.3f\n", g, mean(z[g, grp == "FECD_REP"])))
}
cat("骨架轴各基因的 z 均值（FECD_REP 组）：\n")
for (g in axis_cyto) {
  cat(sprintf("  %-8s %+.3f\n", g, mean(z[g, grp == "FECD_REP"])))
}

## ---------- 4. 组间描述与检验 ----------
cat("\n========== 4. 组间比较 ==========\n")
cat("\n各轴在四组的均值（SD）：\n")
for (ax in c("energy", "cyto")) {
  cat("---", ax, "---\n")
  for (g in levels(grp)) {
    v <- score[group == g][[ax]]
    cat(sprintf("  %-9s %+0.3f (%.3f)  n=%d\n", g, mean(v), sd(v), length(v)))
  }
}

cat("\nKruskal-Wallis（四组）：\n")
for (ax in c("energy", "cyto")) {
  kt <- kruskal.test(score[[ax]] ~ factor(score$group, levels = levels(grp)))
  cat(sprintf("  %-8s chi2=%.2f  df=%d  p=%.4f\n",
              ax, kt$statistic, kt$parameter, kt$p.value))
}

cat("\nWilcoxon 两两比较：\n")
pairs_list <- list(c("Control", "Pre_S"), c("Control", "FECD_REP"),
                   c("Control", "FECD_NR"), c("FECD_REP", "FECD_NR"))
pair_res <- rbindlist(lapply(c("energy", "cyto"), function(ax) {
  rbindlist(lapply(pairs_list, function(pr) {
    x <- score[group == pr[1]][[ax]]
    y <- score[group == pr[2]][[ax]]
    wt <- wilcox.test(x, y)
    data.table(axis = ax, group1 = pr[1], group2 = pr[2],
               mean1 = mean(x), mean2 = mean(y), diff = mean(y) - mean(x),
               p = wt$p.value)
  }))
}))
print(pair_res, row.names = FALSE)

## ---------- 5. 保存表格 ----------
fwrite(score, file.path(out_tab, "T21_axis_score_perSample.csv"))

summary_tab <- rbindlist(lapply(c("energy", "cyto"), function(ax) {
  rbindlist(lapply(levels(grp), function(g) {
    v <- score[group == g][[ax]]
    data.table(axis = ax, group = g, n = length(v), mean = mean(v), sd = sd(v))
  }))
}))
fwrite(summary_tab, file.path(out_tab, "T21b_axis_score_summary.csv"))
fwrite(pair_res,   file.path(out_tab, "T21c_axis_score_pairwise.csv"))
cat("\nT21 / T21b / T21c 已保存\n")

## ---------- 6. 敏感性分析：加入 NUBPL ----------
## NUBPL 在训练集里方向与其余能量基因相反，纳入会抵消信号。
## 单独算一个含 NUBPL 的 5 基因版本，若结论不变则说明主评分稳健。
axis_energy5 <- c(axis_energy, "NUBPL")
score5 <- colMeans(z[axis_energy5, , drop = FALSE])
cat("\n敏感性分析（能量轴加入 NUBPL，5 基因）：\n")
for (g in levels(grp)) cat(sprintf("  %-9s %+0.3f\n", g, mean(score5[grp == g])))
cat("与主评分（4 基因）的 Spearman 相关：",
    round(cor(score$energy, score5, method = "spearman"), 3), "\n")

## ---------- 7. Fig5C ----------
cat("\n========== 7. 绘制 Fig5C ==========\n")
long <- melt(score, id.vars = c("sample", "group"),
             measure.vars = c("energy", "cyto"),
             variable.name = "axis", value.name = "zscore")
long[, axis := factor(axis, levels = c("energy", "cyto"),
                      labels = c("Energy / anabolic axis (4 genes)",
                                 "Cytoskeleton axis (12 genes)"))]
long[, group := factor(group, levels = levels(grp))]

## 组均值与 SEM，用于画趋势线
agg <- long[, .(mean_z = mean(zscore),
                sem    = sd(zscore) / sqrt(.N)), by = .(axis, group)]

grp_cols <- c("Control" = "#2C6FAC", "Pre_S" = "#7FB3D5",
              "FECD_REP" = "#C0392B", "FECD_NR" = "#8C8C8C")

## 每个 panel 的 Kruskal p，标注在图上
kp <- data.table(
  axis = factor(c("Energy / anabolic axis (4 genes)",
                  "Cytoskeleton axis (12 genes)"),
                levels = levels(long$axis)),
  p = c(kruskal.test(score$energy ~ factor(score$group, levels = levels(grp)))$p.value,
        kruskal.test(score$cyto   ~ factor(score$group, levels = levels(grp)))$p.value)
)
kp[, lab := ifelse(p < 0.001,
                   "Kruskal-Wallis p < 0.001",
                   sprintf("Kruskal-Wallis p = %.3f", p))]

p_axis <- ggplot(long, aes(x = group, y = zscore)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(aes(fill = group), width = 0.5, outlier.shape = NA,
               linewidth = 0.35, colour = "grey30", alpha = 0.95) +
  geom_jitter(width = 0.13, size = 0.9, colour = "grey25", alpha = 0.6) +
  ## 趋势线只连扩增阳性进展轴（Control - Pre_S - FECD_REP）
  geom_line(data = agg[group %in% c("Control", "Pre_S", "FECD_REP")],
            aes(x = group, y = mean_z, group = axis),
            colour = "#7B241C", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(data = agg[group %in% c("Control", "Pre_S", "FECD_REP")],
             aes(x = group, y = mean_z), colour = "#7B241C", size = 2.0,
             inherit.aes = FALSE) +
  geom_errorbar(data = agg,
                aes(x = group, ymin = mean_z - sem, ymax = mean_z + sem),
                width = 0.18, colour = "#7B241C", linewidth = 0.35,
                inherit.aes = FALSE) +
  geom_text(data = kp, aes(x = 2.5, y = Inf, label = lab),
            vjust = 1.6, size = 2.6, colour = "grey25", inherit.aes = FALSE) +
  facet_wrap(~ axis, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, name = NULL) +
  labs(title = "Pathway axis scores across the FECD progression series",
       subtitle = paste0("z-score mean of axis genes across all 25 samples; ",
                         "0 = cohort average; dashed line = no change"),
       x = NULL, y = "Axis score (mean z-score)") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_rect(fill = "grey92", colour = NA),
        strip.text         = element_text(size = 7.8, face = "bold"),
        axis.text.x        = element_text(size = 7.2),
        axis.text.y        = element_text(size = 7),
        axis.title.y       = element_text(size = 8),
        legend.position    = "bottom",
        legend.text        = element_text(size = 7),
        legend.key.size    = unit(0.32, "cm"),
        plot.title         = element_text(size = 9.5, face = "bold"),
        plot.subtitle      = element_text(size = 6.8, colour = "grey35"),
        plot.margin        = margin(10, 10, 6, 6))

ggsave(file.path(out_fig, "Fig5C_pathway_axis_score.png"), p_axis,
       width = 7.4, height = 4.0, dpi = 300)
ggsave(file.path(out_fig, "Fig5C_pathway_axis_score.pdf"), p_axis,
       width = 7.4, height = 4.0)
cat("Fig5C 已输出\n")

cat("\n========== 完成 ==========\n")
