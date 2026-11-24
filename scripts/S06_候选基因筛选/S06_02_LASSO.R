## =====================================================================
## S06_02  LASSO 特征筛选（9 个候选基因 -> 收缩）
## ---------------------------------------------------------------------
## 背景与取舍：
##   参考文章用了 LASSO + 随机森林 + SVM-RFE 三算法交叉。本课题训练集
##   只有 15 例（9 Control + 6 FECD_REP），特征 9 个，特征样本比 1.67:1。
##   随机森林在 15 例上每棵树只看到极少样本，特征重要性排序基本是噪声；
##   SVM-RFE 每轮递归都要重训决策边界，结果高度依赖随机种子。
##   因此本步骤只做 LASSO —— L1 正则化本来就是为 n<p 场景设计的，
##   在 15 例、9 特征上仍能给出可解释的收缩路径。
##
## 方法细节：
##   · 表达量：log2(CPM + 1)，消除文库大小差异
##   · 标准化：glmnet 内部 standardize=TRUE，让效应量不同的基因公平竞争
##   · CV：nfolds = 15（留一法），每折留 1 例，用尽样本量
##   · 评估指标：deviance（二分类默认；样本太少时 AUC 每折无法计算，
##     因为单折内可能只有正例或只有负例）
##   · 输出两个解：lambda.min（CV 误差最小）与 lambda.1se（一倍标准误内最简）
##
## 预期与解读原则（重要）：
##   若 lambda.1se 下系数全为 0，属正常现象——说明正则化强度已经超过
##   这批基因的信息量，不应强行取 lambda.min 来"凑"几个基因；
##   此时候选基因以 S06_01 的三层证据为准，LASSO 结果仅作旁证。
##
## 输入：T16_candidateGene_intersection.csv
##       GSE142538_counts_gene.csv / GSE142538_sampleInfo.csv
## 输出：T17_LASSO_coefficients.csv
##       FigS2A_LASSO_cv.png/.pdf
##       FigS2B_LASSO_coefPath.png/.pdf
## 随机种子：20261121
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
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
proc    <- file.path(root, "02_原始数据/03_处理后数据")

set.seed(20261121)

## ---------- 1. 载入候选基因 ----------
cat("\n========== 1. 候选基因 ==========\n")
t16  <- fread(file.path(out_tab, "T16_candidateGene_intersection.csv"))
cand <- t16[pass_all == TRUE, gene]
cat("三层证据全通过：", length(cand), "个\n")
cat(paste(cand, collapse = ", "), "\n")

## ---------- 2. 载入表达矩阵与样本分组 ----------
cat("\n========== 2. 表达矩阵 ==========\n")
cnt  <- fread(file.path(proc, "GSE142538_counts_gene.csv"))
info <- fread(file.path(proc, "GSE142538_sampleInfo.csv"))
setnames(cnt, 1, "gene")
cat("表达矩阵：", nrow(cnt), "基因 × ", ncol(cnt) - 1, "样本\n")

## 只用 Control 与 FECD_REP 两组（与 S03 差异表达的口径完全一致）
## Pre_S（前症状）与 FECD_NR（无重复扩增）不纳入，避免标签混杂
sel <- info[group %in% c("Control", "FECD_REP"), sample]
sel <- sel[sel %in% names(cnt)]
cat("建模样本：", length(sel), "例\n")
print(table(info[match(sel, sample), group]))

## 取 9 个候选基因的原始 counts
sub <- cnt[gene %in% cand]
cat("命中基因：", nrow(sub), " / ", length(cand), "\n")

mat <- as.matrix(sub[, ..sel])
rownames(mat) <- sub$gene

## ---------- 3. log2(CPM+1) 归一化 ----------
## 用全部基因的文库总量做分母（不是只用 9 个基因），才符合 CPM 定义
libsize <- colSums(cnt[, ..sel])
cpm     <- sweep(mat, 2, libsize, "/") * 1e6
logcpm  <- log2(cpm + 1)

cat("\nlog2CPM 范围：", round(min(logcpm), 2), "~", round(max(logcpm), 2), "\n")
cat("每个基因的均值（按组）：\n")
grp <- info[match(sel, sample), group]
for (g in rownames(logcpm)) {
  cat(sprintf("  %-8s Control=%.2f  FECD_REP=%.2f\n",
              g,
              mean(logcpm[g, grp == "Control"]),
              mean(logcpm[g, grp == "FECD_REP"])))
}

## glmnet 要求 样本 × 特征 的矩阵
x <- t(logcpm)
y <- as.numeric(grp == "FECD_REP")
cat("\n设计矩阵：", nrow(x), "样本 × ", ncol(x), "特征\n")
cat("标签分布：0 =", sum(y == 0), "  1 =", sum(y == 1), "\n")

## ---------- 4. LASSO + 留一交叉验证 ----------
cat("\n========== 4. LASSO（alpha=1，LOOCV） ==========\n")
cvfit <- cv.glmnet(x, y,
                   family      = "binomial",
                   alpha       = 1,
                   nfolds      = nrow(x),      # 留一法
                   type.measure = "deviance",
                   standardize = TRUE,
                   intercept   = TRUE)

cat("lambda 序列：", length(cvfit$lambda), "个，范围",
    signif(min(cvfit$lambda), 3), "~", signif(max(cvfit$lambda), 3), "\n")
cat("lambda.min =", signif(cvfit$lambda.min, 5),
    "  CV deviance =", round(cvfit$cvm[cvfit$lambda == cvfit$lambda.min], 4), "\n")
cat("lambda.1se =", signif(cvfit$lambda.1se, 5),
    "  CV deviance =", round(cvfit$cvm[cvfit$lambda == cvfit$lambda.1se], 4), "\n")

## ---------- 5. 两个解下的系数 ----------
cat("\n========== 5. 系数 ==========\n")
get_coef <- function(s) {
  cf <- as.matrix(coef(cvfit, s = s))
  data.frame(gene      = rownames(cf),
             coef      = as.numeric(cf[, 1]),
             stringsAsFactors = FALSE)
}
coef_min <- get_coef("lambda.min")
coef_1se <- get_coef("lambda.1se")

cat("\n--- lambda.min ---\n")
print(coef_min[coef_min$coef != 0, ], row.names = FALSE)
cat("非零特征数：", sum(coef_min$coef[-1] != 0), " / ", ncol(x), "\n")

cat("\n--- lambda.1se ---\n")
print(coef_1se[coef_1se$coef != 0, ], row.names = FALSE)
cat("非零特征数：", sum(coef_1se$coef[-1] != 0), " / ", ncol(x), "\n")

## 合并成一张表：每个基因在两个解下的系数
out <- merge(coef_min, coef_1se, by = "gene", suffixes = c("_min", "_1se"))
out <- out[out$gene != "(Intercept)", ]
out <- out[order(-abs(out$coef_min)), ]
out$selected_min <- out$coef_min != 0
out$selected_1se <- out$coef_1se != 0
write.csv(out, file.path(out_tab, "T17_LASSO_coefficients.csv"), row.names = FALSE)
cat("\nT17 已保存：", nrow(out), "个基因\n")

## ---------- 6. FigS2A：CV 曲线 ----------
cat("\n========== 6. 绘制 FigS2A ==========\n")

cv_df <- data.frame(log_lambda = log(cvfit$lambda),
                    cvm        = cvfit$cvm,
                    cvsd       = cvfit$cvsd,
                    nzero      = cvfit$nzero)

p_cv <- ggplot(cv_df, aes(x = log_lambda)) +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd),
                width = 0.03, colour = "grey65", linewidth = 0.35) +
  geom_point(aes(y = cvm), size = 0.85, colour = "#C0392B") +
  geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed",
             colour = "#2C6FAC", linewidth = 0.5) +
  geom_vline(xintercept = log(cvfit$lambda.1se), linetype = "dotted",
             colour = "grey30", linewidth = 0.5) +
  annotate("text", x = log(cvfit$lambda.min), y = max(cv_df$cvm + cv_df$cvsd),
           label = "lambda.min", hjust = 1.2, vjust = 0, size = 2.6, colour = "#2C6FAC") +
  annotate("text", x = log(cvfit$lambda.1se), y = max(cv_df$cvm + cv_df$cvsd),
           label = "lambda.1se", hjust = -0.2, vjust = 0, size = 2.6, colour = "grey30") +
  labs(title = "LASSO binomial deviance (leave-one-out CV)",
       subtitle = paste0("9 candidate genes, ", nrow(x), " samples (",
                         sum(y == 0), " control / ", sum(y == 1), " FECD)"),
       x = expression(log(lambda)), y = "Binomial deviance") +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        axis.title    = element_text(size = 8),
        axis.text     = element_text(size = 7.3, colour = "black"),
        axis.line     = element_line(linewidth = 0.35),
        plot.margin   = margin(8, 10, 6, 6))

ggsave(file.path(out_fig, "FigS2A_LASSO_cv.png"), p_cv, width = 5.2, height = 4.2, dpi = 300)
ggsave(file.path(out_fig, "FigS2A_LASSO_cv.pdf"), p_cv, width = 5.2, height = 4.2)
cat("FigS2A 已输出\n")

## ---------- 7. FigS2B：系数收缩路径 ----------
cat("\n========== 7. 绘制 FigS2B ==========\n")

## cvfit$glmnet.fit 是同一 lambda 序列上的完整拟合，直接取其 beta
beta <- as.matrix(cvfit$glmnet.fit$beta)
beta_df <- data.frame(
  gene       = rep(rownames(beta), times = ncol(beta)),
  log_lambda = rep(log(cvfit$glmnet.fit$lambda), each = nrow(beta)),
  coef       = as.numeric(beta)
)

## 给每条线配色：按最终系数绝对值从大到小依次分配固定调色板
## （不用随机数，保证每次运行颜色完全一致，图片可复现）
PAL <- c("#C0392B", "#E08214", "#B8860B", "#6B8E23", "#2E7D32",
         "#1F6F8B", "#2C6FAC", "#5B4B8A", "#8E44AD")
ord <- rownames(beta)[order(-abs(beta[, ncol(beta)]))]
pal <- setNames(PAL[seq_along(ord)], ord)

p_path <- ggplot(beta_df, aes(x = log_lambda, y = coef, colour = gene)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(title = "LASSO coefficient shrinkage path",
       subtitle = "Each line is one candidate gene; dashed vertical line = lambda.min",
       x = expression(log(lambda)), y = "Coefficient") +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        axis.title    = element_text(size = 8),
        axis.text     = element_text(size = 7.3, colour = "black"),
        axis.line     = element_line(linewidth = 0.35),
        legend.position = "right",
        legend.text   = element_text(size = 6.6),
        legend.key.height = unit(0.3, "cm"),
        plot.margin   = margin(8, 6, 6, 6))

ggsave(file.path(out_fig, "FigS2B_LASSO_coefPath.png"), p_path,
       width = 6.0, height = 4.2, dpi = 300)
ggsave(file.path(out_fig, "FigS2B_LASSO_coefPath.pdf"), p_path,
       width = 6.0, height = 4.2)
cat("FigS2B 已输出\n")

cat("\n========== 完成 ==========\n")
