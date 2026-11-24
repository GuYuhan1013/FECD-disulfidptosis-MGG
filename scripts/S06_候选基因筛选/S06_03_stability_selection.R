## =====================================================================
## S06_03  稳定性选择（stability selection）
## ---------------------------------------------------------------------
## 为什么必须做这一步：
##   S06_02 的单次 LASSO 结果不可用。
##   · CV 曲线自左至右单调上升，没有 U 形，lambda.min 落在 lambda 序列
##     最左端（0.000415），等于"完全不加正则化"；
##   · 该点系数达到 OXSM -6.98、IQGAP1 +4.71、LRPPRC -4.69，
##     换算成 odds 是 100~1000 倍量级，在 n=15 的角膜内皮数据上不可能成立；
##   · 9 个基因只收缩掉 1 个（SLC7A11），等于没筛选。
##   本质是 n=15、p=9 且 9 个基因都经过三层筛选（每个都带信息），
##   L1 惩罚找不到稀疏解，模型只能靠记忆样本把训练集误差压到最低。
##
## 本步骤的思路转换：
##   不再问"这一次 LASSO 选了哪几个基因"（答案对随机种子极敏感），
##   改问"反复抽样下，哪些基因总是被选中"。
##   这就是 Meinshausen & Bühlmann (2010) 的稳定性选择。
##
## 关键设计：
##   ① 分层子抽样：每次从 9 例 Control 抽 5 例、6 例 FECD 抽 3 例（共 8 例）
##      —— 分层是必须的，否则子样本可能只含单一类别，glmnet 无法拟合
##   ② 重复 B = 300 次
##   ③ 同时跑两种惩罚：
##        alpha = 1.0（纯 LASSO，L1）
##        alpha = 0.5（弹性网，L1+L2 各半）
##      本课题 9 个基因中有 4 个 actin 骨架蛋白（FLNA/FLNB/IQGAP1/NCKAP1）
##      和 3 个能量代谢基因（OXSM/LRPPRC/GYS1），组内表达高度同步。
##      纯 L1 在共线特征上的行为是"随机留一个、其余归零"，因此入选名单
##      会随种子剧烈摆动；L2 项则让共线基因成组进入或成组退出。
##      两个 alpha 的结果并列，可以直观看出这一差别。
##   ④ lambda 固定：不取 lambda.min（那是过拟合区），而是取全样本拟合中
##      使非零特征数接近 4 个的 lambda，即要求模型在"稀疏"区间工作
##   ⑤ 入选频率 ≥ 0.8 定义为稳定特征（Meinshausen 原文建议 0.6~0.9，
##      取偏严的 0.8，宁可少留几个）
##
## 输入：T16_candidateGene_intersection.csv
##       GSE142538_counts_gene.csv / GSE142538_sampleInfo.csv
## 输出：T18_stability_selection.csv
##       Fig4B_stability_frequency.png/.pdf
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
out_fig <- file.path(root, "04_分析结果/02_图片/主图")
proc    <- file.path(root, "02_原始数据/03_处理后数据")

## ---------- 1. 重建建模矩阵（与 S06_02 完全同一口径） ----------
cat("\n========== 1. 建模矩阵 ==========\n")
t16  <- fread(file.path(out_tab, "T16_candidateGene_intersection.csv"))
cand <- t16[pass_all == TRUE, gene]

cnt  <- fread(file.path(proc, "GSE142538_counts_gene.csv"))
info <- fread(file.path(proc, "GSE142538_sampleInfo.csv"))
setnames(cnt, 1, "gene")

sel <- info[group %in% c("Control", "FECD_REP"), sample]
sel <- sel[sel %in% names(cnt)]

sub <- cnt[gene %in% cand]
mat <- as.matrix(sub[, ..sel])
rownames(mat) <- sub$gene

libsize <- colSums(cnt[, ..sel])
logcpm  <- log2(sweep(mat, 2, libsize, "/") * 1e6 + 1)

x <- t(logcpm)
y <- as.numeric(info[match(sel, sample), group] == "FECD_REP")
cat("设计矩阵：", nrow(x), "样本 × ", ncol(x), "特征\n")
cat("标签：0 =", sum(y == 0), "  1 =", sum(y == 1), "\n")

## ---------- 2. 为两种 alpha 各选一个固定 lambda ----------
## 目标：让全样本拟合下非零特征数接近 4 个（9 个的一半），迫使模型在稀疏区工作
cat("\n========== 2. 固定 lambda 的选择 ==========\n")

pick_lambda <- function(alpha, target_nz = 4) {
  fit <- glmnet(x, y, family = "binomial", alpha = alpha, standardize = TRUE)
  ## 每个 lambda 对应的非零特征数
  nz <- vapply(seq_along(fit$lambda), function(k) {
    cf <- as.matrix(fit$beta[, k])
    sum(cf != 0)
  }, integer(1))
  ## 取最接近 target 的那个（lambda 从大到小，nz 从小到大）
  k <- which.min(abs(nz - target_nz))
  list(lambda = fit$lambda[k], nz = nz[k], fit = fit, nz_seq = nz)
}

lam_l1 <- pick_lambda(1.0, target_nz = 4)
lam_en <- pick_lambda(0.5, target_nz = 4)
cat("alpha=1.0 : lambda =", signif(lam_l1$lambda, 4), " 非零 =", lam_l1$nz, "\n")
cat("alpha=0.5 : lambda =", signif(lam_en$lambda, 4), " 非零 =", lam_en$nz, "\n")

## ---------- 3. 稳定性选择主循环 ----------
cat("\n========== 3. 稳定性选择（B = 300） ==========\n")

B          <- 300
n_ctrl_sub <- 5     # 每次从 9 例对照抽 5 例
n_case_sub <- 3     # 每次从 6 例 FECD 抽 3 例

run_stability <- function(alpha, lam, B, seed) {
  set.seed(seed)
  freq <- setNames(rep(0, ncol(x)), colnames(x))
  n_fail <- 0
  for (b in seq_len(B)) {
    ## 分层抽样：保证子样本里两类都在
    idx <- c(sample(which(y == 0), n_ctrl_sub),
             sample(which(y == 1), n_case_sub))
    fit <- tryCatch(
      glmnet(x[idx, , drop = FALSE], y[idx],
             family = "binomial", alpha = alpha,
             standardize = TRUE, lambda = lam),
      error = function(e) NULL)
    if (is.null(fit)) { n_fail <- n_fail + 1; next }
    cf <- as.matrix(coef(fit))[-1, 1]     # 去掉截距
    freq <- freq + as.numeric(cf != 0)
  }
  if (n_fail > 0) cat("  拟合失败次数：", n_fail, "\n")
  freq / B
}

## 两次用相同种子，保证抽到的子样本序列一致，比较才公平
freq_l1 <- run_stability(alpha = 1.0, lam = lam_l1$lambda, B = B, seed = 20261121)
freq_en <- run_stability(alpha = 0.5, lam = lam_en$lambda, B = B, seed = 20261121)

res <- data.frame(
  gene        = names(freq_l1),
  freq_lasso  = as.numeric(freq_l1),
  freq_enet   = as.numeric(freq_en),
  stringsAsFactors = FALSE
)
res <- res[order(-res$freq_enet, -res$freq_lasso), ]

cat("\n========== 4. 入选频率 ==========\n")
print(res, row.names = FALSE)

cat("\n--- 按 alpha=0.5（弹性网）结果排序 ---\n")
for (i in seq_len(nrow(res))) {
  cat(sprintf("  %-8s  弹性网 %.2f   纯LASSO %.2f\n",
              res$gene[i], res$freq_enet[i], res$freq_lasso[i]))
}

THR <- 0.8
stable_en <- res$gene[res$freq_enet  >= THR]
stable_l1 <- res$gene[res$freq_lasso >= THR]
cat("\n频率 >=", THR, " 的稳定特征：\n")
cat("  弹性网 alpha=0.5：", length(stable_en), "个 -> ",
    paste(stable_en, collapse = ", "), "\n", sep = "")
cat("  纯LASSO alpha=1.0：", length(stable_l1), "个 -> ",
    paste(stable_l1, collapse = ", "), "\n", sep = "")
cat("  两者共有：", length(intersect(stable_en, stable_l1)), "个 -> ",
    paste(intersect(stable_en, stable_l1), collapse = ", "), "\n", sep = "")

## ---------- 5. 保存 ----------
res$stable_enet <- res$freq_enet >= THR
res$stable_lasso <- res$freq_lasso >= THR
write.csv(res, file.path(out_tab, "T18_stability_selection.csv"), row.names = FALSE)
cat("\nT18 已保存\n")

## ---------- 6. Fig4B：入选频率条形图 ----------
cat("\n========== 6. 绘制 Fig4B ==========\n")

long <- rbind(
  data.frame(gene = res$gene, method = "Elastic net (alpha=0.5)",
             freq = res$freq_enet, stringsAsFactors = FALSE),
  data.frame(gene = res$gene, method = "LASSO (alpha=1.0)",
             freq = res$freq_lasso, stringsAsFactors = FALSE)
)
## y 轴按弹性网频率排序，条形从高到低
long$gene <- factor(long$gene, levels = rev(res$gene))
long$method <- factor(long$method,
                      levels = c("Elastic net (alpha=0.5)", "LASSO (alpha=1.0)"))

p_freq <- ggplot(long, aes(x = freq, y = gene, fill = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           colour = "grey25", linewidth = 0.2) +
  geom_vline(xintercept = THR, linetype = "dashed", colour = "grey30", linewidth = 0.45) +
  annotate("text", x = THR, y = 0.6, label = "0.8", hjust = -0.25,
           size = 2.6, colour = "grey30") +
  scale_fill_manual(values = c("Elastic net (alpha=0.5)" = "#C0392B",
                               "LASSO (alpha=1.0)"      = "#8FA8C8"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Stability selection of disulfidptosis candidate genes",
       subtitle = paste0(B, " stratified subsamples (5 control + 3 FECD each);",
                         " selection frequency shown"),
       x = "Selection frequency", y = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 9.8, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        axis.text.y   = element_text(size = 7.6, colour = "black"),
        axis.text.x   = element_text(size = 7.4, colour = "black"),
        axis.title.x  = element_text(size = 8),
        axis.line     = element_line(linewidth = 0.35),
        axis.ticks    = element_line(linewidth = 0.3),
        legend.position = "bottom",
        legend.text   = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"),
        plot.margin   = margin(8, 12, 6, 6))

ggsave(file.path(out_fig, "Fig4B_stability_frequency.png"), p_freq,
       width = 5.4, height = 4.6, dpi = 300)
ggsave(file.path(out_fig, "Fig4B_stability_frequency.pdf"), p_freq,
       width = 5.4, height = 4.6)
cat("Fig4B 已输出\n")

cat("\n========== 完成 ==========\n")
