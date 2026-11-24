## =====================================================================
## S06_01  候选基因交集筛选（三层证据交叠）
## ---------------------------------------------------------------------
## 目的：把 19 个可检测的双硫死亡基因，按三层独立证据逐级收窄，
##       得到一个"站得住"的候选池，再决定要不要进机器学习。
##       先看规模，再决定方法——避免 15 个样本硬上三算法导致过拟合。
##
## 三层证据：
##   ① 差异表达层：limma padj < 0.05
##      （本项目中双硫死亡基因效应量普遍温和，未叠加 |log2FC|>1）
##   ② 网络归属层：|kME| > 0.5
##      kME 衡量基因与该模块特征向量的相关，绝对值越高说明它越"属于"这个模块；
##      0.5 是 WGCNA 常用的中等归属阈值
##   ③ 模块性状层：模块与 FECD 状态 |r| > 0.5 且 p < 0.05
##      基因再"属于"某个模块，若该模块本身与疾病无关，也没有意义
##      （lightcyan r=0.45 p=0.090、black r=0.088 p=0.756 因此被排除）
##
## 附加自洽性检验：
##   基因对疾病的推断方向 = sign(kME × 模块-性状r)
##   把它与 logFC 的符号对比，两者一致说明网络推断与差异表达同向
##
## 输入：T03c_disulfidptosisGene_DEG_status.csv
##       T08b_disulfidptosisGene_bestModule.csv
##       T05_WGCNA_module_trait.csv
## 输出：T16_candidateGene_intersection.csv
##       Fig4A_candidateEvidence_matrix.png/.pdf
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
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

## ---------- 1. 读入并合并三张表 ----------
cat("\n========== 1. 读入三张证据表 ==========\n")
deg <- fread(file.path(out_tab, "T03c_disulfidptosisGene_DEG_status.csv"))
mod <- fread(file.path(out_tab, "T08b_disulfidptosisGene_bestModule.csv"))
tr  <- fread(file.path(out_tab, "T05_WGCNA_module_trait.csv"))
cat("DEG 表：", nrow(deg), "行\n")
cat("模块归属表：", nrow(mod), "行\n")
cat("模块-性状表：", nrow(tr), "行\n")

## 以基因为键，把 DEG 结果、最佳模块、模块-性状关联合到一张表
m <- merge(deg[, .(gene, logFC, adj.P.Val, category)],
           mod[, .(gene, best_module, best_kME)],
           by = "gene", all.x = TRUE)
m <- merge(m,
           tr[, .(module, module_r = cor, module_p = p_value, module_n = n_genes)],
           by.x = "best_module", by.y = "module", all.x = TRUE)

cat("\n合并后：", nrow(m), "行 × ", ncol(m), "列\n")
cat("模块缺失：", sum(is.na(m$module_r)), "个\n")

## ---------- 2. 三层的通过标记 ----------
m[, pass_deg := adj.P.Val < 0.05]
m[, pass_kme := abs(best_kME) > 0.5]
m[, pass_mod := abs(module_r) > 0.5 & module_p < 0.05]
m[, pass_all := pass_deg & pass_kme & pass_mod]

## 方向自洽性：网络推断方向 vs 差异表达方向
m[, infer_dir := sign(best_kME * module_r)]
m[, deg_dir   := sign(logFC)]
m[, consistent := infer_dir == deg_dir]

## 按"是否入选 -> |kME| 大小"排序，便于阅读
m <- m[order(-pass_all, -abs(best_kME))]

cat("\n========== 2. 三层筛选逐步收窄 ==========\n")
cat("起点（可检测的双硫死亡基因）：", nrow(m), "个\n")
cat("  通过①DEG        ：", sum(m$pass_deg), "个\n")
cat("  通过②|kME|>0.5  ：", sum(m$pass_kme), "个\n")
cat("  通过③模块显著   ：", sum(m$pass_mod), "个\n")
cat("  三层全通过      ：", sum(m$pass_all), "个\n")

## ---------- 3. 不同口径的候选名单 ----------
cat("\n========== 3. 不同口径的候选名单 ==========\n")

## 口径 I：只看差异表达（最宽）
setI <- m[pass_deg == TRUE, ]
cat("\n[口径 I] 仅 DEG padj<0.05：", nrow(setI), "个\n")
print(setI[, .(gene, logFC, adj.P.Val, best_module, best_kME,
               module_r, consistent)], row.names = FALSE)

## 口径 II：DEG + 模块显著（不看 kME）
setII <- m[pass_deg == TRUE & pass_mod == TRUE, ]
cat("\n[口径 II] DEG + 模块显著：", nrow(setII), "个\n")
print(setII[, .(gene, best_module, best_kME, module_r)], row.names = FALSE)

## 口径 III：三层全通过（最严）
setIII <- m[pass_all == TRUE, ]
cat("\n[口径 III] 三层全通过：", nrow(setIII), "个\n")
print(setIII[, .(gene, logFC, adj.P.Val, best_module, best_kME,
                 module_r, module_p, consistent)], row.names = FALSE)

## 参考口径：仅限 turquoise / blue 两大模块（上一轮讨论中的方案）
setTB <- m[pass_deg == TRUE & best_module %in% c("MEturquoise", "MEblue"), ]
cat("\n[参考] 限定 turquoise/blue：", nrow(setTB), "个 -> ",
    paste(setTB$gene, collapse = ", "), "\n", sep = "")

## 被三层筛掉的基因，逐个说明原因，便于方法学交代
cat("\n========== 被筛掉的基因及原因 ==========\n")
out <- m[pass_all == FALSE, ]
for (i in seq_len(nrow(out))) {
  rs <- c()
  if (!out$pass_deg[i])  rs <- c(rs, sprintf("DEG不显著(padj=%.3f)", out$adj.P.Val[i]))
  if (!out$pass_kme[i])  rs <- c(rs, sprintf("|kME|=%.2f不足0.5", abs(out$best_kME[i])))
  if (!out$pass_mod[i])  rs <- c(rs, sprintf("模块%s不显著(r=%.2f,p=%.3f)",
                                             out$best_module[i], out$module_r[i], out$module_p[i]))
  cat(sprintf("  %-8s %s\n", out$gene[i], paste(rs, collapse = "；")))
}

## ---------- 4. 方向自洽性汇总 ----------
cat("\n========== 4. 方向自洽性（kME×r 与 logFC 符号） ==========\n")
n_con <- sum(m$consistent, na.rm = TRUE)
cat("19 个基因中方向一致：", n_con, " / ", nrow(m), "\n", sep = "")
if (any(!m$consistent, na.rm = TRUE)) {
  cat("不一致的基因：\n")
  print(m[consistent == FALSE, .(gene, logFC, best_module, best_kME,
                                module_r, infer_dir)], row.names = FALSE)
}

## ---------- 5. 保存表格 ----------
cat("\n========== 5. 保存 ==========\n")
save_cols <- c("gene", "category", "logFC", "adj.P.Val", "best_module",
               "best_kME", "module_r", "module_p", "module_n",
               "pass_deg", "pass_kme", "pass_mod", "pass_all",
               "infer_dir", "consistent")
write.csv(m[, ..save_cols], file.path(out_tab, "T16_candidateGene_intersection.csv"),
          row.names = FALSE)
cat("T16 已保存：", nrow(m), "行\n")

## ---------- 6. Fig4A 证据矩阵 ----------
cat("\n========== 6. 绘制 Fig4A ==========\n")

## 长格式：每个基因 × 三条标准
g_order <- m$gene                       # 已按 pass_all 与 |kME| 排好
long <- rbind(
  data.frame(gene = g_order,
             criterion = "DEG (padj<0.05)",
             pass = m$pass_deg, stringsAsFactors = FALSE),
  data.frame(gene = g_order,
             criterion = "|kME| > 0.5",
             pass = m$pass_kme, stringsAsFactors = FALSE),
  data.frame(gene = g_order,
             criterion = "module |r| > 0.5\n& p < 0.05",
             pass = m$pass_mod, stringsAsFactors = FALSE)
)
## y 轴从下到上按 g_order 倒序，让第一个基因在最上方
long$gene <- factor(long$gene, levels = rev(g_order))
long$criterion <- factor(long$criterion,
                         levels = c("DEG (padj<0.05)", "|kME| > 0.5",
                                    "module |r| > 0.5\n& p < 0.05"))

## 右侧附加一列：direction（上/下调），用文字标注
dir_df <- data.frame(
  gene = factor(m$gene, levels = rev(g_order)),
  criterion = factor("log2FC", levels = c(levels(long$criterion), "log2FC")),
  lab = sprintf("%+.2f", m$logFC)
)

p_ev <- ggplot(long, aes(x = criterion, y = gene, fill = pass)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  ## 注意：colour 必须写在 aes 之外。写进 aes 会被当成一个映射变量，
  ## ggplot 会额外生成一条 "colour" 图例，并把文字涂成默认调色板的第一个颜色
  geom_text(aes(label = ifelse(pass, "PASS", "")), colour = "white",
            size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#1B5E20", `FALSE` = "#DCDCDC"), name = NULL) +
  labs(title = "Three-layer evidence for disulfidptosis candidate genes",
       subtitle = "19 detected genes; criteria shown separately; candidate set is their intersection",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid   = element_blank(),
        axis.text.y  = element_text(size = 7.2, colour = "black"),
        axis.text.x  = element_text(size = 7.0, colour = "black"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        legend.position = "none",
        legend.text  = element_text(size = 7),
        plot.margin  = margin(8, 8, 6, 6))

ggsave(file.path(out_fig, "Fig4A_candidateEvidence_matrix.png"), p_ev,
       width = 4.6, height = 5.6, dpi = 600)
ggsave(file.path(out_fig, "Fig4A_candidateEvidence_matrix.pdf"), p_ev,
       width = 4.6, height = 5.6)
cat("Fig4A 已输出\n")

cat("\n========== 完成 ==========\n")
