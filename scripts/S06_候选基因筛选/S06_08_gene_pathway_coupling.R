## =====================================================================
## S06_08  基因-通路耦合分析（Figure 6）
## ---------------------------------------------------------------------
## 目的：
##   Fig4 用三层证据 + 稳定性选择筛出 4 个核心基因（OXSM/IQGAP1/GYS1/FLNB），
##   但"这些基因归属哪条功能轴"此前只是基于模块隶属的推断。
##   本步骤用两种独立方法验证这一归属是否真实：
##
##   方法一（Panel A）：通路评分相关性
##     6 条显著通路各自计算 z-score 均值评分，再与 4 个核心基因的表达量
##     做 Spearman 相关。若 OXSM/GYS1 与能量通路正相关、IQGAP1/FLNB 与
##     骨架/ECM 通路正相关，则基因-通路归属成立。
##
##   方法二（Panel B/C）：单基因 GSEA
##     对目标基因，计算全转录组每个基因与它的 Spearman 相关，得到全基因排序，
##     再对该排序做 GSEA。若 OXSM 的正相关基因集中在 OXPHOS、IQGAP1 的正相关
##     基因集中在 actin/黏附，则说明这两个基因确实是相应通路的"共变核心"。
##
##   Panel D：9 个候选基因之间的共表达相关性，检验"两个阵营"是否在表达层面
##     也彼此分离。
##
## 性能优化说明（重要）：
##   第一版试图把 2 万个 Entrez 全量映射为 Symbol，再对 14,209 个基因集逐个
##   转 symbol，耗时超过 10 分钟仍未完成。现改为：
##     · GSEA 全程用 Entrez ID 作为 key（org.Hs.egGO2ALLEGS 原生就是 Entrez）
##     · 只把表达矩阵的 18,357 个 Symbol 映射为 Entrez（一次）
##     · 只有 6 条目标通路需要 Entrez → Symbol（几百个基因，瞬间完成）
##   这是本脚本能在数分钟内跑完的关键。
##
## 输入：GSE142538_counts_gene.csv、GSE142538_sampleInfo.csv、
##       T14_GSEA_significant.csv、T16_candidateGene_intersection.csv、
##       T18_stability_selection.csv
## 输出：T22_gene_pathway_correlation.csv
##       T23_singleGeneGSEA_OXSM.csv / T23_singleGeneGSEA_IQGAP1.csv
##       T24_candidate_coexpression.csv
##       Fig6_gene_pathway_coupling.png/.pdf（4 panel 拼图）
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(fgsea)
  library(ggplot2)
  library(patchwork)
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

## ---------- 1. 重建 25 例表达矩阵 ----------
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

## ---------- 2. 定位 6 条目标通路 ----------
cat("\n========== 2. 定位目标通路 ==========\n")
t14 <- fread(file.path(out_tab, "T14_GSEA_significant.csv"))
t14[, go_id := sub("^(GO:[0-9]+):.*", "\\1", pathway)]

targets <- data.table(
  go_id = c("GO:0042773", "GO:0006119", "GO:0006099",
            "GO:0043062", "GO:0098609", "GO:0061572"),
  label = c("ATP synthesis coupled ET",
            "Oxidative phosphorylation",
            "TCA cycle",
            "ECM organization",
            "Cell-cell adhesion",
            "Actin filament bundle organization"),
  axis  = c("Energy", "Energy", "Energy", "Cytoskeleton", "Cytoskeleton", "Cytoskeleton")
)

hit <- merge(targets,
             unique(t14[, .(go_id, size, NES, padj)]),
             by = "go_id", all.x = TRUE)

cat("\n目标通路（NES 来自全转录组 GSEA）：\n")
print(hit, row.names = FALSE)

## ---------- 3. 构建基因集（全程 Entrez，避免 2 万次 symbol 映射） ----------
cat("\n========== 3. 构建 GO 基因集 ==========\n")

## 3.1 表达基因 Symbol → Entrez（只映射一次，18k 个）
## 注意：org.Hs.egSYMBOL 是 Entrez→Symbol 的正向表，反向查询必须用
## mapIds(keytype = "SYMBOL", column = "ENTREZID")，否则会得到全 NA。
cat("Symbol → Entrez 映射（", nrow(logcpm), " 个表达基因）...\n", sep = "")
sym2eg <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = rownames(logcpm),
         keytype = "SYMBOL", column = "ENTREZID")
)
cat("映射成功：", sum(!is.na(sym2eg)), " / ", length(sym2eg), "\n")

## 3.2 6 条目标通路的基因集（Entrez → Symbol，仅几百个基因）
go2all <- as.list(org.Hs.egGO2ALLEGS)
cat("org.Hs.egGO2ALLEGS 条目数：", length(go2all), "\n")

pw_sets <- lapply(seq_len(nrow(hit)), function(i) {
  ids <- go2all[[hit$go_id[i]]]
  if (is.null(ids)) return(character(0))
  s <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = as.character(ids),
           keytype = "ENTREZID", column = "SYMBOL")
  )
  unique(unname(na.omit(s)))
})
names(pw_sets) <- hit$label
cat("\n各通路基因数（映射到 Symbol 后）：\n")
for (nm in names(pw_sets)) cat(sprintf("  %-30s %d\n", nm, length(pw_sets[[nm]])))

## 3.3 全 GO-BP 基因集（保持 Entrez，只保留 15–500）
sizes  <- lengths(go2all)
go_f   <- go2all[sizes >= 15 & sizes <= 500]
cat("\n可用于 GSEA 的基因集（15–500）：", length(go_f), "\n")

## ---------- 4. Panel A：通路评分 vs 核心基因 ----------
cat("\n========== 4. Panel A：通路评分相关性 ==========\n")
z <- t(scale(t(logcpm)))

pw_score <- sapply(pw_sets, function(s) {
  g <- intersect(s, rownames(z))
  if (length(g) < 5) return(rep(NA_real_, ncol(z)))
  colMeans(z[g, , drop = FALSE])
})
cat("通路评分矩阵：", nrow(pw_score), " 例 × ", ncol(pw_score), " 通路\n")

t16  <- fread(file.path(out_tab, "T16_candidateGene_intersection.csv"))
t18  <- fread(file.path(out_tab, "T18_stability_selection.csv"))
core4 <- t18[freq_enet >= 0.5, gene]
cat("核心基因：", paste(core4, collapse = ", "), "\n")

cor_long <- rbindlist(lapply(core4, function(g) {
  rbindlist(lapply(colnames(pw_score), function(p) {
    data.table(gene = g, pathway = p,
               rho = cor(logcpm[g, ], pw_score[, p], method = "spearman"))
  }))
}))
cor_long[, pathway := factor(pathway, levels = rev(targets$label))]
cor_long[, gene := factor(gene, levels = core4)]

print(dcast(cor_long, pathway ~ gene, value.var = "rho"), row.names = FALSE)
fwrite(cor_long, file.path(out_tab, "T22_gene_pathway_correlation.csv"))

pA <- ggplot(cor_long, aes(x = gene, y = pathway, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%+.2f", rho)), size = 2.8, colour = "grey10") +
  scale_fill_gradient2(low = "#2C6FAC", mid = "#FAFAFA", high = "#C0392B",
                       midpoint = 0, name = "Spearman rho", limits = c(-0.8, 0.8)) +
  labs(title = "Candidate gene vs pathway score", subtitle = "Spearman, n = 25",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, face = "italic", colour = "black"),
        axis.text.y = element_text(size = 7.2, colour = "black"),
        legend.title = element_text(size = 6.8),
        legend.text  = element_text(size = 6.4),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.24, "cm"),
        plot.title    = element_text(size = 8.6, face = "bold"),
        plot.subtitle = element_text(size = 6.6, colour = "grey35"),
        plot.margin   = margin(6, 8, 4, 6))

## ---------- 5. Panel B/C：单基因 GSEA（Entrez key） ----------
cat("\n========== 5. 单基因 GSEA ==========\n")
ok <- !is.na(sym2eg)
logcpm_eg <- logcpm[ok, , drop = FALSE]
eg_id     <- unname(sym2eg[ok])
cat("进入 GSEA 的基因：", nrow(logcpm_eg), "\n")

cat("计算全基因 rank 矩阵...\n")
rk <- t(apply(logcpm_eg, 1, rank))   ## 行名仍为 Symbol，便于按基因名索引

run_single_gsea <- function(gene) {
  cat("\n>>> 单基因 GSEA：", gene, "\n")
  y <- rk[gene, ]
  rho <- as.numeric(cor(t(rk), y))
  names(rho) <- eg_id              ## 长度与行数一致，按行序赋值
  rho <- rho[!is.na(rho)]
  stats_sorted <- sort(rho, decreasing = TRUE)
  cat("  排序值范围：", round(min(stats_sorted), 3), "~", round(max(stats_sorted), 3), "\n")

  res <- fgsea(pathways = go_f, stats = stats_sorted,
               minSize = 15, maxSize = 500,
               nPermSimple = 10000, eps = 1e-10, nproc = 1)
  res <- res[order(-abs(NES))]
  res[, gene := gene]
  ## 把 GO ID 换成可读的 term 名（只对结果行做，几十次查询，很快）
  res[, term := vapply(pathway, function(g) {
    t <- tryCatch(AnnotationDbi::Term(g), error = function(e) NA_character_)
    if (length(t) == 1 && !is.na(t)) t else g
  }, character(1))]
  cat("  padj<0.05 通路数：", sum(res$padj < 0.05, na.rm = TRUE), "\n")
  cat("  top 5：\n")
  print(head(res[, .(term, NES, padj, size)], 5), row.names = FALSE)
  fwrite(res, file.path(out_tab, paste0("T23_singleGeneGSEA_", gene, ".csv")))
  res
}

gsea_oxsm  <- run_single_gsea("OXSM")
gsea_iqgap <- run_single_gsea("IQGAP1")

## ---------- 6. Panel D：9 个候选基因共表达相关 ----------
cat("\n========== 6. Panel D：候选基因共表达 ==========\n")
cand9 <- t16[pass_all == TRUE, gene]
cm9 <- cor(t(logcpm[cand9, ]), method = "spearman")
print(round(cm9, 2))

cm_long <- as.data.table(as.table(cm9))
setnames(cm_long, c("gene1", "gene2", "rho"))
ord9 <- c(intersect(c("OXSM", "GYS1", "LRPPRC", "RPN1"), cand9),
          intersect(c("IQGAP1", "FLNB", "FLNA", "NCKAP1", "SLC7A11"), cand9))
cm_long[, gene1 := factor(gene1, levels = rev(ord9))]
cm_long[, gene2 := factor(gene2, levels = ord9)]
fwrite(cm_long, file.path(out_tab, "T24_candidate_coexpression.csv"))

pD <- ggplot(cm_long, aes(x = gene2, y = gene1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.3, colour = "grey10") +
  scale_fill_gradient2(low = "#2C6FAC", mid = "#FAFAFA", high = "#C0392B",
                       midpoint = 0, name = "Spearman rho", limits = c(-0.7, 1)) +
  labs(title = "Candidate gene co-expression", subtitle = "Spearman, n = 25",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7, angle = 35, hjust = 1,
                                   face = "italic", colour = "black"),
        axis.text.y = element_text(size = 7, face = "italic", colour = "black"),
        legend.title = element_text(size = 6.8),
        legend.text  = element_text(size = 6.4),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.24, "cm"),
        plot.title    = element_text(size = 8.6, face = "bold"),
        plot.subtitle = element_text(size = 6.6, colour = "grey35"),
        plot.margin   = margin(6, 8, 4, 6))

## ---------- 7. Panel B/C 绘图 ----------
cat("\n========== 7. 绘制 Panel B/C ==========\n")

noise_pat <- paste(c("cardiac", "smooth muscle", "skeletal muscle", "neuron",
                     "neuronal", "brain", "synap", "sperm", "oocyte", "myotube",
                     "osteoblast", "keratinocyte", "hepatocyte", "adipocyte",
                     "chondrocyte"), collapse = "|")

make_nes_bar <- function(res, gene, n_top = 8, title) {
  r <- copy(res)
  r <- r[!grepl(noise_pat, term, ignore.case = TRUE)]
  r <- r[!is.na(padj)]
  r <- r[order(-abs(NES))][1:min(n_top, .N)]
  r[, label := ifelse(nchar(term) > 44, paste0(substr(term, 1, 42), "..."), term)]
  r <- r[order(NES)]
  r[, label := factor(label, levels = label)]
  ggplot(r, aes(x = NES, y = label, fill = NES > 0)) +
    geom_col(width = 0.62, colour = "grey30", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2C6FAC"),
                      guide = "none") +
    labs(title = title,
         subtitle = paste0("Ranked by Spearman rho with ", gene, " expression"),
         x = "NES", y = NULL) +
    theme_classic(base_size = 7.5) +
    theme(axis.text.y  = element_text(size = 6.6, colour = "black"),
          axis.text.x  = element_text(size = 6.6, colour = "black"),
          axis.title.x = element_text(size = 7),
          axis.line    = element_line(linewidth = 0.3),
          plot.title   = element_text(size = 8.4, face = "bold"),
          plot.subtitle = element_text(size = 6.4, colour = "grey35"),
          plot.margin  = margin(6, 8, 4, 6))
}

pB <- make_nes_bar(gsea_oxsm,  "OXSM",   8, "OXSM-correlated pathways")
pC <- make_nes_bar(gsea_iqgap, "IQGAP1", 8, "IQGAP1-correlated pathways")

## ---------- 8. 拼图 ----------
cat("\n========== 8. 拼图 ==========\n")
layout <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A",
                  title = "Gene-pathway coupling of disulfidptosis candidates",
                  theme = theme(plot.title = element_text(size = 10, face = "bold")))

ggsave(file.path(out_fig, "Fig6_gene_pathway_coupling.png"), layout,
       width = 10.5, height = 9.0, dpi = 300)
ggsave(file.path(out_fig, "Fig6_gene_pathway_coupling.pdf"), layout,
       width = 10.5, height = 9.0)
cat("Fig6 已输出\n")

## ---------- 9. 关键数字汇总 ----------
cat("\n========== 9. 关键数字汇总 ==========\n")
cat("通路评分相关性 top 8：\n")
print(cor_long[order(-abs(rho))][1:8], row.names = FALSE)

cat("\nOXSM 单基因 GSEA —— 能量相关通路：\n")
print(gsea_oxsm[grepl("oxidat|electron transport|tricarbox|respirat", term,
                      ignore.case = TRUE)][order(-abs(NES))][1:5,
                      .(term, NES, padj)], row.names = FALSE)

cat("\nIQGAP1 单基因 GSEA —— 骨架相关通路：\n")
print(gsea_iqgap[grepl("actin|adhes|cytoske|extracellular", term,
                       ignore.case = TRUE)][order(-abs(NES))][1:5,
                       .(term, NES, padj)], row.names = FALSE)

cat("\n========== 完成 ==========\n")
