## =====================================================================
## S03_01 GSE142538 训练集差异表达：FECD_REP vs Control
## 设计：9 Control + 6 FECD_REP = 15 例；limma-voom + TMM + filterByExpr + BH
## 输出：T03_DEG_full.csv（全基因结果）
##       T03b_DEG_significant.csv（padj<0.05 & |log2FC|>1）
##       T03c_disulfidptosisGene_DEG_status.csv（19 个基因的 DEG 状态）
##       Fig1A_volcano.png / Fig1B_heatmap_top50DEG.png
## 参数：随机种子 20261121；limma 3.x；edgeR TMM
## =====================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
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
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

set.seed(20261121)

## ---------- 1. 读数据 ----------
cat("=== 读入数据 ===\n")
counts <- as.matrix(read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_counts_gene.csv"),
                             row.names = 1, check.names = FALSE))
info <- read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_sampleInfo.csv"),
                 stringsAsFactors = FALSE)
cat("全矩阵：", nrow(counts), "genes x", ncol(counts), "samples\n")

## 只用 Control + FECD_REP
sel <- info$group %in% c("Control", "FECD_REP")
samples_sel <- info$sample[sel]
group <- factor(info$group[sel], levels = c("Control", "FECD_REP"))
cat("分析样本：", length(samples_sel), "例 (", paste(table(group), collapse = " vs "), ")\n")

cnt <- counts[, samples_sel, drop = FALSE]
storage.mode(cnt) <- "integer"

## ---------- 2. 过滤与归一化 ----------
dge <- DGEList(counts = cnt, group = group)
keep <- filterByExpr(dge, group = group)
cat("\nfilterByExpr：", sum(keep), "/", nrow(dge), " 基因保留 (",
    round(sum(keep)/nrow(dge)*100, 1), "%)\n")

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat("TMM 归一化因子：\n"); print(round(dge$samples$norm.factors, 3))

## ---------- 3. voom + limma ----------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
contr <- makeContrasts(FECD_REP - Control, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

res <- topTable(fit2, number = Inf, sort.by = "P")
res$gene <- rownames(res)

cat("\n=== 差异表达结果 ===\n")
cat("总基因：", nrow(res), "\n")
cat("padj < 0.05：", sum(res$adj.P.Val < 0.05), "\n")
cat("padj < 0.05 & |log2FC| > 1：", sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1), "\n")
cat("上调（log2FC > 1）：", sum(res$adj.P.Val < 0.05 & res$logFC > 1), "\n")
cat("下调（log2FC < -1）：", sum(res$adj.P.Val < 0.05 & res$logFC < -1), "\n")

## ---------- 4. 双硫死亡基因的 DEG 状态 ----------
dsg <- read.csv(file.path(root, "02_原始数据/02_基因集/双硫死亡基因集/DisulfidptosisGenes_Liu2023_v01.csv"),
                stringsAsFactors = FALSE)
dsg_genes <- unique(dsg$gene)
hit <- intersect(dsg_genes, rownames(res))

dsg_res <- res[hit, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
dsg_res$category <- dsg$category[match(dsg_res$gene, dsg$gene)]
dsg_res <- dsg_res[order(dsg_res$P.Value), ]

cat("\n=== 双硫死亡基因 DEG 状态（", length(hit), " 个命中）===\n")
print(dsg_res, row.names = FALSE)

cat("\n达 padj<0.05：", sum(dsg_res$adj.P.Val < 0.05), "/", nrow(dsg_res), "\n")
cat("达 padj<0.05 且 |log2FC|>1：", sum(dsg_res$adj.P.Val < 0.05 & abs(dsg_res$logFC) > 1),
    "/", nrow(dsg_res), "\n")

## ---------- 5. Fig1A 火山图 ----------
res$sig <- ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) > 1,
           ifelse(res$logFC > 0, "Up", "Down"), "NS")
res$sig <- factor(res$sig, levels = c("Up", "Down", "NS"))

## 标注：双硫死亡基因 + Top10 显著
label_genes <- union(dsg_res$gene[dsg_res$adj.P.Val < 0.05],
                     head(res$gene[res$sig != "NS"], 10))
res$label <- ifelse(res$gene %in% label_genes, res$gene, NA)

p1 <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(size = 1.4, alpha = 0.75) +
  scale_color_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 30,
                  box.padding = 0.35, min.segment.length = 0,
                  color = "black", seed = 20261121) +
  labs(x = expression(log[2]*"FC (FECD_REP vs Control)"),
       y = expression(-log[10]*"(adjusted P)"),
       title = "Differential expression in FECD corneal endothelium",
       subtitle = sprintf("GSE142538 | 9 Control vs 6 FECD_REP | %d DEGs (padj<0.05, |log2FC|>1)",
                          sum(res$sig != "NS"))) +
  theme_classic(base_size = 10) +
  theme(legend.title = element_blank(), legend.position = "top",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(out_fig, "Fig1A_volcano.png"), p1, width = 5.6, height = 4.8, dpi = 300)
ggsave(file.path(out_fig, "Fig1A_volcano.pdf"), p1, width = 5.6, height = 4.8)

## ---------- 6. Fig1B Top50 DEG 热图 ----------
sig_genes <- res$gene[res$sig != "NS"]
top50 <- head(sig_genes[order(res$adj.P.Val[res$sig != "NS"])], 50)

lg <- log2(sweep(cnt, 2, colSums(cnt), "/") * 1e6 + 1)
lg_top <- lg[top50, , drop = FALSE]

## 按表达标准化
z <- t(scale(t(lg_top)))
z[z > 2] <- 2; z[z < -2] <- -2

ann <- data.frame(Group = group, row.names = samples_sel)
ann_colors <- list(Group = c(Control = "#2C7BB6", FECD_REP = "#D7191C"))

ph <- pheatmap(z,
         annotation_col = ann,
         annotation_colors = ann_colors,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 5.5,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         width = 6.2, height = 8.0,
         main = "Top 50 DEGs (FECD_REP vs Control)",
         silence = TRUE)

## 同一 gtable 分别落 png 与 pdf，保证两版视觉完全一致
png(file.path(out_fig, "Fig1B_heatmap_top50DEG.png"),
    width = 6.2, height = 8.0, units = "in", res = 300)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
dev.off()

pdf(file.path(out_fig, "Fig1B_heatmap_top50DEG.pdf"), width = 6.2, height = 8.0)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
dev.off()

## ---------- 7. 保存 ----------
write.csv(res, file.path(out_tab, "T03_DEG_full.csv"), row.names = FALSE)
write.csv(res[res$sig != "NS", ], file.path(out_tab, "T03b_DEG_significant.csv"), row.names = FALSE)
write.csv(dsg_res, file.path(out_tab, "T03c_disulfidptosisGene_DEG_status.csv"), row.names = FALSE)

cat("\n=== 完成 ===\n")
cat("T03 / T03b / T03c\nFig1A_volcano / Fig1B_heatmap_top50DEG\n")
