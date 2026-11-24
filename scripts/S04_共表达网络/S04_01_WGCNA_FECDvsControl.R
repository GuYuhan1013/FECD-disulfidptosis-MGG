## =====================================================================
## S04_01 WGCNA 共表达网络：GSE142538 训练集（9 Control vs 6 FECD_REP）
## 流程（手动等价版，因 WGCNA 1.74 的 blockwiseModules 在 R 4.6 不可用）：
##   log2(CPM+1) -> MAD top 25% -> signed adjacency -> TOM -> hclust
##   -> cutreeDynamic(deepSplit=2, minClusterSize=30) -> mergeCloseModules(cutHeight=0.25)
## 输出：T04_WGCNA_gene_module.csv（基因-模块归属）
##       T05_WGCNA_module_trait.csv（模块-性状相关）
##       T06_disulfidptosisGene_kME.csv（19 基因 kME 矩阵）
##       Fig2A_softThreshold / Fig2B_dendrogram / Fig2C_moduleTrait / Fig2D_kME
## 参数：随机种子 20261121；power 由 scale-free R2 扫描确定，未达 0.80 时取最高但上限 12
## =====================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(ggrepel)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 4)

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

## ---------- 1. 读数据并构建表达矩阵 ----------
cat("=== 读入数据 ===\n")
counts <- as.matrix(read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_counts_gene.csv"),
                             row.names = 1, check.names = FALSE))
info <- read.csv(file.path(root, "02_原始数据/03_处理后数据/GSE142538_sampleInfo.csv"), stringsAsFactors = FALSE)

sel <- info$group %in% c("Control", "FECD_REP")
samples <- info$sample[sel]
grp <- factor(info$group[sel], levels = c("Control", "FECD_REP"))
cnt <- counts[, samples, drop = FALSE]

## 仅保留至少在 3 个样本中 counts > 0 的基因（与 S03 的 filterByExpr 不同）
dge_ok <- rowSums(cnt > 0) >= 3            # 简化近似：至少 3 样本表达
lib <- colSums(cnt); lib[lib <= 0] <- 1
logcpm <- log2(sweep(cnt, 2, lib, "/") * 1e6 + 1)
logcpm <- logcpm[dge_ok, , drop = FALSE]
cat("过滤后基因：", nrow(logcpm), " x 样本：", ncol(logcpm), "\n")

## ---------- 2. MAD top 25% ----------
mad_val <- apply(logcpm, 1, mad, na.rm = TRUE)
n_top <- round(nrow(logcpm) * 0.25)
top_genes <- names(sort(mad_val, decreasing = TRUE))[1:n_top]
cat("\nMAD top 25%：", n_top, " 个基因\n")

dat <- t(logcpm[top_genes, , drop = FALSE])   # 样本在前，WGCNA 要求
cat("dat 维度：", nrow(dat), " samples x", ncol(dat), " genes\n")

## 检查是否有缺失/低质量样本
gsg <- goodSamplesGenes(dat, verbose = 0)
cat("goodSamplesGenes 通过：", gsg$allOK, "\n")
if (!gsg$allOK) {
  dat <- dat[gsg$goodSamples, gsg$goodGenes]
  cat("修正后：", nrow(dat), " x", ncol(dat), "\n")
}

## ---------- 3. 软阈值选择 ----------
cat("\n=== 软阈值扫描 ===\n")
powers <- c(1:20)
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed",
                         verbose = 0, corFnc = "cor")
sft_tab <- sft$powerEstimate
r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
cat("各 power 的 signed R2：\n")
df_sft <- data.frame(power = sft$fitIndices[, 1], SFT_R2 = r2, slope = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2])

## 选第一个 R2 >= 0.80 的 power；若都不达标（15 样本常见），
## 取 R2 最高的 power，但上限设为 12 以避免网络过度软化
ok <- which(df_sft$SFT_R2 >= 0.80)
if (length(ok) > 0) {
  power_sel <- df_sft$power[min(ok)]
} else {
  best_p <- df_sft$power[which.max(df_sft$SFT_R2)]
  power_sel <- min(best_p, 12)
  cat("[提示] signed R2 未达 0.80（最高", round(max(df_sft$SFT_R2), 3),
      "），采用 power =", power_sel, "（上限 12）\n")
}
r2_sel <- df_sft$SFT_R2[df_sft$power == power_sel]
cat("选定 power =", power_sel, " (signed R2 =", round(r2_sel, 3), ")\n")

p_a <- ggplot(df_sft, aes(x = power, y = SFT_R2)) +
  geom_line(color = "grey40") + geom_point(size = 2.2, color = "#2C7BB6") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  geom_vline(xintercept = power_sel, linetype = "dotted", color = "#B2182B", linewidth = 0.5) +
  annotate("text", x = power_sel + 0.6, y = 0.30, label = paste0("power = ", power_sel),
           hjust = 0, size = 3, color = "#B2182B") +
  labs(x = "Soft-thresholding power", y = expression("Scale-free topology fit "*R^2),
       title = "Soft-threshold selection (signed network)") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 10.5, face = "bold"))
ggsave(file.path(out_fig, "Fig2A_softThreshold.png"), p_a, width = 4.4, height = 3.8, dpi = 300)
ggsave(file.path(out_fig, "Fig2A_softThreshold.pdf"), p_a, width = 4.4, height = 3.8)

## ---------- 4. 网络构建（手动等价流程） ----------
cat("\n=== 构建 adjacency 与 TOM ===\n")
adj <- adjacency(dat, power = power_sel, type = "signed")
cat("adjacency 完成，维度", dim(adj), "\n")

TOM <- TOMsimilarity(adj, TOMType = "signed", verbose = 0)
dissTOM <- 1 - TOM
cat("TOM 完成\n")

geneTree <- hclust(as.dist(dissTOM), method = "average")
cat("层次聚类完成，叶节点数：", length(geneTree$order), "\n")

## ---------- 5. 动态剪切与模块合并 ----------
cat("\n=== 动态剪切 ===\n")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30, verbose = 0)
moduleColors <- labels2colors(dynamicMods)
cat("剪切模块数：", length(unique(moduleColors)), "\n")
print(table(moduleColors))

merges <- mergeCloseModules(dat, moduleColors, cutHeight = 0.25, verbose = 0)
mergedColors <- merges$colors
## 注意：mergeCloseModules 的 $MEs 在部分版本返回 NULL，改用 moduleEigengenes 重新计算
mergedMEs <- moduleEigengenes(dat, mergedColors, verbose = 0)$eigengenes
cat("\n合并后模块数：", length(unique(mergedColors)), "\n")
print(table(mergedColors))
cat("MEs 维度：", dim(mergedMEs), "\n")

## ---------- 6. 模块-性状相关性 ----------
cat("\n=== 模块-性状相关 ===\n")
trait <- data.frame(FECD = as.numeric(grp == "FECD_REP"))
rownames(trait) <- samples

moduleTraitCor <- cor(mergedMEs, trait$FECD, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(dat))
mod_names <- rownames(moduleTraitCor)

mt <- data.frame(
  module = mod_names,
  n_genes = as.integer(table(mergedColors)[sub("^ME", "", mod_names)]),
  cor = as.numeric(moduleTraitCor),
  p_value = as.numeric(moduleTraitP),
  stringsAsFactors = FALSE
)
mt <- mt[order(-abs(mt$cor)), ]
print(mt, row.names = FALSE)

## Fig2C 模块-性状条形图
mt$module <- factor(mt$module, levels = mt$module[order(mt$cor)])
## Fig2C 模块-性状条形图（模块名与基因数合并到 y 轴，避免标签重叠）
mt$axis_label <- paste0(sub("^ME", "", mt$module), " (n=", mt$n_genes, ")")
mt$axis_label <- factor(mt$axis_label, levels = mt$axis_label[order(mt$cor)])

p_c <- ggplot(mt, aes(x = axis_label, y = cor, fill = cor)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = length(unique(mt$axis_label)) - 0.5 + 0.5, color = NA) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  coord_flip() +
  labs(x = NULL, y = "Module-trait correlation (FECD)",
       title = "Module correlation with FECD status",
       subtitle = paste0("GSE142538 | ", nrow(dat), " samples, ", ncol(dat), " genes")) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 10.5, face = "bold"),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "right")
ggsave(file.path(out_fig, "Fig2C_moduleTrait.png"), p_c, width = 5.0, height = 4.4, dpi = 300)
ggsave(file.path(out_fig, "Fig2C_moduleTrait.pdf"), p_c, width = 5.0, height = 4.4)

## Fig2B 树状图 + 模块色带
pdf(file.path(out_fig, "Fig2B_dendrogram.pdf"), width = 8, height = 4.5)
plotDendroAndColors(geneTree, mergedColors, "Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
png(file.path(out_fig, "Fig2B_dendrogram.png"), width = 2400, height = 1350, res = 300)
plotDendroAndColors(geneTree, mergedColors, "Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

## ---------- 7. 双硫死亡基因的模块归属 ----------
cat("\n=== 双硫死亡基因 kME ===\n")
dsg <- read.csv(file.path(root, "02_原始数据/02_基因集/双硫死亡基因集/DisulfidptosisGenes_Liu2023_v01.csv"), stringsAsFactors = FALSE)
dsg_genes <- intersect(unique(dsg$gene), colnames(dat))
cat("在 WGCNA 输入矩阵中的双硫死亡基因：", length(dsg_genes), "/", length(unique(dsg$gene)), "\n")
cat(paste(dsg_genes, collapse = ", "), "\n")

## 全部基因的 kME（对所有模块）
all_kME <- cor(dat, mergedMEs, use = "p")
dsg_kME <- all_kME[dsg_genes, , drop = FALSE]

## 每个基因的最大 kME 模块
dsg_best <- data.frame(
  gene = dsg_genes,
  best_module = colnames(dsg_kME)[apply(abs(dsg_kME), 1, which.max)],
  best_kME = apply(dsg_kME, 1, function(x) x[which.max(abs(x))]),
  module_color = mergedColors[match(dsg_genes, colnames(dat))],
  stringsAsFactors = FALSE
)
dsg_best$category <- dsg$category[match(dsg_best$gene, dsg$gene)]
print(dsg_best[order(-abs(dsg_best$best_kME)), ], row.names = FALSE)

## 哪些模块含双硫死亡基因
cat("\n各模块含双硫死亡基因数：\n")
tab_mod <- table(dsg_best$module_color)
print(tab_mod)

## 模块富集（双硫死亡基因 vs 该模块全部基因）
cat("\n=== 双硫死亡基因在各模块的富集 ===\n")
enr <- data.frame()
for (m in names(table(mergedColors))) {
  in_mod <- sum(mergedColors == m)
  in_dsg <- sum(dsg_best$module_color == m)
  fisher <- fisher.test(matrix(c(in_dsg, in_mod - in_dsg,
                                  length(dsg_genes) - in_dsg, ncol(dat) - in_mod - (length(dsg_genes) - in_dsg)),
                                nrow = 2))
  enr <- rbind(enr, data.frame(module = m, n_module = in_mod, n_dsg = in_dsg,
                               p_value = fisher$p.value, stringsAsFactors = FALSE))
}
enr <- enr[order(enr$p_value), ]
print(enr, row.names = FALSE)

## ---------- 8. 保存 ----------
gene_module <- data.frame(gene = colnames(dat), module = mergedColors, stringsAsFactors = FALSE)
write.csv(gene_module, file.path(out_tab, "T04_WGCNA_gene_module.csv"), row.names = FALSE)
write.csv(mt, file.path(out_tab, "T05_WGCNA_module_trait.csv"), row.names = FALSE)
write.csv(mergedMEs, file.path(out_tab, "T05c_moduleEigengenes.csv"), row.names = TRUE)
write.csv(data.frame(gene = colnames(dat), module = mergedColors, stringsAsFactors = FALSE),
          file.path(out_tab, "T05d_WGCNA_input_gene_module.csv"), row.names = FALSE)
write.csv(dsg_best, file.path(out_tab, "T06_disulfidptosisGene_module.csv"), row.names = FALSE)
write.csv(as.data.frame(dsg_kME), file.path(out_tab, "T06b_disulfidptosisGene_kME.csv"))
write.csv(enr, file.path(out_tab, "T07_module_enrichment_dsg.csv"), row.names = FALSE)
write.csv(data.frame(power = df_sft$power, signed_R2 = df_sft$SFT_R2),
          file.path(out_tab, "T05b_softThreshold.csv"), row.names = FALSE)

cat("\n=== 完成 ===\n")
cat("T04 / T05 / T05b / T06 / T06b / T07\n")
cat("Fig2A / Fig2B / Fig2C\n")
