## =====================================================================
## S05_01 turquoise / blue 模块的 GO-BP 与 KEGG 富集
## 目的：验证故事线 —— turquoise（FECD 上调）是否富集在细胞骨架/应激，
##       blue（FECD 下调）是否富集在能量代谢/氧化还原
## 输入：04_分析结果/01_结果表格/T05d_WGCNA_input_gene_module.csv
## 输出：T09_GO_BP_turquoise.csv / T10_GO_BP_blue.csv
##       T11_KEGG_turquoise.csv / T12_KEGG_blue.csv
##       Fig3A_GO_turquoise / Fig3B_GO_blue / Fig3C_KEGG_turquoise / Fig3D_KEGG_blue
## 参数：BH 校正，padj<0.05；GO 仅 BP；KEGG organism=hsa
## =====================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
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

set.seed(20261121)

## ---------- 1. 读模块基因 ----------
gm <- read.csv(file.path(out_tab, "T05d_WGCNA_input_gene_module.csv"), stringsAsFactors = FALSE)
cat("总基因：", nrow(gm), "\n"); print(table(gm$module))

turq_sym <- unique(gm$gene[gm$module == "turquoise"])
blue_sym <- unique(gm$gene[gm$module == "blue"])
cat("\nturquoise：", length(turq_sym), " blue：", length(blue_sym), "\n")

## 背景基因（WGCNA 输入的全部 6720 基因）
bg_sym <- unique(gm$gene)

## ---------- 2. ID 转换 ----------
conv <- function(sym) {
  bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}
turq_id <- conv(turq_sym)
blue_id <- conv(blue_sym)
bg_id   <- conv(bg_sym)
cat("\nID 转换：turquoise", nrow(turq_id), " blue", nrow(blue_id), " 背景", nrow(bg_id), "\n")

## ---------- 3. GO BP 富集 ----------
run_go <- function(ids, universe, label_tag) {
  cat("\n--- GO BP:", label_tag, "---\n")
  ego <- enrichGO(gene = ids, universe = universe, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                  readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("无显著 GO 结果\n"); return(NULL)
  }
  df <- as.data.frame(ego)
  cat("显著条目：", nrow(df), "\n")
  print(head(df[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "Count")], 12), row.names = FALSE)
  ego
}

ego_t <- run_go(turq_id$ENTREZID, bg_id$ENTREZID, "turquoise")
ego_b <- run_go(blue_id$ENTREZID, bg_id$ENTREZID, "blue")

## ---------- 4. KEGG 富集（尝试联网，失败则标记） ----------
run_kegg <- function(ids, universe, label_tag) {
  cat("\n--- KEGG:", label_tag, "---\n")
  res <- tryCatch({
    ek <- enrichKEGG(gene = ids, universe = universe, organism = "hsa",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05)
    if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
      df <- setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      cat("显著通路：", nrow(as.data.frame(df)), "\n")
      print(head(as.data.frame(df)[, c("Description", "GeneRatio", "p.adjust", "Count")], 10), row.names = FALSE)
      df
    } else { cat("无显著 KEGG 结果\n"); NULL }
  }, error = function(e) {
    cat("[KEGG 失败]", conditionMessage(e), "\n"); NULL
  })
  res
}

ekegg_t <- run_kegg(turq_id$ENTREZID, bg_id$ENTREZID, "turquoise")
ekegg_b <- run_kegg(blue_id$ENTREZID, bg_id$ENTREZID, "blue")

## ---------- 5. 保存表格 ----------
if (!is.null(ego_t)) write.csv(as.data.frame(ego_t), file.path(out_tab, "T09_GO_BP_turquoise.csv"), row.names = FALSE)
if (!is.null(ego_b)) write.csv(as.data.frame(ego_b), file.path(out_tab, "T10_GO_BP_blue.csv"), row.names = FALSE)
if (!is.null(ekegg_t)) write.csv(as.data.frame(ekegg_t), file.path(out_tab, "T11_KEGG_turquoise.csv"), row.names = FALSE)
if (!is.null(ekegg_b)) write.csv(as.data.frame(ekegg_b), file.path(out_tab, "T12_KEGG_blue.csv"), row.names = FALSE)

## ---------- 6. 气泡图 ----------
dotplot_wrap <- function(obj, title, n_show = 15) {
  if (is.null(obj)) return(NULL)
  p <- dotplot(obj, showCategory = n_show, font.size = 8) +
    labs(title = title) +
    theme(plot.title = element_text(size = 9.5, face = "bold"),
          axis.text.y = element_text(size = 7))
  p
}

p_go_t <- dotplot_wrap(ego_t, "GO BP — turquoise (up in FECD)")
p_go_b <- dotplot_wrap(ego_b, "GO BP — blue (down in FECD)")
p_kg_t <- dotplot_wrap(ekegg_t, "KEGG — turquoise (up in FECD)")
p_kg_b <- dotplot_wrap(ekegg_b, "KEGG — blue (down in FECD)")

if (!is.null(p_go_t)) {
  ggsave(file.path(out_fig, "Fig3A_GO_turquoise.png"), p_go_t, width = 6.4, height = 5.4, dpi = 300)
  ggsave(file.path(out_fig, "Fig3A_GO_turquoise.pdf"), p_go_t, width = 6.4, height = 5.4)
}
if (!is.null(p_go_b)) {
  ggsave(file.path(out_fig, "Fig3B_GO_blue.png"), p_go_b, width = 6.4, height = 5.4, dpi = 300)
  ggsave(file.path(out_fig, "Fig3B_GO_blue.pdf"), p_go_b, width = 6.4, height = 5.4)
}
if (!is.null(p_kg_t)) {
  ggsave(file.path(out_fig, "Fig3C_KEGG_turquoise.png"), p_kg_t, width = 6.4, height = 5.0, dpi = 300)
  ggsave(file.path(out_fig, "Fig3C_KEGG_turquoise.pdf"), p_kg_t, width = 6.4, height = 5.0)
}
if (!is.null(p_kg_b)) {
  ggsave(file.path(out_fig, "Fig3D_KEGG_blue.png"), p_kg_b, width = 6.4, height = 5.0, dpi = 300)
  ggsave(file.path(out_fig, "Fig3D_KEGG_blue.pdf"), p_kg_b, width = 6.4, height = 5.0)
}

cat("\n=== 完成 ===\n")
cat("GO/KEGG 表格与气泡图已输出（若有显著结果）\n")
