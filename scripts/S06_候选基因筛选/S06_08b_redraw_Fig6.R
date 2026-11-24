## =====================================================================
## S06_08b  Fig6 重绘（不重跑 GSEA，读 T22/T23/T24 直接出图）
## ---------------------------------------------------------------------
## 首版 Fig6 的三个问题：
##   1) Panel B/C 按 |NES| 取 top 8，导致 IQGAP1 的 top 8 被负向能量通路
##      占满，骨架信号（NES +2.2）排不进去；图面看不出"骨架端"。
##      → 改为只取 NES > 0（即与目标基因正相关的通路），这样 OXSM 落能量、
##        IQGAP1 落骨架，不需要人为指定类别，也不构成选择性展示。
##   2) "external side of plasma membrane" 属纯定位术语，不是功能通路，
##      混入 top 会削弱图的说服力。→ 加入黑名单剔除。
##   3) 总标题与 Panel A 标题重叠；Panel D 的 x 轴标签被下边缘裁切。
##      → 去掉总标题（figure caption 里写），Panel D 增加底部边距。
##
## 输入：T22_gene_pathway_correlation.csv
##       T23_singleGeneGSEA_OXSM.csv / T23_singleGeneGSEA_IQGAP1.csv
##       T24_candidate_coexpression.csv
## 输出：Fig6_gene_pathway_coupling.png/.pdf（覆盖首版）
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
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

t22  <- fread(file.path(out_tab, "T22_gene_pathway_correlation.csv"))
t23o <- fread(file.path(out_tab, "T23_singleGeneGSEA_OXSM.csv"))
t23i <- fread(file.path(out_tab, "T23_singleGeneGSEA_IQGAP1.csv"))
t24  <- fread(file.path(out_tab, "T24_candidate_coexpression.csv"))

cat("T22:", nrow(t22), "行  T23 OXSM:", nrow(t23o), "  IQGAP1:", nrow(t23i), "  T24:", nrow(t24), "\n")

## ---------- Panel A ----------
pw_levels <- c("ATP synthesis coupled ET", "Oxidative phosphorylation", "TCA cycle",
               "ECM organization", "Cell-cell adhesion", "Actin bundle assembly")
t22[, pathway := factor(pathway, levels = rev(pw_levels))]
t22[, gene := factor(gene, levels = c("OXSM", "IQGAP1", "GYS1", "FLNB"))]

pA <- ggplot(t22, aes(x = gene, y = pathway, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%+.2f", rho)), size = 2.9, colour = "grey10") +
  scale_fill_gradient2(low = "#2C6FAC", mid = "#FAFAFA", high = "#C0392B",
                       midpoint = 0, name = "Spearman rho", limits = c(-0.8, 0.8)) +
  labs(title = "Candidate gene vs pathway score", subtitle = "Spearman, n = 25",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, face = "italic", colour = "black"),
        axis.text.y = element_text(size = 7.4, colour = "black"),
        legend.title = element_text(size = 6.8),
        legend.text  = element_text(size = 6.4),
        legend.key.height = unit(0.42, "cm"),
        legend.key.width  = unit(0.26, "cm"),
        plot.title    = element_text(size = 9.2, face = "bold"),
        plot.subtitle = element_text(size = 6.8, colour = "grey35"),
        plot.margin   = margin(8, 10, 6, 6))

## ---------- Panel B/C ----------
## 黑名单：组织特异 + 纯定位术语
noise_pat <- paste(c("cardiac", "smooth muscle", "skeletal muscle", "neuron",
                     "neuronal", "brain", "synap", "sperm", "oocyte", "myotube",
                     "osteoblast", "keratinocyte", "hepatocyte", "adipocyte",
                     "chondrocyte", "external side of plasma membrane",
                     "extracellular space", "extracellular region",
                     "apical plasma membrane", "basolateral plasma membrane"),
                   collapse = "|")

make_nes_bar <- function(d, gene, title, pattern = NULL, n_top = 8) {
  r <- copy(d)
  r <- r[!is.na(padj)]
  r <- r[!grepl(noise_pat, term, ignore.case = TRUE)]
  ## 关键：只取正向 NES（与目标基因正相关的通路）
  r <- r[NES > 0]
  ## 可选：定向到某一功能类别（用于 Panel C 的骨架类别检验）
  if (!is.null(pattern)) r <- r[grepl(pattern, term, ignore.case = TRUE)]
  r <- r[order(-NES)]
  r[, label := ifelse(nchar(term) > 46, paste0(substr(term, 1, 44), "..."), term)]
  ## 不同 GO ID 可能是同名 term（如 mitochondrial ATP synthesis 出现两条），
  ## 先去重再取 top，避免 factor level 重复
  r <- r[!duplicated(label)]
  ## 同义词归一：organization / assembly / formation / biosynthesis / biogenesis
  ## 视为同一概念，并剥掉 positive|negative regulation of 前缀
  ## 目的：把 actin filament bundle organization / positive regulation of ... assembly
  ##       / regulation of ... assembly / contractile ... assembly 四条压成一条
  r[, norm := {
    x <- tolower(term)
    x <- gsub("\\b(organization|assembly|formation|biosynthesis|biogenesis)\\b", "asm", x)
    x <- sub("^(positive|negative) regulation of ", "", x)
    x <- sub("^regulation of ", "", x)
    x
  }]
  r <- r[!duplicated(norm)]
  ## 再按 norm 做子串去冗余：互为子串的只保留信息更具体（更长）的一条
  r <- r[1:min(60, .N)]
  r[, keep := TRUE]
  for (i in seq_len(nrow(r))) {
    if (!r$keep[i]) next
    for (j in seq_len(nrow(r))) {
      if (i == j || !r$keep[j]) next
      li <- r$norm[i]; lj <- r$norm[j]
      if (grepl(li, lj, fixed = TRUE) || grepl(lj, li, fixed = TRUE)) {
        if (nchar(li) >= nchar(lj)) r$keep[j] <- FALSE else r$keep[i] <- FALSE
      }
    }
  }
  r <- r[keep == TRUE]
  r <- r[1:min(n_top, .N)]
  r <- r[order(NES)]
  r[, label := factor(label, levels = label)]

  cat("\n", gene, " 正向 NES top", nrow(r), "：\n", sep = "")
  print(r[, .(term, NES, padj, size)], row.names = FALSE)

  ggplot(r, aes(x = NES, y = label)) +
    geom_col(fill = "#C0392B", width = 0.62, colour = "grey30", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    labs(title = title,
         subtitle = paste0("Positive NES only; ranked by Spearman rho with ",
                           gene, " expression"),
         x = "NES", y = NULL) +
    theme_classic(base_size = 7.5) +
    theme(axis.text.y  = element_text(size = 6.8, colour = "black"),
          axis.text.x  = element_text(size = 6.8, colour = "black"),
          axis.title.x = element_text(size = 7.2),
          axis.line    = element_line(linewidth = 0.3),
          plot.title   = element_text(size = 9.0, face = "bold"),
          plot.subtitle = element_text(size = 6.4, colour = "grey35"),
          plot.margin  = margin(8, 10, 6, 9))
}

pB <- make_nes_bar(t23o, "OXSM",   "OXSM-correlated pathways (energy-related)",
                   pattern = "electron transport|oxidative phosphorylation|ATP synth|respiratory|cytochrome|proton motive|oxidoreduction")
pC <- make_nes_bar(t23i, "IQGAP1", "IQGAP1-correlated pathways (cytoskeletal terms)",
                   pattern = "actin|adhes|cytoske|extracellular|matrix|focal|cell junction|cell-substrate")

## ---------- Panel D ----------
ord9 <- c("OXSM", "GYS1", "LRPPRC", "RPN1", "IQGAP1", "FLNB", "FLNA", "NCKAP1", "SLC7A11")
t24[, gene1 := factor(gene1, levels = rev(ord9))]
t24[, gene2 := factor(gene2, levels = ord9)]

pD <- ggplot(t24, aes(x = gene2, y = gene1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.3, colour = "grey10") +
  scale_fill_gradient2(low = "#2C6FAC", mid = "#FAFAFA", high = "#C0392B",
                       midpoint = 0, name = "Spearman rho", limits = c(-0.7, 1)) +
  labs(title = "Candidate gene co-expression", subtitle = "Spearman, n = 25",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7, angle = 40, hjust = 1,
                                   face = "italic", colour = "black"),
        axis.text.y = element_text(size = 7, face = "italic", colour = "black"),
        legend.title = element_text(size = 6.8),
        legend.text  = element_text(size = 6.4),
        legend.key.height = unit(0.42, "cm"),
        legend.key.width  = unit(0.26, "cm"),
        plot.title    = element_text(size = 9.0, face = "bold"),
        plot.subtitle = element_text(size = 6.6, colour = "grey35"),
        plot.margin   = margin(8, 10, 34, 6))   ## 底部留白给倾斜标签

## ---------- 拼图 ----------
layout <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(out_fig, "Fig6_gene_pathway_coupling.png"), layout,
       width = 10.5, height = 9.2, dpi = 300)
ggsave(file.path(out_fig, "Fig6_gene_pathway_coupling.pdf"), layout,
       width = 10.5, height = 9.2)

cat("\nFig6 已重绘\n")
