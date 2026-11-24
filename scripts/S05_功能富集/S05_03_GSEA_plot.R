## =====================================================================
## S05_03  GSEA 结果出图（Fig3C / Fig3D 重制版）
## ---------------------------------------------------------------------
## 为什么单独拆一个出图脚本：
##   1) GSEA 主计算（S05_02）耗时约 3 分钟，改图时不必每次重跑
##   2) S05_02 首版 Fig3C 有两个问题需要修：
##      · 标题误写为 "Hallmark"，实际基因集是 GO-BP
##      · y 轴术语名过长（GO ID + 全称），把绘图区挤到画幅外
##      · 显著通路里存在大量同义条目（ECM 三条、电子传递链四条、
##        抗原提呈三条），直接取 top 会造成"同一件事画三遍"
##   3) S05_02 的 Fig3D 因 fgsea 1.38 的 plotEnrichment() 返回结构变化，
##      提取 $data 失败而跳过；本脚本改为自行计算 running ES 曲线，
##      不依赖该函数内部结构
##
## 改进点：
##   · 术语名去掉 "GO:xxxxxxx: " 前缀，过长者截断
##   · 用词袋 Jaccard 相似度去冗余（阈值 0.55），同义条目只留 NES 最强的
##   · Fig3D 固定选 6 条覆盖故事线四大要素的代表通路
##   · 曲线自算：calcGseaStat(returnAllExtREM=TRUE) → tops/bottoms → 折线
##
## 输入：T13_GSEA_allPathways.csv / T03_DEG_full.csv
## 输出：Fig3C_GSEA_NES_barplot.png/.pdf（覆盖原版）
##       Fig3D_GSEA_keyPathway_curves.png/.pdf（新增）
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
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
out_tab  <- file.path(root, "04_分析结果/01_结果表格")
out_fig  <- file.path(root, "04_分析结果/02_图片/主图")
out_mid  <- file.path(root, "04_分析结果/03_中间文件")

set.seed(20261121)

COL_UP   <- "#C0392B"
COL_DOWN <- "#2C6FAC"

## ---------- 1. 读 GSEA 结果并清洗术语名 ----------
cat("\n========== 1. 读取 GSEA 结果 ==========\n")
res_df <- fread(file.path(out_tab, "T13_GSEA_allPathways.csv"))
cat("通路总数：", nrow(res_df), "\n")

## 去掉 "GO:0006119: " 这类前缀，只保留人类可读的术语名
res_df$term <- sub("^GO:[0-9]+: ", "", res_df$pathway)
res_df$term <- tolower(res_df$term)

## 术语名过长会撑爆 y 轴，保留到 72 字符，避免关键 GO 术语被截成半句
trunc_term <- function(x, n = 72) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}
res_df$term_short <- trunc_term(res_df$term, 72)

## ---------- 2. 去冗余（词袋 Jaccard） ----------
cat("\n========== 2. 同义条目去冗余 ==========\n")

## 逻辑：按 |NES| 从大到小依次考察，判为同义则跳过。两条规则：
##   规则1 子集关系：若 A 的每个实词都作为子串出现在 B 的词里，
##          说明 A 是 B 的泛化说法（如 taxis ⊂ chemotaxis、
##          actin filament organization ⊂ actin filament bundle organization）
##   规则2 Jaccard 相似度：实词交集/并集达到阈值
## 先剔除虚词，避免 "regulation of X" 与 "regulation of X" 因 of 拉高相似度
STOPWORDS <- c("of", "via", "to", "from", "in", "the", "and", "by", "with",
               "as", "on", "for", "a", "an", "at", "or", "into", "that", "which")

## 同义词归一：organization / assembly / formation / biosynthesis 在 GO 里
## 都表示"组装、形成"，不归一的话 "nucleosome assembly" 与
## "nucleosome organization" 会被判成两个不同的条目并列出现
SYNONYMS <- c(assembly = "x_asm", organization = "x_asm", formation = "x_asm",
              biosynthesis = "x_asm", biosynthetic = "x_asm", biogenesis = "x_asm")

content_words <- function(x) {
  w <- strsplit(x, "[ ,]+")[[1]]
  w <- w[!w %in% STOPWORDS & nchar(w) > 1]
  w <- ifelse(w %in% names(SYNONYMS), SYNONYMS[w], w)
  unique(w)
}

is_subword_of <- function(a, b) {
  length(a) > 0 && length(b) > 0 &&
    all(vapply(a, function(wa) any(grepl(wa, b, fixed = TRUE)), logical(1)))
}

dedup_terms <- function(df, thr = 0.50) {
  if (nrow(df) == 0) return(df)
  keep  <- integer(0)
  words <- lapply(df$term, content_words)
  for (i in seq_len(nrow(df))) {
    dup <- FALSE
    for (j in keep) {
      if (is_subword_of(words[[i]], words[[j]]) ||
          is_subword_of(words[[j]], words[[i]])) {
        dup <- TRUE; break
      }
      inter <- length(intersect(words[[i]], words[[j]]))
      uni   <- length(union(words[[i]], words[[j]]))
      if (uni > 0 && inter / uni >= thr) { dup <- TRUE; break }

      ## 规则3 近似同义：共享 >= 2 个实词、且 NES 相差不足 0.05
      ##   例：extracellular structure organization 与
      ##       external encapsulating structure organization 是父子术语，
      ##       共享 structure/organization，NES 只差 0.02，属于同一件事
      if (inter >= 2 && abs(df$NES[i] - df$NES[j]) < 0.05) { dup <- TRUE; break }
    }
    if (!dup) keep <- c(keep, i)
  }
  df[keep, ]
}

sig <- res_df[res_df$padj < 0.05, ]
cat("显著通路：", nrow(sig), "\n")

## 剔除与角膜内皮无关的组织/细胞类型条目
## 原因：GO 注释会把通用离子通道、收缩蛋白基因一并归到
##   "cardiac muscle cell action potential"、"smooth muscle contraction"
##   这类组织特异条目上；在内皮转录组里它们没有生物学解释力，
##   属于注释体系带来的噪声，作图时剔除（方法学需注明本条筛选规则）
EXCLUDE_PAT <- paste(c("cardiac", "smooth muscle", "skeletal muscle",
                       "neuron", "neuronal", "brain", "synap",
                       "sperm", "oocyte", "myotube", "osteoblast",
                       "keratinocyte", "hepatocyte", "adipocyte",
                       "chondrocyte"), collapse = "|")
n_before <- nrow(sig)
sig <- sig[!grepl(EXCLUDE_PAT, sig$term, ignore.case = TRUE), ]
cat("剔除组织特异噪声条目：", n_before, "->", nrow(sig), "\n")

up_raw   <- sig[sig$NES > 0, ]
up_raw   <- up_raw[order(-up_raw$NES), ]
down_raw <- sig[sig$NES < 0, ]
down_raw <- down_raw[order(down_raw$NES), ]

up_dd   <- dedup_terms(up_raw,   thr = 0.50)
down_dd <- dedup_terms(down_raw, thr = 0.50)
cat("去冗余后：上调", nrow(up_raw), "->", nrow(up_dd),
    "  下调", nrow(down_raw), "->", nrow(down_dd), "\n")

## 各取 10 条
up_show   <- head(up_dd, 10)
down_show <- head(down_dd, 10)
cat("\n--- 展示（上调）---\n"); print(up_show[, c("term_short", "size", "NES", "padj")], row.names = FALSE)
cat("\n--- 展示（下调）---\n"); print(down_show[, c("term_short", "size", "NES", "padj")], row.names = FALSE)

## ---------- 3. Fig3C：双向 NES 条形图 ----------
cat("\n========== 3. 绘制 Fig3C ==========\n")

plot_df <- rbind(up_show, down_show)
plot_df <- plot_df[order(plot_df$NES, decreasing = TRUE), ]
## 截断后的术语名可能重名（两条 rna splicing 条目截断后完全一样），
## 因此 y 轴改用唯一数字 id，显示文本通过 scale_y_discrete(labels=) 映射
plot_df$uid <- factor(as.character(seq_len(nrow(plot_df))),
                      levels = rev(as.character(seq_len(nrow(plot_df)))))
lab_map <- setNames(plot_df$term_short, as.character(seq_len(nrow(plot_df))))
plot_df$dir <- ifelse(plot_df$NES > 0, "FECD-up", "FECD-down")

p_bar <- ggplot(plot_df, aes(x = NES, y = uid, fill = dir)) +
  geom_col(width = 0.70, colour = "grey25", linewidth = 0.22) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.35) +
  scale_fill_manual(values = c("FECD-up" = COL_UP, "FECD-down" = COL_DOWN), name = NULL) +
  scale_y_discrete(labels = lab_map) +
  scale_x_continuous(limits = c(-3.05, 3.05), breaks = seq(-3, 3, 1)) +
  labs(title = "GSEA — GO biological process",
       subtitle = "FECD_REP vs Control | all terms padj < 0.05 | redundant terms collapsed",
       x = "Normalized enrichment score (NES)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title      = element_text(size = 10.5, face = "bold"),
        plot.subtitle   = element_text(size = 7.2, colour = "grey35"),
        axis.text.y     = element_text(size = 7.0, colour = "black"),
        axis.text.x     = element_text(size = 7.5, colour = "black"),
        axis.title.x    = element_text(size = 8),
        axis.line       = element_line(linewidth = 0.35),
        axis.ticks      = element_line(linewidth = 0.3),
        legend.position = "bottom",
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"),
        plot.margin     = margin(8, 16, 6, 6))

ggsave(file.path(out_fig, "Fig3C_GSEA_NES_barplot.png"), p_bar,
       width = 8.2, height = 6.2, dpi = 600)
ggsave(file.path(out_fig, "Fig3C_GSEA_NES_barplot.pdf"), p_bar,
       width = 8.2, height = 6.2)
cat("Fig3C 已重制\n")

## ---------- 4. 重建 GO-BP 基因集（画曲线要用） ----------
cat("\n========== 4. 重建 GO-BP 基因集 ==========\n")

all_go <- suppressMessages(
  AnnotationDbi::select(GO.db::GO.db,
                        keys    = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
                        columns = c("TERM", "ONTOLOGY"),
                        keytype = "GOID")
)
bp_go  <- all_go[all_go$ONTOLOGY == "BP" & !is.na(all_go$TERM), ]
go2eg  <- as.list(org.Hs.egGO2ALLEGS)
go2eg  <- go2eg[intersect(bp_go$GOID, names(go2eg))]
eg2sym <- as.list(org.Hs.egSYMBOL)
pathways <- lapply(go2eg, function(eg) {
  s <- unlist(eg2sym[eg], use.names = FALSE)
  unique(s[!is.na(s)])
})
nm_map <- setNames(bp_go$TERM, bp_go$GOID)
names(pathways) <- paste0(names(pathways), ": ", nm_map[names(pathways)])
cat("基因集：", length(pathways), "个\n")

## ---------- 5. 重建排序向量 ----------
cat("\n========== 5. 重建排序向量 ==========\n")
deg <- fread(file.path(out_tab, "T03_DEG_full.csv"))
deg <- deg[!is.na(t) & is.finite(t) & !is.na(gene) & gene != ""]
deg <- deg[order(-abs(t))]
deg <- deg[!duplicated(gene)]
deg <- deg[order(t, decreasing = TRUE)]
stats <- setNames(deg$t, deg$gene)
cat("排序向量：", length(stats), "\n")

## ---------- 6. 选 6 条故事线代表通路 ----------
cat("\n========== 6. 选取代表通路 ==========\n")

## 按故事线四大要素各指定 1-2 个关键词，在显著结果里取 |NES| 最大的一条
##   能量代谢 -> oxidative phosphorylation / proton motive force-driven atp synthesis
##   细胞外基质 -> extracellular matrix organization
##   细胞黏附 -> cell-cell adhesion
##   细胞骨架 -> actin filament bundle organization
##   炎症 -> inflammatory response
REP_KW <- c(
  "oxidative phosphorylation",
  "proton motive force-driven atp synthesis",
  "extracellular matrix organization",
  "cell-cell adhesion",
  "actin filament bundle organization",
  "inflammatory response"
)

pick_one <- function(kw) {
  cand <- res_df[res_df$padj < 0.05 & grepl(kw, res_df$term, fixed = TRUE), ]
  if (nrow(cand) == 0) {
    cat("  [未命中]", kw, "\n"); return(NULL)
  }
  cand[which.max(abs(cand$NES)), ]
}
sel <- do.call(rbind, lapply(REP_KW, pick_one))
if (is.null(sel) || nrow(sel) == 0) stop("没有选到任何代表通路")
cat("选中：\n"); print(sel[, c("term", "size", "NES", "padj")], row.names = FALSE)

## ---------- 7. 自算 running enrichment score 曲线 ----------
cat("\n========== 7. 计算富集曲线 ==========\n")

## 自算 running enrichment score（GSEA 原始定义）
## 背景：fgsea 1.38 的 calcGseaStat() 已不再接受 returnAllExtREM 参数，
##       拿不到全程轨迹，故此处按经典定义自行累加：
##         命中基因贡献  + |t|^p / NR      （NR 为全部命中基因权重之和）
##         未命中基因贡献 - 1 / (N - Nh)
##       ES 曲线即上述增量的累积和；曲线最高点（正）或最低点（负）
##       应当等于 fgsea 输出的 ES，脚本会做一致性核对
run_es <- function(genes, stats, gseaParam = 1) {
  n   <- length(stats)
  hit <- names(stats) %in% genes
  nh  <- sum(hit)
  if (nh < 5 || nh >= n) return(NULL)
  w  <- abs(stats) ^ gseaParam
  nr <- sum(w[hit])
  if (nr <= 0) return(NULL)
  inc <- ifelse(hit, w / nr, -1 / (n - nh))
  es  <- cumsum(inc)
  list(x = c(0, seq_len(n)), y = c(0, es), es = es)
}

es_list <- list()
for (i in seq_len(nrow(sel))) {
  pw  <- sel$pathway[i]
  nm  <- sel$term[i]
  nes <- sel$NES[i]
  if (!pw %in% names(pathways)) { cat("  [缺基因集]", nm, "\n"); next }

  r <- tryCatch(run_es(pathways[[pw]], stats),
                error = function(e) {
                  cat("  [出错]", nm, ":", conditionMessage(e), "\n"); NULL
                })
  if (is.null(r)) next

  ## 一致性核对：自算曲线的极值应等于 fgsea 报告的 ES
  ##   正 NES 取曲线最高点，负 NES 取最低点
  my_es <- if (nes > 0) max(r$es) else min(r$es)
  cat(sprintf("  %-50s 自算ES=%6.3f  fgseaES=%6.3f\n", nm, my_es, sel$ES[i]))

  ## 图例标签：短显示名 + NES；短名只用于图例，完整 GO 术语留在图注
  short_lab <- c(
    "oxidative phosphorylation"                 = "Oxidative phosphorylation",
    "proton motive force-driven atp synthesis"  = "ATP synthesis (PMF)",
    "extracellular matrix organization"         = "ECM organization",
    "cell-cell adhesion"                        = "Cell-cell adhesion",
    "actin filament bundle organization"        = "Actin bundle organization",
    "inflammatory response"                     = "Inflammatory response")
  lab_nm <- nm
  if (nm %in% names(short_lab)) lab_nm <- short_lab[[nm]]
  d <- data.frame(x = r$x, y = r$y)
  d$label <- sprintf("%s  (NES %.2f)", lab_nm, nes)
  d$nes   <- nes
  es_list[[length(es_list) + 1]] <- d
}

if (length(es_list) == 0) stop("所有曲线计算失败")

es_all <- do.call(rbind, es_list)
## 按 NES 正负决定线色（与 Fig3C 语义一致）
labs_order <- unique(es_all$label)
lab_nes <- sapply(labs_order, function(s) {
  as.numeric(sub(".*NES ([−-]?[0-9.]+).*", "\\1", s))
})
## 兜底：若正则未匹配，直接用该 label 的首行 nes
if (any(is.na(lab_nes))) {
  for (k in seq_along(labs_order)) {
    if (is.na(lab_nes[k])) {
      lab_nes[k] <- es_all$nes[match(labs_order[k], es_all$label)][1]
    }
  }
}
col_map <- ifelse(lab_nes > 0, COL_UP, COL_DOWN)
names(col_map) <- labs_order
es_all$label <- factor(es_all$label, levels = labs_order)

## ---------- 8. Fig3D：富集曲线 ----------
cat("\n========== 8. 绘制 Fig3D ==========\n")

n_gene <- length(stats)
p_curve <- ggplot(es_all, aes(x = x, y = y, colour = label)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.28, linetype = "dashed") +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(values = col_map, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     breaks = c(1, round(n_gene / 3), round(2 * n_gene / 3), n_gene),
                     labels = c("1", format(round(n_gene / 3), big.mark = ","),
                                format(round(2 * n_gene / 3), big.mark = ","),
                                format(n_gene, big.mark = ","))) +
  annotate("text", x = n_gene * 0.02, y = max(es_all$y) * 0.92,
           label = "FECD-up", hjust = 0, size = 2.9, colour = COL_UP, fontface = "bold") +
  annotate("text", x = n_gene * 0.98, y = max(es_all$y) * 0.92,
           label = "FECD-down", hjust = 1, size = 2.9, colour = COL_DOWN, fontface = "bold") +
  labs(title = "GSEA running enrichment score — representative GO terms",
       subtitle = "Genes ranked by limma moderated t (FECD_REP vs Control); all terms padj < 0.05",
       x = "Gene rank", y = "Enrichment score (ES)") +
  theme_classic(base_size = 9) +
  theme(plot.title      = element_text(size = 10, face = "bold"),
        plot.subtitle   = element_text(size = 7.2, colour = "grey35"),
        axis.title      = element_text(size = 8),
        axis.text       = element_text(size = 7.3, colour = "black"),
        axis.line       = element_line(linewidth = 0.35),
        axis.ticks      = element_line(linewidth = 0.3),
        legend.position = "bottom",
        legend.text     = element_text(size = 6.5),
        legend.key.width  = unit(0.45, "cm"),
        legend.key.height = unit(0.30, "cm"),
        legend.spacing.y  = unit(0.05, "cm"),
        plot.margin     = margin(8, 10, 6, 6)) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

ggsave(file.path(out_fig, "Fig3D_GSEA_keyPathway_curves.png"), p_curve,
       width = 6.6, height = 5.4, dpi = 300)
ggsave(file.path(out_fig, "Fig3D_GSEA_keyPathway_curves.pdf"), p_curve,
       width = 6.6, height = 5.4)
cat("Fig3D 已输出\n")

## ---------- 9. 摘要 ----------
sink(file.path(out_mid, "S05_03_plot_summary.txt"))
cat("S05_03 GSEA 出图摘要\n")
cat("时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Fig3C 展示通路（去冗余后各 10 条）：\n")
print(rbind(up_show, down_show)[, c("term_short", "size", "NES", "padj")], row.names = FALSE)
cat("\nFig3D 代表通路：\n")
print(sel[, c("term", "size", "NES", "padj")], row.names = FALSE)
sink()

cat("\n========== 完成 ==========\n")
cat("Fig3C / Fig3D 已重制\n")
