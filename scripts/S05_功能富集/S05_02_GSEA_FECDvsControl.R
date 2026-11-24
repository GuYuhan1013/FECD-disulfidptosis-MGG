## =====================================================================
## S05_02  GSEA（全转录组排序富集）—— 替代原计划中的 KEGG 富集
## ---------------------------------------------------------------------
## 为什么做这一步：
##   S05_01 的 GO 富集是 ORA（过表达分析），只能回答
##   "turquoise / blue 这两个模块内部显著富集了什么"。
##   它有个先天局限：输入基因已经被 WGCNA 的模块划分和 MAD 筛选过，
##   相当于"先切块再看块内成分"。
##
##   GSEA 则完全不同：它把全部 16,754 个基因按 limma t 值从高到低排成一列，
##   然后问"某个通路的成员是不是整体偏在列表的顶端或底端"。
##   不设阈值、不预先筛基因，因此可以捕捉到"温和但一致的位移"。
##
##   两者视角不同，构成互补的两层证据：
##     GO ORA   → 模块内看（骨架 vs 能量）
##     GSEA     → 全转录组看（整条转录组是否朝同一方向漂移）
##
## 故事线预期（这是本脚本的验证目标）：
##   ① HALLMARK_OXIDATIVE_PHOSPHORYLATION 应在 t<0 端富集（NES<0）
##      → 与 blue 模块的 OXPHOS 富集互相印证
##   ② EMT / APICAL_JUNCTION / TGF_BETA 应在 t>0 端富集（NES>0）
##      → 对应 turquoise 的细胞黏附/运动主题
##   ③ UNFOLDED_PROTEIN_RESPONSE 应在 t>0 端富集
##      → 对应 Top DEG 里的 DDIT3(CHOP)，即 ER 应激
##   ④ TNFA_SIGNALING_VIA_NFKB 应在 t>0 端富集
##      → 对应 NFKBIA 大幅下调（NFKBIA 是 NF-kB 的抑制子，
##         它下调意味着 NF-kB 活性被释放）
##
## 输入：04_分析结果/01_结果表格/T03_DEG_full.csv
##       （列：logFC / AveExpr / t / P.Value / adj.P.Val / B / gene / sig / label）
## 输出：T13_GSEA_allPathways.csv          全部通路（含不显著的）
##       T14_GSEA_significant.csv          padj<0.05 的子集
##       T15_GSEA_keyPathways.csv          故事线关键通路专项表
##       Fig3C_GSEA_NES_barplot.png/.pdf   双向 NES 条形图
##       Fig3D_GSEA_keyPathway_curves.png/.pdf  关键通路富集曲线
##
## 参数：fgsea 1.38.0；minSize=15，maxSize=500；gseaParam=1
##       nPermSimple=10000（置换次数，决定 p 值分辨率）
##       eps=1e-10（p 值下限，避免 adaptive 模式带来的不确定性）
##       nproc=1（单核，保证结果可复现）
## 随机种子：20261121（与项目其他步骤一致）
## 日期：2026-11-21
## =====================================================================

suppressPackageStartupMessages({
  library(fgsea)       # GSEA 主算法
  library(data.table)  # 快速读写 CSV
  library(ggplot2)     # 画图
})

options(connectionObserver = NULL)

## ---------- 0. 路径与全局参数 ----------
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
out_mid <- file.path(root, "04_分析结果/03_中间文件")

set.seed(20261121)   # 固定种子：fgsea 内部置换抽样可复现

## 配色（沿用项目主图风格：上调红、下调蓝）
COL_UP   <- "#C0392B"   # FECD 上调端
COL_DOWN <- "#2C6FAC"   # FECD 下调端

## ---------- 1. 读取差异表达结果并构造排序向量 ----------
cat("\n========== 1. 读取 DEG 表 ==========\n")
deg <- fread(file.path(out_tab, "T03_DEG_full.csv"))
cat("原始行数：", nrow(deg), "\n")
cat("列名：", paste(names(deg), collapse = " | "), "\n")

## 剔除 t 值为 NA 的行（极少数基因因方差为 0 会出现 NA/Inf）
deg <- deg[!is.na(t) & is.finite(t)]
## 剔除空基因名
deg <- deg[!is.na(gene) & gene != ""]
cat("过滤后行数：", nrow(deg), "\n")

## 去重：同一 symbol 若出现多行（理论上 limma 输出唯一），保留 |t| 最大的一行
## 这一步是保险，避免 fgsea 因重复名报错
deg <- deg[order(-abs(t))]
deg <- deg[!duplicated(gene)]
cat("去重后行数：", nrow(deg), "\n")

## 按 t 值降序排列（t 正 = FECD 高表达，t 负 = FECD 低表达）
deg <- deg[order(t, decreasing = TRUE)]
stats <- setNames(deg$t, deg$gene)

cat("排序向量长度：", length(stats), "\n")
cat("t 值范围：", round(min(stats), 2), "~", round(max(stats), 2), "\n")
cat("顶端 5 个（FECD 上调最强）：", paste(head(names(stats), 5), collapse = ", "), "\n")
cat("底端 5 个（FECD 下调最强）：", paste(tail(names(stats), 5), collapse = ", "), "\n")

## 检查几个故事线关键基因在排序列表中的位置
## 位置越靠两端，说明它们在 GSEA 排序中的贡献越大
probe <- c("SLC7A11", "FLNA", "OXSM", "GYS1", "NFKBIA", "DDIT3")
probe_pos <- match(probe, names(stats))
cat("\n关键基因在排序列表中的位次（总长 ", length(stats), "）：\n", sep = "")
for (i in seq_along(probe)) {
  cat(sprintf("  %-8s t=%7.2f  位次=%6d  (%.1f%%)\n",
              probe[i],
              ifelse(is.na(probe_pos[i]), NA, stats[probe_pos[i]]),
              probe_pos[i],
              100 * probe_pos[i] / length(stats)))
}
## 说明：位次越接近 1 表示越靠 FECD 上调端，越接近总长表示越靠下调端

## ---------- 2. 获取基因集 ----------
cat("\n========== 2. 获取基因集 ==========\n")

## 开关：是否尝试联网获取 MSigDB Hallmark
## 2026-11-21 实测：msigdbr 26.1.1 首次调用需从 zenodo.org 下载约 40 MB 基因集包，
##   本机网络下传输中断，报 "Transferred a partial file ... 39571762 bytes missing"，
##   与先前 KEGG（rest.kegg.jp）超时属同一类网络问题。
## 因此本轮改走完全离线路线：由 org.Hs.eg.db 直接构建 GO-BP 基因集，再交给 fgsea。
## 日后若网络恢复、或改装内置数据版 msigdbr（<= 7.5.1），把本开关改为 TRUE 即可。
TRY_HALLMARK <- FALSE

gs_source <- NA_character_   # 记录实际使用的基因集来源，供图标题使用
pathways  <- NULL

## --- 2a. MSigDB Hallmark（可选，需联网） ---
if (TRY_HALLMARK) {
  hm <- tryCatch({
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  }, error = function(e) {
    cat("  [collection= 失败]", conditionMessage(e), "\n"); NULL
  })
  if (is.null(hm)) {
    hm <- tryCatch({
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }, error = function(e) {
      cat("  [category= 失败]", conditionMessage(e), "\n"); NULL
    })
  }
  if (!is.null(hm) && nrow(hm) > 0) {
    sym_col <- intersect(c("gene_symbol", "db_gene_symbol", "human_gene_symbol"),
                         names(hm))[1]
    pathways  <- split(hm[[sym_col]], hm$gs_name)
    gs_source <- "MSigDB Hallmark (H)"
    cat("  Hallmark 载入成功：", length(pathways), "个基因集\n")
  }
}

## --- 2b. 兜底：由 org.Hs.eg.db 直接构建 GO-BP 基因集（完全离线） ---
## 口径说明：
##   org.Hs.egGO2ALLEGS 给出 "GO 术语 -> 全部 Entrez 基因"（已含后代术语传播），
##   这与 S05_01 里 enrichGO 的默认口径一致，因此下面两层结果可以直接对话。
##   再经 org.Hs.egSYMBOL 把 Entrez 换回基因符号，与 stats 的命名对齐。
if (is.null(pathways)) {
  cat("  使用离线 GO-BP 基因集（org.Hs.eg.db）\n")
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(GO.db)
    library(AnnotationDbi)
  })

  ## (1) 取出全部 GO 术语及其可读名称，只保留 BP 分支
  all_go <- suppressMessages(
    AnnotationDbi::select(GO.db::GO.db,
                          keys    = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
                          columns = c("TERM", "ONTOLOGY"),
                          keytype = "GOID")
  )
  bp_go <- all_go[all_go$ONTOLOGY == "BP" & !is.na(all_go$TERM), ]
  cat("  GO-BP 术语数：", nrow(bp_go), "\n")

  ## (2) BP 术语 -> Entrez 基因列表（含后代传播）
  go2eg <- as.list(org.Hs.egGO2ALLEGS)
  go2eg <- go2eg[intersect(bp_go$GOID, names(go2eg))]

  ## (3) Entrez -> 基因符号
  eg2sym <- as.list(org.Hs.egSYMBOL)

  ## (4) 逐术语换成符号，并拼上可读名
  ##     用 unlist(use.names = FALSE) 而非 sapply，快一个量级
  pathways <- lapply(go2eg, function(eg) {
    s <- unlist(eg2sym[eg], use.names = FALSE)
    unique(s[!is.na(s)])
  })
  nm_map <- setNames(bp_go$TERM, bp_go$GOID)
  names(pathways) <- paste0(names(pathways), ": ", nm_map[names(pathways)])

  gs_source <- "GO Biological Process (org.Hs.eg.db, offline)"
  cat("  构建基因集：", length(pathways), "个\n")
}

## 看一眼基因集规模，确认过滤阈值是否合理
gs_size <- vapply(pathways, length, integer(1))
cat("  基因集大小：min=", min(gs_size), "  median=", median(gs_size),
    "  max=", max(gs_size), "\n", sep = "")
cat("  落在 15-500 区间的基因集：", sum(gs_size >= 15 & gs_size <= 500), "个\n")
cat("  基因集来源：", gs_source, "\n")

## ---------- 3. 运行 fgsea ----------
cat("\n========== 3. GSEA (fgsea) ==========\n")
cat("  最小基因集 15，最大 500，置换 10000 次，单核\n")

## eps=1e-10：把 p 值下限固定住，避免 adaptive 多级模式带来的不确定性
## nproc=1：单核，保证同一随机种子下结果完全可复现
## 注意：fgsea 1.38 的 fgsea() 已不再接受 verbose 参数，故此处不传
res <- fgsea(pathways    = pathways,
             stats       = stats,
             minSize     = 15,
             maxSize     = 500,
             gseaParam   = 1,
             nPermSimple = 10000,
             eps         = 1e-10,
             nproc       = 1)

## ---------- 4. 结果整理 ----------
cat("\n========== 4. 结果整理 ==========\n")
res_df <- as.data.frame(res)
## leadingEdge 是 list 列，写 CSV 前转成逗号分隔字符串
res_df$leadingEdge <- vapply(res_df$leadingEdge,
                             function(x) paste(head(x, 30), collapse = ","),
                             character(1))
res_df <- res_df[order(res_df$NES, decreasing = TRUE), ]
rownames(res_df) <- NULL

## 统一通路显示名：
##   Hallmark 去掉 "HALLMARK_" 前缀，下划线转空格，整体小写
##   GO 条目保留 "GO:xxxxxxx: Description" 格式
res_df$pathway_display <- gsub("HALLMARK_", "", res_df$pathway)
res_df$pathway_display <- gsub("_", " ", res_df$pathway_display)
res_df$pathway_display <- tolower(res_df$pathway_display)

cat("通路总数：", nrow(res_df), "\n")
cat("padj<0.05：", sum(res_df$padj < 0.05, na.rm = TRUE),
    "  其中 NES>0：", sum(res_df$padj < 0.05 & res_df$NES > 0, na.rm = TRUE),
    "  NES<0：", sum(res_df$padj < 0.05 & res_df$NES < 0, na.rm = TRUE), "\n\n")

cat("--- NES 最高的 10 条（FECD 上调端）---\n")
print(head(res_df[, c("pathway_display", "size", "NES", "padj")], 10), row.names = FALSE)

cat("\n--- NES 最低的 10 条（FECD 下调端）---\n")
print(tail(res_df[, c("pathway_display", "size", "NES", "padj")], 10), row.names = FALSE)

## ---------- 5. 故事线关键通路专项核查 ----------
cat("\n========== 5. 故事线关键通路 ==========\n")

## 按故事线的四个要素分组定义待核查通路
##（若走 gseGO 分支，这些 Hallmark 名称不会命中，脚本会自动跳过并在下方兜底）
KEY_LIST <- list(
  "能量代谢" = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                 "HALLMARK_GLYCOLYSIS",
                 "HALLMARK_HYPOXIA",
                 "HALLMARK_FATTY_ACID_METABOLISM"),
  "骨架与黏附" = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                   "HALLMARK_APICAL_JUNCTION",
                   "HALLMARK_TGF_BETA_SIGNALING",
                   "HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
  "应激与死亡" = c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                   "HALLMARK_APOPTOSIS",
                   "HALLMARK_P53_PATHWAY",
                   "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  "炎症" = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
             "HALLMARK_INFLAMMATORY_RESPONSE",
             "HALLMARK_IL6_JAK_STAT3_SIGNALING",
             "HALLMARK_COMPLEMENT")
)

key_hit <- res_df[res_df$pathway %in% unlist(KEY_LIST), ]
if (nrow(key_hit) > 0) {
  ## 标注所属故事线要素
  key_hit$story_layer <- NA_character_
  for (nm in names(KEY_LIST)) {
    key_hit$story_layer[key_hit$pathway %in% KEY_LIST[[nm]]] <- nm
  }
  key_hit <- key_hit[order(key_hit$NES, decreasing = TRUE), ]
  cat("命中的关键通路：", nrow(key_hit), "\n")
  print(key_hit[, c("story_layer", "pathway_display", "size", "NES", "padj")],
        row.names = FALSE)
  write.csv(key_hit, file.path(out_tab, "T15_GSEA_keyPathways.csv"), row.names = FALSE)
} else {
  cat("Hallmark 关键通路未命中（当前为离线 GO-BP 分支），改用关键词匹配\n")
  ## gseGO 分支下按关键词抓取，便于同样得到一份关键通路表
  kw <- c("oxidative phosphorylation", "respiration", "electron transport",
          "extracellular matrix", "cell adhesion", "cell junction", "actin",
          "unfolded protein", "response to oxidative stress", "nf-kappa",
          "inflammatory", "apoptotic")
  pat <- paste(kw, collapse = "|")
  key_hit <- res_df[grepl(pat, res_df$pathway_display, ignore.case = TRUE), ]
  key_hit <- key_hit[order(key_hit$NES, decreasing = TRUE), ]
  key_hit$story_layer <- "关键词匹配"
  cat("关键词命中的关键通路：", nrow(key_hit), "\n")
  print(head(key_hit[, c("story_layer", "pathway_display", "size", "NES", "padj")], 15),
        row.names = FALSE)
  write.csv(key_hit, file.path(out_tab, "T15_GSEA_keyPathways.csv"), row.names = FALSE)
}

## ---------- 6. 保存表格 ----------
cat("\n========== 6. 保存表格 ==========\n")
write.csv(res_df, file.path(out_tab, "T13_GSEA_allPathways.csv"), row.names = FALSE)
sig_df <- res_df[res_df$padj < 0.05, ]
write.csv(sig_df, file.path(out_tab, "T14_GSEA_significant.csv"), row.names = FALSE)
cat("T13 全部：", nrow(res_df), "行\n")
cat("T14 显著：", nrow(sig_df), "行\n")
cat("T15 关键：", nrow(key_hit), "行\n")

## ---------- 7. Fig3C：双向 NES 条形图 ----------
cat("\n========== 7. 绘制 Fig3C ==========\n")

## 选取显著通路；若显著通路超过 24 条，按 |NES| 取上下各 12 条
## 排序：NES 从高到低，保证 y 轴从上到下是"上调端 -> 下调端"
plot_df <- sig_df
if (nrow(plot_df) > 24) {
  up   <- head(plot_df[plot_df$NES > 0, ], 12)
  down <- tail(plot_df[plot_df$NES < 0, ], 12)
  plot_df <- rbind(up, down)
}
plot_df <- plot_df[order(plot_df$NES, decreasing = TRUE), ]
## 因子化 y 轴，避免 ggplot 按字母序重排
plot_df$pathway_display <- factor(plot_df$pathway_display, levels = rev(plot_df$pathway_display))
plot_df$dir <- ifelse(plot_df$NES > 0, "FECD-up", "FECD-down")

p_bar <- ggplot(plot_df, aes(x = NES, y = pathway_display, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  scale_fill_manual(values = c("FECD-up" = COL_UP, "FECD-down" = COL_DOWN), name = NULL) +
  labs(title = paste0("GSEA — Hallmark pathways (", gs_source, ")"),
       subtitle = "FECD_REP vs Control | ranked by limma moderated t",
       x = "Normalized enrichment score (NES)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        axis.text.y   = element_text(size = 7.2, colour = "black"),
        axis.line     = element_line(linewidth = 0.35),
        legend.position = "bottom",
        legend.text   = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"),
        plot.margin   = margin(6, 10, 6, 6))

ggsave(file.path(out_fig, "Fig3C_GSEA_NES_barplot.png"), p_bar,
       width = 6.2, height = 5.6, dpi = 300)
ggsave(file.path(out_fig, "Fig3C_GSEA_NES_barplot.pdf"), p_bar,
       width = 6.2, height = 5.6)
cat("Fig3C 已输出\n")

## ---------- 8. Fig3D：关键通路富集曲线 ----------
cat("\n========== 8. 绘制 Fig3D ==========\n")

## 挑选要画曲线的通路：优先取 T15 关键通路里 padj<0.05 的，
## 若不足 4 条，用全局 |NES| top 4 补齐
cand <- key_hit[key_hit$padj < 0.05, ]
if (nrow(cand) < 4) {
  sig_all <- res_df[res_df$padj < 0.05, ]
  extra   <- head(sig_all[order(-abs(sig_all$NES)), ], 4)
  cand    <- unique(rbind(cand[, intersect(names(cand), names(res_df))],
                          extra[, intersect(names(cand), names(res_df))]))
}
## 最多画 6 条，避免图面拥挤
cand <- head(cand[order(-abs(cand$NES)), ], 6)
cat("将绘制曲线的通路：\n")
print(cand[, c("pathway_display", "NES", "padj")], row.names = FALSE)

## 逐个用 fgsea::plotEnrichment 计算 running ES 曲线，提取 $data 后叠加
##（plotEnrichment 内部即 calcGseaStat，与 fgsea 主算法同源，故曲线与 NES 自洽）
es_list <- list()
for (i in seq_len(nrow(cand))) {
  pw  <- cand$pathway[i]
  lab <- cand$pathway_display[i]
  nes <- round(cand$NES[i], 2)
  pj  <- cand$padj[i]
  if (!pw %in% names(pathways)) next
  p_obj <- tryCatch(fgsea::plotEnrichment(pathways[[pw]], stats),
                    error = function(e) NULL)
  if (is.null(p_obj)) next
  d <- p_obj$data
  if (is.null(d) || !all(c("x", "y") %in% names(d))) next
  d$pathway <- sprintf("%s  (NES=%s, padj=%.1e)", lab, nes, pj)
  es_list[[length(es_list) + 1]] <- d
}

if (length(es_list) > 0) {
  es_all <- do.call(rbind, es_list)
  ## 按最终 NES 给每条曲线上色，保持与 Fig3C 一致的语义
  es_all$pathway <- factor(es_all$pathway, levels = unique(es_all$pathway))
  lev_nes <- sapply(levels(es_all$pathway), function(s) {
    as.numeric(sub(".*NES=([-0-9.]+).*", "\\1", s))
  })
  col_map <- ifelse(lev_nes > 0, COL_UP, COL_DOWN)
  names(col_map) <- levels(es_all$pathway)

  p_curve <- ggplot(es_all, aes(x = x, y = y, colour = pathway)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "dashed") +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = col_map, name = NULL) +
    labs(title = "GSEA running enrichment score — key pathways",
         subtitle = "Left = FECD-up end of ranked list; right = FECD-down end",
         x = "Gene rank (ranked by limma moderated t)", y = "Enrichment score (ES)") +
    theme_classic(base_size = 9) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 7.3, colour = "grey35"),
          axis.line     = element_line(linewidth = 0.35),
          legend.position  = "bottom",
          legend.text   = element_text(size = 6.6),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.3, "cm"),
          plot.margin   = margin(6, 10, 6, 6))
  ggsave(file.path(out_fig, "Fig3D_GSEA_keyPathway_curves.png"), p_curve,
         width = 6.4, height = 5.2, dpi = 300)
  ggsave(file.path(out_fig, "Fig3D_GSEA_keyPathway_curves.pdf"), p_curve,
         width = 6.4, height = 5.2)
  cat("Fig3D 已输出（叠加曲线版）\n")
} else {
  cat("曲线提取失败，跳过 Fig3D\n")
}

## ---------- 9. 写一份运行摘要 ----------
sink(file.path(out_mid, "S05_02_GSEA_summary.txt"))
cat("S05_02 GSEA 运行摘要\n")
cat("运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("基因集来源：", gs_source, "\n")
cat("排序统计量：limma moderated t，长度", length(stats), "\n")
cat("fgsea 参数：minSize=15 maxSize=500 nPermSimple=10000 eps=1e-10 nproc=1\n")
cat("随机种子：20261121\n\n")
cat("显著通路（padj<0.05）总数：", nrow(sig_df), "\n")
cat("  NES>0：", sum(sig_df$NES > 0), "   NES<0：", sum(sig_df$NES < 0), "\n\n")
cat("NES 最高 10 条：\n")
print(head(res_df[, c("pathway_display", "size", "NES", "padj")], 10), row.names = FALSE)
cat("\nNES 最低 10 条：\n")
print(tail(res_df[, c("pathway_display", "size", "NES", "padj")], 10), row.names = FALSE)
cat("\n故事线关键通路：\n")
print(key_hit[, c("story_layer", "pathway_display", "size", "NES", "padj")], row.names = FALSE)
sink()

cat("\n========== 完成 ==========\n")
cat("表格：T13 / T14 / T15\n")
cat("图片：Fig3C_GSEA_NES_barplot、Fig3D_GSEA_keyPathway_curves\n")
