# ============================================================
# 09_ptec_mito_reanalysis.R
# 功能：仅针对 PTEC (EC-PTC) 重新进行 AKI vs CTRL 差异基因与 GO(BP) 富集,
#       并使用当前版本 GO 本体剔除已废弃词条;
#       同时单独开展线粒体相关功能探究。
# 说明：不修改任何既有文件与数据, 所有输出写入 results/ptec_reanalysis/
# 输入：results/processed_data/seu_subset_reclustered.rds
#       results/tables/GO_BP_PTEC_AKI_vs_CTRL.csv(仅用于废弃词条审计)
# 输出：results/ptec_reanalysis/tables/*.csv
#       results/ptec_reanalysis/figures/*.pdf
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
out_dir  <- file.path(proj_dir, "results", "ptec_reanalysis")
fig_dir  <- file.path(out_dir, "figures")
tab_dir  <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. 注释版本与废弃词条审计 ----------
pkg_ver <- data.frame(
  package = c("R", "org.Hs.eg.db", "GO.db", "clusterProfiler", "AnnotationDbi",
              "Seurat", "SeuratObject"),
  version = c(as.character(getRversion()),
              as.character(packageVersion("org.Hs.eg.db")),
              as.character(packageVersion("GO.db")),
              as.character(packageVersion("clusterProfiler")),
              as.character(packageVersion("AnnotationDbi")),
              as.character(packageVersion("Seurat")),
              as.character(packageVersion("SeuratObject")))
)
write.csv(pkg_ver, file.path(tab_dir, "annotation_versions.csv"), row.names = FALSE)

# 旧结果中已不在当前 GO 本体中的词条(即废弃词条)
active_ids <- keys(GO.db, keytype = "GOID")
old_go <- read.csv(file.path(proj_dir, "results", "tables",
                             "GO_BP_PTEC_AKI_vs_CTRL.csv"),
                   stringsAsFactors = FALSE)
obsolete_old <- old_go %>%
  filter(!ID %in% active_ids) %>%
  select(ID, Description, p.adjust, direction) %>%
  arrange(p.adjust)
write.csv(obsolete_old, file.path(tab_dir, "obsolete_GO_terms_in_old_result.csv"),
          row.names = FALSE)
cat("旧 PTEC GO 结果中废弃词条数:", nrow(obsolete_old), "\n")

# ---------- 2. 读取对象并重新计算 PTEC 差异基因 ----------
seu <- readRDS(file.path(proj_dir, "results", "processed_data",
                         "seu_subset_reclustered.rds"))
DefaultAssay(seu) <- "RNA"
ptec <- subset(seu, subset = cell_group == "PTEC")
ptec$group <- factor(ptec$group, levels = c("CTRL", "AKI"))
Idents(ptec) <- "group"
cat("PTEC 细胞数:", ncol(ptec), "; CTRL:", sum(ptec$group == "CTRL"),
    "; AKI:", sum(ptec$group == "AKI"), "\n")

deg <- FindMarkers(ptec, ident.1 = "AKI", ident.2 = "CTRL",
                   test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
deg$gene      <- rownames(deg)
deg$cell_type <- "PTEC"
deg <- deg %>%
  mutate(direction = case_when(
    p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "Up in AKI",
    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down in AKI",
    TRUE ~ "NS")) %>%
  arrange(p_val_adj)
write.csv(deg, file.path(tab_dir, "PTEC_DEG_AKI_vs_CTRL_updated.csv"),
          row.names = FALSE)
cat("PTEC 显著差异基因: Up =", sum(deg$direction == "Up in AKI"),
    ", Down =", sum(deg$direction == "Down in AKI"),
    ", 全部检验基因 =", nrow(deg), "\n")

universe_entrez <- unique(bitr(deg$gene, fromType = "SYMBOL", toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)$ENTREZID)

# ---------- 3. 重新做 GO(BP) 富集(当前注释 + 剔除废弃词条) ----------
run_go <- function(deg_df, direction) {
  if (direction == "Up_in_AKI") {
    genes <- deg_df$gene[deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC > 0.25]
  } else {
    genes <- deg_df$gene[deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC < -0.25]
  }
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)$ENTREZID
  ego <- enrichGO(gene = ids, universe = universe_entrez,
                  OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  df <- as.data.frame(ego)
  df$direction <- direction
  df$cell_type <- "PTEC"
  df
}

go_up   <- run_go(deg, "Up_in_AKI")
go_down <- run_go(deg, "Down_in_AKI")
go_all  <- bind_rows(go_up, go_down)
n_before <- nrow(go_all)
# 双重保障: 仅保留当前 GO 本体中仍存在的条目
go_all <- go_all %>% filter(ID %in% active_ids)
cat("GO 富集条目:", n_before, "-> 剔除废弃词条后:", nrow(go_all), "\n")
go_all$is_mito <- grepl("mitochondri", go_all$Description, ignore.case = TRUE) |
  grepl("oxidative phosphorylation|respiratory chain|ATP synth", go_all$Description,
        ignore.case = TRUE)
write.csv(go_all, file.path(tab_dir, "PTEC_GO_BP_AKI_vs_CTRL_updated.csv"),
          row.names = FALSE)

# ---------- 4. 线粒体相关 GO(BP) 条目的定向检验 ----------
# 取当前 GO 本体中所有名称含 mitochondrial 的 BP 条目,
# 对上调/下调基因分别做超几何检验(不预先限制显著性, 便于完整探究)
bp_map <- AnnotationDbi::select(GO.db,
                                keys = keys(GO.db, keytype = "GOID"),
                                columns = c("TERM", "ONTOLOGY"),
                                keytype = "GOID") %>%
  filter(ONTOLOGY == "BP")
mito_terms <- bp_map %>%
  filter(grepl("mitochondri", TERM, ignore.case = TRUE)) %>%
  distinct(GOID, TERM)

# GO 条目 -> 基因(ENTREZ) 映射, 限制在本次检验的 universe 内
go2gene <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = mito_terms$GOID,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "GO")
go2gene <- go2gene %>%
  rename(GOID = GO) %>%
  filter(ENTREZID %in% universe_entrez) %>%
  distinct(GOID, ENTREZID, SYMBOL)

deg_flag <- deg %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys = gene, keytype = "SYMBOL",
                         column = "ENTREZID", multiVals = "first")) %>%
  filter(!is.na(entrez))

hyper_test <- function(term_id, term_name, dir_label) {
  in_term <- unique(go2gene$ENTREZID[go2gene$GOID == term_id])
  n_univ  <- length(universe_entrez)
  if (length(in_term) < 5) return(NULL)
  is_deg <- deg_flag$entrez %in% (if (dir_label == "Up_in_AKI") {
    deg_flag$entrez[deg_flag$p_val_adj < 0.05 & deg_flag$avg_log2FC > 0.25]
  } else {
    deg_flag$entrez[deg_flag$p_val_adj < 0.05 & deg_flag$avg_log2FC < -0.25]
  })
  a <- sum(is_deg & deg_flag$entrez %in% in_term)
  b <- sum(is_deg & !deg_flag$entrez %in% in_term)
  c <- sum(!is_deg & deg_flag$entrez %in% in_term)
  d <- sum(!is_deg & !deg_flag$entrez %in% in_term)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  overlap <- deg_flag$gene[is_deg & deg_flag$entrez %in% in_term]
  data.frame(GOID = term_id, Description = term_name, direction = dir_label,
             n_term = length(in_term), n_overlap = a,
             gene_ratio = a / length(in_term),
             odds_ratio = unname(ft$estimate), p_value = ft$p.value,
             overlap_genes = paste(sort(unique(overlap)), collapse = "/"))
}

mito_tests <- bind_rows(lapply(c("Up_in_AKI", "Down_in_AKI"), function(d) {
  bind_rows(lapply(seq_len(nrow(mito_terms)), function(i) {
    hyper_test(mito_terms$GOID[i], mito_terms$TERM[i], d)
  }))
})) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  arrange(direction, p_adjust)
write.csv(mito_tests, file.path(tab_dir, "PTEC_GO_BP_mitochondrial_terms.csv"),
          row.names = FALSE)
cat("线粒体相关 BP 条目检验完成:", nrow(mito_tests), "条(双向)\n")

# ---------- 5. 线粒体基因集整体检验 ----------
# 采用经典 GO 基因集: 线粒体(CC)、氧化磷酸化、线粒体 ATP 合成、线粒体翻译
# 取 GO 条目的全部子孙条目(兼容新旧 GO.db: 新版为 Bimap 对象而非函数)
get_offspring <- function(ids, bimap) {
  lst <- as.list(bimap)
  unique(unlist(lst[ids], use.names = FALSE))
}

mito_sets <- list(
  "Mitochondrion (GO:0005739)"            = c("GO:0005739",
                                              get_offspring("GO:0005739", GOCCOFFSPRING)),
  "Oxidative phosphorylation (GO:0006119)" = c("GO:0006119",
                                               get_offspring("GO:0006119", GOBPOFFSPRING)),
  "Mitochondrial ATP synthesis (GO:0042776)" = c("GO:0042776",
                                                 get_offspring("GO:0042776", GOBPOFFSPRING)),
  "Mitochondrial translation (GO:0032543)" = c("GO:0032543",
                                               get_offspring("GO:0032543", GOBPOFFSPRING))
)

set_test <- function(set_name, ids, dir_label) {
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ids),
                                 columns = "ENTREZID", keytype = "GO")$ENTREZID
  genes <- unique(genes[genes %in% universe_entrez])
  if (length(genes) < 5) return(NULL)
  is_deg <- deg_flag$entrez %in% (if (dir_label == "Up_in_AKI") {
    deg_flag$entrez[deg_flag$p_val_adj < 0.05 & deg_flag$avg_log2FC > 0.25]
  } else {
    deg_flag$entrez[deg_flag$p_val_adj < 0.05 & deg_flag$avg_log2FC < -0.25]
  })
  in_set <- deg_flag$entrez %in% genes
  a <- sum(is_deg & in_set); b <- sum(is_deg & !in_set)
  c <- sum(!is_deg & in_set); d <- sum(!is_deg & !in_set)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.frame(gene_set = set_name, direction = dir_label,
             n_set = length(genes), n_degs = sum(is_deg), n_overlap = a,
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}

set_res <- bind_rows(lapply(names(mito_sets), function(nm) {
  bind_rows(lapply(c("Up_in_AKI", "Down_in_AKI"), function(d) {
    set_test(nm, mito_sets[[nm]], d)
  }))
})) %>% arrange(p_value)
write.csv(set_res, file.path(tab_dir, "PTEC_mito_gene_set_tests.csv"),
          row.names = FALSE)

# ---------- 6. 图 1: 更新后的 GO 双向点图 ----------
plot_go_dot <- function(tab, n_up, n_down, top_n = 10) {
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  df <- tab %>%
    mutate(GeneRatioNum = vapply(strsplit(GeneRatio, "/"),
                                 function(x) as.numeric(x[1]) / as.numeric(x[2]),
                                 numeric(1)),
           Cluster = factor(ifelse(direction == "Up_in_AKI", "up", "down"),
                            levels = c("down", "up")))
  top_down <- df %>% filter(Cluster == "down") %>% arrange(p.adjust) %>% head(top_n)
  top_up   <- df %>% filter(Cluster == "up")   %>% arrange(p.adjust) %>% head(top_n)
  ord <- bind_rows(
    top_down %>% mutate(src = 1),
    top_up %>% filter(!Description %in% top_down$Description) %>% mutate(src = 2)
  ) %>% arrange(src, p.adjust)
  term_levels <- rev(unique(ord$Description))
  df_plot <- df %>%
    filter(Description %in% term_levels) %>%
    mutate(Description = factor(Description, levels = term_levels),
           TermLabel = factor(str_wrap(as.character(Description), width = 48),
                              levels = str_wrap(term_levels, width = 48)))
  ggplot(df_plot, aes(x = Cluster, y = TermLabel)) +
    geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust))) +
    scale_x_discrete(labels = c(down = paste0("down\n(", n_down, ")"),
                                up   = paste0("up\n(", n_up, ")"))) +
    scale_color_gradientn(colours = c("#2C7BB6", "#8073A8", "#D6604D"),
                          name = expression(-log[10]~"(p.adjust)")) +
    scale_size_continuous(range = c(1, 7), name = "GeneRatio") +
    labs(x = "Cluster", y = NULL,
         title = "Enriched GO (BP) of PTEC DEGs (updated annotation)") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_line(colour = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

p1 <- plot_go_dot(go_all,
                  n_up = sum(deg$direction == "Up in AKI"),
                  n_down = sum(deg$direction == "Down in AKI"),
                  top_n = 10)
if (!is.null(p1)) {
  ggsave(file.path(fig_dir, "PTEC_GO_BP_dotplot_updated.pdf"), p1,
         width = 9, height = 8)
}

# ---------- 7. 图 2: 线粒体相关 GO 条目条形图 ----------
mito_sig <- mito_tests %>%
  filter(p_adjust < 0.05) %>%
  mutate(direction = factor(direction, levels = c("Up_in_AKI", "Down_in_AKI"),
                            labels = c("Up in AKI", "Down in AKI")),
         signed_logp = ifelse(direction == "Up in AKI", -log10(p_adjust),
                              log10(p_adjust)))
if (nrow(mito_sig) > 0) {
  top_mito <- mito_sig %>%
    arrange(p_adjust) %>%
    head(20) %>%
    # 同一词条可能在上下调方向均显著, 标签中加入方向避免条形重叠
    mutate(label = paste0(str_wrap(Description, 55), " [",
                          ifelse(direction == "Up in AKI", "up", "down"), "]"))
  p2 <- ggplot(top_mito, aes(x = signed_logp, y = reorder(label, signed_logp),
                             fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = c("Up in AKI" = "#CE6C8E",
                                 "Down in AKI" = "#5E8DBC"), name = NULL) +
    labs(x = expression(signed~-log[10]~"(p.adjust)"),
         y = NULL, title = "Mitochondria-related GO (BP) in PTEC") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9, colour = "black"))
  ggsave(file.path(fig_dir, "PTEC_GO_BP_mitochondrial_barplot.pdf"), p2,
         width = 9.5, height = 6.5)
} else {
  message("线粒体相关 GO 条目无显著结果, 跳过条形图")
}

# ---------- 8. 图 3: 线粒体相关差异基因表达点图 ----------
mito_gene_ids <- unique(c("GO:0005739", get_offspring("GO:0005739", GOCCOFFSPRING)))
mito_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = mito_gene_ids,
                                      columns = "SYMBOL", keytype = "GOALL")$SYMBOL
mito_deg <- deg %>%
  filter(gene %in% mito_symbols, direction != "NS", !grepl("^MT-", gene))

top_up_mito   <- mito_deg %>% filter(direction == "Up in AKI") %>%
  arrange(p_val_adj) %>% head(10)
top_down_mito <- mito_deg %>% filter(direction == "Down in AKI") %>%
  arrange(p_val_adj) %>% head(10)
mito_show <- c(top_up_mito$gene, top_down_mito$gene)
if (length(mito_show) >= 2) {
  p3 <- DotPlot(ptec, features = mito_show, group.by = "group",
                cols = c("CTRL" = "#7FC4C9", "AKI" = "#E17B6E")) +
    RotatedAxis() +
    labs(x = NULL, y = "Condition",
         title = "Mitochondria-related DEGs in PTEC (AKI vs CTRL)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text.x = element_text(face = "italic", size = 8))
  ggsave(file.path(fig_dir, "PTEC_mito_DEG_dotplot.pdf"), p3,
         width = 11, height = 3.6)
}

# ---------- 9. 运行汇总与可重复性记录 ----------
sink(file.path(proj_dir, "logs", "sessionInfo_09_ptec_mito_reanalysis.txt"))
cat("=== 09_ptec_mito_reanalysis.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 注释版本 ===\n"); print(pkg_ver)
cat("\n=== 旧结果中废弃 GO 词条 ===\n"); print(obsolete_old)
cat("\n=== PTEC 差异基因 ===\n")
cat("Up in AKI:", sum(deg$direction == "Up in AKI"),
    "; Down in AKI:", sum(deg$direction == "Down in AKI"), "\n")
cat("\n=== 更新后 GO 富集 ===\n")
cat("条目数:", nrow(go_all), "; 其中线粒体相关:",
    sum(go_all$is_mito), "\n")
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 09_ptec_mito_reanalysis.R 完成")
