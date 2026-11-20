# ============================================================
# rerun_20261119 / 06_nets_scoring.R
# 功能：以 NETs 基因集 (PADI4, ELANE, MPO, CYBB, S100A9, HMGB1, GSDMD)
#       对中性粒细胞 (Neutrophil) 在 AKI 与 CTRL 组之间做打分评价与可视化。
#
# 评价层次:
#   1) 单细胞水平: AddModuleScore 模块评分 + Wilcoxon 检验 (CTRL 35 vs AKI 344)
#   2) 样本水平  : 每样本平均评分 + Wilcoxon 检验 (CTRL 15 vs AKI 7)
#                  (单细胞水平会把同一供体的细胞当作独立样本, 存在伪重复风险,
#                   故同时给出样本水平的稳健性检验)
#   3) 敏感性分析: 去掉中性粒细胞中最丰度的 S100A9 后重新打分, 检验结论是否依赖单基因
#   4) 单基因层面: 7 个基因各自的表达比例与差异倍数
#
# 输入：rerun_20261119/results/seu_selected.rds
#       rerun_20261119/results/tables/DEG_Neutrophil_AKI_vs_CTRL.csv
# 输出：results/figures/NETs_score_violin.pdf
#       results/figures/NETs_panel_dotplot.pdf
#       results/figures/NETs_score_sample_level.pdf
#       results/figures/NETs_panel_pseudobulk_heatmap.pdf
#       results/tables/NETs_score_per_cell.csv
#       results/tables/NETs_score_sample_summary.csv
#       results/tables/NETs_score_tests.csv
#       results/tables/NETs_panel_gene_stats.csv
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

proj_dir  <- "/home/zhouxin/liudi_AKI_project"
rerun_dir <- file.path(proj_dir, "rerun_20261119")
fig_dir   <- file.path(rerun_dir, "results", "figures")
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# 与既有图一致的配色
col_up   <- "#CE6C8E"   # Up in AKI
col_down <- "#5E8DBC"   # Down in AKI
col_aki  <- "#E17B6E"   # AKI
col_ctrl <- "#7FC4C9"   # CTRL

# ---------- 1. 用户指定的 NETs 基因集 ----------
NETs_panel <- c("PADI4", "ELANE", "MPO", "CYBB", "S100A9", "HMGB1", "GSDMD")

# ---------- 2. 加载数据并提取中性粒细胞 ----------
seu <- readRDS(file.path(rerun_dir, "results", "seu_selected.rds"))
DefaultAssay(seu) <- "RNA"
neut <- subset(seu, subset = cell_group == "Neutrophil")
neut$group <- factor(neut$group, levels = c("CTRL", "AKI"))
Idents(neut) <- "group"
neut$specimen <- as.character(neut$specimen)

cat("中性粒细胞数:", ncol(neut),
    "; CTRL:", sum(neut$group == "CTRL"),
    "; AKI:", sum(neut$group == "AKI"), "\n")
cat("样本数: CTRL =", length(unique(neut$specimen[neut$group == "CTRL"])),
    "; AKI =", length(unique(neut$specimen[neut$group == "AKI"])), "\n")

present <- intersect(NETs_panel, rownames(neut))
missing <- setdiff(NETs_panel, present)
if (length(missing) > 0) cat("未检出的基因:", paste(missing, collapse = ", "), "\n")
cat("参与打分的基因:", paste(present, collapse = ", "), "\n")

# ---------- 3. 模块评分(全 7 基因 + 去掉 S100A9 的敏感性版本) ----------
# AddModuleScore 会对基因按表达量分箱, 用表达水平相近的对照基因做背景校正,
# 因此可直接比较组间评分; 列名会自动追加序号 "1"。
neut <- AddModuleScore(neut, features = list(present), name = "NETs")
neut <- AddModuleScore(neut, features = list(setdiff(present, "S100A9")),
                       name = "NETsnoS100")
col_main <- grep("^NETs[0-9]+$", colnames(neut@meta.data), value = TRUE)[1]
col_wo   <- grep("^NETsnoS100", colnames(neut@meta.data), value = TRUE)[1]

# 细胞层面按基因 z-score 后再平均的评分(对高表达基因不敏感, 作为第三种稳健性指标)
data_mat <- tryCatch(
  GetAssayData(neut, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(neut, assay = "RNA", slot = "data"))
z_mat <- t(scale(t(as.matrix(data_mat[present, , drop = FALSE]))))
z_mat[is.na(z_mat)] <- 0
neut$NETs_zmean <- colMeans(z_mat)

meta <- neut@meta.data
meta$cell <- rownames(meta)

# ---------- 4. 统计检验 ----------
# 4.1 单细胞水平 Wilcoxon
cell_p <- wilcox.test(meta[[col_main]] ~ meta$group, exact = FALSE)
cell_p_wo <- wilcox.test(meta[[col_wo]] ~ meta$group, exact = FALSE)
cell_p_z <- wilcox.test(meta$NETs_zmean ~ meta$group, exact = FALSE)

# 4.2 样本水平: 每样本平均评分, 避免把同一供体的细胞当独立样本
sample_summary <- meta %>%
  group_by(specimen, group) %>%
  summarise(n_cells = n(),
            mean_score = mean(.data[[col_main]]),
            median_score = median(.data[[col_main]]),
            mean_score_noS100A9 = mean(.data[[col_wo]]),
            .groups = "drop") %>%
  arrange(group, specimen)
write.csv(sample_summary,
          file.path(tab_dir, "NETs_score_sample_summary.csv"),
          row.names = FALSE)

samp_p <- wilcox.test(mean_score ~ group, data = sample_summary, exact = FALSE)
samp_p_wo <- wilcox.test(mean_score_noS100A9 ~ group, data = sample_summary,
                         exact = FALSE)

fmt_p <- function(p) {
  ifelse(p < 1e-4,
         paste0("P = ", formatC(p, format = "e", digits = 2)),
         paste0("P = ", formatC(p, format = "f", digits = 4)))
}

tests <- data.frame(
  comparison = c("Cell level: NETs score (7 genes)",
                 "Cell level: NETs score (6 genes, -S100A9)",
                 "Cell level: mean gene z-score (7 genes)",
                 "Sample level: mean NETs score (7 genes)",
                 "Sample level: mean NETs score (6 genes, -S100A9)"),
  unit = c("cell", "cell", "cell", "sample", "sample"),
  n_CTRL = c(sum(meta$group == "CTRL"), sum(meta$group == "CTRL"),
             sum(meta$group == "CTRL"), sum(sample_summary$group == "CTRL"),
             sum(sample_summary$group == "CTRL")),
  n_AKI  = c(sum(meta$group == "AKI"), sum(meta$group == "AKI"),
             sum(meta$group == "AKI"), sum(sample_summary$group == "AKI"),
             sum(sample_summary$group == "AKI")),
  mean_CTRL = c(mean(meta[[col_main]][meta$group == "CTRL"]),
                mean(meta[[col_wo]][meta$group == "CTRL"]),
                mean(meta$NETs_zmean[meta$group == "CTRL"]),
                mean(sample_summary$mean_score[sample_summary$group == "CTRL"]),
                mean(sample_summary$mean_score_noS100A9[sample_summary$group == "CTRL"])),
  mean_AKI = c(mean(meta[[col_main]][meta$group == "AKI"]),
               mean(meta[[col_wo]][meta$group == "AKI"]),
               mean(meta$NETs_zmean[meta$group == "AKI"]),
               mean(sample_summary$mean_score[sample_summary$group == "AKI"]),
               mean(sample_summary$mean_score_noS100A9[sample_summary$group == "AKI"])),
  p_value = c(cell_p$p.value, cell_p_wo$p.value, cell_p_z$p.value,
              samp_p$p.value, samp_p_wo$p.value),
  stringsAsFactors = FALSE
)
tests$direction <- ifelse(tests$mean_AKI > tests$mean_CTRL, "Higher in AKI",
                          "Higher in CTRL")
write.csv(tests, file.path(tab_dir, "NETs_score_tests.csv"), row.names = FALSE)
print(tests)

# 4.3 逐细胞评分明细
per_cell <- meta %>%
  transmute(cell, specimen, group,
            NETs_score = .data[[col_main]],
            NETs_score_noS100A9 = .data[[col_wo]],
            NETs_zmean = NETs_zmean)
write.csv(per_cell, file.path(tab_dir, "NETs_score_per_cell.csv"), row.names = FALSE)

# ---------- 5. 单基因层面统计 ----------
deg <- read.csv(file.path(tab_dir, "DEG_Neutrophil_AKI_vs_CTRL.csv"),
                stringsAsFactors = FALSE)
gene_stats <- data.frame(gene = NETs_panel, in_data = NETs_panel %in% present)
gene_stats <- gene_stats %>%
  left_join(deg %>% select(gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2),
            by = "gene")
# pct.1 = AKI 组检出比例, pct.2 = CTRL 组检出比例(FindMarkers ident.1=AKI)
gene_stats <- gene_stats %>%
  rename(pct_AKI = pct.1, pct_CTRL = pct.2) %>%
  mutate(significant = ifelse(!is.na(p_val_adj) & p_val_adj < 0.05, "yes", "no"))
write.csv(gene_stats, file.path(tab_dir, "NETs_panel_gene_stats.csv"),
          row.names = FALSE)
print(gene_stats)

# ============================================================
# 图 1: NETs 模块评分小提琴图(含去 S100A9 的敏感性面板)
# ============================================================
sc <- meta %>%
  select(group, all_of(c(col_main, col_wo))) %>%
  rename(`NETs score (7 genes)` = all_of(col_main),
         `NETs score (6 genes, -S100A9)` = all_of(col_wo)) %>%
  pivot_longer(cols = -group, names_to = "panel", values_to = "score") %>%
  mutate(panel = factor(panel, levels = c("NETs score (7 genes)",
                                          "NETs score (6 genes, -S100A9)")))

# 每个分面各自的范围用于放置显著性括号
y_stats <- sc %>%
  group_by(panel) %>%
  summarise(y_max = max(score, na.rm = TRUE),
            y_min = min(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_rng = pmax(y_max - y_min, 0.05),
         lab = c(fmt_p(cell_p$p.value), fmt_p(cell_p_wo$p.value)),
         y_seg = y_max + y_rng * 0.12,
         y_txt = y_max + y_rng * 0.20)

p1 <- ggplot(sc, aes(x = group, y = score, fill = group)) +
  geom_violin(scale = "width", width = 0.85, colour = "grey35", linewidth = 0.5,
              alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "grey25",
               linewidth = 0.5, fill = "white", alpha = 0.9) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  geom_segment(data = y_stats,
               aes(x = 1, xend = 2, y = y_seg, yend = y_seg),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(data = y_stats,
               aes(x = 1, xend = 1, y = y_seg - y_rng * 0.05, yend = y_seg),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(data = y_stats,
               aes(x = 2, xend = 2, y = y_seg - y_rng * 0.05, yend = y_seg),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = y_stats, aes(x = 1.5, y = y_txt, label = lab),
            size = 3.6, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  labs(x = NULL, y = "Module score (mean expression)",
       title = "NETs programme score in neutrophils (AKI vs CTRL)",
       subtitle = paste0("AddModuleScore on PADI4/ELANE/MPO/CYBB/S100A9/HMGB1/GSDMD; ",
                         "right panel excludes the highly abundant S100A9")) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey30"),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.text = element_text(colour = "black"))
ggsave(file.path(fig_dir, "NETs_score_violin.pdf"), p1,
       width = 8.8, height = 5.0)

# ============================================================
# 图 2: 7 个 NETs 基因的表达比例与表达强度点图
# ============================================================
p2 <- DotPlot(neut, features = rev(present), group.by = "group", dot.scale = 9) +
  scale_colour_gradient2(low = col_down, mid = "white", high = col_up,
                         midpoint = 0) +
  coord_flip() +
  labs(x = NULL, y = NULL,
       colour = "Scaled mean\nexpression",
       title = "NETs panel genes in neutrophils (AKI vs CTRL)",
       subtitle = "Dot size = fraction of expressing cells; colour = scaled mean expression") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 8, 8))
ggsave(file.path(fig_dir, "NETs_panel_dotplot.pdf"), p2,
       width = 5.8, height = 5.2)

# ============================================================
# 图 3: 样本水平评分(每样本一个点, 避免伪重复)
# ============================================================
p3 <- ggplot(sample_summary, aes(x = group, y = mean_score, fill = group)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.45,
               colour = "grey35", linewidth = 0.5) +
  geom_point(aes(size = n_cells), shape = 21, colour = "grey20",
             position = position_jitter(width = 0.12, height = 0, seed = 2026)) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  scale_size_continuous(range = c(1.6, 4.2), guide = "none") +
  annotate("segment", x = 1, xend = 2,
           y = max(sample_summary$mean_score) * 1.10,
           yend = max(sample_summary$mean_score) * 1.10,
           colour = "grey25", linewidth = 0.5) +
  annotate("text", x = 1.5, y = max(sample_summary$mean_score) * 1.16,
           label = fmt_p(samp_p$p.value), size = 3.6) +
  annotate("text", x = 1, y = min(sample_summary$mean_score) - 0.02,
           label = paste0("n = ", sum(sample_summary$group == "CTRL")),
           size = 3.2, colour = "grey35") +
  annotate("text", x = 2, y = min(sample_summary$mean_score) - 0.02,
           label = paste0("n = ", sum(sample_summary$group == "AKI")),
           size = 3.2, colour = "grey35") +
  expand_limits(y = c(min(sample_summary$mean_score) - 0.04, NA)) +
  labs(x = NULL, y = "Mean NETs score per sample",
       title = "Sample-level NETs score in neutrophils",
       subtitle = "Each point is one specimen; symbol size reflects the number of neutrophils") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETs_score_sample_level.pdf"), p3,
       width = 5.2, height = 4.8)

# ============================================================
# 图 4: 各样本 NETs 基因表达的伪批量热图
# ============================================================
pb <- sapply(sample_summary$specimen, function(s) {
  cells <- rownames(meta)[meta$specimen == s]
  Matrix::rowMeans(data_mat[present, cells, drop = FALSE])
})
pb <- pb[, sample_summary$specimen, drop = FALSE]
z <- t(scale(t(pb)))
z[is.na(z)] <- 0
z <- pmax(pmin(z, 2), -2)
colnames(z) <- sample_summary$specimen

ha <- HeatmapAnnotation(
  Group = sample_summary$group,
  col = list(Group = c(CTRL = col_ctrl, AKI = col_aki)),
  annotation_name_side = "left",
  annotation_legend_param = list(Group = list(title = "Group"))
)

ht <- Heatmap(
  z,
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c(col_down, "white", col_up)),
  row_labels = rownames(z),
  row_names_gp = gpar(fontface = "italic", fontsize = 9),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(sample_summary$group, levels = c("CTRL", "AKI")),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 7),
  top_annotation = ha,
  column_gap = unit(1.2, "mm"),
  border = TRUE,
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2),
                              title = "Scaled mean\nexpression (z)")
)
pdf(file.path(fig_dir, "NETs_panel_pseudobulk_heatmap.pdf"),
    width = 7.6, height = 4.4)
draw(ht,
     column_title = "NETs panel genes across neutrophil samples",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     column_title_side = "top",
     padding = unit(c(12, 2, 2, 2), "mm"))
invisible(dev.off())

# ---------- 6. 运行记录 ----------
sink(file.path(log_dir, "sessionInfo_06_nets_scoring.txt"))
cat("=== 06_nets_scoring.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 中性粒细胞构成 ===\n")
print(table(neut$group))
cat("\n=== 样本构成 ===\n")
print(table(sample_summary$group))
cat("\n=== 检验结果 ===\n")
print(tests)
cat("\n=== 单基因统计 ===\n")
print(gene_stats)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 06_nets_scoring.R 完成")
