# ============================================================
# rerun_20261119_noSample / 06_netsig_scoring.R
# 功能：以用户提供的 38 基因 NET 基因集, 对中性粒细胞(Neutrophil)
#       在 AKI 与 CTRL 组之间做打分评价与可视化。
#       本目录已剔除 specimen 以 "Sample" 开头的样本(13 样本: AKI 7 / CTRL 6)。
#
# 基因集(38):
#   AKT1 AKT2 ATG7 CLEC6A CSF3 CTSG CYBB DNASE1 ELANE ENTPD4 F3 HMGB1
#   IL17A IL1B IL6 IL8 IRAK4 ITGAM ITGB2 KCNN3 MAPK1 MAPK3 MMP9 MPO MTOR
#   PADI4 PTAFR PIK3CA RIPK1 RIPK3 SELP SELPLG SIGLEC14 TLR2 TLR4 TLR7 TLR8 TNF
#
# 设计:
#   1) 总体评分: 38 基因 AddModuleScore + 细胞水平/样本水平 Wilcoxon
#   2) 子模块评分: 按功能拆成 4 组(颗粒与 NET 效应分子、细胞因子、
#      受体/黏附、上游信号), 避免广泛表达的信号通路基因稀释结论
#   3) 单基因层面: 38 个基因的表达比例与差异倍数(来自既有 DEG 结果)
#
# 输入：rerun_20261119_noSample/results/seu_selected.rds
#       rerun_20261119_noSample/results/tables/DEG_Neutrophil_AKI_vs_CTRL.csv
# 输出：results/figures/NETsig_total_score_violin.pdf
#       results/figures/NETsig_module_scores_boxplot.pdf
#       results/figures/NETsig_panel_dotplot.pdf
#       results/figures/NETsig_score_sample_level.pdf
#       results/figures/NETsig_pseudobulk_heatmap.pdf
#       results/tables/NETsig_score_per_cell.csv
#       results/tables/NETsig_score_sample_summary.csv
#       results/tables/NETsig_score_tests.csv
#       results/tables/NETsig_gene_stats.csv
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
rerun_dir <- file.path(proj_dir, "rerun_20261119_noSample")
fig_dir   <- file.path(rerun_dir, "results", "figures")
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")

# 与既有图一致的配色
col_up   <- "#CE6C8E"
col_down <- "#5E8DBC"
col_aki  <- "#E17B6E"
col_ctrl <- "#7FC4C9"

# ---------- 1. 基因集 ----------
NET_genes <- c("AKT1","AKT2","ATG7","CLEC6A","CSF3","CTSG","CYBB","DNASE1",
               "ELANE","ENTPD4","F3","HMGB1","IL17A","IL1B","IL6","IL8",
               "IRAK4","ITGAM","ITGB2","KCNN3","MAPK1","MAPK3","MMP9","MPO",
               "MTOR","PADI4","PTAFR","PIK3CA","RIPK1","RIPK3","SELP",
               "SELPLG","SIGLEC14","TLR2","TLR4","TLR7","TLR8","TNF")
stopifnot(length(NET_genes) == 38)

# 按功能拆分为 4 个子模块(合计 38 个基因)
net_modules <- list(
  "Granule / NET effectors" =
    c("ELANE","MPO","CTSG","MMP9","PADI4","CYBB","DNASE1"),
  "Cytokines / chemokines" =
    c("IL1B","IL6","IL8","TNF","CSF3","IL17A","HMGB1"),
  "Receptors / adhesion" =
    c("TLR2","TLR4","TLR7","TLR8","CLEC6A","SIGLEC14","ITGAM","ITGB2",
      "SELP","SELPLG","PTAFR","F3","ENTPD4","KCNN3"),
  "Upstream signalling" =
    c("AKT1","AKT2","MAPK1","MAPK3","MTOR","PIK3CA","IRAK4","RIPK1",
      "RIPK3","ATG7")
)
stopifnot(length(unlist(net_modules)) == 38)

# ---------- 2. 数据 ----------
seu <- readRDS(file.path(rerun_dir, "results", "seu_selected.rds"))
DefaultAssay(seu) <- "RNA"
neut <- subset(seu, subset = cell_group == "Neutrophil")
neut$group <- factor(neut$group, levels = c("CTRL", "AKI"))
Idents(neut) <- "group"
neut$specimen <- as.character(neut$specimen)

cat("中性粒细胞数:", ncol(neut),
    "; CTRL:", sum(neut$group == "CTRL"),
    "; AKI:", sum(neut$group == "AKI"), "\n")

present <- intersect(NET_genes, rownames(neut))
cat("参与打分的基因:", length(present), "/", length(NET_genes), "\n")
if (length(setdiff(NET_genes, present)) > 0) {
  cat("未检出:", paste(setdiff(NET_genes, present), collapse = ", "), "\n")
}
net_modules <- lapply(net_modules, function(x) intersect(x, present))

# ---------- 3. 打分 ----------
# 总体评分
neut <- AddModuleScore(neut, features = list(present), name = "NETsig")
col_total <- grep("^NETsig[0-9]+$", colnames(neut@meta.data), value = TRUE)[1]

# 子模块评分(一次传入 4 个基因集, 依次生成 NETmod1..NETmod4)
neut <- AddModuleScore(neut, features = net_modules, name = "NETmod")
mod_cols <- grep("^NETmod[0-9]+$", colnames(neut@meta.data), value = TRUE)
mod_names <- names(net_modules)
cat("子模块评分列:", paste(mod_cols, collapse = ", "), "\n")

meta <- neut@meta.data
meta$cell <- rownames(meta)

# ---------- 4. 统计检验 ----------
fmt_p <- function(p) {
  ifelse(p < 1e-4,
         paste0("P = ", formatC(p, format = "e", digits = 2)),
         paste0("P = ", formatC(p, format = "f", digits = 3)))
}

# 4.1 细胞水平(注意: 同一供体细胞存在伪重复, 仅作描述)
cell_test <- function(x) wilcox.test(x ~ meta$group, exact = FALSE)

# 4.2 样本水平(每样本均值, 更保守)
sample_summary <- meta %>%
  group_by(specimen, group) %>%
  summarise(n_cells = n(),
            mean_total = mean(.data[[col_total]]),
            across(all_of(mod_cols), mean, .names = "mean_{.col}"),
            .groups = "drop") %>%
  arrange(group, specimen)
write.csv(sample_summary,
          file.path(tab_dir, "NETsig_score_sample_summary.csv"),
          row.names = FALSE)

tests <- data.frame(
  score = c("Total NET score (38 genes)", mod_names),
  n_genes = c(length(present), sapply(net_modules, length)),
  mean_CTRL_cell = NA_real_, mean_AKI_cell = NA_real_,
  p_cell = NA_real_,
  mean_CTRL_sample = NA_real_, mean_AKI_sample = NA_real_,
  p_sample = NA_real_,
  stringsAsFactors = FALSE
)
score_cols <- c(col_total, mod_cols)
sample_cols <- c("mean_total", paste0("mean_", mod_cols))
for (i in seq_along(score_cols)) {
  ct <- cell_test(meta[[score_cols[i]]])
  st <- wilcox.test(sample_summary[[sample_cols[i]]] ~ sample_summary$group,
                    exact = FALSE)
  tests$mean_CTRL_cell[i] <- mean(meta[[score_cols[i]]][meta$group == "CTRL"])
  tests$mean_AKI_cell[i]  <- mean(meta[[score_cols[i]]][meta$group == "AKI"])
  tests$p_cell[i]         <- ct$p.value
  tests$mean_CTRL_sample[i] <- mean(sample_summary[[sample_cols[i]]][sample_summary$group == "CTRL"])
  tests$mean_AKI_sample[i]  <- mean(sample_summary[[sample_cols[i]]][sample_summary$group == "AKI"])
  tests$p_sample[i]         <- st$p.value
}
tests$direction_cell <- ifelse(tests$mean_AKI_cell > tests$mean_CTRL_cell,
                               "Higher in AKI", "Higher in CTRL")
tests$direction_sample <- ifelse(tests$mean_AKI_sample > tests$mean_CTRL_sample,
                                 "Higher in AKI", "Higher in CTRL")
write.csv(tests, file.path(tab_dir, "NETsig_score_tests.csv"), row.names = FALSE)
print(tests, digits = 3)

# 4.3 每细胞评分
per_cell <- meta %>%
  transmute(cell, specimen, group,
            NET_score_total = .data[[col_total]],
            across(all_of(mod_cols), .names = "score_{.col}"))
write.csv(per_cell, file.path(tab_dir, "NETsig_score_per_cell.csv"),
          row.names = FALSE)

# 4.4 单基因统计
deg <- read.csv(file.path(tab_dir, "DEG_Neutrophil_AKI_vs_CTRL.csv"),
                stringsAsFactors = FALSE)
gene_stats <- data.frame(gene = NET_genes, stringsAsFactors = FALSE) %>%
  left_join(data.frame(gene = unlist(net_modules),
                       module = rep(names(net_modules),
                                    sapply(net_modules, length))),
            by = "gene") %>%
  left_join(deg %>% select(gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2),
            by = "gene") %>%
  rename(pct_AKI = pct.1, pct_CTRL = pct.2) %>%
  mutate(in_data = gene %in% present,
         significant = ifelse(!is.na(p_val_adj) & p_val_adj < 0.05,
                              "yes", "no"))
write.csv(gene_stats, file.path(tab_dir, "NETsig_gene_stats.csv"),
          row.names = FALSE)
print(gene_stats %>% select(gene, module, avg_log2FC, p_val_adj, pct_AKI, pct_CTRL))

# ============================================================
# 图 1: 38 基因总体 NET 评分小提琴图
# ============================================================
p1 <- ggplot(meta, aes(x = group, y = .data[[col_total]], fill = group)) +
  geom_violin(scale = "width", width = 0.85, colour = "grey35", linewidth = 0.5,
              alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "grey25",
               linewidth = 0.5, fill = "white", alpha = 0.9) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  annotate("segment", x = 1, xend = 2,
           y = max(meta[[col_total]]) * 1.06,
           yend = max(meta[[col_total]]) * 1.06,
           colour = "grey25", linewidth = 0.5) +
  annotate("text", x = 1.5, y = max(meta[[col_total]]) * 1.19,
           label = paste0("Cell-level ", fmt_p(tests$p_cell[1]),
                          "\nSample-level ", fmt_p(tests$p_sample[1])),
           size = 3.4, lineheight = 0.95) +
  expand_limits(y = c(min(meta[[col_total]]) * 1.05,
                      max(meta[[col_total]]) * 1.33)) +
  labs(x = NULL, y = "NET module score (mean expression)",
       title = "NET programme score in neutrophils (38-gene panel)",
       subtitle = "AKI vs CTRL; cell-level test is descriptive, sample-level test accounts for donor clustering") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_total_score_violin.pdf"), p1,
       width = 5.6, height = 4.8)

# ============================================================
# 图 2: 4 个功能子模块评分
# ============================================================
sc <- meta %>%
  select(group, all_of(mod_cols)) %>%
  setNames(c("group", mod_names)) %>%
  pivot_longer(cols = -group, names_to = "module", values_to = "score") %>%
  mutate(module = factor(module, levels = mod_names))

mod_lab <- tests %>%
  filter(score != "Total NET score (38 genes)") %>%
  mutate(module = factor(score, levels = mod_names),
         lab = paste0("Cell ", fmt_p(p_cell), "\nSample ", fmt_p(p_sample)))

y_stats <- sc %>%
  group_by(module) %>%
  summarise(y_max = max(score, na.rm = TRUE),
            y_min = min(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_rng = pmax(y_max - y_min, 0.04)) %>%
  left_join(mod_lab %>% select(module, lab), by = "module") %>%
  mutate(y_seg = y_max + y_rng * 0.10,
         y_txt = y_max + y_rng * 0.24)

p2 <- ggplot(sc, aes(x = group, y = score, fill = group)) +
  geom_violin(scale = "width", width = 0.85, colour = "grey35", linewidth = 0.5,
              alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "grey25",
               linewidth = 0.5, fill = "white", alpha = 0.9) +
  facet_wrap(~module, scales = "free_y", nrow = 1) +
  geom_segment(data = y_stats,
               aes(x = 1, xend = 2, y = y_seg, yend = y_seg),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = y_stats, aes(x = 1.5, y = y_txt, label = lab),
            size = 2.9, lineheight = 0.95, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  labs(x = NULL, y = "Module score (mean expression)",
       title = "Functional sub-module scores of the 38-gene NET panel",
       subtitle = "Cell-level and sample-level Wilcoxon tests are both shown for each sub-module") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9.5),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_module_scores_boxplot.pdf"), p2,
       width = 10.5, height = 5.0)

# ============================================================
# 图 3: 38 个基因的表达比例与表达强度
# ============================================================
p3 <- DotPlot(neut, features = rev(present), group.by = "group",
              dot.scale = 7) +
  scale_colour_gradient2(low = col_down, mid = "white", high = col_up,
                         midpoint = 0) +
  coord_flip() +
  labs(x = NULL, y = NULL, colour = "Scaled mean\nexpression",
       title = "38-gene NET panel in neutrophils (AKI vs CTRL)",
       subtitle = "Dot size = fraction of expressing cells; colour = scaled mean expression") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 8.5),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 8, 8))
ggsave(file.path(fig_dir, "NETsig_panel_dotplot.pdf"), p3,
       width = 5.8, height = 9.2)

# ============================================================
# 图 4: 样本水平总体评分
# ============================================================
p4 <- ggplot(sample_summary, aes(x = group, y = mean_total, fill = group)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.45,
               colour = "grey35", linewidth = 0.5) +
  geom_point(aes(size = n_cells), shape = 21, colour = "grey20",
             position = position_jitter(width = 0.12, height = 0, seed = 2026)) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  scale_size_continuous(range = c(1.6, 4.2), guide = "none") +
  annotate("segment", x = 1, xend = 2,
           y = max(sample_summary$mean_total) * 1.12,
           yend = max(sample_summary$mean_total) * 1.12,
           colour = "grey25", linewidth = 0.5) +
  annotate("text", x = 1.5, y = max(sample_summary$mean_total) * 1.18,
           label = fmt_p(tests$p_sample[1]), size = 3.5) +
  labs(x = NULL, y = "Mean NET score per sample",
       title = "Sample-level NET score (38-gene panel)",
       subtitle = "Each point is one specimen; symbol size reflects the number of neutrophils") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_score_sample_level.pdf"), p4,
       width = 5.2, height = 4.8)

# ============================================================
# 图 5: 各样本 38 基因伪批量热图
# ============================================================
data_mat <- tryCatch(
  GetAssayData(neut, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(neut, assay = "RNA", slot = "data"))

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
  row_names_gp = gpar(fontface = "italic", fontsize = 7.5),
  row_split = factor(gene_stats$module[match(rownames(z), gene_stats$gene)],
                     levels = names(net_modules)),
  row_title_gp = gpar(fontsize = 7), row_title_rot = 0,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(sample_summary$group, levels = c("CTRL", "AKI")),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 7),
  top_annotation = ha,
  column_gap = unit(1.2, "mm"),
  row_gap = unit(1.2, "mm"),
  border = TRUE,
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2),
                              title = "Scaled mean\nexpression (z)")
)
pdf(file.path(fig_dir, "NETsig_pseudobulk_heatmap.pdf"),
    width = 8.2, height = 7.6)
draw(ht,
     column_title = "38-gene NET panel across neutrophil samples",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     column_title_side = "top",
     padding = unit(c(12, 2, 2, 2), "mm"))
invisible(dev.off())

# ---------- 5. 运行记录 ----------
sink(file.path(log_dir, "sessionInfo_07_netsig_scoring.txt"))
cat("=== 07_netsig_scoring.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 中性粒细胞构成 ===\n")
print(table(neut$group))
cat("\n=== 子模块基因数 ===\n")
print(sapply(net_modules, length))
cat("\n=== 检验结果 ===\n")
print(tests, digits = 3)
cat("\n=== 单基因统计 ===\n")
print(gene_stats)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 07_netsig_scoring.R 完成")
