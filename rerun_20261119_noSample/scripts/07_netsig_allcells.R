# ============================================================
# rerun_20261119_noSample / 07_netsig_allcells.R
# 功能：以 38 基因 NET 基因集, 在全部细胞(不区分细胞类型)上
#       评价 AKI 与 CTRL 两组之间的整体评分差异。
#
# 关键混杂: AKI 组中性粒细胞占比显著高于 CTRL, 而中性粒细胞本身高表达
#           CYBB/ITGB2/HMGB1 等基因, 因此"整体评分升高"有可能只是细胞组成
#           差异所致, 而非同一细胞类型内的真实上调。故同时给出:
#   (1) 未校正的细胞水平检验
#   (2) 校正 cell_group 的线性模型检验 (score ~ group + cell_group)
#   (3) 样本水平检验 (每样本评分均值)
#   (4) 分细胞类型的诊断面板, 用于判断差异来源
#
# 输入：rerun_20261119_noSample/results/seu_selected.rds
# 输出：results/figures/NETsig_allcells_total_score_violin.pdf
#       results/figures/NETsig_allcells_module_scores.pdf
#       results/figures/NETsig_allcells_by_celltype.pdf
#       results/figures/NETsig_allcells_sample_level.pdf
#       results/figures/NETsig_allcells_dotplot.pdf
#       results/tables/NETsig_allcells_score_tests.csv
#       results/tables/NETsig_allcells_score_per_cell.csv
#       results/tables/NETsig_allcells_sample_summary.csv
#       results/tables/NETsig_allcells_celltype_composition.csv
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

proj_dir  <- "/home/zhouxin/liudi_AKI_project"
rerun_dir <- file.path(proj_dir, "rerun_20261119_noSample")
fig_dir   <- file.path(rerun_dir, "results", "figures")
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")

col_up   <- "#CE6C8E"
col_down <- "#5E8DBC"
col_aki  <- "#E17B6E"
col_ctrl <- "#7FC4C9"

# ---------- 1. 基因集与分组模块 ----------
NET_genes <- c("AKT1","AKT2","ATG7","CLEC6A","CSF3","CTSG","CYBB","DNASE1",
               "ELANE","ENTPD4","F3","HMGB1","IL17A","IL1B","IL6","IL8",
               "IRAK4","ITGAM","ITGB2","KCNN3","MAPK1","MAPK3","MMP9","MPO",
               "MTOR","PADI4","PTAFR","PIK3CA","RIPK1","RIPK3","SELP",
               "SELPLG","SIGLEC14","TLR2","TLR4","TLR7","TLR8","TNF")

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

# ---------- 2. 全部细胞 ----------
seu <- readRDS(file.path(rerun_dir, "results", "seu_selected.rds"))
DefaultAssay(seu) <- "RNA"
seu$group <- factor(seu$group, levels = c("CTRL", "AKI"))
seu$cell_group <- factor(seu$cell_group,
                         levels = c("Neutrophil", "PTEC", "GEC", "PT"))
seu$specimen <- as.character(seu$specimen)

cat("全部细胞数:", ncol(seu),
    "; CTRL:", sum(seu$group == "CTRL"),
    "; AKI:", sum(seu$group == "AKI"), "\n")
print(table(seu$cell_group, seu$group))

# 细胞类型组成(判断混杂来源)
composition <- as.data.frame(table(seu$specimen, seu$cell_group)) %>%
  rename(specimen = Var1, cell_group = Var2, n = Freq) %>%
  group_by(specimen) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  left_join(seu@meta.data %>% distinct(specimen, group), by = "specimen") %>%
  arrange(group, specimen, cell_group)
write.csv(as.data.frame(composition),
          file.path(tab_dir, "NETsig_allcells_celltype_composition.csv"),
          row.names = FALSE)
print(composition %>% group_by(group, cell_group) %>%
        summarise(mean_fraction = mean(fraction), .groups = "drop") %>%
        pivot_wider(names_from = cell_group, values_from = mean_fraction))

# ---------- 3. 打分 ----------
present <- intersect(NET_genes, rownames(seu))
net_modules <- lapply(net_modules, function(x) intersect(x, present))
cat("参与打分基因:", length(present), "/ 38\n")

seu <- AddModuleScore(seu, features = list(present), name = "NETsig")
col_total <- grep("^NETsig[0-9]+$", colnames(seu@meta.data), value = TRUE)[1]
seu <- AddModuleScore(seu, features = net_modules, name = "NETmod")
mod_cols <- grep("^NETmod[0-9]+$", colnames(seu@meta.data), value = TRUE)
mod_names <- names(net_modules)

meta <- seu@meta.data
meta$cell <- rownames(meta)

# ---------- 4. 统计 ----------
fmt_p <- function(p) {
  ifelse(p < 1e-4,
         paste0("P = ", formatC(p, format = "e", digits = 2)),
         paste0("P = ", formatC(p, format = "f", digits = 3)))
}

# 4.1 样本水平均值
sample_summary <- meta %>%
  group_by(specimen, group) %>%
  summarise(n_cells = n(),
            mean_total = mean(.data[[col_total]]),
            across(all_of(mod_cols), mean, .names = "mean_{.col}"),
            .groups = "drop") %>%
  arrange(group, specimen)
write.csv(sample_summary,
          file.path(tab_dir, "NETsig_allcells_sample_summary.csv"),
          row.names = FALSE)

# 4.2 三种检验: 未校正细胞水平 / 校正细胞类型 / 样本水平
score_cols <- c(col_total, mod_cols)
sample_cols <- c("mean_total", paste0("mean_", mod_cols))

tests <- data.frame(
  score = c("Total NET score (38 genes)", mod_names),
  n_genes = c(length(present), sapply(net_modules, length)),
  mean_CTRL = NA_real_, mean_AKI = NA_real_,
  p_cell_unadjusted = NA_real_,
  p_cell_celltype_adjusted = NA_real_,
  p_sample_level = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(score_cols)) {
  x <- meta[[score_cols[i]]]
  tests$mean_CTRL[i] <- mean(x[meta$group == "CTRL"])
  tests$mean_AKI[i]  <- mean(x[meta$group == "AKI"])
  tests$p_cell_unadjusted[i] <- wilcox.test(x ~ meta$group, exact = FALSE)$p.value
  # 以 cell_group 为协变量校正细胞组成
  fit <- lm(x ~ meta$group + meta$cell_group)
  tests$p_cell_celltype_adjusted[i] <-
    summary(fit)$coefficients["meta$groupAKI", 4]
  tests$p_sample_level[i] <- wilcox.test(
    sample_summary[[sample_cols[i]]] ~ sample_summary$group, exact = FALSE)$p.value
}
tests$direction <- ifelse(tests$mean_AKI > tests$mean_CTRL,
                          "Higher in AKI", "Higher in CTRL")
write.csv(tests, file.path(tab_dir, "NETsig_allcells_score_tests.csv"),
          row.names = FALSE)
print(tests, digits = 3)

# 4.3 每细胞评分
per_cell <- meta %>%
  transmute(cell, specimen, group, cell_group,
            NET_score_total = .data[[col_total]],
            across(all_of(mod_cols), .names = "score_{.col}"))
write.csv(per_cell, file.path(tab_dir, "NETsig_allcells_score_per_cell.csv"),
          row.names = FALSE)

# ============================================================
# 图 1: 整体评分小提琴图
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
           label = paste0("Unadjusted ", fmt_p(tests$p_cell_unadjusted[1]),
                          "\nCell-type adjusted ", fmt_p(tests$p_cell_celltype_adjusted[1]),
                          "\nSample level ", fmt_p(tests$p_sample_level[1])),
           size = 3.3, lineheight = 0.95) +
  expand_limits(y = c(min(meta[[col_total]]) * 1.05,
                      max(meta[[col_total]]) * 1.38)) +
  labs(x = NULL, y = "NET module score (mean expression)",
       title = "NET programme score across all cells (38-gene panel)",
       subtitle = "All cell types pooled; cell-type-adjusted test uses cell_group as a covariate") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_allcells_total_score_violin.pdf"), p1,
       width = 5.8, height = 5.0)

# ============================================================
# 图 2: 4 个功能子模块(全部细胞)
# ============================================================
sc <- meta %>%
  select(group, all_of(mod_cols)) %>%
  setNames(c("group", mod_names)) %>%
  pivot_longer(cols = -group, names_to = "module", values_to = "score") %>%
  mutate(module = factor(module, levels = mod_names))

mod_lab <- tests %>%
  filter(score != "Total NET score (38 genes)") %>%
  mutate(module = factor(score, levels = mod_names),
         lab = paste0("Unadjusted ", fmt_p(p_cell_unadjusted),
                      "\nCell-type adj. ", fmt_p(p_cell_celltype_adjusted)))

y_stats <- sc %>%
  group_by(module) %>%
  summarise(y_max = max(score, na.rm = TRUE),
            y_min = min(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_rng = pmax(y_max - y_min, 0.04)) %>%
  left_join(mod_lab %>% select(module, lab), by = "module") %>%
  mutate(y_seg = y_max + y_rng * 0.10,
         y_txt = y_max + y_rng * 0.26)

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
            size = 2.8, lineheight = 0.95, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  labs(x = NULL, y = "Module score (mean expression)",
       title = "NET sub-module scores across all cells",
       subtitle = "Unadjusted and cell-type-adjusted P values are shown for each sub-module") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9.5),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_allcells_module_scores.pdf"), p2,
       width = 10.5, height = 5.0)

# ============================================================
# 图 3: 分细胞类型诊断面板(说明整体差异来源)
# ============================================================
ct_p <- meta %>%
  group_by(cell_group) %>%
  summarise(p = wilcox.test(.data[[col_total]] ~ group, exact = FALSE)$p.value,
            n_CTRL = sum(group == "CTRL"), n_AKI = sum(group == "AKI"),
            .groups = "drop")

y_ct <- meta %>%
  group_by(cell_group) %>%
  summarise(y_max = max(.data[[col_total]], na.rm = TRUE),
            y_min = min(.data[[col_total]], na.rm = TRUE), .groups = "drop") %>%
  mutate(y_rng = pmax(y_max - y_min, 0.05)) %>%
  left_join(ct_p, by = "cell_group") %>%
  mutate(y_seg = y_max + y_rng * 0.10,
         y_txt = y_max + y_rng * 0.16,
         lab = paste0(fmt_p(p), "\n(n=", n_CTRL, " vs ", n_AKI, ")"))

p3 <- ggplot(meta, aes(x = group, y = .data[[col_total]], fill = group)) +
  geom_violin(scale = "width", width = 0.85, colour = "grey35", linewidth = 0.5,
              alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "grey25",
               linewidth = 0.5, fill = "white", alpha = 0.9) +
  facet_wrap(~cell_group, scales = "free_y", nrow = 1) +
  geom_segment(data = y_ct,
               aes(x = 1, xend = 2, y = y_seg, yend = y_seg),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = y_ct, aes(x = 1.5, y = y_txt, label = lab),
            size = 2.9, lineheight = 0.95, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  labs(x = NULL, y = "NET module score (mean expression)",
       title = "NET score by cell type (diagnostic for composition confounding)",
       subtitle = "If AKI enrichment disappears within cell types, the pooled difference is driven by cell composition") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_allcells_by_celltype.pdf"), p3,
       width = 10.5, height = 5.0)

# ============================================================
# 图 4: 样本水平整体评分
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
  annotate("text", x = 1.5, y = max(sample_summary$mean_total) * 1.19,
           label = fmt_p(tests$p_sample_level[1]), size = 3.5) +
  expand_limits(y = max(sample_summary$mean_total) * 1.28) +
  labs(x = NULL, y = "Mean NET score per sample",
       title = "Sample-level NET score across all cells",
       subtitle = "Each point is one specimen; symbol size reflects the number of cells") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6))
ggsave(file.path(fig_dir, "NETsig_allcells_sample_level.pdf"), p4,
       width = 5.2, height = 4.8)

# ============================================================
# 图 5: 38 个基因在全部细胞中的表达点图
# ============================================================
p5 <- DotPlot(seu, features = rev(present), group.by = "group", dot.scale = 7) +
  scale_colour_gradient2(low = col_down, mid = "white", high = col_up,
                         midpoint = 0) +
  coord_flip() +
  labs(x = NULL, y = NULL, colour = "Scaled mean\nexpression",
       title = "38-gene NET panel across all cells (AKI vs CTRL)",
       subtitle = "All cell types pooled; dot size = fraction of expressing cells") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 8.5),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 8, 8))
ggsave(file.path(fig_dir, "NETsig_allcells_dotplot.pdf"), p5,
       width = 5.8, height = 9.2)

# ---------- 5. 运行记录 ----------
sink(file.path(log_dir, "sessionInfo_07_netsig_allcells.txt"))
cat("=== 07_netsig_allcells.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 细胞构成 ===\n")
print(table(seu$cell_group, seu$group))
cat("\n=== 细胞类型平均占比 ===\n")
print(composition %>% group_by(group, cell_group) %>%
        summarise(mean_fraction = mean(fraction), .groups = "drop") %>%
        pivot_wider(names_from = cell_group, values_from = mean_fraction))
cat("\n=== 检验结果 ===\n")
print(tests, digits = 3)
cat("\n=== 分细胞类型检验 ===\n")
print(ct_p)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 07_netsig_allcells.R 完成")
