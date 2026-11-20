# ============================================================
# 11_ptec_mito_heatmap.R
# 功能：绘制 PTEC 线粒体相关 TOP50 差异基因的表达热图
#       1) 基因集合: GO:0005739(mitchondrion)及其子孙条目注释到的基因
#          (排除 MT- 线粒体基因组基因, 避免技术性/计量偏差主导配色)
#       2) TOP50 选取: 上调、下调各取 25 个(按 p_val_adj 排序),
#          不足时由另一方向补足
#       3) 热图数值: 各样本(剔除细胞数过少的样本)的伪批量平均表达,
#          按基因做 z-score 标准化
# 说明：不修改任何既有文件与数据, 仅新增 outputs 到 results/ptec_reanalysis/
# 输入：results/processed_data/seu_subset_reclustered.rds
#       results/ptec_reanalysis/tables/PTEC_DEG_AKI_vs_CTRL_updated.csv
# 输出：results/ptec_reanalysis/figures/PTEC_mito_top50_heatmap.pdf
#       results/ptec_reanalysis/tables/PTEC_mito_top50_heatmap_genes.csv
#       results/ptec_reanalysis/tables/PTEC_mito_top50_pseudobulk_zscore.csv
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
out_dir  <- file.path(proj_dir, "results", "ptec_reanalysis")
fig_dir  <- file.path(out_dir, "figures")
tab_dir  <- file.path(out_dir, "tables")

# 与既有图一致的配色: 上调(粉) / 下调(蓝) / 对照组(青绿)
col_up   <- "#CE6C8E"
col_down <- "#5E8DBC"
col_aki  <- "#E17B6E"
col_ctrl <- "#7FC4C9"

# ---------- 1. 读取对象与差异结果 ----------
seu <- readRDS(file.path(proj_dir, "results", "processed_data",
                         "seu_subset_reclustered.rds"))
DefaultAssay(seu) <- "RNA"
ptec <- subset(seu, subset = cell_group == "PTEC")
ptec$group <- factor(ptec$group, levels = c("CTRL", "AKI"))

deg <- read.csv(file.path(tab_dir, "PTEC_DEG_AKI_vs_CTRL_updated.csv"),
                stringsAsFactors = FALSE)

# ---------- 2. 线粒体相关基因集合 ----------
# GO:0005739(mitchondrion) 及其子孙条目 -> 基因(经 GOALL 传递注释)
get_offspring <- function(ids, bimap) {
  lst <- as.list(bimap)
  unique(unlist(lst[ids], use.names = FALSE))
}
mito_cc_ids <- unique(c("GO:0005739",
                        get_offspring("GO:0005739", GOCCOFFSPRING)))
mito_symbols <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = mito_cc_ids,
                        columns = "SYMBOL", keytype = "GOALL")$SYMBOL)
mito_symbols <- unique(mito_symbols[!is.na(mito_symbols)])

# 在本次 PTEC 检验的基因中筛选线粒体相关基因
mito_deg <- deg %>%
  filter(gene %in% mito_symbols, !grepl("^MT-", gene))

# ---------- 3. 选取 TOP50(上调/下调各 25) ----------
n_each <- 25
top_up <- mito_deg %>%
  filter(direction == "Up in AKI") %>% arrange(p_val_adj) %>% head(n_each)
top_down <- mito_deg %>%
  filter(direction == "Down in AKI") %>% arrange(p_val_adj) %>% head(n_each)

# 某一方向不足 25 个时, 由另一方向补足到 50
if (nrow(top_up) + nrow(top_down) < 2 * n_each) {
  rest <- mito_deg %>%
    filter(!gene %in% c(top_up$gene, top_down$gene)) %>%
    arrange(p_val_adj) %>%
    head(2 * n_each - nrow(top_up) - nrow(top_down))
  top_up <- bind_rows(top_up, rest %>% filter(direction == "Up in AKI"))
  top_down <- bind_rows(top_down, rest %>% filter(direction == "Down in AKI"))
}

top50 <- bind_rows(top_up, top_down) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(direction = factor(direction,
                            levels = c("Up in AKI", "Down in AKI")))
genes <- top50$gene
cat("TOP50 基因数:", length(genes),
    "; Up:", sum(top50$direction == "Up in AKI"),
    "; Down:", sum(top50$direction == "Down in AKI"), "\n")

write.csv(top50 %>%
            select(gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, direction),
          file.path(tab_dir, "PTEC_mito_top50_heatmap_genes.csv"),
          row.names = FALSE)

# ---------- 4. 各样本伪批量平均表达 ----------
# 基因表达取 log-normalized data; 每样本计算基因均值
data_mat <- tryCatch(
  GetAssayData(ptec, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(ptec, assay = "RNA", slot = "data"))
data_mat <- data_mat[genes, , drop = FALSE]

meta <- ptec@meta.data
meta$specimen <- as.character(meta$specimen)
# 剔除细胞数过少的样本, 避免伪批量均值不稳定
min_cells <- 10
keep_samples <- names(which(table(meta$specimen) >= min_cells))
meta_keep <- meta %>% filter(specimen %in% keep_samples)

pb <- sapply(keep_samples, function(s) {
  cells <- rownames(meta_keep)[meta_keep$specimen == s]
  Matrix::rowMeans(data_mat[, cells, drop = FALSE])
})
# 样本按 CTRL -> AKI 排序, 组内按样本名排序
sample_group <- meta_keep %>%
  distinct(specimen, group) %>%
  mutate(group = as.character(group)) %>%
  arrange(group, specimen)
pb <- pb[, sample_group$specimen, drop = FALSE]

# 基因方向注释(用于热图行分割)
gene_dir <- top50$direction[match(rownames(pb), top50$gene)]

# ---------- 5. 按基因 z-score 标准化 ----------
z <- t(scale(t(pb)))
z[is.na(z)] <- 0
z <- pmax(pmin(z, 2), -2)   # 截断极端值, 保证配色稳定
colnames(z) <- sample_group$specimen

write.csv(data.frame(gene = rownames(z), z, check.names = FALSE),
          file.path(tab_dir, "PTEC_mito_top50_pseudobulk_zscore.csv"),
          row.names = FALSE)

# ---------- 6. 绘制热图 ----------
# 行顺序: 上调基因(按 log2FC 由高到低) 在前, 下调基因(由低到高) 在后
ord_up   <- top50 %>% filter(direction == "Up in AKI") %>%
  arrange(desc(avg_log2FC)) %>% pull(gene)
ord_down <- top50 %>% filter(direction == "Down in AKI") %>%
  arrange(avg_log2FC) %>% pull(gene)
row_order <- c(ord_up, ord_down)
row_order <- row_order[row_order %in% rownames(z)]
z <- z[row_order, , drop = FALSE]
gene_dir <- gene_dir[match(row_order, rownames(z))]

col_fun <- colorRamp2(c(-2, 0, 2), c(col_down, "white", col_up))

ha <- HeatmapAnnotation(
  Group = sample_group$group,
  col = list(Group = c(CTRL = col_ctrl, AKI = col_aki)),
  annotation_name_side = "left",
  annotation_legend_param = list(Group = list(title = "Group"))
)

ht <- Heatmap(
  z,
  name = "z-score",
  col = col_fun,
  row_labels = rownames(z),
  row_names_gp = gpar(fontface = "italic", fontsize = 7.5),
  row_names_side = "right",
  row_split = factor(gene_dir, levels = c("Up in AKI", "Down in AKI")),
  row_title_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(sample_group$group, levels = c("CTRL", "AKI")),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 7),
  top_annotation = ha,
  column_gap = unit(1.2, "mm"),
  row_gap = unit(1.2, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Scaled mean\nexpression (z)"
  )
)

pdf(file.path(fig_dir, "PTEC_mito_top50_heatmap.pdf"),
    width = 9.5, height = 9.8)
draw(ht,
     column_title = "Mitochondria-related TOP50 genes in PTEC (AKI vs CTRL)",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     column_title_side = "top",
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     # 底部留出空间, 避免旋转后的样本名被页面裁切
     padding = unit(c(14, 2, 2, 2), "mm"))
invisible(dev.off())

# ---------- 7. 运行记录(可重复性) ----------
sink(file.path(proj_dir, "logs", "sessionInfo_11_ptec_mito_heatmap.txt"))
cat("=== 11_ptec_mito_heatmap.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== PTEC 细胞数 ===\n")
print(table(ptec$group))
cat("\n=== TOP50 基因 ===\n"); print(top50[, c("gene", "avg_log2FC", "p_val_adj", "direction")])
cat("\n=== 参与伪批量的样本数 ===\n")
print(table(sample_group$group))
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 11_ptec_mito_heatmap.R 完成")
