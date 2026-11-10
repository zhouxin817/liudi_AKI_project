# ============================================================
# 02_umap_markers.R
# 功能：绘制提取后细胞的 UMAP 图(按细胞群体 / 分组 / 聚类),
#       并绘制各群体特征性 marker 基因的 DotPlot
# 输入：results/processed_data/seu_subset_reclustered.rds
# 输出：results/figures/UMAP_*.pdf, DotPlot_canonical_markers.pdf
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 读取重新降维聚类后的子集对象 ----------
seu_sub <- readRDS(file.path(proj_dir, "results", "processed_data",
                             "seu_subset_reclustered.rds"))
cat("可用 reductions:", paste(names(seu_sub@reductions), collapse = ", "), "\n")

# 统一使用重新计算的 umap 嵌入
stopifnot("umap" %in% names(seu_sub@reductions))

# 目标群体配色
group_cols <- c("Neutrophil" = "#E64B35",
                "PTEC"       = "#4DBBD5",
                "GEC"        = "#00A087",
                "PT"         = "#3C5488")

# ---------- 2. UMAP: 按目标细胞群体着色 ----------
p1 <- DimPlot(seu_sub, reduction = "umap", group.by = "cell_group", pt.size = 0.5) +
  scale_color_manual(values = group_cols) +
  ggtitle("UMAP of Extracted Cell Populations (Reclustered)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_cell_groups.pdf"), p1, width = 8, height = 6)

# ---------- 3. UMAP: 按聚类结果着色 ----------
p2 <- DimPlot(seu_sub, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 4, pt.size = 0.5) +
  ggtitle("UMAP of Unsupervised Clusters (resolution = 0.5)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_clusters.pdf"), p2, width = 8, height = 6)

# ---------- 4. UMAP: 按 AKI / CTRL 分组着色 ----------
p3 <- DimPlot(seu_sub, reduction = "umap", group.by = "group", pt.size = 0.5,
              cols = c("CTRL" = "#7F7F7F", "AKI" = "#DC0000")) +
  ggtitle("UMAP by Condition (AKI vs CTRL)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_AKI_vs_CTRL.pdf"), p3, width = 8, height = 6)

# ---------- 5. UMAP: 按分组分面展示各群体分布 ----------
p4 <- DimPlot(seu_sub, reduction = "umap", group.by = "cell_group",
              split.by = "group", pt.size = 0.5, cols = group_cols) +
  ggtitle("UMAP Split by Condition (Colored by Cell Population)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_split_by_condition.pdf"), p4, width = 13, height = 6)

# ---------- 6. 特征性 marker DotPlot ----------
# 各群体经典 marker(依据肾脏单细胞文献, 并经本数据表达验证):
#   Neutrophil: S100A8, S100A9, CSF3R
#   PTEC(肾小管周毛细血管内皮): PLVAP, VWF, KDR
#   GEC(肾小球毛细血管内皮)   : EHD3, SOST, EMCN
#   PT(健康近端小管)          : LRP2, CUBN, SLC22A6
markers <- c("S100A8", "S100A9", "CSF3R",
             "PLVAP", "VWF", "KDR",
             "EHD3", "SOST", "EMCN",
             "LRP2", "CUBN", "SLC22A6")

markers_found   <- markers[markers %in% rownames(seu_sub)]
markers_missing <- setdiff(markers, markers_found)
if (length(markers_missing) > 0) {
  cat("以下 marker 在数据中缺失, 已跳过:", paste(markers_missing, collapse = ", "), "\n")
}

Idents(seu_sub) <- "cell_group"
# 注: DotPlot 默认对平均表达做 z-score, 故颜色 midpoint 设为 0
p5 <- DotPlot(seu_sub, features = markers_found, group.by = "cell_group") +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                        midpoint = 0, limits = c(-2.5, 2.5),
                        oob = scales::squish) +
  labs(x = "Marker genes", y = "Cell populations",
       title = "Canonical Markers of Target Cell Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "DotPlot_canonical_markers.pdf"), p5, width = 10, height = 5)

# ---------- 7. 记录运行环境 ----------
sink(file.path(proj_dir, "logs", "sessionInfo_02_umap_markers.txt"))
cat("=== 02_umap_markers.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 02_umap_markers.R 完成")
