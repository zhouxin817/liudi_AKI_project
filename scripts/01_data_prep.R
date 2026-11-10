# ============================================================
# 01_data_prep.R
# 功能：加载 GSE183276 h5Seurat 数据, 提取目标细胞群体
#       (Neutrophil / GEC / PTEC / PT) 及 AKI 与 Ref(CTRL) 样本,
#       对提取后的细胞重新进行标准化、降维与聚类,
#       保存子集 Seurat 对象供下游分析
# 数据集：GSE183276 Kidney Healthy-Injury Cell Atlas (scCv3)
# ============================================================

# ---------- 0. 环境准备 ----------
# 固定随机种子, 保证全流程结果可重复
set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
})

# 项目路径
proj_dir <- "/home/zhouxin/liudi_AKI_project"
h5_file  <- file.path(proj_dir, "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
out_dir  <- file.path(proj_dir, "results", "processed_data")
tab_dir  <- file.path(proj_dir, "results", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 加载 h5Seurat 原始数据 ----------
# 该对象包含原始 counts、meta.data 与图谱原始 UMAP 嵌入
message(">>> 加载 h5Seurat 数据 ...")
seu <- LoadH5Seurat(h5_file)
message("全部细胞数: ", ncol(seu))

# ---------- 2. 定义目标细胞群体 ----------
# 依据 subclass.l2 注释定义四个目标群体:
#   Neutrophil : N                —— 中性粒细胞
#   GEC        : EC-GC            —— 肾小球毛细血管内皮细胞
#   PT         : PT-S1/2, PT-S3   —— 健康近端肾小管上皮细胞
#   PTEC       : EC-PTC           —— 肾小管周毛细血管内皮细胞
group_map <- list(
  Neutrophil = "N",
  GEC        = "EC-GC",
  PT         = c("PT-S1/2", "PT-S3"),
  PTEC       = "EC-PTC"
)
target_subclasses <- unlist(group_map, use.names = FALSE)

# ---------- 3. 提取目标细胞与目标分组 ----------
# 分组依据 condition.l1: AKI(疾病组) 与 Ref(正常对照, 记为 CTRL)
# CKD 相关样本不参与本次分析
cells_keep <- WhichCells(
  seu,
  expression = condition.l1 %in% c("AKI", "Ref") & subclass.l2 %in% target_subclasses
)
message("符合条件的细胞数: ", length(cells_keep))

seu_sub <- subset(seu, cells = cells_keep)

# ---------- 4. 添加分析用标签 ----------
# group: AKI / CTRL(二分组)
seu_sub$group <- factor(ifelse(seu_sub$condition.l1 == "AKI", "AKI", "CTRL"),
                        levels = c("CTRL", "AKI"))

# cell_group: 目标细胞群体名称
seu_sub$cell_group <- dplyr::case_when(
  seu_sub$subclass.l2 == "N"                     ~ "Neutrophil",
  seu_sub$subclass.l2 == "EC-GC"                 ~ "GEC",
  seu_sub$subclass.l2 %in% c("PT-S1/2", "PT-S3") ~ "PT",
  seu_sub$subclass.l2 == "EC-PTC"                ~ "PTEC"
)
seu_sub$cell_group <- factor(seu_sub$cell_group,
                             levels = c("Neutrophil", "PTEC", "GEC", "PT"))

# 各群体 × 分组细胞数核查并导出
cnt <- as.data.frame(table(CellGroup = seu_sub$cell_group, Group = seu_sub$group))
write.csv(cnt, file.path(tab_dir, "Cell_counts_by_group.csv"), row.names = FALSE)
cat("\n=== 各细胞群体 × 分组细胞数 ===\n")
print(table(seu_sub$cell_group, seu_sub$group))

# ---------- 5. 重新降维与聚类 ----------
# 说明: 不使用图谱原有嵌入, 对提取后的 4 类细胞重新走标准流程,
#       以获得仅针对本分析子集的 UMAP 与聚类结果
seu_sub <- NormalizeData(seu_sub, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)

seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)

# 仅对高变基因做标准化, 降低内存占用
seu_sub <- ScaleData(seu_sub, features = VariableFeatures(seu_sub), verbose = FALSE)

seu_sub <- RunPCA(seu_sub, npcs = 30, features = VariableFeatures(seu_sub),
                  seed.use = 2026, verbose = FALSE)

seu_sub <- FindNeighbors(seu_sub, dims = 1:20, verbose = FALSE)

# 聚类: resolution = 0.5, 固定随机种子
seu_sub <- FindClusters(seu_sub, resolution = 0.5, algorithm = 1,
                        random.seed = 2026, verbose = FALSE)

seu_sub <- RunUMAP(seu_sub, dims = 1:20, seed.use = 2026, verbose = FALSE)

# 聚类结果与人工注释的一致性核查
cat("\n=== 聚类结果 × 目标细胞群体 ===\n")
print(table(seu_sub$seurat_clusters, seu_sub$cell_group))

# 保存聚类-注释对应表
cluster_tab <- as.data.frame(table(Cluster = seu_sub$seurat_clusters,
                                   CellGroup = seu_sub$cell_group))
write.csv(cluster_tab, file.path(tab_dir, "Cluster_vs_cell_group.csv"), row.names = FALSE)

# ---------- 6. 默认 assay 与 identity ----------
DefaultAssay(seu_sub) <- "RNA"
Idents(seu_sub) <- "cell_group"

# ---------- 7. 保存子集对象 ----------
saveRDS(seu_sub, file.path(out_dir, "seu_subset_reclustered.rds"))
message(">>> 子集对象已保存 (含新 PCA / UMAP / 聚类结果)")

# ---------- 8. 记录运行环境(可重复性) ----------
sink(file.path(proj_dir, "logs", "sessionInfo_01_data_prep.txt"))
cat("=== 01_data_prep.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 01_data_prep.R 完成")
