# ============================================================
# rerun_20261119_noSample / 01_prep_subset.R
# 功能：从全量原始 metadata + 4 细胞子集 Seurat 对象出发,
#       按样本质量与目标细胞数筛选合适样本, 完成 cell_group / group 标注
#       本目录版本额外剔除 specimen 以 "Sample" 开头的样本(共 9 个, 均为 Ref/CTRL)
#       说明: 前序脚本已经产出了 4 类细胞的 Seurat 子集, 这里直接从该子集上
#             取筛选结果, 避免重复加载 5.8GB 的原始 h5Seurat
# 目标细胞: Neutrophil(N), PTEC(EC-PTC), GEC(EC-GC), PT(PT-S1/2 + PT-S3)
# 筛选原则(保守, 兼顾细胞数与样本质量):
#   - 全量总细胞数 >= 800(样本质量下限)
#   - PTEC >= 10 / GEC >= 5(内皮差异分析最低细胞数)
#   - 剔除 CKD(condition.l1 != AKI / Ref)
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

proj_dir  <- "/home/zhouxin/liudi_AKI_project"
rerun_dir <- file.path(proj_dir, "rerun_20261119_noSample")
rerun_fig <- file.path(rerun_dir, "results", "figures")
rerun_tab <- file.path(rerun_dir, "results", "tables")
rerun_log <- file.path(rerun_dir, "logs")
dir.create(rerun_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(rerun_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(rerun_log, recursive = TRUE, showWarnings = FALSE)

# ---------- 0. 筛选阈值 ----------
min_total  <- 800
min_ptec   <- 10
min_gec    <- 5

# ---------- 1. 从原始 metadata 算每样本全量总细胞 ----------
meta_file <- file.path(proj_dir,
  "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")
header_line <- readLines(gzfile(meta_file), n = 1)
col_names <- c("cell_id", strsplit(gsub('"', "", header_line), "\t")[[1]])
keep_idx <- which(col_names %in% c("cell_id", "specimen", "condition.l1"))
meta_full <- fread(meta_file, sep = "\t", header = FALSE,
                   skip = 1, select = keep_idx)
setnames(meta_full, col_names[keep_idx])

full_total <- meta_full %>%
  group_by(specimen) %>%
  summarise(condition_l1 = first(condition.l1),
            total_all   = n(),
            .groups     = "drop") %>%
  filter(condition_l1 %in% c("AKI", "Ref"))

# ---------- 1b. 按要求剔除样本名以 "Sample" 开头的样本 ----------
excluded_sample <- full_total %>%
  filter(grepl("^Sample", specimen)) %>%
  arrange(specimen)
full_total <- full_total %>% filter(!grepl("^Sample", specimen))
cat("\n按要求剔除 specimen 以 Sample 开头的样本:", nrow(excluded_sample), "个\n")
print(excluded_sample %>% count(condition_l1))
write.csv(excluded_sample,
          file.path(rerun_tab, "excluded_Sample_prefixed_specimens.csv"),
          row.names = FALSE)

# ---------- 2. 4 细胞子集 Seurat 对象: 每样本目标细胞计数 ----------
seu <- readRDS(file.path(proj_dir, "results", "processed_data",
                         "seu_subset_reclustered.rds"))
DefaultAssay(seu) <- "RNA"

seu$cell_group <- dplyr::case_when(
  seu$subclass.l2 == "N"                 ~ "Neutrophil",
  seu$subclass.l2 == "EC-PTC"            ~ "PTEC",
  seu$subclass.l2 == "EC-GC"             ~ "GEC",
  seu$subclass.l2 %in% c("PT-S1/2", "PT-S3") ~ "PT",
  TRUE ~ NA_character_)

rds_counts <- seu@meta.data %>%
  filter(!is.na(cell_group)) %>%
  group_by(specimen) %>%
  summarise(Neutrophil = sum(cell_group == "Neutrophil"),
            PTEC       = sum(cell_group == "PTEC"),
            GEC        = sum(cell_group == "GEC"),
            PT         = sum(cell_group == "PT"),
            .groups    = "drop")

# ---------- 3. 合并并筛选 ----------
per_sample <- full_total %>%
  left_join(rds_counts, by = "specimen") %>%
  mutate(across(c(Neutrophil, PTEC, GEC, PT), ~ replace_na(.x, 0L)),
         group = recode(condition_l1, "AKI" = "AKI", "Ref" = "CTRL"),
         pass  = total_all >= min_total & PTEC >= min_ptec & GEC >= min_gec,
         note  = case_when(
           total_all < min_total ~ paste0("total_all<", min_total),
           PTEC      < min_ptec  ~ paste0("PTEC<", min_ptec),
           GEC       < min_gec   ~ paste0("GEC<", min_gec),
           TRUE ~ "PASS")) %>%
  arrange(group, desc(total_all))

cat("=== 筛选结果(按质量从高到低) ===\n")
print(as.data.frame(per_sample), right = FALSE)
cat("\n通过样本:", sum(per_sample$pass),
    "  剔除:", sum(!per_sample$pass), "\n")
print(per_sample %>% filter(pass) %>% count(group))

# ---------- 4. 提取保留样本的细胞 ----------
keep_specimens <- per_sample$specimen[per_sample$pass]
seu <- subset(seu, subset = specimen %in% keep_specimens)

# group 标注
seu$group <- factor(
  recode(as.character(seu$condition.l1),
         "AKI" = "AKI", "Ref" = "CTRL"),
  levels = c("CTRL", "AKI"))
seu$cell_group <- factor(seu$cell_group,
                         levels = c("Neutrophil", "PTEC", "GEC", "PT"))

cat("\n=== 保留细胞数(cell_group x group) ===\n")
print(table(seu$cell_group, seu$group))

# ---------- 5. 保存 ----------
saveRDS(seu, file.path(rerun_dir, "results", "seu_selected.rds"))
write.csv(per_sample, file.path(rerun_tab, "per_sample_before_filter.csv"),
          row.names = FALSE)
write.csv(per_sample %>% filter(pass) %>% select(-pass),
          file.path(rerun_tab, "per_sample_selected.csv"),
          row.names = FALSE)

cell_count <- seu@meta.data %>%
  count(cell_group, group, name = "n_cells") %>%
  arrange(group, cell_group)
write.csv(cell_count,
          file.path(rerun_tab, "cell_counts_by_group.csv"),
          row.names = FALSE)

# ---------- 6. 运行记录 ----------
sink(file.path(rerun_log, "sessionInfo_01_prep_subset.txt"))
cat("=== 01_prep_subset.R 运行时间 ===\n"); print(Sys.time())
cat("\n=== 筛选阈值 ===\n")
cat("total_all >=", min_total, "  PTEC >=", min_ptec, "  GEC >=", min_gec, "\n")
cat("\n=== 筛选结果 ===\n"); print(as.data.frame(per_sample))
cat("\n=== 最终细胞计数 ===\n"); print(cell_count)
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 01_prep_subset.R 完成")
