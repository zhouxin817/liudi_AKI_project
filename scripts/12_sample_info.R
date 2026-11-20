# ============================================================
# 12_sample_info.R
# 功能：整理数据集中各样本(specimen)的详细信息
#       - 分组(condition.l1/l2/l3)、患者编号、性别、种族、组织类型
#       - 各样本总细胞数及目标细胞群体细胞数
#       - 是否纳入本次 AKI vs CTRL 分析
# 说明：公开元数据中不含年龄字段, 因此年龄列标注为未提供
# 输入：GSE183276_..._Metadata_03282022.txt.gz
# 输出：results/tables/Sample_information.csv
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
tab_dir  <- file.path(proj_dir, "results", "tables")
meta_file <- file.path(proj_dir,
  "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")

# ---------- 1. 读取所需元数据列 ----------
header_line <- readLines(gzfile(meta_file), n = 1)
col_names <- c("cell_id", strsplit(gsub('"', "", header_line), "\t")[[1]])
keep_cols <- c("cell_id", "specimen", "patient", "condition.l1", "condition.l2",
               "condition.l3", "sex", "race", "tissue_type", "experiment.short",
               "subclass.l2")
keep_idx <- which(col_names %in% keep_cols)
meta <- fread(meta_file, sep = "\t", header = FALSE, skip = 1, select = keep_idx)
setnames(meta, col_names[keep_idx])

# ---------- 2. 每个样本的患者/分组/人口学信息 ----------
# 同一 specimen 的信息可能因个别细胞缺失而出现多值, 取首个非缺失值
first_non_na <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) NA_character_ else x[1]
}

sample_info <- meta %>%
  group_by(specimen) %>%
  summarise(
    patient          = first_non_na(patient),
    condition_l1     = first_non_na(condition.l1),
    condition_l2     = first_non_na(condition.l2),
    condition_l3     = first_non_na(condition.l3),
    sex              = first_non_na(sex),
    race             = first_non_na(race),
    tissue_type      = first_non_na(tissue_type),
    experiment       = first_non_na(experiment.short),
    total_cells      = n(),
    .groups = "drop"
  )

# ---------- 3. 各目标细胞群体细胞数 ----------
count_subclass <- function(ct) {
  meta %>%
    filter(subclass.l2 == ct) %>%
    count(specimen, name = ct)
}
ct_tab <- Reduce(function(a, b) full_join(a, b, by = "specimen"),
                 lapply(c("N", "EC-PTC", "EC-GC", "PT-S1/2", "PT-S3"),
                        count_subclass))
ct_tab[is.na(ct_tab)] <- 0

sample_info <- sample_info %>%
  left_join(ct_tab, by = "specimen") %>%
  mutate(
    Neutrophil = N,
    PTEC       = `EC-PTC`,
    GEC        = `EC-GC`,
    PT         = `PT-S1/2` + `PT-S3`
  ) %>%
  select(-N, -`EC-PTC`, -`EC-GC`, -`PT-S1/2`, -`PT-S3`)

# ---------- 4. 分组与纳入标记 ----------
sample_info <- sample_info %>%
  mutate(
    # 本次分析只纳入 AKI 与 Ref(对照)
    analysis_group = case_when(
      condition_l1 == "AKI" ~ "AKI",
      condition_l1 == "Ref" ~ "CTRL",
      TRUE                  ~ "Excluded (CKD)"),
    included       = analysis_group %in% c("AKI", "CTRL"),
    # 公开元数据无年龄字段
    age            = NA_character_,
    age_note       = "Not provided in public metadata"
  ) %>%
  arrange(factor(analysis_group,
                 levels = c("AKI", "CTRL", "Excluded (CKD)")),
          specimen)

# 输出列顺序
sample_info <- sample_info %>%
  select(specimen, analysis_group, included, patient, condition_l1,
         condition_l2, condition_l3, sex, race, tissue_type, experiment,
         age, age_note, total_cells, Neutrophil, PTEC, GEC, PT,
         everything())

write.csv(sample_info, file.path(tab_dir, "Sample_information.csv"),
          row.names = FALSE)

# ---------- 5. 汇总输出 ----------
cat("=== 样本数 ===\n"); print(table(sample_info$analysis_group))
cat("\n=== 患者数 ===\n")
cat("AKI:", n_distinct(sample_info$patient[sample_info$analysis_group == "AKI"]),
    "; CTRL(Ref):",
    n_distinct(sample_info$patient[sample_info$analysis_group == "CTRL"]),
    "; CKD:", n_distinct(sample_info$patient[sample_info$analysis_group == "Excluded (CKD)"]), "\n")
cat("\n=== 性别分布(样本数) ===\n")
print(table(sample_info$analysis_group, sample_info$sex))
cat("\n=== 本次分析样本(AKI + CTRL) ===\n")
print(as.data.frame(
  sample_info %>% filter(included) %>%
    select(specimen, analysis_group, patient, condition_l2, sex, race,
           total_cells, Neutrophil, PTEC, GEC, PT)), right = FALSE)

# ---------- 6. 运行记录 ----------
sink(file.path(proj_dir, "logs", "sessionInfo_12_sample_info.txt"))
cat("=== 12_sample_info.R 运行时间 ===\n"); print(Sys.time())
cat("\n=== 样本信息 ===\n"); print(as.data.frame(sample_info), right = FALSE)
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 12_sample_info.R 完成")
