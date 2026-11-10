# ============================================================
# 03_neutrophil_stats.R
# 功能：统计 AKI 与 CTRL 各样本中中性粒细胞(Neutrophil)数量及
#       占样本总细胞数的比例, 进行 Wilcoxon 秩和检验并绘图
# 输入：原始 metadata(全量细胞, 用于计算样本总细胞数作为分母)
# 输出：results/figures/Neutrophil_*_stats.pdf
#       results/tables/Neutrophil_per_sample_stats.csv
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
tab_dir  <- file.path(proj_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 读取全量 metadata ----------
# 注意: 该文件第一列为细胞行名(无列名), header 列数比数据列少 1,
#       需先读取 header 再手动对齐列名, 避免列错位
meta_file <- file.path(proj_dir, "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")

header_line <- readLines(gzfile(meta_file), n = 1)
col_names <- c("cell_id", strsplit(gsub('"', '', header_line), "\t")[[1]])

# 仅提取分析所需列(按列位置索引, 防止列名错位)
keep_idx <- which(col_names %in% c("cell_id", "specimen", "patient",
                                   "condition.l1", "subclass.l2"))
meta <- fread(meta_file, sep = "\t", header = FALSE, skip = 1, select = keep_idx)
setnames(meta, col_names[keep_idx])

# ---------- 2. 仅保留 AKI 与 Ref(CTRL)样本 ----------
meta <- meta %>% filter(condition.l1 %in% c("AKI", "Ref"))
meta$group <- factor(ifelse(meta$condition.l1 == "AKI", "AKI", "CTRL"),
                     levels = c("CTRL", "AKI"))

# ---------- 3. 按样本统计中性粒细胞数量与占比 ----------
# 分子: subclass.l2 == "N" 的细胞数; 分母: 该样本全部细胞数
sample_total <- meta %>%
  group_by(specimen, group) %>%
  summarise(total_cells = n(), .groups = "drop")

sample_neu <- meta %>%
  filter(subclass.l2 == "N") %>%
  group_by(specimen, group) %>%
  summarise(neutrophils = n(), .groups = "drop")

stat_df <- sample_total %>%
  left_join(sample_neu, by = c("specimen", "group")) %>%
  mutate(neutrophils = ifelse(is.na(neutrophils), 0L, neutrophils),
         proportion  = neutrophils / total_cells * 100)

cat("=== 各样本中性粒细胞统计 ===\n")
print(as.data.frame(stat_df))

fwrite(stat_df, file.path(tab_dir, "Neutrophil_per_sample_stats.csv"))

# ---------- 4. 统计学检验: Wilcoxon 秩和检验 ----------
wt <- wilcox.test(proportion ~ group, data = stat_df, exact = FALSE)
cat("\n中性粒细胞占比 Wilcoxon 检验 p =", format.pval(wt$p.value, digits = 4), "\n")

wt_cnt <- wilcox.test(neutrophils ~ group, data = stat_df, exact = FALSE)
cat("中性粒细胞数量 Wilcoxon 检验 p =", format.pval(wt_cnt$p.value, digits = 4), "\n")

# ---------- 5. 绘图: 各样本中性粒细胞占比(带统计学标注) ----------
p1 <- ggplot(stat_df, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("CTRL" = "#7F7F7F", "AKI" = "#DC0000")) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("CTRL", "AKI")),
                     label = "p.format", size = 5) +
  labs(x = NULL, y = "Neutrophil proportion (%)",
       title = "Neutrophil Proportion per Sample") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "Neutrophil_proportion_stats.pdf"), p1, width = 5, height = 6)

# ---------- 6. 绘图: 各样本中性粒细胞绝对数量 ----------
p2 <- ggplot(stat_df, aes(x = group, y = neutrophils, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("CTRL" = "#7F7F7F", "AKI" = "#DC0000")) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("CTRL", "AKI")),
                     label = "p.format", size = 5) +
  labs(x = NULL, y = "Neutrophil count",
       title = "Neutrophil Count per Sample") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "Neutrophil_count_stats.pdf"), p2, width = 5, height = 6)

# ---------- 7. 记录运行环境 ----------
sink(file.path(proj_dir, "logs", "sessionInfo_03_neutrophil_stats.txt"))
cat("=== 03_neutrophil_stats.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 占比 Wilcoxon 检验 ===\n"); print(wt)
cat("\n=== 数量 Wilcoxon 检验 ===\n"); print(wt_cnt)
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 03_neutrophil_stats.R 完成")
