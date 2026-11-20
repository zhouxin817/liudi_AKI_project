# ============================================================
# rerun_20261119 / 04_visualizations.R
# 功能:
#   1) 中性粒细胞占比箱线图: 每个样本占比 + Wilcoxon 检验
#      分母 = 原始 metadata 中该样本的全量总细胞数
#   2) PTEC GO(BP) 富集双向条形图: 复用指定样式
#      (去冗余 + 配色 + 侧边竖排标注 + theme_classic)
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

proj_dir  <- "/home/zhouxin/liudi_AKI_project"
rerun_dir <- file.path(proj_dir, "rerun_20261119")
fig_dir   <- file.path(rerun_dir, "results", "figures")
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")

col_up   <- "#CE6C8E"   # PTEC GO Up in AKI
col_down <- "#5E8DBC"   # PTEC GO Down in AKI
col_aki  <- "#E17B6E"   # AKI 组
col_ctrl <- "#7FC4C9"   # CTRL 组

# ---------- 1. 中性粒细胞占比箱线图 ----------
meta_file <- file.path(proj_dir,
  "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")
header_line <- readLines(gzfile(meta_file), n = 1)
col_names <- c("cell_id", strsplit(gsub('"', "", header_line), "\t")[[1]])
keep_idx <- which(col_names %in%
  c("cell_id", "specimen", "condition.l1", "subclass.l2"))
meta_full <- fread(meta_file, sep = "\t", header = FALSE,
                   skip = 1, select = keep_idx)
setnames(meta_full, col_names[keep_idx])

# 4 细胞子集里保留的样本
kept <- read.csv(file.path(tab_dir, "per_sample_selected.csv"),
                 stringsAsFactors = FALSE)

# 每个样本: N 细胞数 + 全量总细胞数
per_specimen <- meta_full %>%
  filter(specimen %in% kept$specimen) %>%
  group_by(specimen) %>%
  summarise(Neutrophil = sum(subclass.l2 == "N"),
            total_all  = n(),
            group      = first(
              recode(as.character(condition.l1),
                     "AKI" = "AKI", "Ref" = "CTRL")),
            .groups = "drop") %>%
  mutate(group = factor(group, levels = c("CTRL", "AKI")),
         pct   = Neutrophil / total_all * 100)

write.csv(per_specimen,
          file.path(tab_dir, "Neutrophil_proportion_per_sample.csv"),
          row.names = FALSE)

# 箱线图: 散点叠箱线, 组间括号 + Wilcoxon p 值
wt <- wilcox.test(pct ~ group, data = per_specimen, exact = FALSE)
p_lab <- ifelse(wt$p.value < 1e-4,
                paste0("P = ", formatC(wt$p.value, format = "e", digits = 2)),
                paste0("P = ", formatC(wt$p.value, format = "f", digits = 4)))

smry <- per_specimen %>% group_by(group) %>%
  summarise(mean = mean(pct), sd = sd(pct), n = n(), .groups = "drop")

y_top  <- max(per_specimen$pct, na.rm = TRUE)
y_brk  <- y_top * 1.18
y_tick <- y_top * 0.03

p_neu <- ggplot(per_specimen, aes(x = group, y = pct, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "grey35", linewidth = 0.6) +
  geom_jitter(shape = 21, size = 3.8, colour = "grey20",
              stroke = 0.9, width = 0.08, alpha = 0.95) +
  geom_errorbar(data = smry, aes(x = group, y = mean,
                                  ymin = mean - sd, ymax = mean + sd),
                width = 0.24, colour = "grey35", linewidth = 0.5,
                inherit.aes = FALSE) +
  geom_segment(data = smry, aes(x = as.numeric(group) - 0.08,
                                 xend = as.numeric(group) + 0.08,
                                 y = mean, yend = mean),
               colour = "grey35", linewidth = 0.9,
               inherit.aes = FALSE) +
  annotate("segment", x = 1, xend = 2,
           y = y_brk, yend = y_brk,
           colour = "grey25", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1,
           y = y_brk - y_tick, yend = y_brk,
           colour = "grey25", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2,
           y = y_brk - y_tick, yend = y_brk,
           colour = "grey25", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_brk * 1.04,
           label = p_lab, size = 3.8, colour = "black") +
  annotate("text", x = 1, y = -Inf, vjust = -0.5,
           label = paste0("n = ", smry$n[smry$group == "CTRL"]),
           size = 3.2, colour = "grey35") +
  annotate("text", x = 2, y = -Inf, vjust = -0.5,
           label = paste0("n = ", smry$n[smry$group == "AKI"]),
           size = 3.2, colour = "grey35") +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki),
                    guide = "none") +
  labs(x = NULL,
       y = "Neutrophil proportion (%)",
       title = "Neutrophil proportion per sample") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = margin(10, 10, 12, 6))

ggsave(file.path(fig_dir, "Neutrophil_proportion_boxplot.pdf"),
       p_neu, width = 4.6, height = 5.4)

# ---------- 2. PTEC GO(BP) 双向条形图 ----------
go_tab <- read.csv(file.path(tab_dir, "GO_BP_PTEC_AKI_vs_CTRL.csv"),
                   stringsAsFactors = FALSE)

dedup_terms <- function(df, top_n = 8, jaccard_cut = 0.7) {
  df <- df[order(df$p.adjust), ]
  kept <- character(0)
  kept_genes <- list()
  for (i in seq_len(nrow(df))) {
    g <- strsplit(df$geneID[i], "/")[[1]]
    redundant <- FALSE
    for (kg in kept_genes) {
      inter <- length(intersect(g, kg))
      uni   <- length(union(g, kg))
      if (uni > 0 && inter / uni > jaccard_cut) {
        redundant <- TRUE; break
      }
    }
    if (!redundant) {
      kept <- c(kept, df$Description[i])
      kept_genes[[length(kept_genes) + 1]] <- g
    }
    if (length(kept) >= top_n) break
  }
  df[df$Description %in% kept, ]
}

up_df   <- dedup_terms(go_tab[go_tab$direction == "Up_in_AKI",   ], top_n = 8)
down_df <- dedup_terms(go_tab[go_tab$direction == "Down_in_AKI", ], top_n = 10)

bar_df <- bind_rows(
  up_df   %>% mutate(signed_logp = -log10(p.adjust), dir = "up"),
  down_df %>% mutate(signed_logp =  log10(p.adjust), dir = "down"))

# y 轴顺序: 上调由高到低在上部, 下调由高到低在下部, 中线在中间
ord_up   <- up_df$Description[order(-(-log10(up_df$p.adjust)))]
ord_down <- down_df$Description[order(-(-log10(down_df$p.adjust)))]
term_levels <- c(rev(ord_down), rev(ord_up))

bar_df$TermLabel <- factor(
  str_wrap(bar_df$Description, width = 55),
  levels = str_wrap(term_levels, width = 55))
bar_df$dir <- factor(bar_df$dir, levels = c("down", "up"))

x_max  <- max(abs(bar_df$signed_logp))
x_lim  <- x_max * 1.38
lab_up <- x_max * 1.20
lab_dn <- -x_max * 1.20

n_down <- length(ord_down)
n_up   <- length(ord_up)
y_down <- c(0.6, n_down + 0.4)
y_up   <- c(n_down + 0.6, n_down + n_up + 0.4)

p_go <- ggplot(bar_df, aes(x = signed_logp, y = TermLabel, fill = dir)) +
  geom_col(width = 0.68) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  # 右侧竖线 + "enriched in AKI PTEC"
  annotate("segment",
           x = lab_up + x_max * 0.05,
           xend = lab_up + x_max * 0.05,
           y = y_up[1], yend = y_up[2],
           colour = col_up, linewidth = 0.5) +
  annotate("text", x = lab_up, y = mean(y_up),
           label = "enriched in\nAKI PTEC",
           colour = col_up, angle = 90,
           size = 3.4, lineheight = 0.95) +
  # 左侧竖线 + "enriched in Control PTEC"
  annotate("segment",
           x = lab_dn - x_max * 0.05,
           xend = lab_dn - x_max * 0.05,
           y = y_down[1], yend = y_down[2],
           colour = col_down, linewidth = 0.5) +
  annotate("text", x = lab_dn, y = mean(y_down),
           label = "enriched in\nControl PTEC",
           colour = col_down, angle = 90,
           size = 3.4, lineheight = 0.95) +
  scale_fill_manual(values = c(down = col_down, up = col_up),
                    guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = pretty(c(-x_lim, x_lim), n = 7)) +
  labs(x = expression(-log[10]~"(P.adjust)"), y = NULL,
       title = "PTEC AKI vs Control geneset enrichment analysis") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(
          colour = ifelse(
            levels(bar_df$TermLabel) %in% str_wrap(ord_up, 55),
            col_up, col_down),
          size = 8.3, face = "italic"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black"),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = margin(10, 12, 8, 6))

ggsave(file.path(fig_dir, "PTEC_GO_BP_barplot.pdf"),
       p_go, width = 9.0, height = 6.6)

# ---------- 3. 运行记录 ----------
sink(file.path(log_dir, "sessionInfo_04_visualizations.txt"))
cat("=== 04_visualizations.R ===\n"); print(Sys.time())
cat("\n=== Neutrophil 占比 ===\n")
cat("Wilcoxon p = ", wt$p.value, "\n")
print(as.data.frame(per_specimen))
cat("\n=== PTEC GO bar data ===\n"); print(bar_df[, c("Description","dir","signed_logp")])
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()
message(">>> 04_visualizations.R 完成")
