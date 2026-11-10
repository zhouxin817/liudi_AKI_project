# ============================================================
# 06_extra_figures.R
# 功能：按指定图形样式补充三类结果图
#   1) PTEC 差异基因 GO(BP) 富集双向条形图
#      (x 轴为带方向的 -log10(P.adjust), 左蓝 = CTRL 富集, 右粉 = AKI 富集)
#   2) GEC 差异基因火山图(标注 top10 基因)
#   3) PTEC / GEC 各样本细胞比例统计图
#      (散点 + 均值±标准差 + Wilcoxon 检验 p 值括号标注)
# 输入：results/tables/DEG_GEC_AKI_vs_CTRL.csv
#       results/tables/GO_BP_PTEC_AKI_vs_CTRL.csv
#       原始 metadata(计算各样本细胞比例的分母)
# 输出：results/figures/GO_BP_PTEC_enrichment_barplot.pdf
#       results/figures/Volcano_GEC_AKI_vs_CTRL.pdf
#       results/figures/Cell_proportion_PTEC.pdf
#       results/figures/Cell_proportion_GEC.pdf
#       results/tables/Cell_proportion_PTEC_GEC.csv
# ============================================================

# ---------- 0. 环境准备 ----------
set.seed(2026)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
tab_dir  <- file.path(proj_dir, "results", "tables")

# 统一配色: 疾病组(AKI, 鲑鱼粉) 与 对照组(CTRL, 青绿)
col_aki   <- "#E17B6E"
col_ctrl  <- "#7FC4C9"
# GO 条形图配色: 上调(粉) / 下调(蓝)
col_up    <- "#CE6C8E"
col_down  <- "#5E8DBC"

# ============================================================
# 一、PTEC GO(BP) 富集双向条形图
# ============================================================
go_tab <- read.csv(file.path(tab_dir, "GO_BP_PTEC_AKI_vs_CTRL.csv"),
                   stringsAsFactors = FALSE)

# 简化冗余条目: GO 条目间存在父子/近义关系, 若与已保留条目的基因重叠
# (Jaccard 系数) 超过阈值, 则视为冗余并跳过, 保证图中条目互不重复
dedup_terms <- function(df, top_n = 8, jaccard_cut = 0.7) {
  df <- df %>% arrange(p.adjust)
  kept <- character(0)   # 已保留条目名称
  kept_genes <- list()   # 已保留条目的基因集合
  for (i in seq_len(nrow(df))) {
    g <- strsplit(df$geneID[i], "/")[[1]]
    redundant <- FALSE
    for (kg in kept_genes) {
      inter <- length(intersect(g, kg))
      union <- length(union(g, kg))
      if (union > 0 && inter / union > jaccard_cut) { redundant <- TRUE; break }
    }
    if (!redundant) {
      kept <- c(kept, df$Description[i])
      kept_genes[[length(kept_genes) + 1]] <- g
    }
    if (length(kept) >= top_n) break
  }
  df %>% filter(Description %in% kept) %>% arrange(p.adjust)
}

up_df   <- dedup_terms(go_tab %>% filter(direction == "Up_in_AKI"),   top_n = 8)
down_df <- dedup_terms(go_tab %>% filter(direction == "Down_in_AKI"), top_n = 10)

# 组装绘图数据: 上调记正值, 下调记负值
bar_df <- bind_rows(
  up_df   %>% mutate(signed_logp =  -log10(p.adjust), dir = "up"),
  down_df %>% mutate(signed_logp = log10(p.adjust),   dir = "down")
)

# y 轴顺序(自下而上): 先排下调块(最不显著在下, 最显著贴近中线),
# 再排上调块(最显著置于图最上方), 与参考图一致
ord_up   <- up_df$Description[order(-(-log10(up_df$p.adjust)))]
ord_down <- down_df$Description[order(-(-log10(down_df$p.adjust)))]
term_levels <- c(rev(ord_down), rev(ord_up))

bar_df <- bar_df %>%
  mutate(TermLabel = factor(str_wrap(Description, width = 55),
                            levels = str_wrap(term_levels, width = 55)),
         dir = factor(dir, levels = c("down", "up")))

# x 轴范围与两侧文字标注位置
x_max  <- max(abs(bar_df$signed_logp))
x_lim  <- x_max * 1.35
lab_aki <- x_max * 1.16
lab_ctl <- -x_max * 1.16

# 上/下两个条目块在 y 轴上的位置范围(用于绘制侧边竖线)
n_down  <- length(ord_down)
n_up    <- length(ord_up)
y_down  <- c(0.6, n_down + 0.4)               # 下调块
y_up    <- c(n_down + 0.6, n_down + n_up + 0.4)  # 上调块

p_go <- ggplot(bar_df, aes(x = signed_logp, y = TermLabel, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  # 右侧竖线与竖排文字: 标注 AKI 富集(上调块)
  annotate("segment", x = lab_aki + x_max * 0.045, xend = lab_aki + x_max * 0.045,
           y = y_up[1], yend = y_up[2],
           colour = col_up, linewidth = 0.5) +
  annotate("text", x = lab_aki, y = mean(y_up),
           label = "enriched in\nAKI PTEC", colour = col_up,
           angle = 90, size = 3.4, lineheight = 0.95) +
  # 左侧竖线与竖排文字: 标注 CTRL 富集(下调块)
  annotate("segment", x = lab_ctl - x_max * 0.045, xend = lab_ctl - x_max * 0.045,
           y = y_down[1], yend = y_down[2],
           colour = col_down, linewidth = 0.5) +
  annotate("text", x = lab_ctl, y = mean(y_down),
           label = "enriched in\nControl PTEC", colour = col_down,
           angle = 90, size = 3.4, lineheight = 0.95) +
  scale_fill_manual(values = c(down = col_down, up = col_up), guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = seq(-8, 8, 2)) +
  labs(x = expression(-log[10]~"(P.adjust)"), y = NULL,
       title = "PTEC AKI vs Control geneset enrichment analysis") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(colour = ifelse(
          levels(bar_df$TermLabel) %in% str_wrap(up_df$Description, 55),
          col_up, col_down),
          size = 8.5, face = "italic"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.line.x = element_line(colour = "grey35"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = margin(10, 10, 8, 6))

ggsave(file.path(fig_dir, "GO_BP_PTEC_enrichment_barplot.pdf"), p_go,
       width = 8.6, height = 6.4)

# ============================================================
# 二、GEC 差异基因火山图(标注 top10)
# ============================================================
gec_deg <- read.csv(file.path(tab_dir, "DEG_GEC_AKI_vs_CTRL.csv"),
                    stringsAsFactors = FALSE)

# top10 标注: 排除线粒体(MT-)与核糖体(RPL/RPS)基因, 上/下调各取 top5
df <- gec_deg %>% mutate(neg_log10_p = -log10(p_val + 1e-300))
df_filt <- df %>% filter(!grepl("^MT-|^RPL|^RPS", gene))

top_up <- df_filt %>% filter(direction == "Up in AKI") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>% head(5)
top_down <- df_filt %>% filter(direction == "Down in AKI") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>% head(5)
top10 <- bind_rows(top_up, top_down)

# 若某方向显著基因不足, 从另一方向补足
if (nrow(top10) < 10) {
  extra <- df_filt %>%
    filter(p_val_adj < 0.05, !gene %in% top10$gene) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>% head(10 - nrow(top10))
  top10 <- bind_rows(top10, extra)
}

p_volc <- ggplot(df, aes(x = avg_log2FC, y = neg_log10_p, color = direction)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("Up in AKI" = "#DC0000",
                                "Down in AKI" = "#2166AC",
                                "NS" = "grey70")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top10, aes(label = gene),
                  size = 4, fontface = "italic", color = "black",
                  max.overlaps = Inf, box.padding = 1.0,
                  force = 3, force_pull = 0.5,
                  segment.size = 0.3, seed = 2026) +
  labs(x = expression(log[2]~"Fold Change (AKI vs CTRL)"),
       y = expression(-log[10]~"(p-value)"),
       title = "Differentially Expressed Genes in GEC (AKI vs CTRL)",
       color = "Regulation") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "Volcano_GEC_AKI_vs_CTRL.pdf"), p_volc,
       width = 8, height = 7)

# ============================================================
# 三、PTEC / GEC 各样本细胞比例统计图
# ============================================================
# 读取全量 metadata(第一列为无列名的细胞行名, 需手动补齐列名)
meta_file <- file.path(proj_dir, "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")
header_line <- readLines(gzfile(meta_file), n = 1)
col_names <- c("cell_id", strsplit(gsub('"', '', header_line), "\t")[[1]])
keep_idx <- which(col_names %in% c("cell_id", "specimen", "condition.l1", "subclass.l2"))
meta <- fread(meta_file, sep = "\t", header = FALSE, skip = 1, select = keep_idx)
setnames(meta, col_names[keep_idx])

meta <- meta %>% filter(condition.l1 %in% c("AKI", "Ref")) %>%
  mutate(group = factor(ifelse(condition.l1 == "AKI", "AKI", "CTRL"),
                        levels = c("AKI", "CTRL")))

# 分母为该样本全部细胞数, 分子为目标细胞数
sample_total <- meta %>% group_by(specimen, group) %>%
  summarise(total_cells = n(), .groups = "drop")

calc_prop <- function(ct) {
  sample_total %>%
    left_join(meta %>% filter(subclass.l2 == ct) %>%
                group_by(specimen, group) %>%
                summarise(n_cells = n(), .groups = "drop"),
              by = c("specimen", "group")) %>%
    mutate(n_cells    = ifelse(is.na(n_cells), 0L, n_cells),
           percentage = n_cells / total_cells * 100,
           cell_type  = ct)
}

prop_df <- bind_rows(calc_prop("EC-PTC"), calc_prop("EC-GC"))
write.csv(prop_df, file.path(tab_dir, "Cell_proportion_PTEC_GEC.csv"), row.names = FALSE)

# 绘制单个细胞类型的比例图(散点 + 均值±标准差 + Wilcoxon p 值)
plot_prop <- function(dat, title, file_name) {
  # Wilcoxon 秩和检验(AKI vs CTRL)
  wt <- wilcox.test(percentage ~ group, data = dat, exact = FALSE)
  p_lab <- formatC(wt$p.value, format = "e", digits = 1)

  # 均值与标准差(用于误差棒)
  smry <- dat %>% group_by(group) %>%
    summarise(mean = mean(percentage), sd = sd(percentage), .groups = "drop")

  # 括号标注位置
  y_top  <- max(dat$percentage + smry$sd[1])
  y_brk  <- max(dat$percentage) * 1.16
  y_tick <- max(dat$percentage) * 0.03

  ggplot(dat, aes(x = group, y = percentage)) +
    # 散点: 大圆点 + 深色描边
    geom_point(aes(fill = group), shape = 21, size = 4.2,
               colour = "grey20", stroke = 0.9, alpha = 0.95) +
    # 均值与标准差误差棒
    geom_errorbar(data = smry, aes(x = group, y = mean, ymin = mean - sd, ymax = mean + sd),
                  width = 0.28, colour = "grey35", linewidth = 0.7, inherit.aes = FALSE) +
    geom_segment(data = smry, aes(x = as.numeric(group) - 0.09, xend = as.numeric(group) + 0.09,
                                  y = mean, yend = mean),
                 colour = "grey35", linewidth = 0.9, inherit.aes = FALSE) +
    # 组间比较括号与 p 值
    annotate("segment", x = 1, xend = 2, y = y_brk, yend = y_brk,
             colour = "grey25", linewidth = 0.5) +
    annotate("segment", x = 1, xend = 1, y = y_brk - y_tick, yend = y_brk,
             colour = "grey25", linewidth = 0.5) +
    annotate("segment", x = 2, xend = 2, y = y_brk - y_tick, yend = y_brk,
             colour = "grey25", linewidth = 0.5) +
    annotate("text", x = 1.5, y = y_brk * 1.03, label = p_lab,
             size = 4, colour = "black") +
    scale_fill_manual(values = c(AKI = col_aki, CTRL = col_ctrl), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = "Percentage", title = title) +
    theme_classic(base_size = 12) +
    theme(axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.text   = element_text(colour = "black"),
          axis.ticks  = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 11, face = "plain"),
          plot.margin = margin(10, 14, 8, 8))
}

p_ptec <- plot_prop(prop_df %>% filter(cell_type == "EC-PTC"),
                    "Peritubular Capillary Endothelial Cell (PTEC)",
                    "Cell_proportion_PTEC.pdf")
ggsave(file.path(fig_dir, "Cell_proportion_PTEC.pdf"), p_ptec, width = 3.6, height = 5)

p_gec <- plot_prop(prop_df %>% filter(cell_type == "EC-GC"),
                   "Glomerular Endothelial Cell (GEC)",
                   "Cell_proportion_GEC.pdf")
ggsave(file.path(fig_dir, "Cell_proportion_GEC.pdf"), p_gec, width = 3.6, height = 5)

# 检验结果输出
cat("PTEC 比例 Wilcoxon p =",
    formatC(wilcox.test(percentage ~ group,
                        data = filter(prop_df, cell_type == "EC-PTC"),
                        exact = FALSE)$p.value, format = "e", digits = 3), "\n")
cat("GEC  比例 Wilcoxon p =",
    formatC(wilcox.test(percentage ~ group,
                        data = filter(prop_df, cell_type == "EC-GC"),
                        exact = FALSE)$p.value, format = "e", digits = 3), "\n")

# ---------- 记录运行环境(可重复性) ----------
sink(file.path(proj_dir, "logs", "sessionInfo_06_extra_figures.txt"))
cat("=== 06_extra_figures.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== PTEC GO 图中展示条目 ===\n")
print(bar_df[, c("Description", "signed_logp", "dir")])
cat("\n=== PTEC / GEC 比例 Wilcoxon 检验 ===\n")
print(wilcox.test(percentage ~ group, data = prop_df %>% filter(cell_type == "EC-PTC"), exact = FALSE))
print(wilcox.test(percentage ~ group, data = prop_df %>% filter(cell_type == "EC-GC"),  exact = FALSE))
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 06_extra_figures.R 完成")
