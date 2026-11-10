# ============================================================
# 08_boxplot_style_stats.R
# 功能：按用户提供的参考图样式(箱线图 + 散点 + 括号 p 值)重新绘制
#       PTEC / GEC 的细胞数量与比例差异
# 说明：不改变任何既有结果, 仅读取 06 生成的逐样本统计表重新绘图
# 输入：results/tables/Cell_proportion_PTEC_GEC.csv
# 输出：results/figures/Cell_count_PTEC_GEC_boxplot.pdf
#       results/figures/Cell_proportion_PTEC_GEC_boxplot.pdf
# ============================================================

# ---------- 0. 环境准备 ----------
set.seed(2026)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
tab_dir  <- file.path(proj_dir, "results", "tables")

# 参考图配色: CTRL 灰, AKI 砖红
col_ctrl <- "#9E9E9E"
col_aki  <- "#C0504D"

# ---------- 1. 读取逐样本统计表 ----------
endo_df <- read.csv(file.path(tab_dir, "Cell_proportion_PTEC_GEC.csv"),
                    stringsAsFactors = FALSE) %>%
  mutate(group     = factor(group, levels = c("CTRL", "AKI")),
         cell_type = recode(cell_type, "EC-PTC" = "PTEC", "EC-GC" = "GEC"),
         cell_type = factor(cell_type, levels = c("PTEC", "GEC")))

# ---------- 2. 参考图样式绘图函数 ----------
# dat      : 含 group / cell_type / value 三列
# y_label  : y 轴标题
# title    : 图标题
# 说明：每个细胞类型单独做 Wilcoxon 秩和检验, 括号标于最高点上方的留白区
plot_box <- function(dat, y_label, title) {
  # 各分面的检验结果与括号位置
  stat_df <- dat %>%
    group_by(cell_type) %>%
    summarise(y_max = max(value),
              p_val = wilcox.test(value ~ group, exact = FALSE)$p.value,
              .groups = "drop") %>%
    mutate(p_label = ifelse(p_val < 0.001,
                            formatC(p_val, format = "e", digits = 1),
                            formatC(p_val, format = "f", digits = 3)),
           y_brk  = y_max * 1.10,
           y_tick = y_max * 0.030,
           y_text = y_max * 1.145)

  ggplot(dat, aes(x = group, y = value)) +
    # 箱线图: 灰色/砖红填充, 深色边框, 隐藏离群点(由散点层展示)
    geom_boxplot(aes(fill = group), width = 0.62, outlier.shape = NA,
                 colour = "grey15", linewidth = 0.7) +
    # 逐样本散点(黑色实心点, 水平抖动)
    geom_point(shape = 16, size = 1.9, colour = "grey10", alpha = 0.9,
               position = position_jitter(width = 0.13, height = 0, seed = 2026)) +
    # 组间比较括号
    geom_segment(data = stat_df,
                 aes(x = 1, xend = 2, y = y_brk, yend = y_brk),
                 colour = "grey15", linewidth = 0.55, inherit.aes = FALSE) +
    geom_segment(data = stat_df,
                 aes(x = 1, xend = 1, y = y_brk - y_tick, yend = y_brk),
                 colour = "grey15", linewidth = 0.55, inherit.aes = FALSE) +
    geom_segment(data = stat_df,
                 aes(x = 2, xend = 2, y = y_brk - y_tick, yend = y_brk),
                 colour = "grey15", linewidth = 0.55, inherit.aes = FALSE) +
    geom_text(data = stat_df, aes(x = 1.5, y = y_text, label = p_label),
              size = 4.2, colour = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = NULL, y = y_label, title = title) +
    facet_wrap(~ cell_type, nrow = 1, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(axis.line    = element_line(colour = "black", linewidth = 0.9),
          axis.ticks   = element_line(colour = "black", linewidth = 0.9),
          axis.text    = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 12),
          strip.background = element_blank(),
          strip.text   = element_text(size = 12, face = "bold"),
          plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
          panel.spacing = unit(1.2, "lines"),
          plot.margin  = margin(10, 12, 8, 8))
}

# ---------- 3. 细胞数量图 ----------
count_df <- endo_df %>%
  transmute(specimen, group, cell_type, value = n_cells)

p_count <- plot_box(count_df, "Cell count", "PTEC and GEC Count per Sample")
ggsave(file.path(fig_dir, "Cell_count_PTEC_GEC_boxplot.pdf"), p_count,
       width = 7.6, height = 5.2)

# ---------- 4. 细胞比例图 ----------
pct_df <- endo_df %>%
  transmute(specimen, group, cell_type, value = percentage)

p_pct <- plot_box(pct_df, "Cell proportion (%)", "PTEC and GEC Proportion per Sample")
ggsave(file.path(fig_dir, "Cell_proportion_PTEC_GEC_boxplot.pdf"), p_pct,
       width = 7.6, height = 5.2)

# ---------- 5. 控制台核对统计量 ----------
cat("=== PTEC / GEC 箱线图统计核对(Wilcoxon, AKI vs CTRL) ===\n")
for (m in c("count", "pct")) {
  d <- if (m == "count") count_df else pct_df
  for (ct in c("PTEC", "GEC")) {
    sub <- d[d$cell_type == ct, ]
    p <- wilcox.test(value ~ group, data = sub, exact = FALSE)$p.value
    cat(sprintf("%-5s %-2s : CTRL mean = %.3f, AKI mean = %.3f, p = %.4g\n",
                ct, m, mean(sub$value[sub$group == "CTRL"]),
                mean(sub$value[sub$group == "AKI"]), p))
  }
}

# ---------- 6. 记录运行环境(可重复性) ----------
sink(file.path(proj_dir, "logs", "sessionInfo_08_boxplot_style_stats.txt"))
cat("=== 08_boxplot_style_stats.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 08_boxplot_style_stats.R 完成")
