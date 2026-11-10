# ============================================================
# 07_supplementary_visuals.R
# 功能：在不改动既有结果的前提下, 补充可视化与统计汇总
#   1) 中性粒细胞: 各样本数量(count)与比例(percentage)双面板图
#   2) GEC / PTEC : 各样本数量与比例 2x2 分面图
#   3) 上述三群体数量/比例的统计汇总表
# 图形样式: 散点(深描边大圆点) + 均值±标准差 + 组间括号与 p 值
# 输入：results/tables/Neutrophil_per_sample_stats.csv
#       results/tables/Cell_proportion_PTEC_GEC.csv
# 输出：results/figures/Neutrophil_count_proportion.pdf
#       results/figures/Endothelial_count_proportion_stats.pdf
#       results/tables/Supplementary_stats_summary.csv
# ============================================================

# ---------- 0. 环境准备 ----------
set.seed(2026)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
tab_dir  <- file.path(proj_dir, "results", "tables")

# 与既有图保持一致的配色: AKI 鲑鱼粉, CTRL 青绿
col_aki  <- "#E17B6E"
col_ctrl <- "#7FC4C9"

# ---------- 1. 读取既有统计表(不重算, 保证与已有结果一致) ----------
neu_df <- read.csv(file.path(tab_dir, "Neutrophil_per_sample_stats.csv"),
                   stringsAsFactors = FALSE) %>%
  mutate(group = factor(group, levels = c("AKI", "CTRL")))

# 内皮细胞表: 统一细胞类型缩写
endo_df <- read.csv(file.path(tab_dir, "Cell_proportion_PTEC_GEC.csv"),
                    stringsAsFactors = FALSE) %>%
  mutate(group = factor(group, levels = c("AKI", "CTRL")),
         cell_type = recode(cell_type, "EC-PTC" = "PTEC", "EC-GC" = "GEC"))

# ---------- 2. 通用绘图函数 ----------
# dat : 含 group / value / panel 三列(panel 为分面变量)
# 说明: 每个 panel 单独做 Wilcoxon 秩和检验, 并在图上绘制括号与 p 值
plot_scatter <- function(dat, panel_col = "panel", nrow = 1) {
  # 每个 panel 的统计量(均值/标准差/p 值)
  stat_df <- dat %>%
    group_by(.data[[panel_col]]) %>%
    summarise(mean_aki  = mean(value[group == "AKI"]),
              sd_aki    = sd(value[group == "AKI"]),
              mean_ctrl = mean(value[group == "CTRL"]),
              sd_ctrl   = sd(value[group == "CTRL"]),
              y_max     = max(value),
              p_val     = wilcox.test(value ~ group, exact = FALSE)$p.value,
              .groups   = "drop") %>%
    mutate(p_label = formatC(p_val, format = "e", digits = 1),
           y_brk   = y_max * 1.10,
           y_tick  = y_max * 0.025,
           y_text  = y_max * 1.15)

  # 均值与标准差的误差棒数据(宽表转长表)
  err_df <- stat_df %>%
    select(all_of(panel_col), mean_aki, sd_aki, mean_ctrl, sd_ctrl) %>%
    pivot_longer(-all_of(panel_col), names_to = c(".value", "group"),
                 names_pattern = "(mean|sd)_(aki|ctrl)") %>%
    mutate(group = factor(ifelse(group == "aki", "AKI", "CTRL"),
                          levels = c("AKI", "CTRL")))

  p <- ggplot(dat, aes(x = group, y = value)) +
    # 样本散点: 大圆点 + 深色描边; 轻微水平抖动避免重叠点完全遮盖
    geom_point(aes(fill = group), shape = 21, size = 3.6,
               colour = "grey20", stroke = 0.8, alpha = 0.95,
               position = position_jitter(width = 0.06, height = 0, seed = 2026)) +
    # 均值 ± 标准差
    geom_errorbar(data = err_df,
                  aes(x = group, y = mean, ymin = mean - sd, ymax = mean + sd),
                  width = 0.26, colour = "grey35", linewidth = 0.6,
                  inherit.aes = FALSE) +
    geom_segment(data = err_df,
                 aes(x = as.numeric(group) - 0.08, xend = as.numeric(group) + 0.08,
                     y = mean, yend = mean),
                 colour = "grey35", linewidth = 0.9, inherit.aes = FALSE) +
    # 组间比较括号与 p 值
    geom_segment(data = stat_df,
                 aes(x = 1, xend = 2, y = y_brk, yend = y_brk),
                 colour = "grey25", linewidth = 0.45, inherit.aes = FALSE) +
    geom_segment(data = stat_df,
                 aes(x = 1, xend = 1, y = y_brk - y_tick, yend = y_brk),
                 colour = "grey25", linewidth = 0.45, inherit.aes = FALSE) +
    geom_segment(data = stat_df,
                 aes(x = 2, xend = 2, y = y_brk - y_tick, yend = y_brk),
                 colour = "grey25", linewidth = 0.45, inherit.aes = FALSE) +
    geom_text(data = stat_df, aes(x = 1.5, y = y_text, label = p_label),
              size = 3.4, colour = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = c(AKI = col_aki, CTRL = col_ctrl), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.14))) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.text   = element_text(colour = "black"),
          axis.ticks  = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text  = element_text(size = 10.5, colour = "black"),
          panel.spacing = unit(0.9, "lines"),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.margin = margin(8, 10, 6, 6)) +
    facet_wrap(as.formula(paste("~", panel_col)), nrow = nrow, scales = "free_y")

  p
}

# ---------- 3. 中性粒细胞: 数量 + 比例 双面板 ----------
neu_long <- bind_rows(
  neu_df %>% transmute(specimen, group, panel = "Count (cells)", value = neutrophils),
  neu_df %>% transmute(specimen, group, panel = "Percentage (%)", value = proportion)
) %>%
  mutate(panel = factor(panel, levels = c("Count (cells)", "Percentage (%)")))

p_neu <- plot_scatter(neu_long) +
  labs(title = "Neutrophil Count and Proportion per Sample") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))

ggsave(file.path(fig_dir, "Neutrophil_count_proportion.pdf"), p_neu,
       width = 7.2, height = 4.8)

# ---------- 4. GEC / PTEC: 数量 + 比例 2x2 分面 ----------
endo_long <- bind_rows(
  endo_df %>% transmute(specimen, group, cell_type, metric = "Count (cells)", value = n_cells),
  endo_df %>% transmute(specimen, group, cell_type, metric = "Percentage (%)", value = percentage)
) %>%
  mutate(cell_type = factor(cell_type, levels = c("PTEC", "GEC")),
         metric    = factor(metric, levels = c("Count (cells)", "Percentage (%)")),
         panel     = factor(paste(cell_type, metric, sep = " | "),
                            levels = c("PTEC | Count (cells)", "PTEC | Percentage (%)",
                                       "GEC | Count (cells)",  "GEC | Percentage (%)")))

# 逐面板统计量(用于导出汇总表)
endo_stat <- endo_long %>%
  group_by(cell_type, metric) %>%
  summarise(n_AKI     = sum(group == "AKI"),
            n_CTRL    = sum(group == "CTRL"),
            median_AKI  = median(value[group == "AKI"]),
            median_CTRL = median(value[group == "CTRL"]),
            mean_AKI    = mean(value[group == "AKI"]),
            mean_CTRL   = mean(value[group == "CTRL"]),
            sd_AKI      = sd(value[group == "AKI"]),
            sd_CTRL     = sd(value[group == "CTRL"]),
            p_value   = wilcox.test(value ~ group, exact = FALSE)$p.value,
            .groups   = "drop")

p_endo <- plot_scatter(endo_long, panel_col = "panel", nrow = 2) +
  labs(title = "PTEC and GEC Count and Proportion per Sample") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))

ggsave(file.path(fig_dir, "Endothelial_count_proportion_stats.pdf"), p_endo,
       width = 7.2, height = 7.2)

# ---------- 5. 统计汇总表 ----------
# 中性粒细胞两部分(数量/比例)的统计
neu_stat <- neu_long %>%
  group_by(panel) %>%
  summarise(n_AKI       = sum(group == "AKI"),
            n_CTRL      = sum(group == "CTRL"),
            median_AKI  = median(value[group == "AKI"]),
            median_CTRL = median(value[group == "CTRL"]),
            mean_AKI    = mean(value[group == "AKI"]),
            mean_CTRL   = mean(value[group == "CTRL"]),
            sd_AKI      = sd(value[group == "AKI"]),
            sd_CTRL     = sd(value[group == "CTRL"]),
            p_value     = wilcox.test(value ~ group, exact = FALSE)$p.value,
            .groups     = "drop") %>%
  mutate(population = "Neutrophil", metric = as.character(panel))

summary_df <- bind_rows(
  neu_stat,
  endo_stat %>% mutate(population = as.character(cell_type),
                       metric = as.character(metric))
) %>%
  select(population, metric, n_AKI, n_CTRL, mean_AKI, mean_CTRL,
         median_AKI, median_CTRL, sd_AKI, sd_CTRL, p_value)

write.csv(summary_df, file.path(tab_dir, "Supplementary_stats_summary.csv"),
          row.names = FALSE)

cat("=== 补充分析统计汇总 ===\n")
print(as.data.frame(summary_df), digits = 4)

# ---------- 6. 记录运行环境(可重复性) ----------
sink(file.path(proj_dir, "logs", "sessionInfo_07_supplementary_visuals.txt"))
cat("=== 07_supplementary_visuals.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 统计汇总 ===\n")
print(summary_df)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

message(">>> 07_supplementary_visuals.R 完成")
