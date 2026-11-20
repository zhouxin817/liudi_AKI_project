# ============================================================
# 10_ptec_go_barplot.R
# 功能：将 PTEC 更新后的 GO(BP) 富集结果绘制为双向条形图
#       (x 轴为带方向的 -log10(P.adjust); 左蓝 = Control 富集,
#        右粉 = AKI 富集; 条目文字按方向着色并斜体)
# 输入：results/ptec_reanalysis/tables/PTEC_GO_BP_AKI_vs_CTRL_updated.csv
# 输出：results/ptec_reanalysis/figures/PTEC_GO_BP_barplot_updated.pdf
#       results/ptec_reanalysis/tables/PTEC_GO_BP_barplot_data.csv
# 说明：不修改既有图与表, 仅新增文件
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
in_tab   <- file.path(proj_dir, "results", "ptec_reanalysis", "tables",
                      "PTEC_GO_BP_AKI_vs_CTRL_updated.csv")
fig_dir  <- file.path(proj_dir, "results", "ptec_reanalysis", "figures")
tab_dir  <- file.path(proj_dir, "results", "ptec_reanalysis", "tables")

# 与既有 GO 图一致的配色: 上调(粉) / 下调(蓝)
col_up   <- "#CE6C8E"
col_down <- "#5E8DBC"

go_tab <- read.csv(in_tab, stringsAsFactors = FALSE)

# ---------- 1. 去冗余并挑选条目 ----------
# GO 条目间存在父子/近义关系, 基因重叠(Jaccard)超过阈值视为冗余,
# 保证图中条目互不重复
dedup_terms <- function(df, top_n = 8, jaccard_cut = 0.7) {
  df <- df %>% arrange(p.adjust)
  kept <- character(0)
  kept_genes <- list()
  for (i in seq_len(nrow(df))) {
    g <- strsplit(df$geneID[i], "/")[[1]]
    redundant <- FALSE
    for (kg in kept_genes) {
      inter <- length(intersect(g, kg))
      uni   <- length(union(g, kg))
      if (uni > 0 && inter / uni > jaccard_cut) { redundant <- TRUE; break }
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

# 上调记正值, 下调记负值
bar_df <- bind_rows(
  up_df   %>% mutate(signed_logp =  -log10(p.adjust), dir = "up"),
  down_df %>% mutate(signed_logp = log10(p.adjust),   dir = "down")
)
write.csv(bar_df, file.path(tab_dir, "PTEC_GO_BP_barplot_data.csv"),
          row.names = FALSE)

# y 轴顺序(自下而上): 先下调块(最显著者贴近中线), 再上调块(最显著者置顶)
ord_up   <- up_df$Description[order(log10(up_df$p.adjust))]
ord_down <- down_df$Description[order(log10(down_df$p.adjust))]
term_levels <- c(rev(ord_down), rev(ord_up))

bar_df <- bar_df %>%
  mutate(TermLabel = factor(str_wrap(Description, width = 55),
                            levels = str_wrap(term_levels, width = 55)),
         dir = factor(dir, levels = c("down", "up")))

# ---------- 2. 坐标范围与侧边标注位置 ----------
# 完整展示全部条形的真实长度(不截断), 按数据实际范围等比留出标注空间
x_max  <- max(abs(bar_df$signed_logp))
x_lim  <- x_max * 1.42
lab_up <- x_max * 1.18
lab_dn <- -x_max * 1.18

n_down <- length(ord_down)
n_up   <- length(ord_up)
y_down <- c(0.6, n_down + 0.4)
y_up   <- c(n_down + 0.6, n_down + n_up + 0.4)

# ---------- 3. 绘图 ----------
p <- ggplot(bar_df, aes(x = signed_logp, y = TermLabel, fill = dir)) +
  geom_col(width = 0.68) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  # 右侧: AKI 富集标注
  annotate("segment", x = lab_up + x_max * 0.05, xend = lab_up + x_max * 0.05,
           y = y_up[1], yend = y_up[2], colour = col_up, linewidth = 0.5) +
  annotate("text", x = lab_up, y = mean(y_up),
           label = "enriched in\nAKI PTEC", colour = col_up,
           angle = 90, size = 3.4, lineheight = 0.95) +
  # 左侧: Control 富集标注
  annotate("segment", x = lab_dn - x_max * 0.05, xend = lab_dn - x_max * 0.05,
           y = y_down[1], yend = y_down[2], colour = col_down, linewidth = 0.5) +
  annotate("text", x = lab_dn, y = mean(y_down),
           label = "enriched in\nControl PTEC", colour = col_down,
           angle = 90, size = 3.4, lineheight = 0.95) +
  scale_fill_manual(values = c(down = col_down, up = col_up), guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = pretty(c(-x_lim, x_lim), n = 7)) +
  labs(x = expression(-log[10]~"(P.adjust)"), y = NULL,
       title = "PTEC AKI vs Control geneset enrichment analysis") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(
          colour = ifelse(levels(bar_df$TermLabel) %in%
                            str_wrap(up_df$Description, 55), col_up, col_down),
          size = 8.5, face = "italic"),
        axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(colour = "black"),
        axis.line.x  = element_line(colour = "grey35"),
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin  = margin(10, 10, 8, 6))

ggsave(file.path(fig_dir, "PTEC_GO_BP_barplot_updated.pdf"), p,
       width = 8.6, height = 6.4)

message(">>> 10_ptec_go_barplot.R 完成")
