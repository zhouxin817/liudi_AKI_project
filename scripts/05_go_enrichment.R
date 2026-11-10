# ============================================================
# 05_go_enrichment.R
# 功能：对 PTEC / GEC / PT 差异基因(AKI vs CTRL)进行 GO(BP)富集分析,
#       导出富集结果表格, 并绘制 GEC 的 GO 富集双向点图
# 输入：results/tables/DEG_*.csv
# 输出：results/tables/GO_BP_*.csv
#       results/figures/GO_BP_GEC_dotplot.pdf
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
tab_dir  <- file.path(proj_dir, "results", "tables")

# ---------- 1. 读取差异基因表格 ----------
deg_files <- c("PTEC" = "DEG_PTEC_AKI_vs_CTRL.csv",
               "GEC"  = "DEG_GEC_AKI_vs_CTRL.csv",
               "PT"   = "DEG_PT_AKI_vs_CTRL.csv")
deg_list <- lapply(deg_files, function(f) {
  read.csv(file.path(tab_dir, f), stringsAsFactors = FALSE)
})

# ---------- 2. GO 富集分析函数 ----------
# 显著差异基因: p_val_adj < 0.05 且 |avg_log2FC| > 0.25
# 上调/下调基因分别富集, 结果合并并保留 direction 列
run_go <- function(deg_df, cell_type) {
  sig_up   <- deg_df$gene[deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC > 0.25]
  sig_down <- deg_df$gene[deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC < -0.25]

  map_up   <- bitr(sig_up,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  map_down <- bitr(sig_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  res_all <- NULL
  for (dir in c("Up_in_AKI", "Down_in_AKI")) {
    ids <- if (dir == "Up_in_AKI") map_up$ENTREZID else map_down$ENTREZID
    if (length(ids) < 10) {
      message(cell_type, " ", dir, " 基因数不足(<10), 跳过富集分析")
      next
    }
    ego <- enrichGO(gene          = ids,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      tmp <- as.data.frame(ego)
      tmp$direction <- dir
      tmp$cell_type <- cell_type
      res_all <- rbind(res_all, tmp)
    }
  }
  res_all
}

# ---------- 3. 对三个群体分别做 GO 富集 ----------
go_results <- list()
for (ct in c("PTEC", "GEC", "PT")) {
  message(">>> GO 富集分析: ", ct)
  go_results[[ct]] <- run_go(deg_list[[ct]], ct)
  if (!is.null(go_results[[ct]]) && nrow(go_results[[ct]]) > 0) {
    write.csv(go_results[[ct]],
              file.path(tab_dir, paste0("GO_BP_", ct, "_AKI_vs_CTRL.csv")),
              row.names = FALSE)
    cat(ct, "富集到显著 GO BP 条目数:", nrow(go_results[[ct]]), "\n")
  } else {
    cat(ct, "未富集到显著 GO 条目\n")
  }
}

# ---------- 4. GEC 的 GO 富集结果图(双向点图) ----------
# 图形样式: y 轴为 GO 条目, x 轴为 Cluster(down / up 两个方向),
#           点大小 = GeneRatio, 点颜色 = p.adjust(越显著越红)
gec_deg <- deg_list[["GEC"]]
n_up    <- sum(gec_deg$p_val_adj < 0.05 & gec_deg$avg_log2FC > 0.25)
n_down  <- sum(gec_deg$p_val_adj < 0.05 & gec_deg$avg_log2FC < -0.25)
cat("GEC 差异基因数: up =", n_up, ", down =", n_down, "\n")

plot_go_dot <- function(tab, n_up, n_down, title, top_n = 10) {
  if (is.null(tab) || nrow(tab) == 0) return(NULL)

  # GeneRatio 由 "a/b" 字符串转换为数值
  df <- tab %>%
    mutate(GeneRatioNum = vapply(strsplit(GeneRatio, "/"),
                                 function(x) as.numeric(x[1]) / as.numeric(x[2]),
                                 numeric(1)),
           Cluster = factor(ifelse(direction == "Up_in_AKI", "up", "down"),
                            levels = c("down", "up")))

  # 分别取各方向 top_n 条目(按 p.adjust), 再取并集
  top_down <- df %>% filter(Cluster == "down") %>% arrange(p.adjust) %>% head(top_n)
  top_up   <- df %>% filter(Cluster == "up")   %>% arrange(p.adjust) %>% head(top_n)
  if (nrow(top_down) == 0 && nrow(top_up) == 0) return(NULL)

  # 条目排序: down 方向条目排图上端, up 方向特有条目排下端
  ord <- bind_rows(
    top_down %>% mutate(src = 1),
    top_up %>% filter(!Description %in% top_down$Description) %>% mutate(src = 2)
  ) %>% arrange(src, p.adjust)
  term_levels <- rev(unique(ord$Description))

  df_plot <- df %>%
    filter(Description %in% term_levels) %>%
    mutate(Description = factor(Description, levels = term_levels),
           # 过长条目自动换行, 避免图左侧被裁切
           TermLabel = factor(str_wrap(as.character(Description), width = 45),
                              levels = str_wrap(term_levels, width = 45)),
           neg_log10_padj = -log10(p.adjust))

  # 颜色断点(以 p.adjust 科学计数法标注), 超出数据范围则剔除
  brks <- c(1e-2, 1e-4, 1e-6, 1e-8)
  brks <- brks[brks <= max(df_plot$p.adjust) & brks >= min(df_plot$p.adjust)]
  if (length(brks) == 0) brks <- signif(max(df_plot$p.adjust), 1)

  ggplot(df_plot, aes(x = Cluster, y = TermLabel)) +
    geom_point(aes(size = GeneRatioNum, color = neg_log10_padj)) +
    scale_x_discrete(labels = c(down = paste0("down\n(", n_down, ")"),
                                up   = paste0("up\n(", n_up, ")"))) +
    scale_color_gradientn(
      colours = c("#2C7BB6", "#8073A8", "#D6604D"),
      limits  = c(min(df_plot$neg_log10_padj), max(df_plot$neg_log10_padj)),
      breaks  = -log10(brks), labels = format(brks, scientific = TRUE),
      name    = "p.adjust", oob = scales::squish) +
    scale_size_continuous(range = c(1, 7), name = "GeneRatio") +
    labs(x = "Cluster", y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_line(colour = "grey90"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey30", fill = NA),
          axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
          axis.text.x = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.key.height = grid::unit(1.1, "cm"))
}

p_go <- plot_go_dot(go_results[["GEC"]], n_up, n_down, "Enriched GO of DEGs (GEC)")
if (!is.null(p_go)) {
  ggsave(file.path(fig_dir, "GO_BP_GEC_dotplot.pdf"), p_go, width = 8, height = 8)
} else {
  message("GEC 无显著富集条目, 跳过绘图")
}

# ---------- 5. 记录运行环境 ----------
sink(file.path(proj_dir, "logs", "sessionInfo_05_go_enrichment.txt"))
cat("=== 05_go_enrichment.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== GO 富集条目统计 ===\n")
for (ct in names(go_results)) {
  n <- ifelse(is.null(go_results[[ct]]), 0, nrow(go_results[[ct]]))
  cat(ct, ":", n, "条显著 GO BP 条目\n")
}
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 05_go_enrichment.R 完成")
