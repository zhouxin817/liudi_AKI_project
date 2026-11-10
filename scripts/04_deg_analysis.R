# ============================================================
# 04_deg_analysis.R
# 功能：对 PTEC / GEC / PT 三个群体分别进行 AKI vs CTRL 差异基因分析
#       (Wilcoxon 秩和检验), 导出差异基因表格,
#       并绘制 PTEC 差异基因火山图(标注 top10 基因)
# 输入：results/processed_data/seu_subset_reclustered.rds
# 输出：results/tables/DEG_*.csv
#       results/figures/Volcano_PTEC_AKI_vs_CTRL.pdf
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

proj_dir <- "/home/zhouxin/liudi_AKI_project"
fig_dir  <- file.path(proj_dir, "results", "figures")
tab_dir  <- file.path(proj_dir, "results", "tables")

# ---------- 1. 读取子集对象 ----------
seu_sub <- readRDS(file.path(proj_dir, "results", "processed_data",
                             "seu_subset_reclustered.rds"))
DefaultAssay(seu_sub) <- "RNA"

# ---------- 2. 差异基因分析函数 ----------
# 对每个细胞群体执行 FindMarkers(AKI vs CTRL, Wilcoxon 秩和检验)
run_deg <- function(obj, cell_type) {
  message(">>> 差异分析: ", cell_type)
  sub <- subset(obj, subset = cell_group == cell_type)
  Idents(sub) <- "group"
  sub$group <- droplevels(sub$group)
  # 不设 logfc/pct 过滤, 保留全部结果便于导出完整表格
  deg <- FindMarkers(sub, ident.1 = "AKI", ident.2 = "CTRL",
                     test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  deg$gene      <- rownames(deg)
  deg$cell_type <- cell_type
  # 按校正后 p 值与 log2FC 标注差异方向
  deg <- deg %>%
    mutate(direction = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Up in AKI",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down in AKI",
      TRUE ~ "NS"
    )) %>%
    arrange(p_val_adj)
  deg
}

# ---------- 3. 分别对 PTEC / GEC / PT 做差异分析 ----------
deg_list <- list()
for (ct in c("PTEC", "GEC", "PT")) {
  deg_list[[ct]] <- run_deg(seu_sub, ct)
  write.csv(deg_list[[ct]],
            file.path(tab_dir, paste0("DEG_", ct, "_AKI_vs_CTRL.csv")),
            row.names = FALSE)
  cat(ct, "显著差异基因数(p_val_adj<0.05 & |log2FC|>0.25):",
      sum(deg_list[[ct]]$direction != "NS"), "\n")
}

# ---------- 4. PTEC 火山图(标注 top10 基因) ----------
plot_volcano <- function(deg_df, title) {
  df <- deg_df %>% mutate(neg_log10_p = -log10(p_val + 1e-300))

  # top10 标注策略: 排除线粒体(MT-)与核糖体(RPL/RPS)基因,
  # 上调与下调方向各取 top5(按校正 p 值排序, 兼顾 |log2FC|)
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
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      head(10 - nrow(top10))
    top10 <- bind_rows(top10, extra)
  }

  ggplot(df, aes(x = avg_log2FC, y = neg_log10_p, color = direction)) +
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
         title = title, color = "Regulation") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_volc <- plot_volcano(deg_list[["PTEC"]],
                       "Differentially Expressed Genes in PTEC (AKI vs CTRL)")
ggsave(file.path(fig_dir, "Volcano_PTEC_AKI_vs_CTRL.pdf"), p_volc,
       width = 8, height = 7)

# ---------- 5. 记录运行环境 ----------
sink(file.path(proj_dir, "logs", "sessionInfo_04_deg_analysis.txt"))
cat("=== 04_deg_analysis.R 运行时间 ===\n")
print(Sys.time())
cat("\n=== 各群体显著差异基因统计 ===\n")
for (ct in names(deg_list)) {
  cat(ct, ": 总基因数", nrow(deg_list[[ct]]),
      "; Up in AKI:", sum(deg_list[[ct]]$direction == "Up in AKI"),
      "; Down in AKI:", sum(deg_list[[ct]]$direction == "Down in AKI"), "\n")
}
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 04_deg_analysis.R 完成")
