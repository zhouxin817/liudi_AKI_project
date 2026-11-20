# ============================================================
# rerun_20261119 / 02_deg_analysis.R
# 功能：四类细胞 AKI vs CTRL 差异基因分析 (Wilcoxon 秩和检验)
# 备注：Seurat 6 + dplyr 存在列名遮蔽, 全部 direction 判定采用 base R
# ============================================================
set.seed(2026)
suppressPackageStartupMessages({library(Seurat);library(dplyr)})

rerun_dir <- "/home/zhouxin/liudi_AKI_project/rerun_20261119"
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")

seu <- readRDS(file.path(rerun_dir, "results", "seu_selected.rds"))
DefaultAssay(seu) <- "RNA"; Idents(seu) <- "group"

cell_groups <- c("PTEC", "GEC", "PT", "Neutrophil")
for (cg in cell_groups) {
  message(">>> ", cg)
  sub <- subset(seu, subset = cell_group == cg)
  sub$group <- droplevels(sub$group)
  Idents(sub) <- "group"

  deg <- FindMarkers(sub, ident.1 = "AKI", ident.2 = "CTRL",
                     test.use = "wilcox",
                     logfc.threshold = 0, min.pct = 0)
  deg$gene <- rownames(deg)
  # base R direction 判定, 避开 dplyr 闭包变量遮蔽
  p_adj   <- deg$p_val_adj
  fc      <- deg$avg_log2FC
  dir_vec <- rep("NS", nrow(deg))
  dir_vec[p_adj < 0.05 & fc >  0.25] <- "Up in AKI"
  dir_vec[p_adj < 0.05 & fc < -0.25] <- "Down in AKI"
  deg$direction <- dir_vec
  deg <- deg[order(deg$p_val_adj), ]

  write.csv(deg, file.path(tab_dir,
            paste0("DEG_", cg, "_AKI_vs_CTRL.csv")), row.names = FALSE)
  cat(cg, ": total", nrow(deg),
      "  Up:", sum(deg$direction == "Up in AKI"),
      "  Down:", sum(deg$direction == "Down in AKI"), "\n")
}

sink(file.path(log_dir, "sessionInfo_02_deg_analysis.txt"))
cat("=== 02_deg_analysis.R ===\n"); print(Sys.time())
for (cg in cell_groups) {
  d <- read.csv(file.path(tab_dir,
               paste0("DEG_", cg, "_AKI_vs_CTRL.csv")),
               stringsAsFactors = FALSE)
  cat(cg, ": total", nrow(d),
      "  Up:", sum(d$direction == "Up in AKI"),
      "  Down:", sum(d$direction == "Down in AKI"), "\n")
}
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()
message(">>> 02_deg_analysis.R 完成")
