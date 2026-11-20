# ============================================================
# rerun_20261119 / 05_mitochondrial_visualizations.R
# 功能: 解释 PTEC GO 条形图中两个 ATP 合成词条方向相反的现象,
#       并补充更多线粒体调控结果的可视化。
#
# 核心现象:
#   GO:0015986 "proton motive force-driven ATP synthesis"  -> 富集于 AKI(上调)
#   GO:0042776 "proton motive force-driven mitochondrial ATP synthesis" -> 富集于 CTRL(下调)
# 根因: 两个词条虽然功能相关, 但在 PTEC DEG 中命中的是完全不同的基因子集:
#   - 上调(AKI)侧: 几乎全是核编码 ATP 合酶 F1/FO 亚基(ATP5F1A/B/C/D/E,
#     ATP5PB/MC1/2/3, ATP5ME/MF/MG/PF/PO, ATP5IF1) -> 核基因组代偿性升高
#   - 下调(AKI)侧: 几乎全是线粒体基因组编码的呼吸链亚基(MT-ND1/2/3/4L/5/6,
#     MT-CO2/3, MT-CYB, MT-ATP8) -> mtDNA 编码的 OXPHOS 受抑
#   即 PTEC 在 AKI 中呈现 "核编码 ATP 合酶升高 / mtDNA 呼吸链降低" 的
#   核-线粒体解耦(nuclear-mitochondrial decoupling)。
#
# 本脚本输出(均存 rerun_20261119/results/):
#   tables/PTEC_ATP_synthesis_gene_breakdown.csv   ATP 合成基因拆解表(证据)
#   tables/PTEC_mito_module_scores.csv             各细胞线粒体模块评分
#   tables/PTEC_mito_top50_heatmap_genes.csv       热图基因清单
#   tables/PTEC_mito_top50_pseudobulk_zscore.csv   热图 z-score 矩阵
#   figures/PTEC_ATP_nuclear_vs_mt_violin.pdf      图1: 核编码 vs mt 编码解耦小提琴图
#   figures/PTEC_mito_module_scores_boxplot.pdf    图2: 线粒体功能模块评分箱线图
#   figures/PTEC_mito_top50_heatmap.pdf            图3: 线粒体 DEG top50 热图(按模块分组)
#   figures/PTEC_GO_BP_mitochondrial_barplot.pdf   图4: 线粒体专用 GO 双向条形图
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
})

proj_dir  <- "/home/zhouxin/liudi_AKI_project"
rerun_dir <- file.path(proj_dir, "rerun_20261119_noSample")
fig_dir   <- file.path(rerun_dir, "results", "figures")
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# 与既有图一致的配色
col_up   <- "#CE6C8E"   # Up in AKI
col_down <- "#5E8DBC"   # Down in AKI
col_aki  <- "#E17B6E"   # AKI
col_ctrl <- "#7FC4C9"   # CTRL

# ---------- 1. 加载数据 ----------
seu <- readRDS(file.path(rerun_dir, "results", "seu_selected.rds"))
DefaultAssay(seu) <- "RNA"
ptec <- subset(seu, subset = cell_group == "PTEC")
ptec$group <- factor(ptec$group, levels = c("CTRL", "AKI"))
Idents(ptec) <- "group"
cat("PTEC 细胞数:", ncol(ptec),
    "; CTRL:", sum(ptec$group == "CTRL"),
    "; AKI:", sum(ptec$group == "AKI"), "\n")

deg    <- read.csv(file.path(tab_dir, "DEG_PTEC_AKI_vs_CTRL.csv"),
                   stringsAsFactors = FALSE)
go_tab <- read.csv(file.path(tab_dir, "GO_BP_PTEC_AKI_vs_CTRL.csv"),
                   stringsAsFactors = FALSE)

# Seurat 6 兼容取表达矩阵
get_data_layer <- function(obj, assay = "RNA") {
  tryCatch(GetAssayData(obj, assay = assay, layer = "data"),
           error = function(e) GetAssayData(obj, assay = assay, slot = "data"))
}
data_mat <- get_data_layer(ptec)

# ---------- 2. 辅助函数 ----------
# GO ID(可含 CC/BP) -> gene symbol, 经 GOALL 传递注释
# 注意: 部分基础 GO ID(如某些 CC 词条)不是 org.Hs.eg.db 的有效 GOALL 键,
#       但其子孙条目通常有效; 这里先过滤无效键, 靠 GOALL 传递注释从子孙恢复基因。
go_genes <- function(ids) {
  ids <- unique(ids[!is.na(ids)])
  valid_keys <- keys(org.Hs.eg.db, keytype = "GOALL")
  ids <- intersect(ids, valid_keys)
  if (length(ids) == 0) return(character(0))
  suppressWarnings(
    AnnotationDbi::select(org.Hs.eg.db, keys = ids,
                          columns = "SYMBOL", keytype = "GOALL")$SYMBOL
  ) %>% na.omit() %>% unique() %>% as.character()
}
# GO 子孙
get_offspring <- function(ids, bimap) {
  lst <- as.list(bimap)
  unique(unlist(lst[ids], use.names = FALSE))
}
# 基因是否线粒体基因组编码
is_mt_encoded <- function(g) grepl("^MT-", g)

# ============================================================
# 3. ATP 合成基因拆解表(回答用户疑问的证据)
# ============================================================
# 从 GO 富集表中抓取 GO:0015986 与 GO:0042776 在 Up/Down 方向的基因
atp_go_ids <- c("GO:0015986", "GO:0042776")
atp_rows <- go_tab %>% filter(ID %in% atp_go_ids) %>%
  select(ID, Description, direction, p.adjust, geneID, Count)

breakdown <- bind_rows(lapply(seq_len(nrow(atp_rows)), function(i) {
  g <- strsplit(atp_rows$geneID[i], "/")[[1]]
  data.frame(
    GO_ID        = atp_rows$ID[i],
    GO_Term      = atp_rows$Description[i],
    direction    = atp_rows$direction[i],
    go_p_adj     = atp_rows$p.adjust[i],
    gene         = g,
    stringsAsFactors = FALSE)
}))
# 合并 DEG 统计 + 编码来源标注
breakdown <- breakdown %>%
  left_join(deg %>% select(gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2),
            by = "gene") %>%
  mutate(encoding = ifelse(is_mt_encoded(gene), "Mitochondrial genome",
                           "Nuclear genome"),
         deg_significant = ifelse(p_val_adj < 0.05, "yes", "no")) %>%
  arrange(GO_ID, direction, encoding, desc(abs(avg_log2FC)))
write.csv(breakdown,
          file.path(tab_dir, "PTEC_ATP_synthesis_gene_breakdown.csv"),
          row.names = FALSE)

cat("\n=== ATP 合成词条基因拆解 ===\n")
print(breakdown %>% count(GO_ID, GO_Term, direction, encoding))

# ============================================================
# 图 1: 核编码 ATP 合酶 vs 线粒体编码 OXPHOS 的解耦小提琴图
# ============================================================
# 核编码 ATP 合酶亚基(F1 + FO + 附属)
nuclear_atp <- c("ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
                 "ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3",
                 "ATP5ME","ATP5MF","ATP5MG","ATP5PF","ATP5PO","ATP5IF1")
# 线粒体基因组编码 OXPHOS 亚基(复合体 I/III/IV/V)
mt_oxphos <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
               "MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ATP6","MT-ATP8")
nuclear_atp <- intersect(nuclear_atp, rownames(data_mat))
mt_oxphos   <- intersect(mt_oxphos,   rownames(data_mat))
cat("\n核编码 ATP 合酶基因(在数据中):", length(nuclear_atp), "\n")
cat("线粒体编码 OXPHOS 基因(在数据中):", length(mt_oxphos), "\n")

# 用 AddModuleScore 给每个细胞打分(即基因集平均表达, 经对照基因归一)
ptec <- AddModuleScore(ptec, features = list(nuclear_atp), name = "NucATP")
ptec <- AddModuleScore(ptec, features = list(mt_oxphos),   name = "MtOxphos")
# AddModuleScore 自动在列名尾追加序号
col_nuc <- grep("^NucATP", colnames(ptec@meta.data), value = TRUE)[1]
col_mt  <- grep("^MtOxphos", colnames(ptec@meta.data), value = TRUE)[1]

sc <- ptec@meta.data %>%
  select(group, all_of(c(col_nuc, col_mt))) %>%
  rename(Nuclear_ATP_synthase = all_of(col_nuc),
         Mt_encoded_OXPHOS   = all_of(col_mt)) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = c(Nuclear_ATP_synthase, Mt_encoded_OXPHOS),
               names_to = "module", values_to = "score") %>%
  mutate(module = factor(module,
            levels = c("Nuclear_ATP_synthase", "Mt_encoded_OXPHOS"),
            labels = c("Nuclear-encoded ATP synthase",
                       "Mitochondria-encoded OXPHOS")))

# 每个模块的 Wilcoxon 检验
wilcox_p <- sc %>%
  group_by(module) %>%
  summarise(p = wilcox.test(score ~ group, exact = FALSE)$p.value,
            .groups = "drop")
wilcox_p$lab <- ifelse(wilcox_p$p < 1e-4,
                       paste0("P = ", formatC(wilcox_p$p, format = "e", digits = 2)),
                       paste0("P = ", formatC(wilcox_p$p, format = "f", digits = 3)))
print(wilcox_p)

# 每个分面各自的数据范围(用于放置显著性括号), 避免 free_y 下轴范围失衡
y_stats <- sc %>% group_by(module) %>%
  summarise(y_max = max(score, na.rm = TRUE),
            y_min = min(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_rng = y_max - y_min)
wilcox_p <- wilcox_p %>% left_join(y_stats, by = "module")

p1 <- ggplot(sc, aes(x = group, y = score, fill = group)) +
  geom_violin(scale = "width", width = 0.85, colour = "grey35", linewidth = 0.5,
              alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "grey25",
               linewidth = 0.5, fill = "white", alpha = 0.9) +
  facet_wrap(~module, scales = "free_y", nrow = 1) +
  geom_segment(data = wilcox_p,
               aes(x = 1, xend = 2,
                   y = y_max + y_rng * 0.12, yend = y_max + y_rng * 0.12),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(data = wilcox_p,
               aes(x = 1, xend = 1,
                   y = y_max + y_rng * 0.07, yend = y_max + y_rng * 0.12),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(data = wilcox_p,
               aes(x = 2, xend = 2,
                   y = y_max + y_rng * 0.07, yend = y_max + y_rng * 0.12),
               colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = wilcox_p,
            aes(x = 1.5, y = y_max + y_rng * 0.20, label = lab),
            size = 3.6, colour = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), guide = "none") +
  labs(x = NULL, y = "Module score (mean expression)",
       title = "Nuclear-mitochondrial decoupling of ATP synthesis in PTEC (AKI vs CTRL)",
       subtitle = "Nuclear-encoded ATP synthase subunits are induced, whereas mitochondria-encoded OXPHOS subunits are repressed in AKI") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10.5),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.8, colour = "grey30"),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.text = element_text(colour = "black"))

ggsave(file.path(fig_dir, "PTEC_ATP_nuclear_vs_mt_violin.pdf"),
       p1, width = 8.6, height = 5.0)

# ============================================================
# 图 2: 多个线粒体功能模块评分(AKI vs CTRL 箱线图)
# ============================================================
# 用 GO 定义各功能模块基因集, 再与数据取交集
module_defs <- list(
  "ATP synthase (Complex V)" =
    c("GO:0042776", get_offspring("GO:0042776", GOBPOFFSPRING),
      "GO:0005753", get_offspring("GO:0005753", GOCCOFFSPRING)),
  "Complex I" =
    c("GO:0005747", get_offspring("GO:0005747", GOCCOFFSPRING)),
  "Complex III" =
    c("GO:0005750", get_offspring("GO:0005750", GOCCOFFSPRING)),
  "Complex IV" =
    c("GO:0005751", get_offspring("GO:0005751", GOCCOFFSPRING)),
  "Mitochondrial translation" =
    c("GO:0032543", get_offspring("GO:0032543", GOBPOFFSPRING)),
  "Mitochondrial fusion" =
    c("GO:0008053", get_offspring("GO:0008053", GOBPOFFSPRING)),
  "Mitochondrial fission" =
    c("GO:0000266", get_offspring("GO:0000266", GOBPOFFSPRING)),
  "Mitophagy" =
    c("GO:0000423", get_offspring("GO:0000423", GOBPOFFSPRING)),
  "Response to oxidative stress" =
    c("GO:0006979", get_offspring("GO:0006979", GOBPOFFSPRING))
)
# curated 标志基因(数据驱动, 按基因名模式从表达矩阵抓取), 作为 GO 注释的补充:
# 部分基础 CC 词条(如呼吸链复合体)未被 org.Hs.eg.db 直接注释, 靠此保证核心模块不丢失
all_genes <- rownames(data_mat)
curated_genes <- list(
  "ATP synthase (Complex V)" =
    c(grep("^ATP5|^ATP6V|^ATP5I|^ATP5J", all_genes, value = TRUE),
      grep("^MT-ATP", all_genes, value = TRUE)),
  "Complex I"  = c(grep("^NDUF", all_genes, value = TRUE),
                   grep("^MT-ND", all_genes, value = TRUE)),
  "Complex III"= c(grep("^UQCR|^CYC1", all_genes, value = TRUE),
                   grep("^MT-CYB", all_genes, value = TRUE)),
  "Complex IV" = c(grep("^COX", all_genes, value = TRUE),
                   grep("^MT-CO", all_genes, value = TRUE)),
  "Mitochondrial translation" = grep("^MRP|^MTERF", all_genes, value = TRUE),
  "Mitochondrial fusion"  = intersect(c("MFN1","MFN2","OPA1"), all_genes),
  "Mitochondrial fission" = intersect(c("DNM1L","FIS1","MFF"), all_genes),
  "Mitophagy" = intersect(c("PINK1","PRKN","OPTN","CALCOCO2","NBR1",
                            "SQSTM1","BNIP3","BNIP3L","FUNDC1"), all_genes),
  "Response to oxidative stress" = character(0)
)
module_genes <- lapply(names(module_defs), function(nm) {
  unique(c(go_genes(module_defs[[nm]]), curated_genes[[nm]]))
})
names(module_genes) <- names(module_defs)
module_genes <- lapply(module_genes, function(g) intersect(g, all_genes))
module_genes <- module_genes[lengths(module_genes) >= 5]
cat("\n各线粒体模块基因数:\n"); print(lengths(module_genes))

# 批量 AddModuleScore
score_names <- paste0("M", seq_along(module_genes), "_")
for (i in seq_along(module_genes)) {
  ptec <- AddModuleScore(ptec, features = list(module_genes[[i]]),
                         name = score_names[i])
}
score_cols <- paste0(score_names, "1")
score_cols <- intersect(score_cols, colnames(ptec@meta.data))

sc_mod <- ptec@meta.data %>%
  select(group, all_of(score_cols))
colnames(sc_mod) <- c("group", names(module_genes))
sc_mod$cell <- rownames(sc_mod)
sc_long <- sc_mod %>%
  pivot_longer(cols = -c(group, cell), names_to = "module", values_to = "score") %>%
  mutate(module = factor(module, levels = names(module_genes)))

# 各模块 Wilcoxon
mod_p <- sc_long %>%
  group_by(module) %>%
  summarise(p = wilcox.test(score ~ group, exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(lab = ifelse(p < 1e-4,
                      paste0("P=", formatC(p, format = "e", digits = 1)),
                      paste0("P=", formatC(p, format = "f", digits = 2))),
         sig = ifelse(p < 0.001, "***",
               ifelse(p < 0.01, "**",
               ifelse(p < 0.05, "*", "ns"))))

# 保存模块评分表(便于后续复核)
write.csv(sc_mod,
          file.path(tab_dir, "PTEC_mito_module_scores.csv"),
          row.names = FALSE)
write.csv(mod_p,
          file.path(tab_dir, "PTEC_mito_module_Wilcoxon.csv"),
          row.names = FALSE)
cat("\n各线粒体模块 Wilcoxon P:\n"); print(mod_p)

# 箱线图: 每模块一个箱线, AKI vs CTRL 并列, 标注显著性
y_top_mod <- max(sc_long$score, na.rm = TRUE)
y_span_mod <- diff(range(sc_long$score, na.rm = TRUE))

p2 <- ggplot(sc_long, aes(x = module, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, position = position_dodge(0.7),
               colour = "grey30", linewidth = 0.45, alpha = 0.9) +
  geom_point(aes(fill = group), shape = 21, size = 0.9, stroke = 0.25,
             colour = "grey25", alpha = 0.35,
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.7)) +
  geom_text(data = mod_p,
            aes(x = module, y = y_top_mod + y_span_mod * 0.08, label = sig),
            size = 3.8, colour = "grey15", inherit.aes = FALSE) +
  scale_fill_manual(values = c(CTRL = col_ctrl, AKI = col_aki), name = NULL) +
  labs(x = NULL, y = "Module score",
       title = "Mitochondrial functional module activity in PTEC (AKI vs CTRL)",
       subtitle = "Module score = mean expression of module genes (AddModuleScore)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = "black",
                                   size = 9),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8.8, colour = "grey30"),
        plot.margin = margin(10, 12, 8, 32),
        legend.position = "top")

ggsave(file.path(fig_dir, "PTEC_mito_module_scores_boxplot.pdf"),
       p2, width = 9.2, height = 5.6)

# ============================================================
# 图 3: 线粒体 DEG top50 热图(按功能模块分组行)
# ============================================================
# 线粒体基因全集(GO:0005739 + 子孙, GOALL)
mito_cc_ids <- unique(c("GO:0005739",
                        get_offspring("GO:0005739", GOCCOFFSPRING)))
# 基因池 = GO:0005739 注释 ∪ curated 模块基因(弥补 GO 注释缺口, 如 MRP63)
mito_symbols <- unique(c(go_genes(mito_cc_ids), unlist(curated_genes)))
# 热图基因池: 取线粒体基因中 p_val_adj < 0.05 者(不强制 |log2FC|>0.25),
# 这样融合/分裂/线粒体自噬等模块中变化幅度较小但方向协调的领头基因也能入图,
# 更全面地呈现线粒体调控; 行内再按模块定额选取最显著者。
mito_deg <- deg %>%
  filter(gene %in% mito_symbols, gene %in% rownames(data_mat),
         p_val_adj < 0.05)

# 基因 -> 功能模块分类(用于热图行分组)
classify_module <- function(g) {
  case_when(
    grepl("^MT-ND|^NDUF", g)   ~ "Complex I",
    grepl("^MT-CO|^COX|^COA", g) ~ "Complex IV",
    grepl("^MT-CYB|^UQCR|^CYC1|^UQCRC", g) ~ "Complex III",
    grepl("^ATP5|^MT-ATP|^ATP6V|^ATP5I|^ATP5J", g) ~ "ATP synthase",
    grepl("^MRPS|^MRPL|^MRP|^MTERF", g) ~ "Mito translation",
    grepl("^MFN|^OPA1", g) ~ "Fusion",
    grepl("^DNM1L|^FIS1|^MFF", g) ~ "Fission",
    grepl("^PINK1|^PRKN|^OPTN|^NDP52|^CALCOCO|^BNIP|^FUNDC", g) ~ "Mitophagy",
    TRUE ~ "Other mitochondrial")
}
mito_deg$module <- vapply(mito_deg$gene, classify_module, character(1))

# 选取 top50: 先给每个模块定额名额(保证各功能模块都有代表),
# 剩余名额再按全局 p_val_adj 从其余线粒体 DEG 中补足
set.seed(2026)
modules_order <- c("ATP synthase","Complex I","Complex III","Complex IV",
                   "Mito translation","Fusion","Fission","Mitophagy",
                   "Other mitochondrial")
quota <- 6
selected <- character(0)
for (m in modules_order) {
  sub <- mito_deg %>% filter(module == m) %>% arrange(p_val_adj)
  take <- head(sub$gene, quota)
  selected <- c(selected, take)
}
# 补足到 50: 从剩余线粒体 DEG 按 p_val_adj 取
rest <- mito_deg %>%
  filter(!gene %in% selected) %>%
  arrange(p_val_adj) %>%
  head(max(0, 50 - length(selected)))
selected <- c(selected, rest$gene)
top50 <- mito_deg %>% filter(gene %in% selected) %>%
  distinct(gene, .keep_all = TRUE)
# 行顺序: 先按模块顺序, 模块内按 log2FC 降序
top50$module <- factor(top50$module, levels = modules_order)
top50 <- top50 %>% arrange(module, desc(avg_log2FC))
genes_hm <- top50$gene
cat("\n热图基因数:", length(genes_hm),
    "; 模块分布:\n"); print(table(top50$module))

write.csv(top50 %>% select(gene, module, avg_log2FC, p_val, p_val_adj),
          file.path(tab_dir, "PTEC_mito_top50_heatmap_genes.csv"),
          row.names = FALSE)

# 伪批量: 每样本(细胞数 >= 10)基因均值
meta_hm <- ptec@meta.data
meta_hm$specimen <- as.character(meta_hm$specimen)
keep_samples <- names(which(table(meta_hm$specimen) >= 10))
meta_keep <- meta_hm %>% filter(specimen %in% keep_samples)
pb <- sapply(keep_samples, function(s) {
  cells <- rownames(meta_keep)[meta_keep$specimen == s]
  Matrix::rowMeans(data_mat[genes_hm, cells, drop = FALSE])
})
sample_group <- meta_keep %>%
  distinct(specimen, group) %>%
  mutate(group = as.character(group)) %>%
  arrange(group, specimen)
pb <- pb[, sample_group$specimen, drop = FALSE]

# z-score by gene
z <- t(scale(t(pb)))
z[is.na(z)] <- 0
z <- pmax(pmin(z, 2), -2)
write.csv(data.frame(gene = rownames(z), module = top50$module, z,
                     check.names = FALSE),
          file.path(tab_dir, "PTEC_mito_top50_pseudobulk_zscore.csv"),
          row.names = FALSE)

row_mod <- factor(top50$module, levels = modules_order)
col_fun <- colorRamp2(c(-2, 0, 2), c(col_down, "white", col_up))
# 模块配色
mod_pal <- c("ATP synthase" = "#C77CFF","Complex I" = "#7CAE00",
             "Complex III" = "#00BFC4","Complex IV" = "#F8766D",
             "Mito translation" = "#E76BF3","Fusion" = "#619CFF",
             "Fission" = "#FF7C00","Mitophagy" = "#00BA38",
             "Other mitochondrial" = "grey60")

ha <- HeatmapAnnotation(
  Group = sample_group$group,
  col = list(Group = c(CTRL = col_ctrl, AKI = col_aki)),
  annotation_name_side = "left",
  annotation_legend_param = list(Group = list(title = "Group"))
)
ra <- rowAnnotation(
  Module = row_mod,
  col = list(Module = mod_pal),
  annotation_name_side = "top",
  annotation_legend_param = list(Module = list(title = "Module", ncol = 2))
)

ht <- Heatmap(
  z, name = "z-score", col = col_fun,
  row_labels = rownames(z),
  row_names_gp = gpar(fontface = "italic", fontsize = 7),
  row_names_side = "right",
  row_split = row_mod,
  row_title_gp = gpar(fontsize = 7), row_title_rot = 0,
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = factor(sample_group$group, levels = c("CTRL", "AKI")),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 90, column_names_gp = gpar(fontsize = 6.5),
  top_annotation = ha, left_annotation = ra,
  column_gap = unit(1.2, "mm"), row_gap = unit(1.0, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Scaled mean\nexpression (z)")
)

pdf(file.path(fig_dir, "PTEC_mito_top50_heatmap.pdf"),
    width = 9.5, height = 10.2)
draw(ht,
     column_title = "Mitochondria-related TOP50 DEGs in PTEC (AKI vs CTRL)",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(14, 2, 2, 2), "mm"))
invisible(dev.off())

# ============================================================
# 图 4: 线粒体专用 GO(BP) 双向条形图
# ============================================================
mito_bp_ids <- unique(get_offspring("GO:0005739", GOBPOFFSPRING))
is_mito_term <- go_tab$ID %in% mito_bp_ids |
  grepl("mitochond|oxidative phosphorylation|electron transport|oxphos|ATP synth|respiratory chain",
        go_tab$Description, ignore.case = TRUE)
go_mito <- go_tab %>% filter(is_mito_term)
cat("\n线粒体相关 GO(BP) 条目(双向):", nrow(go_mito), "\n")

dedup_terms <- function(df, top_n = 8, jaccard_cut = 0.7) {
  df <- df[order(df$p.adjust), ]
  kept <- character(0); kept_genes <- list()
  for (i in seq_len(nrow(df))) {
    g <- strsplit(df$geneID[i], "/")[[1]]
    redundant <- FALSE
    for (kg in kept_genes) {
      inter <- length(intersect(g, kg)); uni <- length(union(g, kg))
      if (uni > 0 && inter / uni > jaccard_cut) { redundant <- TRUE; break }
    }
    if (!redundant) {
      kept <- c(kept, df$Description[i])
      kept_genes[[length(kept_genes) + 1]] <- g
    }
    if (length(kept) >= top_n) break
  }
  df[df$Description %in% kept, ]
}

up_df   <- dedup_terms(go_mito[go_mito$direction == "Up_in_AKI",   ], top_n = 6)
down_df <- dedup_terms(go_mito[go_mito$direction == "Down_in_AKI", ], top_n = 8)

bar_df <- bind_rows(
  up_df   %>% mutate(signed_logp = -log10(p.adjust), dir = "up"),
  down_df %>% mutate(signed_logp =  log10(p.adjust), dir = "down"))
ord_up   <- up_df$Description[order(-(-log10(up_df$p.adjust)))]
ord_down <- down_df$Description[order(-(-log10(down_df$p.adjust)))]
term_levels <- c(rev(ord_down), rev(ord_up))
bar_df$TermLabel <- factor(str_wrap(bar_df$Description, width = 55),
                           levels = str_wrap(term_levels, width = 55))
bar_df$dir <- factor(bar_df$dir, levels = c("down", "up"))

x_max  <- max(abs(bar_df$signed_logp))
x_lim  <- x_max * 1.38
lab_up <- x_max * 1.20
lab_dn <- -x_max * 1.20
n_down <- length(ord_down); n_up <- length(ord_up)
y_down <- c(0.6, n_down + 0.4)
y_up   <- c(n_down + 0.6, n_down + n_up + 0.4)

p4 <- ggplot(bar_df, aes(x = signed_logp, y = TermLabel, fill = dir)) +
  geom_col(width = 0.68) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  annotate("segment", x = lab_up + x_max * 0.05, xend = lab_up + x_max * 0.05,
           y = y_up[1], yend = y_up[2], colour = col_up, linewidth = 0.5) +
  annotate("text", x = lab_up, y = mean(y_up), label = "enriched in\nAKI PTEC",
           colour = col_up, angle = 90, size = 3.4, lineheight = 0.95) +
  annotate("segment", x = lab_dn - x_max * 0.05, xend = lab_dn - x_max * 0.05,
           y = y_down[1], yend = y_down[2], colour = col_down, linewidth = 0.5) +
  annotate("text", x = lab_dn, y = mean(y_down),
           label = "enriched in\nControl PTEC",
           colour = col_down, angle = 90, size = 3.4, lineheight = 0.95) +
  scale_fill_manual(values = c(down = col_down, up = col_up), guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = pretty(c(-x_lim, x_lim), n = 7)) +
  labs(x = expression(-log[10]~"(P.adjust)"), y = NULL,
       title = "Mitochondria-related GO (BP) in PTEC: AKI vs Control") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(
          colour = ifelse(
            levels(bar_df$TermLabel) %in% str_wrap(ord_up, 55),
            col_up, col_down),
          size = 8.5, face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = margin(10, 12, 8, 6))

ggsave(file.path(fig_dir, "PTEC_GO_BP_mitochondrial_barplot.pdf"),
       p4, width = 9.0, height = 6.0)

# ---------- 运行记录 ----------
sink(file.path(log_dir, "sessionInfo_05_mitochondrial_visualizations.txt"))
cat("=== 05_mitochondrial_visualizations.R ===\n"); print(Sys.time())
cat("\n=== PTEC 细胞数 ===\n"); print(table(ptec$group))
cat("\n=== ATP 词条基因拆解 ===\n")
print(breakdown %>% count(GO_ID, GO_Term, direction, encoding))
cat("\n=== 模块基因数 ===\n"); print(lengths(module_genes))
cat("\n=== 模块 Wilcoxon P ===\n"); print(mod_p)
cat("\n=== 热图基因模块分布 ===\n"); print(table(top50$module))
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 05_mitochondrial_visualizations.R 完成")
