# ============================================================
# rerun_20261119 / 03_go_enrichment.R
# 功能：四类细胞 AKI vs CTRL 差异基因的 GO(BP) 功能富集
#       - 每类细胞分别富集上调基因、下调基因, 保留显著条目(p.adjust < 0.05)
#       - PTEC 额外: 筛选线粒体 CC(BP) 条目 + Fisher 精确检验评估
#         上调/下调集分别富集线粒体基因集的程度
#       - 废弃词条审计(剔除已废弃 GO ID)
# 输入：tables/DEG_{Neutrophil,PTEC,GEC,PT}_AKI_vs_CTRL.csv
# 输出：tables/GO_BP_*.csv
#       tables/annotation_versions.csv
#       tables/obsolete_GO_terms_audit.csv
#       tables/PTEC_mito_term_enrichment_Fisher.csv
# ============================================================

set.seed(2026)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
})

rerun_dir <- "/home/zhouxin/liudi_AKI_project/rerun_20261119"
tab_dir   <- file.path(rerun_dir, "results", "tables")
log_dir   <- file.path(rerun_dir, "logs")

# ---------- 0. 记录注释包版本 ----------
ann_ver <- data.frame(
  package  = c("org.Hs.eg.db", "GO.db", "clusterProfiler"),
  version  = c(as.character(packageVersion("org.Hs.eg.db")),
               as.character(packageVersion("GO.db")),
               as.character(packageVersion("clusterProfiler"))))
write.csv(ann_ver, file.path(tab_dir, "annotation_versions.csv"),
          row.names = FALSE)

# ---------- 1. 辅助函数 ----------
# 废弃 GO ID 审计: 只保留当前 GO.db 中有效的 ID
active_go <- keys(GO.db, keytype = "GOID")

# symbol -> entrez 映射
sym2entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys  = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL") %>%
  distinct(SYMBOL, .keep_all = TRUE)

get_entrez <- function(syms) {
  sym2entrez %>% filter(SYMBOL %in% syms) %>% pull(ENTREZID) %>%
    unique() %>% as.character()
}

# 有效 GO ID 过滤
filter_active_go <- function(res) {
  res %>% mutate(is_active = ID %in% active_go)
}

# BP 富集: 分别富集上/下调基因, 合并并加 direction 列
run_go_bp <- function(up_syms, down_syms, universe_syms,
                      cell_group, ont = "BP",
                      pval_cutoff = 0.05, qval_cutoff = 0.05) {
  up_ids   <- get_entrez(up_syms)
  down_ids <- get_entrez(down_syms)
  uni_ids  <- get_entrez(universe_syms)

  do_enrich <- function(ids, direction_label) {
    if (length(ids) < 5) {
      warning(cell_group, ": ", direction_label,
              "基因 entrez 不足, 跳过富集")
      return(NULL)
    }
    tryCatch({
      e <- enrichGO(
        gene          = ids,
        universe      = uni_ids,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cutoff,
        qvalueCutoff  = qval_cutoff,
        readable      = TRUE)
      as.data.frame(e) %>%
        mutate(direction = direction_label,
               cell_group = cell_group) %>%
        filter(ID %in% active_go)
    }, error = function(e) {
      warning("enrichGO error (", direction_label, "): ", conditionMessage(e))
      NULL
    })
  }

  up_df   <- do_enrich(up_ids,   "Up_in_AKI")
  down_df <- do_enrich(down_ids, "Down_in_AKI")

  res <- bind_rows(up_df, down_df)
  if (nrow(res) > 0) {
    res %>% arrange(direction, p.adjust) %>%
      select(cell_group, direction, ID, Description,
             GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
  } else NULL
}

# ---------- 2. 四类细胞 GO 富集 ----------
cell_groups <- c("Neutrophil", "PTEC", "GEC", "PT")
all_go <- list()
all_obsolete <- list()

for (cg in cell_groups) {
  deg <- read.csv(file.path(tab_dir, paste0("DEG_", cg, "_AKI_vs_CTRL.csv")),
                  stringsAsFactors = FALSE)
  # 显著差异基因(p_val_adj < 0.05 & |log2FC| > 0.25)
  sig <- deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
  up   <- sig %>% filter(avg_log2FC > 0) %>% pull(gene)
  down <- sig %>% filter(avg_log2FC < 0) %>% pull(gene)
  uni  <- deg %>% pull(gene)        # universe = 本次检测到的全部基因

  cat(cg, ": sig =", nrow(sig),
      "  Up:", length(up), "  Down:", length(down), "\n")

  go <- run_go_bp(up, down, uni, cell_group = cg)
  if (!is.null(go)) {
    all_go[[cg]] <- go
    write.csv(go,
              file.path(tab_dir, paste0("GO_BP_", cg, "_AKI_vs_CTRL.csv")),
              row.names = FALSE)
    cat("  GO(BP) entries:", nrow(go),
        "(Up:", sum(go$direction == "Up_in_AKI"),
        " Down:", sum(go$direction == "Down_in_AKI"), ")\n")
    # 废弃 ID 审计
    obso <- go %>% filter(!ID %in% active_go) %>%
      transmute(cell_group, direction, ID, Description, p.adjust)
    if (nrow(obso) > 0) all_obsolete[[cg]] <- obso
  }
}

obso_df <- bind_rows(all_obsolete)
write.csv(obso_df, file.path(tab_dir, "obsolete_GO_terms_audit.csv"),
          row.names = FALSE)
cat("\n废弃 GO 条目数:", nrow(obso_df),
    "(已在富集输出中自动剔除)\n")

# ---------- 3. PTEC 线粒体定向 Fisher 检验 ----------
# GO:0005739(mitchondrion, CC) 及其子孙 BP 条目
get_offspring <- function(ids, bimap) {
  lst <- as.list(bimap)
  unique(unlist(lst[ids], use.names = FALSE))
}
mito_bp_ids <- unique(c(get_offspring("GO:0005739", GOBPOFFSPRING)))
mito_cc_ids <- unique(c("GO:0005739",
                        get_offspring("GO:0005739", GOCCOFFSPRING)))

# PTEC DEG
ptec_deg <- read.csv(file.path(tab_dir, "DEG_PTEC_AKI_vs_CTRL.csv"),
                     stringsAsFactors = FALSE)
ptec_sig <- ptec_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
ptec_up   <- ptec_sig %>% filter(avg_log2FC > 0) %>% pull(gene)
ptec_down <- ptec_sig %>% filter(avg_log2FC < 0) %>% pull(gene)
ptec_uni  <- ptec_deg %>% pull(gene)

# CC 线粒体基因集(GOALL 传递注释)
mito_syms_cc <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = mito_cc_ids,
                        columns = "SYMBOL", keytype = "GOALL")$SYMBOL) %>%
  na.omit() %>% unique()

# Fisher 检验(2x2: 线粒体 vs 非线粒体 x DEG vs 非 DEG)
run_fisher <- function(sym_subset, sym_universe, sym_mito, label) {
  mat <- matrix(c(
    length(intersect(sym_subset, sym_mito)),
    length(setdiff(sym_subset,   sym_mito)),
    length(intersect(sym_universe, sym_mito)) - length(intersect(sym_subset, sym_mito)),
    length(setdiff(sym_universe,   sym_mito)) - length(setdiff(sym_subset,   sym_mito))
  ), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  data.frame(label             = label,
             n_deg_mito        = length(intersect(sym_subset, sym_mito)),
             n_deg_non_mito    = length(setdiff(sym_subset,   sym_mito)),
             n_universe_mito   = length(intersect(sym_universe, sym_mito)),
             n_universe_non_mito= length(setdiff(sym_universe, sym_mito)),
             odds_ratio        = ft$estimate,
             p_value           = ft$p.value,
             stringsAsFactors  = FALSE)
}

ft_df <- bind_rows(
  run_fisher(ptec_up,   ptec_uni, mito_syms_cc, "PTEC Up in AKI  vs mito"),
  run_fisher(ptec_down, ptec_uni, mito_syms_cc, "PTEC Down in AKI vs mito"))
ft_df$p_adj_fisher <- p.adjust(ft_df$p_value, method = "BH")
write.csv(ft_df, file.path(tab_dir, "PTEC_mito_term_enrichment_Fisher.csv"),
          row.names = FALSE)
cat("\n=== PTEC 线粒体 Fisher 检验 ===\n")
print(ft_df %>% select(label, odds_ratio, p_value, p_adj_fisher))

# PTEC BP 富集中的线粒体条目审计
ptec_go <- all_go[["PTEC"]]
if (!is.null(ptec_go)) {
  mito_in_bp <- ptec_go %>%
    filter(ID %in% mito_bp_ids | grepl("mitochond|oxidative phosphorylation|electron transport|oxphos",
                                       Description, ignore.case = TRUE))
  write.csv(mito_in_bp,
            file.path(tab_dir, "PTEC_GO_BP_mitochondrial_terms.csv"),
            row.names = FALSE)
  cat("\nPTEC BP 富集中的线粒体相关条目:", nrow(mito_in_bp),
      "\n"); print(mito_in_bp %>% select(direction, Description, p.adjust))
}

# ---------- 4. 运行记录 ----------
sink(file.path(log_dir, "sessionInfo_03_go_enrichment.txt"))
cat("=== 03_go_enrichment.R 运行时间 ===\n"); print(Sys.time())
cat("\n=== 注释包版本 ===\n"); print(ann_ver)
cat("\n=== Fisher 检验 ===\n"); print(ft_df)
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()

message(">>> 03_go_enrichment.R 完成")
