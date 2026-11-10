# GSE183276 AKI vs CTRL Analysis Record (Reanalysis)

## 1. 数据来源
- GEO 数据集: [GSE183276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276)
- Kidney Healthy-Injury Cell Atlas (scCv3, 10X scRNA-seq)
- 输入文件: `..._Seurat_03282022.h5Seurat`(counts + metadata)、`..._Metadata_03282022.txt.gz`(全量注释)

## 2. 分析设计

### 分组
- **AKI**: Acute Kidney Injury
- **CTRL**: Normal Reference (`condition.l1 == "Ref"`)
- CKD(DKD / H-CKD)样本不参与分析

### 目标细胞群体(`subclass.l2`)
| 群体 | 取值 | 说明 | CTRL | AKI |
|---|---|---|---|---|
| Neutrophil | `N` | 中性粒细胞 | 50 | 407 |
| GEC | `EC-GC` | 肾小球毛细血管内皮细胞 | 895 | 206 |
| PTEC | `EC-PTC` | 肾小管周毛细血管内皮细胞 | 1489 | 644 |
| PT | `PT-S1/2`, `PT-S3` | 健康近端肾小管上皮细胞 | 3087 | 1079 |

合计 7857 个细胞。其余细胞类型与样本不参与任何分析。

### 降维聚类流程
提取 4 类细胞后重新执行: NormalizeData(LogNormalize) → FindVariableFeatures(vst, 2000) →
ScaleData → RunPCA(30 PCs) → FindNeighbors(dims 1:20) → FindClusters(resolution 0.5) →
RunUMAP(dims 1:20)。得到 15 个无监督 cluster，与人工注释高度一致(见 `Cluster_vs_cell_group.csv`)。

## 3. 脚本(按执行顺序)

| 脚本 | 功能 |
|---|---|
| `scripts/01_data_prep.R` | 提取 4 群体 + AKI/CTRL; 标准化、PCA、聚类、UMAP |
| `scripts/02_umap_markers.R` | 重聚类 UMAP(群体/聚类/分组/分面) + marker DotPlot |
| `scripts/03_neutrophil_stats.R` | 各样本中性粒细胞数量与占比 + Wilcoxon 检验 |
| `scripts/04_deg_analysis.R` | PTEC/GEC/PT 差异基因 + PTEC 火山图 |
| `scripts/05_go_enrichment.R` | 三群体 GO(BP) 富集 + GEC 双向点图 |
| `scripts/06_extra_figures.R` | PTEC GO 双向条形图 + GEC 火山图 + PTEC/GEC 比例统计图 |
| `scripts/07_supplementary_visuals.R` | 补充: 中性粒细胞数量/比例双面板图 + PTEC/GEC 数量与比例 2x2 图 + 统计汇总表 |
| `scripts/08_boxplot_style_stats.R` | 按参考图样式(箱线图+散点+括号 p 值)重绘 PTEC/GEC 数量与比例 |

### Marker 基因(经数据表达验证)
- Neutrophil: S100A8, S100A9, CSF3R
- PTEC: PLVAP, VWF, KDR
- GEC: EHD3, SOST, EMCN
- PT: LRP2, CUBN, SLC22A6

## 4. 主要结果

### 中性粒细胞
- 占比: AKI 显著高于 CTRL, Wilcoxon **p = 0.0053**
- 数量: AKI 显著高于 CTRL, Wilcoxon **p = 0.0014**

### 血管内皮细胞比例
| 群体 | CTRL 均值 | AKI 均值 | Wilcoxon p |
|---|---|---|---|
| PTEC (EC-PTC) | 更高 | 更低 | 4.94e-04 |
| GEC (EC-GC) | 更高 | 更低 | 1.36e-04 |

比例分母为该样本全部细胞数(全量 metadata)。

### 补充: 数量与比例的完整统计(Wilcoxon, AKI vs CTRL)
| 群体 | 指标 | AKI 均值 | CTRL 均值 | p 值 | 结论 |
|---|---|---|---|---|---|
| Neutrophil | 细胞数 | 33.9 | 2.5 | 1.44e-03 | AKI 升高 |
| Neutrophil | 比例 (%) | 0.979 | 0.262 | 5.30e-03 | AKI 升高 |
| PTEC | 细胞数 | 53.7 | 74.5 | 1.15e-01 | 无显著差异 |
| PTEC | 比例 (%) | 2.52 | 7.05 | 4.94e-04 | CTRL 更高 |
| GEC | 细胞数 | 17.2 | 44.8 | 2.27e-02 | CTRL 更高 |
| GEC | 比例 (%) | 0.727 | 3.82 | 1.36e-04 | CTRL 更高 |

注: PTEC 绝对细胞数差异未达显著(p = 0.115), 但其占样本细胞总数的比例在 AKI 显著下降;
GEC 的绝对数量与比例在 AKI 均显著降低。n = 12 (AKI) / 20 (CTRL)。

### 差异基因(AKI vs CTRL, Wilcoxon, p_adj<0.05 & |log2FC|>0.25)
| 群体 | 显著 DEG |
|---|---|
| PTEC | 1006 |
| GEC | 342 |
| PT | 3096 |

- PTEC top DEG: 上调 TGM2, MARCKS, ZFP36L1, SELENOM, PECAM1; 下调 MT2A, CEBPD, MT1X, MT1M, APOLD1

### GO(BP)富集(BH 校正, p<0.05 & q<0.05)
| 群体 | 显著条目 |
|---|---|
| PTEC | 526 |
| GEC | 36 |
| PT | 864 |

- GEC 上调基因富集于 ATP 合成/质子跨膜转运等能量代谢通路; 下调富集于金属离子应答、糖皮质激素应答、细胞黏附负调控

## 5. 输出文件

### 图(`results/figures/`, PDF, 全英文)
- `UMAP_cell_groups.pdf` — 重聚类 UMAP(按细胞群体着色)
- `UMAP_clusters.pdf` — 无监督聚类 UMAP
- `UMAP_AKI_vs_CTRL.pdf` — 按分组着色
- `UMAP_split_by_condition.pdf` — 分组分面
- `DotPlot_canonical_markers.pdf` — 特征性 marker 点图
- `Neutrophil_proportion_stats.pdf` / `Neutrophil_count_stats.pdf` — 中性粒细胞统计(含 Wilcoxon p)
- `Volcano_PTEC_AKI_vs_CTRL.pdf` — PTEC 火山图(top10 标注)
- `GO_BP_GEC_dotplot.pdf` — GEC GO 双向点图(cluster=down/up, 点大小=GeneRatio, 颜色=p.adjust)
- `GO_BP_PTEC_enrichment_barplot.pdf` — PTEC GO 双向条形图(x 为带方向 -log10(P.adjust), 粉=AKI 富集, 蓝=CTRL 富集)
- `Volcano_GEC_AKI_vs_CTRL.pdf` — GEC 火山图(top10 标注)
- `Cell_proportion_PTEC.pdf` / `Cell_proportion_GEC.pdf` — 各样本细胞比例(散点+均值±SD+Wilcoxon p 值)
- `Neutrophil_count_proportion.pdf` — 中性粒细胞数量与比例双面板(散点+均值±SD+p 值)
- `Endothelial_count_proportion_stats.pdf` — PTEC / GEC 数量与比例 2x2 分面统计图
- `Cell_count_PTEC_GEC_boxplot.pdf` — PTEC / GEC 细胞数量箱线图(散点+括号 p 值, 参考图样式)
- `Cell_proportion_PTEC_GEC_boxplot.pdf` — PTEC / GEC 细胞比例箱线图(散点+括号 p 值, 参考图样式)

### 表(`results/tables/`)
- `DEG_{PTEC,GEC,PT}_AKI_vs_CTRL.csv` — 差异基因完整表
- `GO_BP_{PTEC,GEC,PT}_AKI_vs_CTRL.csv` — GO 富集表
- `Neutrophil_per_sample_stats.csv` — 各样本中性粒细胞统计
- `Cell_proportion_PTEC_GEC.csv` — 各样本 PTEC/GEC 细胞数与比例
- `Supplementary_stats_summary.csv` — 三群体数量/比例的 n、均值、中位数、SD 与 Wilcoxon p 值汇总
- `Cell_counts_by_group.csv`, `Cluster_vs_cell_group.csv` — 细胞数与聚类核查

### 中间数据(`results/processed_data/`)
- `seu_subset_reclustered.rds` — 提取并重新降维聚类的 Seurat 对象

### 运行记录(`logs/`)
- `run_0X_*.log`(运行日志) 与 `sessionInfo_0X_*.txt`(R 环境与包版本)

## 6. 可重复性
- 所有脚本 `set.seed(2026)`, 降维/聚类/UMAP 均显式指定 seed
- 每个脚本输出 `sessionInfo()`
- 重现: `Rscript scripts/01_data_prep.R` → `02` → `03` → `04` → `05` → `06` → `07` → `08`
  (07/08 只读取 03/06 已生成的统计表, 不重算数据, 不影响前面任何结果)
- 环境: R 4.3.1, Seurat 5.3.0, clusterProfiler 4.8.3
