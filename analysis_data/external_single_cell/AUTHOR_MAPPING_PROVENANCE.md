# Voigt 2019 作者 cluster 映射来源

## 原始来源

- 论文：Voigt AP et al. *Molecular characterization of foveal versus peripheral human retina by single-cell RNA sequencing*. Experimental Eye Research. 2019;184:234--242.
- DOI：`10.1016/j.exer.2019.05.001`
- 核验位置：accepted manuscript PDF 第 38 页（文档共 42 页）的 Figure 1F。
- 核验日期：2026-08-13。源 PDF 因版权不随复现包分发。

## 映射证据

`voigt2019_author_cluster_mapping.csv` 的 `author_label` 逐项抄录 Figure 1F：

- clusters 1--2：Rods
- clusters 3--4：Cones
- clusters 5--6：Bipolar cells
- cluster 7：Retinal ganglion cells
- cluster 8A：Horizontal cells
- cluster 8B：Amacrine cells
- cluster 9：Unknown
- cluster 10：Pericytes
- cluster 11：Endothelial cells
- cluster 12：Microglia
- clusters 13--17：Glial cells

原文 Results 3.2 提供文字交叉核验：cluster 9 未表达所选细胞特异基因并被定义为 unknown；cluster 10 为 pericytes/smooth-muscle-like mural cells；cluster 11 为 endothelial cells；cluster 12 为 microglia；clusters 13--17 为 Müller cells and/or astrocytes 的 glial cells。

## 分析标签边界

- Cluster 9 从定位汇总中排除，未重新分类。
- Clusters 13--17 的原作者标签始终保留为 `Glial cells`。
- 分析列使用 `Müller-enriched glia`，依据原文说明这五个 cluster 均高表达 `RLBP1`，而 `ALDH1L1` 与 `GFAP` 较低；这是来源约束下的分析标签，不声称作者把每个细胞明确标为 Müller cell。
- 未进行 marker-based reclustering、自动注释或自行推断作者 cluster 映射。
