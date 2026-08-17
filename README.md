# Discover Applied Sciences Research Article 复现包

本目录对应 Research Article：**Donor-aware cross-cohort transcriptomics prioritizes P2RX4 in human diabetic retinopathy**。公共代码仓库地址为 <https://github.com/xumedlab/DRnet>。

## 研究定位

本研究在 GSE160306 的 26 名糖尿病供体中进行供体级黄斑转录组分析，分别估计 DR severity 总关联和 DME 条件关联。P2RX4 在两个估计目标中均排名第一。外部眼内数据仅用于单目标压力测试，不参与候选筛选。

主要结果与证据边界如下：

- P2RX4 总关联系数为 `0.688`，HC3+t 双侧 `P=0.00149`；158 基因校正后没有主估计目标达到 FDR<0.05。
- P2RX4 的供体 bootstrap top-5 频率为 `61.65%`，LODO top-5 频率为 `100%`；因此结论是候选优先级，不是稳定多基因签名。
- GSE276892 的 31 个 FASTQ lanes 已按 17 个患者样本聚合，并以 STAR 2.7.8a、featureCounts 2.0.1 和 GENCODE GRCh38.p13 release 42 重建计数。疾病单变量负二项模型得到 log2 fold change `0.295`、双侧 `P=0.615`；来源、年龄/性别、去除 PDR_S10 和双参考敏感性均未形成精确疾病效应。
- GSE179568 的未调整膜区室比较较强，但 PDR 与 macular-pucker 对照年龄完全不重叠。年龄/性别 HC3 模型的疾病系数双侧 `P=0.773`，不能从年龄和组织组成中分离疾病效应。
- GSE94019 仅为组织不匹配的方向检查；GSE130636 仅提供正常视网膜定位背景。不同 GEO accession 不等于患者身份独立。

因此，稿件支持 P2RX4 作为后续独立视网膜队列和实验研究的聚焦候选，不声称复制、诊断性能、确定细胞来源、受体活性或因果机制。

## 目录结构

- `analysis_scripts/`：发现分析、原始计数负二项模型、外部敏感性、正常视网膜定位、验证器和发布树生成器。
- `analysis_data/`：外部数据、来源记录、协议、FASTQ manifests 和远程重建结果。
- `analysis_results/`：正文与补充材料引用的 CSV/JSON 结果。
- `figures/`：Figure 1–5 及两个新增补充图的 PDF/PNG。
- `project_inputs/`：冻结的发现队列输入。
- `manuscript/`：主文、Supplementary Information、cover letter 和参考文献源文件。
- `submission_files/`：投稿 PDF、S1–S28 工作簿、源码包和校验清单。
- `supplementary_tables/`：与投稿文件一致的 S1–S28 工作簿副本。

## 正式复现

在复现包根目录运行：

```powershell
D:\testproject\env_code\.venv\Scripts\python.exe analysis_scripts\run_final_research_article_pipeline.py
D:\testproject\env_code\.venv\Scripts\python.exe -m pytest -q
D:\testproject\env_code\.venv\Scripts\python.exe -m ruff check analysis_scripts tests
```

正式入口固定随机种子 `20260813`，执行 2,000 次供体 bootstrap、4,999 次 studentized wild-bootstrap-t 和 10,000 次外部队列 bootstrap。`--quick` 仅用于烟雾测试，其输出不能替代正式结果。

GSE276892 原始计数重建使用的归档文件为：

```text
analysis_data/independent_validation/remote_results/DRnet_GSE276892_remote_results.tar.gz
```

其 SHA-256 为：

```text
CF6784D6852D30A7A2A67FDFC0C7FB93E0A3ABBF3207CE775F860709795D65A6
```

该 444,319,883-byte 归档适合作为 GitHub Release 附件，不进入普通 Git 树。公开树中的脚本会从上述固定路径读取该附件。
