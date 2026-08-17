# 复现说明

## 单一正式入口

在复现包根目录运行：

```powershell
D:\testproject\env_code\.venv\Scripts\python.exe analysis_scripts\run_final_research_article_pipeline.py
D:\testproject\env_code\.venv\Scripts\python.exe -m pytest -q
D:\testproject\env_code\.venv\Scripts\python.exe -m ruff check analysis_scripts tests
```

正式入口依次执行：

1. `33_final_discovery_statistics.py`：2,000 次 stage-stratified donor bootstrap、4,999 次 HC3-scaled studentized wild-bootstrap-t、LODO、影响诊断、稳健回归、严重度编码和预测审计；
2. `38_raw_count_p2rx4_validation.py`：核验远程归档与内部 manifests，读取两套 STAR/featureCounts 重建计数，拟合 P2RX4 负二项模型并生成来源、QC、PDR_S10 和参考定义敏感性；
3. `32_independent_p2rx4_validation.py`：对 GSE276892 processed values、GSE179568 临床混杂和 GSE94019 方向检查执行患者级统计；
4. `25_voigt_single_cell_localization.py`：对六个 GSE130636 donor-region libraries 执行供体聚合和三折 LODO 定位；
5. `26_updated_study_design.py`：生成最终 Figure 1；
6. `34_validate_final_submission.py`：核对稿件措辞、关键数值、原始重建完整性、S1–S28 工作簿和全部最终图。

随机种子固定为 `20260813`。正式参数不能由 `--quick` 结果替代。

## 输入与来源

- 发现输入位于 `project_inputs/`。
- 外部 processed matrices、临床来源、FASTQ manifests 与协议位于 `analysis_data/independent_validation/`。
- Voigt et al. 作者 cluster 映射及其 PDF 页码证据位于 `analysis_data/external_single_cell/AUTHOR_MAPPING_PROVENANCE.md`。
- P2RX4 本地协议 SHA-256 为 `FA0EBEEFE45709178E197B6F0819F7AA1E1692E50AA2C7F47E2601D7CA3C8EED`。该哈希验证文件内容，不构成创建时间证明或公开预注册。

## GSE276892 原始计数归档

脚本期望下列文件存在：

```text
analysis_data/independent_validation/remote_results/DRnet_GSE276892_remote_results.tar.gz
analysis_data/independent_validation/remote_results/DRnet_GSE276892_remote_results.tar.gz.sha256
```

归档大小为 `444,319,883` bytes，SHA-256 为：

```text
CF6784D6852D30A7A2A67FDFC0C7FB93E0A3ABBF3207CE775F860709795D65A6
```

归档含 521 个 tar members；两套工作流各有 221 个 manifest entries，内部校验失败数均为 0。该文件作为 GitHub Release 附件分发，普通 Git 树仅保留脚本、manifests、派生结果和归档哈希。

## 环境与输出

- 精确 Python 依赖见 `requirements.txt`，运行环境见 `COMPUTATIONAL_ENVIRONMENT.md`。
- Linux 原始重建使用 STAR 2.7.8a、featureCounts 2.0.1 和 GENCODE GRCh38.p13 release 42。
- 最终机器可读结果位于 `analysis_results/`；正文图与补充图位于 `figures/`；补充表位于 `supplementary_tables/Discover_Applied_Sciences_supplementary_tables.xlsx`。
- 公共代码地址为 <https://github.com/xumedlab/DRnet>。
