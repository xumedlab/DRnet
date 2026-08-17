# 计算环境

- 开发和验证操作系统：Windows
- 运行接口：Python 3.10.11
- 随机种子：`20260813`
- LaTeX：Tectonic 0.16.9
- 工作簿生成：Codex bundled Node.js 和 `@oai/artifact-tool`
- 精确 Python 包版本：`requirements.txt`

## 正式参数

```text
discovery_stage_stratified_donor_bootstrap=2000
studentized_wild_bootstrap_t=4999
wild_bootstrap_weights=Rademacher
wild_bootstrap_null_residual_scaling=HC3
independent_validation_bootstrap=10000
random_seed=20260813
```

Tectonic 编译和 PDF 渲染只改变排版产物，不参与统计计算。工作簿由结果 CSV 机械生成；S22 依次包含六个 donor-region library 的 dominant cell type 和供体聚合定位。
