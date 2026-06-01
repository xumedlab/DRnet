# 论文生信结果整理

## 数据集与分组
- GSE160306 中最终纳入黄斑区样本 39 例。
- 分组构成：healthy control 10 例，diabetic 10 例，NPDR 9 例，NPDR/PDR + DME 10 例。

## 关键生信结果
- 以 healthy control vs NPDR/PDR + DME 为主比较，DESeq2 共识别显著差异基因 468 个。
- 与 HALLMARK_INFLAMMATORY_RESPONSE 取交后得到炎症核心基因 11 个：MSR1, TIMP1, OPRK1, LYN, FZD5, TLR3, NDP, NMI, CLEC5A, CYBB, CMKLR1。
- 其中满足正向严重度趋势的 progressive inflammatory genes 共 8 个：MSR1, TIMP1, OPRK1, LYN, FZD5, TLR3, NMI, CMKLR1。
- LASSO 最终保留的诊断候选基因为：MSR1, NMI, FZD5, TIMP1, CMKLR1, LYN, TLR3。
- 组合签名的 OOF AUC 为 0.870。
- 炎症 ssGSEA 主比较 P 值为 0.0963，严重度相关 rho 为 0.340。
- 免疫 ssGSEA 中最值得关注的细胞类型：Immature_B_cell, Macrophage。
- 复现度最高的 Hallmark 通路：HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, HALLMARK_INTERFERON_ALPHA_RESPONSE, HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_COAGULATION。

## 稿件撰写时应直接修正的地方
- 旧稿中 `limma / GSVA / glmnet / clusterProfiler` 的方法描述已经与当前 Python 管线不符，必须整体替换。
- 旧稿中“7 个诊断基因”“45 个 DE-INFGs”“39 条通路”“3 类免疫细胞”等数字不是当前结果，继续保留就是错稿。
- 当前仓库只覆盖生信部分，不包含临床验证和细胞实验的原始数据，Methods 2.7 以后及 Results 3.6 以后不能凭空补数字。

## 重点人工核查
- `progressive inflammatory genes` 的判定阈值使用 Spearman 正相关且 FDR < 0.1，属于探索性筛选，需要在文中明确写清楚。
- 免疫浸润来自 bulk RNA 的 ssGSEA 推断，只能写 supportive evidence，不能写成真实浸润比例。
- 组合签名 AUC 来自 OOF 预测而非独立外部验证，不能夸大成已验证的临床诊断模型。
