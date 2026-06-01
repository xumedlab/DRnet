import pandas as pd
import re

def build_filtered_metadata(matrix_file, counts_file, output_file):
    print(f"正在解析 GEO 临床数据文件: {matrix_file} ...")
    
    sample_data = {}
    with open(matrix_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # 提取样本 title
            if line.startswith('!Sample_title'):
                titles = line.split('\t')[1:]
                titles = [t.strip('"') for t in titles]
                for i, title in enumerate(titles):
                    sample_data[i] = {'SampleID_GEO': title}
            
            # 提取样本临床特征
            elif line.startswith('!Sample_characteristics_ch1'):
                chars = line.split('\t')[1:]
                chars = [c.strip('"') for c in chars]
                for i, char in enumerate(chars):
                    if i in sample_data and ':' in char:
                        key, value = char.split(':', 1)
                        key = key.strip().lower()
                        value = value.strip()
                        sample_data[i][key] = value

    df_meta_full = pd.DataFrame.from_dict(sample_data, orient='index')
    
    # 明确指定我们从报错日志中看到的真实列名
    region_key = 'sample_site'
    disease_key = 'disease_group'
    
    if region_key not in df_meta_full.columns or disease_key not in df_meta_full.columns:
        raise KeyError(f"列名匹配依然失败。现有列名: {df_meta_full.columns.tolist()}")

    # 映射疾病状态
    def map_group(disease_state):
        state = str(disease_state).lower()
        if 'control' in state or 'normal' in state or 'healthy' in state:
            return 'Control'
        elif 'dr' in state or 'diabetic' in state or 'pdr' in state or 'npdr' in state:
            return 'DR'
        return 'Other'

    df_meta_full['Group'] = df_meta_full[disease_key].apply(map_group)
    df_meta_full['Region_Clean'] = df_meta_full[region_key].apply(lambda x: str(x).lower())

    # 执行过滤：仅保留黄斑区 (Macula) 且标签为 Control 或 DR 的样本
    print("正在执行过滤：仅保留黄斑区 (Macula) 样本...")
    df_filtered = df_meta_full[
        (df_meta_full['Region_Clean'].str.contains('macula')) & 
        (df_meta_full['Group'].isin(['Control', 'DR']))
    ].copy()

    print(f"过滤后符合条件的 Macula 样本数: {len(df_filtered)}")

    # 对齐 counts 矩阵
    try:
        counts_df = pd.read_csv(counts_file, sep='\t', nrows=0)
        counts_columns = set(counts_df.columns)
    except Exception as e:
        raise FileNotFoundError(f"无法读取 counts 文件以验证列名: {e}")

    valid_samples = []
    for _, row in df_filtered.iterrows():
        matched = False
        # 1. 尝试匹配作者自己定义的 sampleid 列 (如 'sample_09')
        if 'sampleid' in row and row['sampleid'] in counts_columns:
            valid_samples.append(row['sampleid'])
            matched = True
        # 2. 尝试匹配 GEO 的 SampleID_GEO (如 'Macula_DR_09')
        elif not matched and row['SampleID_GEO'] in counts_columns:
            valid_samples.append(row['SampleID_GEO'])
            matched = True
        # 3. 尝试使用正则硬抓取 'sample_xx'
        elif not matched:
            match = re.search(r'(sample_\d+)', row['SampleID_GEO'], re.IGNORECASE)
            if match and match.group(1).lower() in counts_columns:
                valid_samples.append(match.group(1).lower())
                matched = True
            elif 'sampleid' in row:
                match2 = re.search(r'(sample_\d+)', row['sampleid'], re.IGNORECASE)
                if match2 and match2.group(1).lower() in counts_columns:
                    valid_samples.append(match2.group(1).lower())
                    matched = True
        
        if not matched:
            print(f"警告: 样本 {row.to_dict()} 无法在 counts 矩阵中找到对应列！")
            valid_samples.append(None) # 占位符，随后丢弃

    df_filtered['SampleID_Matched'] = valid_samples
    # 丢弃未匹配上的样本
    df_filtered = df_filtered.dropna(subset=['SampleID_Matched'])

    # 组装最终可用的 metadata
    final_metadata = df_filtered[['SampleID_Matched', 'Group']].rename(columns={'SampleID_Matched': 'SampleID'})
    
    final_metadata.to_csv(output_file, index=False)
    print(f"\n成功生成标准化临床信息文件: {output_file}")
    print(f"最终纳入差异分析队列: Control 组 {len(final_metadata[final_metadata['Group']=='Control'])} 例，DR 组 {len(final_metadata[final_metadata['Group']=='DR'])} 例。")

if __name__ == "__main__":
    matrix_file = '/home/wensheng/gjq_workspace/DRnet/GSE160306_series_matrix.txt'
    counts_file = '/home/wensheng/gjq_workspace/DRnet/GSE160306_human_retina_DR_totalRNA_counts.txt'
    output_file = '/home/wensheng/gjq_workspace/DRnet/metadata.csv'
    
    build_filtered_metadata(matrix_file, counts_file, output_file)