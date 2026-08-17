# GSE276892 原始 reads 的远程 Linux 重建说明

## 为什么正式分析必须改用 Linux

STAR 2.7.8a 官方发布物仅提供 Linux/macOS 二进制，官方支持平台不包含原生 Windows。本地 MSYS2 兼容编译已经完成绝大多数目标文件，但最终链接仍受 `vfork`/`kill` POSIX 接口影响。即使继续打补丁得到 Windows 可执行文件，也会引入无法与作者/审稿人标准环境等同的非官方平台差异。因此，论文正式结果采用官方 Linux x86_64 static STAR 2.7.8a 和 Linux Subread/featureCounts 2.0.1；Windows 尝试仅保留为环境诊断记录，不进入主分析。

本包执行两个流程：

1. `all_regions_source_reconstruction`：GENCODE GRCh38.p13 release 42 all-regions FASTA/GTF，作为最接近原研究描述的 raw-read reconstruction。
2. `primary_assembly_sensitivity`：GENCODE release 42 primary assembly FASTA/GTF，作为标准参考敏感性分析。

两套流程均执行 STAR 2.7.8a 比对、featureCounts 2.0.1 exon/gene_id 非链特异计数。21 个 GSE147657 文件是 7 名复用对照的三个技术 lanes；脚本先逐 lane 比对和计数，再按供体求和，最终形成 17 个生物样本（8 PDR、9 control），不会把技术 lane 当作独立患者。

## 服务器要求

- Linux x86_64；不支持 ARM，也不建议 WSL 作为正式稿环境。
- RAM：最低 32 GiB，建议 64 GiB。
- CPU：建议 16–32 threads。
- 空闲磁盘：输入与计算工作区位于同一文件系统时，`KEEP_BAM=0` 至少 180 GiB；如果输入放在 NFS、计算中间件放本地，`CLEAN_GENERATED_REFERENCE_AND_INDEX=1` 时本地工作区建议至少 100 GiB。若保留 BAM，建议至少 400 GiB，并设置 `KEEP_BAM=1`。
- 软件：`bash`、`python3`、`tar`、`gzip`、GNU `coreutils`。不需要 root 权限，也不需要 Conda。
- 建议将工作目录置于本地 SSD/NVMe 文件系统，不要直接在高延迟网络盘上建 STAR index。

## 需要从 Windows 上传的内容

以远程目录 `/data/$USER/DRnet_GSE276892_reconstruction` 为例，必须保持下列相对路径：

```text
DRnet_GSE276892_reconstruction/
├── remote_linux_reconstruction/
│   ├── README_zh.md
│   ├── preflight.sh
│   ├── run_reconstruction.sh
│   ├── collect_results.sh
│   ├── config/
│   │   ├── sample_sheet.tsv
│   │   ├── input_files.tsv
│   │   └── provenance_sources.tsv
│   └── scripts/
│       ├── verify_remote_inputs.py
│       └── aggregate_featurecounts.py
├── analysis_data/independent_validation/raw_reads/
│   ├── GSE276892/                 # 10 FASTQ
│   └── GSE147657/                 # 21 FASTQ lanes
└── tools/
    ├── downloads/
    │   ├── STAR-2.7.8a.tar.gz
    │   └── subread-2.0.1-Linux-x86_64.tar.gz
    └── reference/
        ├── gencode_v42_all_regions/
        │   ├── GRCh38.p13.genome.fa.gz
        │   └── gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz
        └── gencode_v42_primary_assembly/
            ├── GRCh38.primary_assembly.genome.fa.gz
            └── gencode.v42.primary_assembly.annotation.gtf.gz
```

Windows 本地源目录为：

```text
D:\testproject\DRnet\DRnet\show\Discover_Applied_Sciences_submission_package
```

不要上传以下内容：`tools/portable/msys64`、Windows 版 `featureCounts.exe`、STAR 的 MSYS2 编译中间文件、`quarantine`、投稿主文或内部审稿记录。

### 推荐上传方式

大文件总量为 31.62 GiB。优先使用 WinSCP 的 SFTP 目录上传并开启断点续传；也可从 PowerShell 使用 `scp -r`，但连接中断后的恢复能力不如 WinSCP/rsync。上传完成后不要手工判断完整性，必须运行下面的 `preflight.sh`，它会对所有 31 个 FASTQ、4 个参考文件和 2 个工具归档同时核验文件大小、MD5 和 SHA-256。

对于 `ai2406` 当前磁盘布局，推荐将上述整个输入包放在 `/data/NFS_DeepLearningResults` 下有写权限的专用目录，将 `WORK_ROOT` 放在本地 `/home/wensheng/gjq_workspace`。NFS 只承载固定输入，STAR index、解压参考和临时 BAM 在本地完成，避免在 NFS 上建立高 I/O index。NFS 当前已使用 95%，上传前仍需确认个人目录写权限和配额。

## 在服务器执行

登录服务器后：

```bash
cd /data/NFS_DeepLearningResults/DRnet_GSE276892_wensheng
chmod +x remote_linux_reconstruction/*.sh
env \
  WORK_ROOT=/home/wensheng/gjq_workspace/DRnet_GSE276892_work \
  MIN_FREE_GB=100 \
  bash remote_linux_reconstruction/preflight.sh
```

只有出现 `Preflight PASS` 才开始计算。建议使用 `nohup`，避免 SSH 断开终止任务：

```bash
cd /data/NFS_DeepLearningResults/DRnet_GSE276892_wensheng
nohup env \
  WORK_ROOT=/home/wensheng/gjq_workspace/DRnet_GSE276892_work \
  MIN_FREE_GB=100 \
  THREADS=16 \
  KEEP_BAM=0 \
  CLEAN_GENERATED_REFERENCE_AND_INDEX=1 \
  RUN_PRIMARY_ASSEMBLY_SENSITIVITY=1 \
  bash remote_linux_reconstruction/run_reconstruction.sh \
  > remote_linux_reconstruction/nohup_reconstruction.log 2>&1 &
echo $! > remote_linux_reconstruction/reconstruction.pid
```

查看进度：

```bash
tail -n 80 -f remote_linux_reconstruction/nohup_reconstruction.log
```

如需确认进程：

```bash
ps -fp "$(cat remote_linux_reconstruction/reconstruction.pid)"
```

脚本可从已完成的 index、lane counts 和 STAR log 恢复。不要同时启动第二个实例；发生中断时，先确认旧 PID 已退出，再用同一条命令重启。

## 需要回传的结果

成功结束后会生成：

```text
remote_work/DRnet_GSE276892_remote_results.tar.gz
remote_work/DRnet_GSE276892_remote_results.tar.gz.sha256
```

只需下载这两个文件，不必下载 STAR index、解压参考或 BAM。结果归档包含：

- 两套参考流程的 17-sample gene-count matrix；
- `P2RX4_counts.tsv`；
- 31 个 lane 的 featureCounts 输出及 summary；
- STAR `Log.final.out` 等 QC 文件；
- sample-level library size、mapping rate、assigned rate、detected-gene count；
- 服务器环境、工具版本、输入核验报告和所有输出 SHA-256。

将它们放回 Windows：

```text
D:\testproject\DRnet\DRnet\show\Discover_Applied_Sciences_submission_package\analysis_data\independent_validation\remote_results\
```

随后本地流水线将完成 P2RX4 count-level negative-binomial 模型、来源批次可识别性、仅 2 个新对照方向敏感性、mapping/library-size/detected-gene QC 相关性，以及两套参考的一致性分析。raw-read 重建是在查看处理矩阵后的 protocol deviation，不能追溯性宣称为完全 confirmatory 或有外部时间戳的 target lock。

## 固定参数与可解释性

- `SJDB_OVERHANG=100`：STAR 默认值，符合“STAR default”重建思路；脚本记录实际值。
- `-t exon -g gene_id -s 0`：featureCounts 的 exon/gene-level、unstranded 设定。
- `KEEP_BAM=0`：计数和 summary 验证完成后删除可再生 BAM，保留全部日志与 count；这不改变分析结果。如服务器磁盘足够可改为 `KEEP_BAM=1`。
- `CLEAN_GENERATED_REFERENCE_AND_INDEX=1`：每套 workflow 的 count、QC 和输出 SHA-256 完成后，删除可从固定输入重建的解压参考和 STAR index；适合 `ai2406` 当前 176 GiB 本地空闲空间。
- 不执行 ComBat。疾病状态与样本来源高度共线，批次校正不能制造缺失的设计重叠；后续模型会直接报告该限制。

权威来源：

- STAR 2.7.8a 官方仓库与平台说明：https://github.com/alexdobin/STAR
- Subread/featureCounts 官方项目：https://subread.sourceforge.net/
- GENCODE human release 42：https://www.gencodegenes.org/human/release_42.html
- ENA Portal API：https://www.ebi.ac.uk/ena/portal/api/
