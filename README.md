# SICER2_workflow

这是一个只负责运行 SICER2 的简洁 Snakemake workflow：

```text
processed SE BAM
  -> 按 group pooling
  -> sicer domain calling
  -> 按 contrast 运行 sicer_df differential calling
```

本仓库不负责 FASTQ、比对、MACS3、IDR、bigWig 或下游注释。

当前阶段的目标是先跑通流程并熟悉 `sicer` / `sicer_df` 的原生结果。因此：

- SICER2 生成什么就保留什么。
- 原生输出文件名不修改。
- 不另外生成标准化结果表。
- 不拆分 increased/decreased 表。
- 不追加 fold-change 列，也不做二次 fold-change 筛选。
- 下游筛选、注释和表格整理以后单独实现。

## 输入表

`config/samples.txt` 每行是一个 biological replicate，使用 tab 分隔：

```text
group  sample  treatment_bam  control_bam  layout
WT_WL  WT_WL_rep1  /path/IP1.bam  /path/IgG_A.bam  SE
WT_WL  WT_WL_rep2  /path/IP2.bam  /path/IgG_A.bam  SE
```

- `sample` 必须唯一。
- `control_bam` 可以是 Input、IgG 或留空。
- 同一 group 必须全部有 control 或全部没有 control。
- 第一版只接受 single-end BAM。除检查 `layout=SE` 外，pooling 前也会检查 BAM
  flag；检测到 paired reads 会停止。
- 多个 treatment BAM 通过内存管道 merge 后转换为一个 read-level BED6，不保存
  pooled BAM。
- pooled BED 只保留 `config/TAIR10.nuclear.chrom.sizes` 中的 Chr1–Chr5。

同一个 control 文件在多行重复填写时，只 pool 一次。文件身份按真实文件判断，所以
普通路径、指向同一文件的 symlink 和 hardlink 会被视作同一个文件。

> Deduplication of shared controls is performed by input-file identity, not genomic read coordinates.

两个独立 BAM 中坐标相同的 reads 都会先进入 pooled BED。随后 SICER 根据
`redundancy_threshold` 执行自身的原生坐标去冗余。

`config/contrasts.txt` 定义 `sicer_df` 的输入顺序：

```text
contrast  test  reference
lowRFR_vs_WL  WT_lowRFR  WT_WL
```

空表（只保留表头）表示只运行 `sicer`，不运行 `sicer_df`。

## 参数

主要参数在 `config/config.toml` 中：

```toml
[sicer2_call]
chrom_sizes = "config/TAIR10.nuclear.chrom.sizes"
redundancy_threshold = 1
window_size = 200
fragment_size = 150
effective_genome_fraction = 0.8
gap_size = 600
false_discovery_rate = 0.01
extra = ""

[sicer2_diff]
false_discovery_rate_df = 0.01
extra = ""
```

两个 FDR 是 SICER2 原生运行参数：

- `false_discovery_rate`：每个 library call enriched islands 时使用。
- `false_discovery_rate_df`：`sicer_df` 比较 union islands 时使用。

`sicer_df` 会同时收到这两个参数。workflow 不在 SICER2 运行结束后再次按 FDR 或
fold change 筛选。

`extra` 原样附加到命令末尾，可用于 `--e_value`、`--significant_reads` 或
`--verbose` 等未单列参数；不要在其中重复已有参数。

当前参数是 **Arabidopsis H2A.Z-informed hybrid preset**：

```text
公开 TAIR10 H2A.Z：1 / 200 / 100 / 0.8 / 600 / 0.01
当前第一版：       1 / 200 / 150 / 0.8 / 600 / 0.01
```

即 fragment size 暂用 SICER2 默认的 150 bp，其余参考 Arabidopsis H2A.Z 公开
设置。以后根据 fragment-size QC 再调整 TOML。

## 输出

workflow 只负责按 group 或 contrast 创建输出目录：

```text
results/
  sicer2/<group>/
  sicer2_diff/<contrast>/
  logs/
```

`sicer2/<group>/` 中是该 group 的 `sicer` 原生输出；
`sicer2_diff/<contrast>/` 中是该 contrast 的 `sicer_df` 原生输出。包括软件生成的
score islands、islands summary、FDR island BED、normalized WIG、union islands 和
increased/decreased islands 等文件。具体有哪些文件由是否提供 control、参数以及
SICER2 版本决定。

pooling 后的 BED basename 是：

```text
<group>.treatment.bed
<group>.control.bed
```

所以 SICER2 原生输出文件名会自然包含 group 名。workflow 不复制、不重命名，也不
解析这些结果文件。

`logs/` 只保存各步骤的标准输出和错误信息，不属于 SICER2 结果整理。

## 运行

从仓库根目录执行：

```bash
source /etc/profile.d/lmod.sh
module load Snakemake/9.27.0

snakemake \
  --snakefile workflow/Snakefile \
  --config toml=config/config.toml \
  --cores 10 \
  --dry-run

snakemake \
  --snakefile workflow/Snakefile \
  --config toml=config/config.toml \
  --cores 10
```

每个 SICER task 默认使用 5 个线程，对应 Chr1–Chr5。

## 实现说明

服务器 module 名为 `SICER2/2.1.1`，实际 Python package metadata 为
`SICER 2.1.0`。该 package 的 `sicer_df` 仍调用新版 SciPy 已移除的
`scipy.array`；wrapper 只在当前 Python 进程内将其映射为 `numpy.array`，不会修改
`/opt` 安装，也不会改变 SICER2 输出。

为了使用自定义 TAIR10 Chr1–Chr5 chromosome sizes，wrapper 在进程内注册该 genome
后直接调用 SICER2 的 `run_SICER` / `run_SICER_df`。除这一注册和 SciPy 兼容处理外，
结果生成仍由 SICER2 原生代码完成。

第一版不支持 PE。若以后支持 PE，需要单独实现和验证真实 fragment midpoint，不能
将完整 fragment BED 与 SE 的 `fragment_size/2` shift 混用。
