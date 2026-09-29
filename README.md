# SICER2_workflow

一个只负责 SICER2 domain calling 和 differential islands 的简洁 Snakemake
workflow。输入必须是已经完成比对和必要过滤的 single-end BAM；本仓库不负责
FASTQ 质控、比对、MACS3、IDR 或 bigWig。

```text
processed SE BAM
  -> pooled treatment/control BED6
  -> SICER2 enriched domains
  -> SICER2 differential islands
  -> optional fold-change filtering
```

## 为什么单独建仓库

MACS3/IDR 主要围绕 peak、summit 和 replicate reproducibility；SICER2 使用
window/gap 将分散信号连接为 domains，并有自己的 differential-island 算法。两类
caller 共享的稳定接口是 processed BAM，没有必要把两套参数和输出揉进同一 DAG。

## 输入

`config/samples.tsv` 每行表示一个 biological replicate：

```text
group  sample  treatment_bam  control_bam  layout
WT_WL  WT_WL_rep1  /path/IP1.bam  /path/IgG_A.bam  SE
WT_WL  WT_WL_rep2  /path/IP2.bam  /path/IgG_A.bam  SE
```

- `control_bam` 可以是 Input DNA、IgG 或留空。
- 同一 group 必须全部有 control 或全部没有 control。
- 第一版只接受 `layout=SE`。除检查表格外，pooling 前还会检查 BAM flag；检测到
  paired reads 会立即终止。
- biological replicates 在内存管道中 merge 后转换成 read-level BED6，不保存大型
  pooled BAM。
- 只保留 `chrom_sizes` 中列出的 contigs；默认是 TAIR10 的 `Chr1`–`Chr5`。
- BED6 保留原始 strand；SICER 的 SE tag position 为正链 `start + 75`、负链
  `end - 1 - 75`（`fragment_size=150`）。

### Shared control 的去重边界

同一个 control BAM 可能在多个 sample 行中重复填写。例如两个 IP replicates 共用
同一个 IgG library。workflow 会根据解析后的真实文件身份去重：相同文件路径、指向
同一文件的 symlink 或 hardlink 只 pool 一次。

> Deduplication of shared controls is performed by input-file identity, not genomic read coordinates.

两个不同 control BAM 即使包含坐标相同的 reads，也都会进入 pooled BED。随后
SICER2 自身仍会按 `redundancy_threshold` 在 pooled tags 上处理相同
`(strand, start, end)` 的 reads；默认阈值 1 因而可能合并来自不同 libraries 的重合
tags。这是保留的 SICER 原生行为，不是 workflow 的文件去重。

`config/contrasts.tsv` 定义 differential 方向：

```text
contrast  test  reference
lowRFR_vs_WL  WT_lowRFR  WT_WL
```

空表（只保留表头）表示不运行 differential calling。

## 参数

TOML 字段沿用 SICER2 的原始参数名称。两个 FDR 必须区分：

- `[sicer2_call].false_discovery_rate` 对应 `--false_discovery_rate/-fdr`，用于每个
  condition 独立 call enriched islands。
- `[sicer2_diff].false_discovery_rate_df` 对应
  `--false_discovery_rate_df/-fdr_df`，用于两个 libraries 的 differential change。
- 运行 differential 时两个参数会同时传入。
- `min_fold_change` 是 SICER statistical significance 之后由 workflow 添加的
  effect-size filter，不是 SICER 原生统计参数。

当前配置是 **Arabidopsis H2A.Z-informed hybrid preset**：

```text
公开 TAIR10 H2A.Z：1 / 200 / 100 / 0.8 / 600 / 0.01
当前第一版：       1 / 200 / 150 / 0.8 / 600 / 0.01
                                 ^ fragment size 暂用 SICER2 默认值
```

这里依次为 redundancy、window、fragment、effective genome fraction、gap 和
call FDR。公开数据记录见
[GSM3908304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3908304)。
目标 PIF 论文只报告使用 TAIR10、WT IgG、合并 biological replicates 和 SICER，
没有公开这些数值参数：
[Willige et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9169284/)。

`extra` 只用于没有单列在 TOML 中的 SICER 原生选项，例如 `--e_value`、
`--significant_reads` 或 `--verbose`。workflow 管理的 treatment/control、chrom sizes、
window、gap、fragment、EGF、两个 FDR 和 CPU 等参数若在 `extra` 中重复，会在建 DAG
时直接报错，避免命令行出现两套相互冲突的值。

### TAIR10 genome size

`119667750` 是 UCSC/Ensembl TAIR10 七条序列的 golden-path length，包含五条核
染色体、叶绿体和线粒体。MACS3 的 `-g` 是背景模型的 genome size，而不是 contig
白名单。本流程不使用该 MACS3 数值；SICER2 从
`config/TAIR10.nuclear.chrom.sizes` 读取 Chr1–Chr5 的实际总长 `119146348`，再
应用 `effective_genome_fraction=0.8`。

服务器当前参考的 `ChrC=154478`、`ChrM=367808`，连同核染色体合计
`119668634`；这两个 organellar contigs 不进入本 workflow。

## Differential 输出定义

若 contrast 是 `test` 对 `reference`：

- `increased.fdr.tsv`：`FDR_test_vs_reference <= false_discovery_rate_df`
- `decreased.fdr.tsv`：`FDR_reference_vs_test <= false_discovery_rate_df`
- `increased.filtered.tsv`：再要求 `test/reference >= min_fold_change`
- `decreased.filtered.tsv`：再要求 `test/reference <= 1/min_fold_change`

`min_fold_change=1.0` 不增加 effect-size 门槛；设为 `1.3` 时，decreased 要求
`test/reference <= 0.7692`。FDR-only 与 FDR+FC 文件始终分别保存。

SICER2 `sicer_df` 比较 pooled libraries，不建模 biological-replicate variance；它
不能替代 DiffBind、csaw 或 DESeq2 类 replicate-aware differential analysis。

## 输出

```text
results/
  sicer2/<group>/
    <group>.islands.bed
    <group>.islands.summary.tsv
    raw/
  sicer2_diff/<contrast>/
    <contrast>.all_islands.tsv
    <contrast>.increased.fdr.tsv
    <contrast>.decreased.fdr.tsv
    <contrast>.increased.filtered.tsv
    <contrast>.decreased.filtered.tsv
    raw/
  run_manifest.tsv
  logs/
```

`raw/` 保留 SICER2 tool-native 文件；顶层文件提供稳定下游接口。manifest 记录每个
group 的真实输入文件、参数、module/package version 和结果数量。
有 control 时，`*.islands.bed` 是通过 call FDR 的 domains；
`*.islands.summary.tsv` 保留全部候选 domains 及其 read counts、p-value、fold
enrichment 和 FDR，便于日后改阈值而不重跑 SICER。无 control 时二者都来自
SICER 的 score-island 输出，summary 只增加稳定列名。
文件身份以 `device:inode:resolved_path` 记录；因此 manifest 同时展示用户提供文件最终
指向的对象，以及每个 pooled group 实际使用的唯一文件列表。differential 行明确记录
`test` 和 `reference`，并分别统计 all/FDR-only/FDR+FC 输出数量。

服务器当前 package metadata 是 `SICER 2.1.0`，module 名是 `SICER2/2.1.1`。
SICER 2.1.0 的 differential 源码仍调用已从新版 SciPy 移除的 `scipy.array`；wrapper
只在自己的 Python 进程内将它兼容映射到 `numpy.array`，不会修改 `/opt` 安装。

## 运行

先复制并修改示例 TOML、samples 和 contrasts，然后从仓库根目录执行：

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

每个 SICER task 默认最多使用 5 个线程，因为 TAIR10 nuclear 配置只有 5 条染色体。

开发测试可运行：

```bash
/opt/conda-envs/SICER2/2.1.1/bin/python -m unittest discover -s tests -v
```

## 可选 fragment-size sanity check

该检查不是 workflow 依赖，也不会自动修改 TOML：

```bash
source /etc/profile.d/lmod.sh
module load MACS3/3.0.4

macs3 predictd \
  -i treatment_rep1.bam treatment_rep2.bam \
  -f BAM \
  -g 119146348 \
  --outdir fragment_size_qc
```

H2A.Z 信号较宽，MACS3 model 可能失败或不稳定，所以这里只用于比较当前文库估计值
与候选的 100/150 bp，最终值必须人工确认后显式写回 TOML。

## PE 数据

第一版不支持 PE。以后若增加支持，需要独立验证 proper-pair 过滤、真实 fragment
midpoint、1-bp midpoint BED，以及普通 call 与 differential call 的一致性，不能把
完整 fragment BED 与 SE 的 `fragment_size/2` shift 混用。
