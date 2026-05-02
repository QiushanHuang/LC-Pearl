# LC-Pearl v1.0.0

[![English](https://img.shields.io/badge/Language-English-24292f)](README.en.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Release](https://img.shields.io/badge/release-v1.0.0-blue.svg)](RELEASE_NOTES_v1.0.0.md)

LC-Pearl 是一个可复现的分析 pipeline，用于定量分析 LAMMPS dump 轨迹中的液晶 mesogen 聚集行为。它面向单链或类链状液晶弹性体模拟体系：在这类体系中，Gay-Berne 吸引、局部取向有序性、链连接关系和外力拉伸相互竞争，并形成类似 pearl-necklace 的结构。

`v1.0.0` 发行版定义了三个层级：

- `mesogen contact`：type-1 ellipsoid 之间具有吸引性、且取向相关的 pair interaction。
- `domain`：由 contact graph 推断出的局部 mesogen 束，并分类为 weak 或 robust。
- `pearl`：由一个或多个 robust domain 组成的 3D 紧凑珠状 bead-like assembly。

pipeline 同时输出轴向力学量（`L_parallel`、`Rg_parallel`、`Rg_perp`、`S2_force`），使聚集行为可以和力、拉伸、温度、时间进行对照。

## 发布内容

本仓库包含完整的 LC-Pearl v1 pipeline：

- `scripts/` 中的 Python 分析脚本。
- 当前目录启动器 `lc_pearl_here.py`。
- TOML 配置驱动启动器 `lc_pearl_cli.py`。
- `templates/lc_pearl_preflight/` 中的 preflight 模板。
- `configs/` 中的快速配置文件。
- `docs/` 中的学术用户手册和算法手册。
- `tests/` 中的 2D lobe threshold prior 回归测试。

软件发行版本为 `v1.0.0`。当前 threshold-prior artifact schema 为 `schema_version = 5`，名称为 `LC Domain-Pearl V2 2D lobe streaming threshold prior`；这是内部算法 artifact schema 名称，不是 GitHub release 编号。

## 安装

使用 Python 3.11 或更高版本。最小运行依赖为 `numpy` 和 `matplotlib`。

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 -m pip install -e .
```

如果不安装 package，也可以在同一个 Python 环境中直接运行脚本：

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

## 推荐的当前目录工作流

进入包含 LAMMPS dump 文件的文件夹并运行：

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

首次使用时，LC-Pearl 会创建 `lc_pearl_preflight/` 并提示缺少哪些文件。请把需要的源文件放到这里：

- `lc_pearl_preflight/lammps/*.in` 或 `*.lmp`：原始 LAMMPS input，需包含 `pair_style gayberne` 和 type 1-1 `pair_coeff`。
- `lc_pearl_preflight/lammps/lammps_executable.txt`：一行 LAMMPS 可执行文件路径，例如 `/Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp`。
- 可选 `lc_pearl_preflight/topology/*.data` 或 `*.dat`：包含 `Atoms` 和 `Bonds` section 的 LAMMPS data 文件。LC-Pearl 会将其转换为 local-pair 和 excluded-pair 表。

然后运行：

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

`auto` 会按顺序执行：

1. 复用或创建 `lc_pearl_preflight/validation/verified_potential.json`。
2. 复用或创建 `lc_pearl_preflight/thresholds/global_thresholds.json`。
3. 在主分析前应用 2D lobe threshold prior。
4. 运行 domain、pearl、OVITO label 和 mechanics 分析。

## 配置文件驱动工作流

编辑 `configs/quick_run.toml`，然后预览命令：

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

运行 pipeline：

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py run --config configs/quick_run.toml
```

只做 validation：

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py validate --config configs/quick_validate.toml
```

## 主要输出

默认输出目录为 `lc_domain_pearl_v2_output/`。

- `aggregation_timeseries.tsv`：逐帧定量总表，包含 aggregation、domain、pearl 和 mechanics。
- `per_file/*_aggregation.tsv`：每个输入文件对应的 time series。
- `per_file/*_summary.json`：紧凑的每文件 summary。
- `per_file/*_lc_labels.dump`：OVITO 可读的 particle labels，启用 OVITO label 输出时生成。
- `per_file/*_lc_cluster_envelopes.dump`：visual cluster envelope particles。
- `per_file/*_lc_contact_edges.dump`：吸引 contact 的 local-style edge records。
- `per_file/*_lc_contact_segments.vtk`：OVITO 可加载的 contact line segments。
- `diagnostics/diagnostic_summary.json`：阈值、计数和 diagnostic provenance。
- `diagnostics/gb_strength_vs_p2.png`：当前主分析 gate 下 accepted-edge diagnostic。
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hist.png`：用于 threshold prior 的全 streaming candidate-pair distribution。
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`：同一全量 2D distribution 的 dot-grid 图，并带有选定阈值。
- `lc_pearl_preflight/thresholds/global_thresholds.json`：可复用的 threshold-prior artifact。

默认情况下，LC-Pearl 不写出完整的 `edge_diagnostics.tsv` 或完整逐帧 JSONL debug records，因为这些文件可能变成多 GB 输出。只有需要时才开启：

```toml
[analysis]
edge_diagnostics_table = "sample"  # off, sample, or full
edge_diagnostics_sample_size = 200000
write_frame_jsonl = true
```

## OVITO 可视化

在 OVITO 中打开 `per_file/*_lc_labels.dump`。常用 scalar fields 包括：

- `lc_cluster`：visual center-distance cluster id。
- `lc_cluster_size`：visual cluster 的大小。
- `lc_contact_degree`：与该 mesogen 相连的 accepted contact edge 数。
- `lc_min_pair_energy`：该 mesogen 参与的最强吸引 pair energy。
- `lc_mean_pair_energy`：该 mesogen 参与的平均吸引 pair energy。
- `lc_max_gb_strength`：该 mesogen 参与的最强 normalized attraction。
- `lc_mean_gb_strength`：该 mesogen 参与的平均 normalized attraction。
- `lc_domain`：robust 或 weak domain id。
- `lc_pearl`：3D pearl id。
- `lc_state`：紧凑状态码，用于表示 unassigned、weak-domain、robust-domain 和 pearl-supported particles。

对于 contact geometry，可把 `*_lc_contact_segments.vtk` 作为额外 OVITO pipeline 加载。如果没有 OVITO Pro，先检查 labels 和 envelopes，再对代表帧单独加载 line segments。

## 算法概述

LC-Pearl 使用 type-1 ellipsoids 作为 mesogen 成员。Type-2 spheres 和 type-3 anchors 不作为聚集成员，但可用于支持链序重建、端点定位、拉伸方向和 topology reconstruction。

对于每个候选 E-E pair，Gay-Berne pair energy 由 dump coordinates、quaternions、shape axes 和 LAMMPS `pair_style/pair_coeff` 源文件重建。归一化吸引强度定义为：

```text
gb_strength = max(0, -U_GB / U_well)
```

取向评分定义为：

```text
P2 = (3 |u_i dot u_j|^2 - 1) / 2
```

当前 threshold prior 从完整的 `gb_strength x P2` 二维分布中选取，而不是从已经 accepted 的 edges 中选取。它识别 high-P2 weak-contact lobe、high-P2 strong-attraction lobe，以及两者之间的 valley。该过程给出 `gb_on`；`gb_off` 是用于 gray/support contacts 的保守 left-lobe shoulder；`p2_cut` 是在 2D split 中一致使用的 orientation gate。

随后从 contact graph 构建 domains。Weak local domains 会被保留。Robust domains 需要尺寸以及取向、persistence、nonlocal support 和 parameter stability 等证据。Pearls 定义在更高一层：只有当 robust domains 的 3D gap、cross contacts、boundary support 和 bead-like aspect ratio 满足 pearl criteria 时，才会合并为同一个 pearl。

## 关键参数

- `contact_mode = "gayberne"`：使用重建的 Gay-Berne pair energy 作为 contact strength。
- `gb_on_strength`：threshold prior 给出的 strong edge 阈值。
- `gb_off_strength`：threshold prior 给出的 gray/support edge 阈值。
- `p2_cut`：pair orientation 阈值。
- `n_min`：robust-domain 尺寸证据所需的最小 mesogen 数。
- `s_excl`：被视为 local support 的 chain-sequence separation。
- `domain_min_lifetime`：作为 persistence evidence 的 processed-frame age。
- `pearl_gap_cut`：pearl merging 允许的最大 3D domain-domain gap。
- `pearl_min_cross_contacts`：pearl merging 所需的最小 cross-domain contacts。
- `pearl_min_boundary_particles`：两个 domain 上所需的最小 supported boundary particles。
- `pearl_max_aspect_ratio`：bead-like compactness 的 aspect ratio 上限。
- `threshold_prior.global_frame_stride`：threshold prior 使用全帧或按 stride 抽帧；`1` 表示所有 candidate frames。
- `workers = "auto"`：根据 CPU 数、输入规模以及环境变量上限 `LC_PEARL_MAX_WORKERS` 或 `LC_PEARL_MAX_AUTO_WORKERS` 自动使用并行 worker。

## Debug 与可靠性检查

建议按以下顺序检查：

1. 确认 potential cache：`lc_pearl_preflight/validation/verified_potential.json`。
2. 检查 `lc_pearl_preflight/thresholds/global_thresholds.json`。
3. 对比 `gb_strength_vs_p2_stream_hist.png` 和 `gb_strength_vs_p2_stream_lobe_split_dotgrid.png`。
4. 检查 `diagnostics/diagnostic_summary.json`，确认主分析实际应用的阈值。
5. 在 OVITO 中检查代表帧 labels。
6. 跨力或温度比较 `N_domain`、`N_pearl`、largest domain fraction、largest pearl fraction、`L_parallel` 和 `Rg_parallel/Rg_perp`。

如果 threshold cuts 看起来不合理，只删除 `lc_pearl_preflight/thresholds/global_thresholds.json` 并重新运行，即可重建 threshold prior。如果 potential 参数改变，则需要重新运行 validation 并重建 threshold prior。

## 文档

完整学术手册为：

- `docs/lc_pearl_user_guide.html`
- `docs/lc_pearl_algorithm_reference.html`

历史过程记录单独保存在 `docs/` 中，不作为当前算法参考。

## 发行版

本仓库已准备为 `LC-Pearl v1.0.0`，即 domain-pearl aggregation analysis pipeline 的第一个稳定发行版。
