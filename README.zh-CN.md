# LC-Pearl v1.0.0 中文说明

[![English](https://img.shields.io/badge/Language-English-24292f)](README.en.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

LC-Pearl 是一个用于分析液晶 mesogen 聚集行为的开源 pipeline。它读取 LAMMPS dump 轨迹，重算或近似重建 type-1 ellipsoid 之间的 Gay-Berne pair attraction，并结合取向序、链拓扑、局部/非局域接触和 3D bead 几何，输出 domain、pearl 和力学响应的定量结果。

## 研究目标

LC-Pearl 解决的核心问题是：在外力拉伸和液晶吸引相互竞争时，如何把肉眼在 OVITO 中看到的“几个中等 cluster 挤在一起、又被部分拉开”的结构转化为可复现、可统计、可审查的数值指标。

它把聚集结构拆成三层：

- `mesogen contact`：两个液晶 ellipsoid 是否存在足够强的吸引和足够一致的取向。
- `domain`：若干 mesogen contact 是否组成局部液晶束，保留 weak local domain，并在满足尺寸、取向、时间、非局域支持和参数稳定性证据时升级为 robust domain。
- `pearl`：一个或多个 robust domain 是否在 3D 空间中组成同一个 bead-like 致密团块。

## 安装

推荐使用 Python 3.11 或更高版本，并安装 `numpy` 和 `matplotlib`。

```bash
cd /path/to/LC-Pearl
python3 -m pip install -e .
```

也可以不安装，直接使用绝对路径运行脚本。

## 最简工作流

进入包含 dump 文件的输出目录：

```bash
cd /path/to/dump_output_folder
python3 /path/to/LC-Pearl/lc_pearl_here.py
```

第一次运行会创建 `lc_pearl_preflight/`，并告诉你缺少哪些文件。通常需要准备：

- `lc_pearl_preflight/lammps/*.in` 或 `*.lmp`：包含 `pair_style gayberne` 和 type 1-1 `pair_coeff` 的 LAMMPS 输入文件。
- `lc_pearl_preflight/lammps/lammps_executable.txt`：一行 LAMMPS 可执行文件路径。
- `lc_pearl_preflight/topology/*.data` 或 `*.dat`：可选 LAMMPS data 文件，用于自动生成 local/exclude pair 表。

准备好后运行：

```bash
python3 /path/to/LC-Pearl/lc_pearl_here.py auto
```

## 输出结果

默认输出目录为 `lc_domain_pearl_v2_output/`。

- `aggregation_timeseries.tsv`：每帧聚集、domain、pearl 和力学量总表。
- `per_file/*_lc_labels.dump`：可在 OVITO 中打开的粒子 label dump。
- `per_file/*_lc_cluster_envelopes.dump`：cluster envelope 可视化辅助粒子。
- `per_file/*_lc_contact_segments.vtk`：contact line segments。
- `diagnostics/diagnostic_summary.json`：主分析实际使用的阈值和诊断摘要。
- `lc_pearl_preflight/validation/verified_potential.json`：GB potential 验证缓存。
- `lc_pearl_preflight/thresholds/global_thresholds.json`：可复用的全局阈值先验。
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`：2D lobe threshold prior 可视化。

## OVITO 字段

- `lc_cluster`：用于直观显示的中心距 cluster id，不等价于最终 robust domain。
- `lc_cluster_size`：visual cluster 中 mesogen 数量。
- `lc_contact_degree`：该 mesogen 参与的 accepted contact 数。
- `lc_min_pair_energy`：该 mesogen 参与的最强吸引 pair energy。
- `lc_mean_pair_energy`：该 mesogen 的平均 contact energy。
- `lc_max_gb_strength`：最大归一化 GB 吸引强度。
- `lc_mean_gb_strength`：平均归一化 GB 吸引强度。
- `lc_domain`：domain id。
- `lc_pearl`：pearl id。
- `lc_state`：未归属、weak、robust、pearl-supported 等状态码。

## 文档

- 用户手册：[docs/lc_pearl_user_guide.html](docs/lc_pearl_user_guide.html)
- 算法参考：[docs/lc_pearl_algorithm_reference.html](docs/lc_pearl_algorithm_reference.html)
- 发行说明：[RELEASE_NOTES_v1.0.0.md](RELEASE_NOTES_v1.0.0.md)

## 贡献者

Qiushan Huang / 黄秋山。

## 开源声明

LC-Pearl 以 MIT License 开源发布。你可以使用、修改、复制和再分发本软件，但必须保留版权和许可证声明。若用于学术论文、组会报告或项目报告，建议引用 `CITATION.cff`，并明确报告使用的 GB potential 验证状态、threshold prior、domain 定义和 pearl 定义。
