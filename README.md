# LC-Pearl v1.0.0

[![English](https://img.shields.io/badge/Language-English-24292f)](#english)
[![简体中文](https://img.shields.io/badge/语言-简体中文-1677ff)](#中文)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Release](https://img.shields.io/badge/release-v1.0.0-blue.svg)](RELEASE_NOTES_v1.0.0.md)

## English

LC-Pearl is a LAMMPS dump analysis pipeline for quantifying liquid-crystal mesogen aggregation under the competition between Gay-Berne attraction, local orientational order, chain connectivity, and external stretching.

This repository is the open-source `v1.0.0` release. See [README.en.md](README.en.md) for the full English overview, [docs/lc_pearl_algorithm_reference.html](docs/lc_pearl_algorithm_reference.html) for the algorithm and parameter reference, and [docs/lc_pearl_user_guide.html](docs/lc_pearl_user_guide.html) for the user manual.

Core hierarchy:

- `mesogen contact`: pairwise attractive and orientationally relevant contact between type-1 ellipsoids.
- `domain`: local mesogen bundle inferred from the contact graph, classified as weak local or robust.
- `pearl`: compact 3D bead-like aggregate composed of one or more robust domains.
- `mechanics`: auxiliary mechanical observables including `L_parallel`, `Rg_parallel`, `Rg_perp`, and `S2_force`.

Minimal workflow:

```bash
cd /path/to/dump_output_folder
python3 /path/to/LC-Pearl/lc_pearl_here.py
```

The first run creates `lc_pearl_preflight/` and reports missing inputs. Add the LAMMPS input, LAMMPS executable path, and optional topology data, then run:

```bash
python3 /path/to/LC-Pearl/lc_pearl_here.py auto
```

Contributor: Qiushan Huang.

Open-source notice: LC-Pearl is released under the MIT License. Please retain the copyright and license notice when using, modifying, or redistributing the code. For academic use, cite `CITATION.cff` and report the threshold prior, domain definition, and pearl definition used in the analysis.

## 中文

LC-Pearl 是一个面向 LAMMPS dump 轨迹的液晶 mesogen 聚集分析 pipeline，用于定量描述 Gay-Berne 吸引、局部取向、链连接和外力拉伸共同作用下形成的 domain-pearl 聚集结构。

本仓库为 `v1.0.0` 开源发行版。完整中文说明见 [README.zh-CN.md](README.zh-CN.md)，算法和参数细节见 [docs/lc_pearl_algorithm_reference.html](docs/lc_pearl_algorithm_reference.html)，使用手册见 [docs/lc_pearl_user_guide.html](docs/lc_pearl_user_guide.html)。

核心层级定义：

- `mesogen contact`：type-1 ellipsoid 之间的吸引和取向相关 pair contact。
- `domain`：由 contact graph 推断的局部 mesogen 束，分为 weak local domain 和 robust domain。
- `pearl`：一个或多个 robust domain 在 3D 空间中组成的 bead-like 致密团块。
- `mechanics`：`L_parallel`、`Rg_parallel`、`Rg_perp`、`S2_force` 等力学辅助量。

最简运行方式：

```bash
cd /path/to/dump_output_folder
python3 /path/to/LC-Pearl/lc_pearl_here.py
```

第一次运行会创建 `lc_pearl_preflight/`，请按提示放入 LAMMPS input、LAMMPS executable 路径和可选 topology data。之后运行：

```bash
python3 /path/to/LC-Pearl/lc_pearl_here.py auto
```

贡献者：Qiushan Huang / 黄秋山。

开源声明：LC-Pearl 以 MIT License 开源发布。使用、修改和再分发时请保留版权与许可证声明；用于学术工作时建议引用 `CITATION.cff` 并说明使用的 threshold prior、domain 和 pearl 定义。
