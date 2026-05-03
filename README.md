# LC-Pearl v2.1.0

[![English](https://img.shields.io/badge/Language-English-24292f)](#english)
[![简体中文](https://img.shields.io/badge/Language-%E7%AE%80%E4%BD%93%E4%B8%AD%E6%96%87-1677ff)](#%E4%B8%AD%E6%96%87)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Release](https://img.shields.io/badge/release-v2.1.0-blue.svg)](RELEASE_NOTES_v2.1.0.md)

## English

LC-Pearl is a reproducible analysis pipeline for quantifying liquid-crystal mesogen aggregation in LAMMPS dump trajectories. It was designed for single-chain or chain-like liquid-crystal simulations where Gay-Berne attraction, local orientational order, chain connectivity, temperature, and external stretching compete to form pearl-necklace-like structures.

The v2.1.0 release defines a three-level physical hierarchy plus two diagnostic core-contact tiers:

- `mesogen contact`: a pairwise attractive and orientationally relevant type-1 ellipsoid interaction.
- `domain`: a local mesogen bundle inferred from the contact graph and classified as weak or robust.
- `pearl`: a compact 3D bead-like assembly of one or more robust domains.
- `core shoulder`: a diagnostic contact tier with high pair alignment, typically `P2 > 0.71`, and `gb_strength >= gb_core_strength`.
- `strict core`: a stricter diagnostic contact tier with `P2 > 0.80` and `gb_strength >= gb_strict_core_strength`.

The main domain and pearl definitions are not replaced by the core tiers. The core tiers exist for visualization, audit, and physical interpretation of the strongest contact regions.

### What Is Released

This repository contains the complete LC-Pearl v2.1.0 pipeline:

- Python analysis scripts in `scripts/`.
- Current-directory launcher `lc_pearl_here.py`.
- TOML-driven launcher `lc_pearl_cli.py`.
- Preflight template in `templates/lc_pearl_preflight/`.
- Quick configuration files in `configs/`.
- Academic manuals in Markdown under `MANUAL.md` and `docs/`.
- Regression tests for the 2D lobe threshold prior and the per-frame domain-size count diagnostic.
- Source release archives under `releases/`.

The software release is `v2.1.0`. The threshold-prior artifact schema is `schema_version = 7`, named `LC-Pearl 2.1.0 core-tier streaming threshold prior`. This schema name identifies the threshold-prior artifact format and algorithm family; it is not a replacement for the GitHub release tag.

### Core V2 Algorithm

The V2 algorithm is centered on a 2D `GB strength x P2` split. It is not a pure distance cutoff, not a pure orientation cutoff, and not a post-run threshold suggestion. The threshold prior is built before the main analysis and then applied to the contact graph.

```text
validated Gay-Berne reconstruction
  -> streaming 2D GB-strength x P2 prior
  -> contact tiers: none / weak-transition / aggregation / core shoulder / strict core
  -> weak and robust domains
  -> compact 3D pearls
  -> mechanics, OVITO labels, and diagnostics
```

The tier logic is:

| Tier | Pair/particle meaning | Rule |
|---:|---|---|
| 0 | no accepted contact | below gray/support gate |
| 1 | weak or transition support | `gb_strength >= gb_off_strength` and `P2 >= p2_cut`, but not strong |
| 2 | aggregation/strong contact | `gb_strength >= gb_on_strength` and `P2 >= p2_cut` |
| 3 | core shoulder | `gb_strength >= gb_core_strength` and `P2 > 0.71` |
| 4 | strict core | `gb_strength >= gb_strict_core_strength` and `P2 > 0.80` |

`lc_aggregation_tier` is the OVITO particle field that reports the maximum incident contact tier for each mesogen. It is the recommended first color-coding field because it shows unaggregated, weak/transition, aggregated, core-shoulder, and strict-core mesogens in one scalar.

### Scientific Motivation

In stretched liquid-crystal chain simulations, a visually large cluster can be a compact aggregate, a set of nearby medium domains, or a partially stretched pearl-necklace state. A simple distance cluster or a pure orientation cluster cannot reliably separate these cases. LC-Pearl therefore separates the analysis into levels:

- Pair contacts measure direct mesogen-mesogen energetic and orientational coupling.
- Domains measure local bundles of mesogens.
- Pearls measure 3D bead-like assemblies of robust domains.
- Mechanics connects these objects to extension and chain shape.

This is intended to support analysis of force-aggregation competition, cooling-induced aggregation, pearl-necklace morphology, and OVITO-assisted inspection.

### Installation

Use Python 3.11 or newer. Runtime dependencies are `numpy` and `matplotlib`.

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 -m pip install -e .
```

If you do not install the package, run scripts directly with the same Python environment:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

### Recommended Current-Directory Workflow

Enter a directory that contains LAMMPS dump files:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

On first use, LC-Pearl creates `lc_pearl_preflight/` and tells you what is missing. Add the source files required by your project:

- `lc_pearl_preflight/lammps/*.in` or `*.lmp`: original LAMMPS input containing `pair_style gayberne` and the type 1-1 `pair_coeff`.
- `lc_pearl_preflight/lammps/lammps_executable.txt`: one line with the LAMMPS executable path, for example `/Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp`.
- Optional `lc_pearl_preflight/topology/*.data` or `*.dat`: LAMMPS data file with `Atoms` and `Bonds` sections. LC-Pearl converts it into local-pair and excluded-pair tables.

Then run the full automatic workflow:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

`auto` performs the following sequence:

1. Reuse or create `lc_pearl_preflight/validation/verified_potential.json`.
2. Reuse or create `lc_pearl_preflight/thresholds/global_thresholds.json`.
3. Apply the threshold prior before the main analysis.
4. Run contact, domain, pearl, OVITO-label, diagnostics, and mechanics analysis.

### Config-Driven Workflow

Edit `configs/quick_run.toml`, then preview the generated command:

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

Run the pipeline:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py run --config configs/quick_run.toml
```

Validate only:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py validate --config configs/quick_validate.toml
```

Use TOML configs when you want reproducible batch analysis across force, temperature, density, chain length, or sampling conditions.

### Preflight Folder

The preflight folder is designed to be portable. If two simulation folders share the same potential and topology logic, you can copy a validated `lc_pearl_preflight/` folder to the new dump folder and rerun LC-Pearl.

| Path | Required | Purpose |
|---|---:|---|
| `lc_pearl_preflight/lammps/*.in` or `*.lmp` | yes | Source for Gay-Berne pair style and coefficients |
| `lc_pearl_preflight/lammps/gb_param_source.in` | generated | Normalized input used by the analysis |
| `lc_pearl_preflight/lammps/lammps_executable.txt` | recommended | LAMMPS executable used by validation |
| `lc_pearl_preflight/topology/*.data` | optional but recommended | Source topology for local/excluded pair tables |
| `lc_pearl_preflight/topology/local_pairs.tsv` | generated | Local chain-neighbor support pairs |
| `lc_pearl_preflight/topology/exclude_pairs.tsv` | generated | Pairs excluded from aggregation/contact consideration |
| `lc_pearl_preflight/validation/verified_potential.json` | generated | Potential-validation cache |
| `lc_pearl_preflight/thresholds/global_thresholds.json` | generated | Reusable threshold-prior artifact |

### Potential Validation

In `contact_mode = "gayberne"`, LC-Pearl reconstructs type-1 ellipsoid pair energies from positions, quaternions, shape axes, and the LAMMPS Gay-Berne parameters. The normalized attraction is

```text
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
```

where `U_GB,ij` is the reconstructed pair energy and `U_well,ij` is the orientation-dependent well depth. Normalization makes attraction strength comparable across different orientations.

The validation cache is:

```text
lc_pearl_preflight/validation/verified_potential.json
```

Revalidate when:

- `pair_style`, `pair_coeff`, units, shape axes, or type mapping changes.
- A new LAMMPS build changes aspherical pair behavior.
- You suspect the cache came from another input deck.

You usually do not need to revalidate only because force, stretch length, or temperature changes, provided the potential settings are identical.

### Threshold Prior

The threshold prior is computed before the main analysis. It uses the candidate-pair 2D distribution of normalized Gay-Berne attraction and pair alignment:

```text
q_ij = |u_i dot u_j|
P2_ij = (3 q_ij^2 - 1) / 2
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
```

The prior estimates:

| Parameter | Source | Role |
|---|---|---|
| `p2_cut` | split between low-orientation background and high-orientation lobe | main pair-orientation gate |
| `gb_off_strength` | conservative weak-contact shoulder | gray/support edge threshold |
| `gb_on_strength` | valley between weak-contact and strong-attraction high-P2 lobes | strong edge threshold |
| `p2_core_cut = 0.71` | fixed high-alignment diagnostic gate | core-shoulder P2 gate |
| `gb_core_strength` | strong-lobe left shoulder | core-shoulder GB gate |
| `p2_strict_core_cut = 0.80` | stricter high-alignment diagnostic gate | strict-core P2 gate |
| `gb_strict_core_strength` | high-strength weighted quantile | strict-core GB gate |

The most important threshold-prior outputs are:

- `lc_pearl_preflight/thresholds/global_thresholds.json`
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hist.png`
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_dotgrid.png`
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hexbin.png`
- `lc_pearl_preflight/thresholds/lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`
- `lc_pearl_preflight/thresholds/lobe_split_preview/gb_core_slice_hist.png`

For final production parameters, prefer a representative trajectory that covers the relevant morphology range. Use `threshold_prior.global_frame_stride = 1` for full frame sampling, or larger strides such as `10` or `100` for quick checks.

### Domain And Pearl Algorithm

LC-Pearl builds a contact graph using two edge classes:

```text
strong edge: gb_strength >= gb_on_strength  and P2 >= p2_cut
gray edge:   gb_strength >= gb_off_strength and P2 >= p2_cut
```

Domains are local mesogen bundles inferred from this graph.

- Weak local domains are retained. They preserve adjacent or weakly supported local bundling.
- Robust domains require evidence beyond simple local adjacency.
- Robust evidence can include size, internal orientation, persistence, nonlocal support, and parameter stability.

Pearls are defined one level above domains. A pearl is a compact 3D bead-like assembly of one or more robust domains. Two robust domains merge into one pearl only when the domain-domain gap, cross contacts, boundary support, and bead-like aspect ratio support the merge.

This distinction is important:

- `domain` answers whether a local liquid-crystal bundle exists.
- `pearl` answers whether those bundles form the same 3D bead.
- axial segmentation describes how beads and connectors arrange along the chain or stretch axis.

### Mechanics Layer

LC-Pearl reports mechanics observables so aggregation can be compared with force, extension, temperature, and time.

| Quantity | Meaning | Use |
|---|---|---|
| `L_parallel` | projected end-to-end length along the stretch axis | primary force-extension observable |
| `Rg_parallel` | parallel radius-of-gyration component | chain-shape response |
| `Rg_perp` | perpendicular radius-of-gyration component | lateral compaction or swelling |
| `S2_force` | mesogen alignment relative to the stretch axis | stretch-induced orientational order |

For multiple constant-force simulations, compare steady-state averages and, when appropriate, estimate

```text
C_parallel = Delta <L_parallel> / Delta F
```

`Rg_parallel` and `Rg_perp` are supporting shape descriptors. They should not replace `L_parallel` as the primary extension measure.

### Main Outputs

Output is written by default to `lc_domain_pearl_v2_output/`.

| Output | Purpose |
|---|---|
| `aggregation_timeseries.tsv` | frame-level quantitative table for contacts, domains, pearls, and mechanics |
| `per_file/*_aggregation.tsv` | per-input-file time series |
| `per_file/*_summary.json` | compact per-file summary |
| `per_file/*_lc_labels.dump` | OVITO-readable particle labels |
| `per_file/*_lc_cluster_envelopes.dump` | visual cluster envelope particles |
| `per_file/*_lc_contact_edges.dump` | contact-edge records |
| `per_file/*_lc_contact_segments.vtk` | OVITO-loadable contact line segments |
| `diagnostics/domain_diagnostics.tsv` | domain evidence, size, S2, and classification |
| `diagnostics/domain_size_vs_s2.png` | scatter of domain size and internal S2 |
| `diagnostics/domain_size_frame_counts.tsv` | per-frame and per-size domain-count table |
| `diagnostics/domain_size_vs_domain_count.png` | scatter where each point is one `(frame, domain size)` group |
| `diagnostics/pearl_candidate_diagnostics.tsv` | domain-domain pearl merge candidate table |
| `diagnostics/diagnostic_summary.json` | applied thresholds, counts, and provenance |

Large debug tables are disabled by default:

```toml
[analysis]
edge_diagnostics_table = "off"  # off, sample, or full
write_frame_jsonl = false
accepted_edge_audit = false
```

Enable them only for targeted debugging.

### OVITO Visualization

Open `per_file/*_lc_labels.dump` in OVITO. Useful scalar fields include:

| Field | Meaning | Recommended use |
|---|---|---|
| `lc_cluster` | center-distance visual cluster id | fast spatial cluster view |
| `lc_cluster_size` | size of visual cluster | identify compact visual groups |
| `lc_aggregation_tier` | `0` none, `1` weak/transition, `2` aggregation, `3` core shoulder, `4` strict core | recommended first color-coding field |
| `lc_contact_degree` | number of accepted contacts incident on a mesogen | identify contact hubs |
| `lc_min_pair_energy` | most attractive pair energy involving a mesogen | locate strong attraction |
| `lc_mean_pair_energy` | mean pair energy around a mesogen | review local energetic environment |
| `lc_max_gb_strength` | strongest normalized attraction involving a mesogen | locate strong GB contacts |
| `lc_mean_gb_strength` | mean normalized attraction | distinguish isolated strong contacts from broader attraction |
| `lc_core_contact_degree` | number of core-shoulder contacts | inspect `P2 > 0.71` core regions |
| `lc_core_nonlocal_degree` | number of nonlocal core-shoulder contacts | inspect nonlocal support |
| `lc_is_core_particle` | has at least one core-shoulder contact | binary core coloring |
| `lc_strict_core_contact_degree` | number of strict-core contacts | inspect strongest core regions |
| `lc_is_strict_core_particle` | has at least one strict-core contact | binary strict-core coloring |
| `lc_domain` | weak or robust domain id | local bundle membership |
| `lc_pearl` | 3D pearl id | bead-like pearl membership |
| `lc_state` | `0` none, `1` weak domain, `2` robust domain | compact domain-state view |

For contact geometry, load `*_lc_contact_segments.vtk` as an additional OVITO pipeline. If OVITO Pro is unavailable, inspect label dumps first and open representative VTK files separately.

### How To Analyze Results

For a single trajectory:

- plot `robust_domain_count`, `weak_domain_count`, `pearl_count`, and `largest_pearl_fraction`;
- compare `local_edge_fraction` and `nonlocal_edge_fraction`;
- inspect `domain_size_vs_s2.png` for domain quality;
- inspect `domain_size_vs_domain_count.png` for how many same-sized domains appear in each frame;
- validate representative frames in OVITO.

For a force series:

- reuse the same potential cache if the potential is unchanged;
- reuse the same threshold prior only if it was built from representative morphology and temperature ranges;
- compare steady-state `L_parallel`, `pearl_count`, `largest_pearl_fraction`, and nonlocal support;
- report uncertainty across windows or replicate simulations.

For a cooling series:

- build the threshold prior from a trajectory covering the full high-to-low temperature morphology range when possible;
- inspect the 2D lobe split figures;
- compare aggregation onset using robust domain count, pearl count, and core-contact fraction.

### Key Parameters

| Parameter | Default or source | Purpose |
|---|---|---|
| `contact_mode` | `gayberne` | use reconstructed Gay-Berne pair energy |
| `mesogen_type` | `1` | aggregation members |
| `anchor_types` | `3` | endpoint and stretch-axis support |
| `gb_off_strength` | threshold prior | gray/support edge threshold |
| `gb_on_strength` | threshold prior | strong edge threshold |
| `p2_cut` | threshold prior | pair-orientation gate |
| `gb_core_strength` | threshold prior | diagnostic core-shoulder GB gate |
| `p2_core_cut` | `0.71` | diagnostic core-shoulder P2 gate |
| `gb_strict_core_strength` | threshold prior | diagnostic strict-core GB gate |
| `p2_strict_core_cut` | `0.80` | diagnostic strict-core P2 gate |
| `s_excl` | `1` | local sequence support range |
| `n_min` | `3` | robust-domain size evidence |
| `robust_min_s2` | `0.70` | robust-domain orientation evidence |
| `domain_min_lifetime` | `2` | persistence evidence |
| `pearl_gap_cut` | `auto` | maximum domain-domain 3D gap for pearl merging |
| `pearl_min_cross_contacts` | `2` | minimum cross-domain contacts |
| `pearl_min_boundary_particles` | `2` | minimum supported boundary particles |
| `pearl_max_aspect_ratio` | `3.0` | bead-like compactness limit |
| `threshold_prior.global_frame_stride` | `1` | threshold-prior frame sampling stride |
| `workers` | `auto` | parallelism capped by CPU, file count, and environment limits |

### Debugging And Reliability Checks

Use this order:

1. Confirm `lc_pearl_preflight/validation/verified_potential.json`.
2. Inspect `lc_pearl_preflight/thresholds/global_thresholds.json`.
3. Compare `gb_strength_vs_p2_stream_hist.png`, `gb_strength_vs_p2_stream_dotgrid.png`, and `gb_strength_vs_p2_stream_hexbin.png`.
4. Check `lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`.
5. Check `diagnostics/diagnostic_summary.json` for actual applied thresholds.
6. Inspect OVITO labels for representative frames.
7. Compare `N_domain`, `N_pearl`, largest domain fraction, largest pearl fraction, `L_parallel`, and `Rg_parallel/Rg_perp`.

If threshold cuts look wrong, rebuild the threshold prior by deleting only `lc_pearl_preflight/thresholds/global_thresholds.json` and rerunning. If potential parameters change, rerun validation and rebuild the threshold prior.

### Documentation

The full manuals are:

- [MANUAL.md](MANUAL.md)
- [docs/lc_pearl_algorithm_details.md](docs/lc_pearl_algorithm_details.md)
- [docs/lc_pearl_user_guide.md](docs/lc_pearl_user_guide.md)
- [docs/lc_pearl_operation_manual.md](docs/lc_pearl_operation_manual.md)
- [docs/LC_Pearl_v2.1_Academic_User_Manual.md](docs/LC_Pearl_v2.1_Academic_User_Manual.md)

Historical process notes and browser HTML copies are kept separately in `docs/`, but the Markdown documents above are the canonical v2.1 documentation.

### Release

This repository is prepared as `LC-Pearl v2.1.0`. Previous source snapshots are preserved under `releases/` for reproducibility.

## 中文

LC-Pearl 是一个可复现的分析 pipeline，用于定量分析 LAMMPS dump 轨迹中的液晶 mesogen 聚集行为。它面向单链或类链状液晶模拟体系：在这类体系中，Gay-Berne 吸引、局部取向有序性、链连接关系、温度变化和外力拉伸相互竞争，并形成类似 pearl-necklace 的结构。

v2.1.0 发行版定义了三个物理层级，并增加两个诊断性 core-contact 层：

- `mesogen contact`：type-1 ellipsoid 之间具有吸引性、且取向相关的 pair interaction。
- `domain`：由 contact graph 推断出的局部 mesogen 束，并分类为 weak 或 robust。
- `pearl`：由一个或多个 robust domain 组成的 3D 紧凑珠状 bead-like assembly。
- `core shoulder`：满足高取向条件的诊断性接触层，通常为 `P2 > 0.71` 且 `gb_strength >= gb_core_strength`。
- `strict core`：更严格的诊断性接触层，满足 `P2 > 0.80` 且 `gb_strength >= gb_strict_core_strength`。

core 层只用于可视化、审查和物理解释，不替代主 domain/pearl 判定。

### 发布内容

本仓库包含完整的 LC-Pearl v2.1.0 pipeline：

- `scripts/` 中的 Python 分析脚本。
- 当前目录启动器 `lc_pearl_here.py`。
- TOML 配置驱动启动器 `lc_pearl_cli.py`。
- `templates/lc_pearl_preflight/` 中的 preflight 模板。
- `configs/` 中的快速配置文件。
- `MANUAL.md` 和 `docs/` 中的 Markdown 学术手册。
- 2D lobe threshold prior 和 per-frame domain-size count 的回归测试。
- `releases/` 中的源码快照。

软件发行版本为 `v2.1.0`。当前 threshold-prior artifact schema 为 `schema_version = 7`，名称为 `LC-Pearl 2.1.0 core-tier streaming threshold prior`。这是阈值先验 artifact 的算法格式名称，不是 GitHub release tag 的替代。

### V2 核心算法

V2 的核心是 `GB strength x P2` 二维切割。它不是纯距离 cutoff，不是纯取向 cutoff，也不是主分析之后才给出的后验建议。阈值 prior 在主分析之前生成，并直接用于 contact graph。

```text
已验证 Gay-Berne 重算
  -> streaming 2D GB-strength x P2 prior
  -> contact tiers: none / weak-transition / aggregation / core shoulder / strict core
  -> weak 和 robust domains
  -> compact 3D pearls
  -> mechanics、OVITO labels 和 diagnostics
```

分层规则为：

| Tier | 含义 | 规则 |
|---:|---|---|
| 0 | 无 accepted contact | 低于 gray/support gate |
| 1 | weak 或 transition support | `gb_strength >= gb_off_strength` 且 `P2 >= p2_cut`，但未达到 strong |
| 2 | aggregation/strong contact | `gb_strength >= gb_on_strength` 且 `P2 >= p2_cut` |
| 3 | core shoulder | `gb_strength >= gb_core_strength` 且 `P2 > 0.71` |
| 4 | strict core | `gb_strength >= gb_strict_core_strength` 且 `P2 > 0.80` |

`lc_aggregation_tier` 是 OVITO 中的粒子字段，表示每个 mesogen 的 incident contacts 中最高的 contact tier。它是最推荐的第一优先上色字段，因为它可以同时显示未聚集、weak/transition、aggregation、core shoulder 和 strict core。

### 科学目的

在拉伸液晶链模拟中，一个肉眼看到的大 cluster 可能是真正紧凑聚集体，也可能是几个中等 domain 彼此贴近，或者是被拉伸后的 pearl-necklace 中间态。单纯距离聚类或单纯取向聚类都很难稳定区分这些情况。LC-Pearl 因此分层定义：

- pair contact 衡量两个 mesogen 是否有直接能量和取向耦合；
- domain 衡量局部 mesogen 束；
- pearl 衡量一个或多个 robust domains 是否构成同一个 3D bead；
- mechanics 将聚集行为与拉伸长度、形状和取向响应联系起来。

### 安装

使用 Python 3.11 或更高版本。运行依赖为 `numpy` 和 `matplotlib`。

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 -m pip install -e .
```

如果不安装 package，也可以在同一个 Python 环境中直接运行脚本：

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

### 推荐的当前目录工作流

进入包含 LAMMPS dump 文件的目录：

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

首次运行时，LC-Pearl 会创建 `lc_pearl_preflight/` 并提示缺少哪些文件。请放入：

- `lc_pearl_preflight/lammps/*.in` 或 `*.lmp`：原始 LAMMPS input，需包含 `pair_style gayberne` 和 type 1-1 `pair_coeff`。
- `lc_pearl_preflight/lammps/lammps_executable.txt`：一行 LAMMPS 可执行文件路径，例如 `/Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp`。
- 可选 `lc_pearl_preflight/topology/*.data` 或 `*.dat`：包含 `Atoms` 和 `Bonds` section 的 LAMMPS data 文件。

然后运行完整自动流程：

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

`auto` 会依次执行：

1. 复用或生成 `lc_pearl_preflight/validation/verified_potential.json`。
2. 复用或生成 `lc_pearl_preflight/thresholds/global_thresholds.json`。
3. 在主分析前应用阈值先验。
4. 运行 contact、domain、pearl、OVITO label、diagnostics 和 mechanics 分析。

### 配置文件工作流

编辑 `configs/quick_run.toml`，然后先预览命令：

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

运行 pipeline：

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py run --config configs/quick_run.toml
```

只做验证：

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py validate --config configs/quick_validate.toml
```

当你需要跨多个力、温度、密度或链长目录重复分析时，推荐用 TOML 配置文件固定参数。

### 前置文件夹

`lc_pearl_preflight/` 被设计为可迁移文件夹。若两个模拟目录使用相同 potential 和 topology 逻辑，可以把已经验证过的 preflight 文件夹复制到新 dump 目录中复用。

| 路径 | 是否必要 | 用途 |
|---|---:|---|
| `lc_pearl_preflight/lammps/*.in` 或 `*.lmp` | 必要 | 提供 Gay-Berne pair style 和 pair coefficients |
| `lc_pearl_preflight/lammps/gb_param_source.in` | 自动生成 | 分析使用的规范化 input |
| `lc_pearl_preflight/lammps/lammps_executable.txt` | 推荐 | 用于 LAMMPS 验证 |
| `lc_pearl_preflight/topology/*.data` | 可选但推荐 | 生成 local/excluded pair 表 |
| `lc_pearl_preflight/topology/local_pairs.tsv` | 自动生成 | 局域链邻接支持 pair |
| `lc_pearl_preflight/topology/exclude_pairs.tsv` | 自动生成 | 排除 pair 表 |
| `lc_pearl_preflight/validation/verified_potential.json` | 自动生成 | potential 验证缓存 |
| `lc_pearl_preflight/thresholds/global_thresholds.json` | 自动生成 | 可复用阈值先验 |

### Potential 验证

在 `contact_mode = "gayberne"` 下，LC-Pearl 从 position、quaternion、shape axes 和 LAMMPS Gay-Berne 参数重算 type-1 ellipsoid pair energy。归一化吸引强度为：

```text
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
```

其中 `U_GB,ij` 是重算 pair energy，`U_well,ij` 是取向相关势阱深度。归一化后，不同取向的吸引强弱可以放到同一尺度比较。

验证缓存为：

```text
lc_pearl_preflight/validation/verified_potential.json
```

以下情况需要重新验证：

- `pair_style`、`pair_coeff`、units、shape axes 或 type mapping 改变；
- LAMMPS asphere pair 行为可能因编译版本改变；
- 怀疑当前 cache 来自其他 input deck。

如果只是外力、拉伸长度或温度改变，而 potential 完全相同，通常不需要重新验证。

### 阈值先验

阈值先验在主分析前计算。它使用 candidate pair 的二维分布：

```text
q_ij = |u_i dot u_j|
P2_ij = (3 q_ij^2 - 1) / 2
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
```

阈值先验估计以下参数：

| 参数 | 来源 | 作用 |
|---|---|---|
| `p2_cut` | 低取向背景与高取向 lobe 的分割 | 主 pair 取向门槛 |
| `gb_off_strength` | 弱接触 lobe 的保守 shoulder | gray/support edge |
| `gb_on_strength` | 高 P2 条件下弱接触 lobe 与强吸引 lobe 的 valley | strong edge |
| `p2_core_cut = 0.71` | 固定高取向诊断门槛 | core-shoulder P2 gate |
| `gb_core_strength` | 强接触 lobe 的 left shoulder | core-shoulder GB gate |
| `p2_strict_core_cut = 0.80` | 更严格高取向诊断门槛 | strict-core P2 gate |
| `gb_strict_core_strength` | 高强度 weighted quantile | strict-core GB gate |

最重要的阈值先验输出为：

- `lc_pearl_preflight/thresholds/global_thresholds.json`
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hist.png`
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_dotgrid.png`
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hexbin.png`
- `lc_pearl_preflight/thresholds/lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`
- `lc_pearl_preflight/thresholds/lobe_split_preview/gb_core_slice_hist.png`

最终生产参数建议使用覆盖主要构型变化的代表性轨迹生成。`threshold_prior.global_frame_stride = 1` 表示全帧采样，`10` 或 `100` 可用于快速检查。

### Domain 和 Pearl 算法

LC-Pearl 用两类边构造 contact graph：

```text
strong edge: gb_strength >= gb_on_strength  and P2 >= p2_cut
gray edge:   gb_strength >= gb_off_strength and P2 >= p2_cut
```

domain 是由 contact graph 推断出的局部 mesogen 束。

- weak local domain 会被保留，用于描述相邻或弱支持的局部束化。
- robust domain 需要额外证据，不能只靠局域相邻接触自动成立。
- robust evidence 包括 size、内部取向、persistence、nonlocal support 和 parameter stability。

pearl 是 domain 之上的层级。pearl 表示一个或多个 robust domains 是否在 3D 空间中组成同一个紧凑 bead-like 珠子。两个 robust domains 只有在 3D gap、cross contacts、boundary support 和 aspect ratio 都满足条件时才合并成一个 pearl。

因此：

- `domain` 回答有没有局部液晶束；
- `pearl` 回答这些束是否构成同一个 3D 珠子；
- axial segmentation 描述珠子和 connector 沿链或拉伸方向如何排列。

### Mechanics 层

LC-Pearl 输出力学相关量，使聚集可以和外力、拉伸、温度、时间对照。

| 量 | 含义 | 用途 |
|---|---|---|
| `L_parallel` | 沿拉伸轴投影端到端长度 | 主 force-extension 量 |
| `Rg_parallel` | 沿拉伸轴回转半径分量 | 链形状响应 |
| `Rg_perp` | 垂直拉伸轴回转半径分量 | 横向聚集或膨胀 |
| `S2_force` | mesogen 相对拉伸轴取向序 | 拉伸诱导取向 |

多个恒力模拟可以比较稳态均值，并在合适时估计：

```text
C_parallel = Delta <L_parallel> / Delta F
```

`Rg_parallel` 和 `Rg_perp` 是辅助形状量，不应替代 `L_parallel` 作为主拉伸长度指标。

### 主要输出

默认输出目录为 `lc_domain_pearl_v2_output/`。

| 输出 | 用途 |
|---|---|
| `aggregation_timeseries.tsv` | 逐帧 contact、domain、pearl、mechanics 主表 |
| `per_file/*_aggregation.tsv` | 每个输入文件的时间序列 |
| `per_file/*_summary.json` | compact 文件摘要 |
| `per_file/*_lc_labels.dump` | OVITO 粒子 label |
| `per_file/*_lc_cluster_envelopes.dump` | visual cluster envelope 粒子 |
| `per_file/*_lc_contact_edges.dump` | contact edge 记录 |
| `per_file/*_lc_contact_segments.vtk` | OVITO contact 线段 |
| `diagnostics/domain_diagnostics.tsv` | domain evidence、size、S2 和 classification |
| `diagnostics/domain_size_vs_s2.png` | domain size 与内部 S2 散点图 |
| `diagnostics/domain_size_frame_counts.tsv` | 每帧、每种 size 的 domain 数量表 |
| `diagnostics/domain_size_vs_domain_count.png` | 每个点代表一个 `(frame, domain size)` 分组的散点图 |
| `diagnostics/pearl_candidate_diagnostics.tsv` | domain-domain pearl merge candidate 表 |
| `diagnostics/diagnostic_summary.json` | 实际阈值、计数和 provenance |

大型 debug 表默认关闭：

```toml
[analysis]
edge_diagnostics_table = "off"  # off, sample, or full
write_frame_jsonl = false
accepted_edge_audit = false
```

只有在定点 debug 时才建议打开。

### OVITO 可视化

在 OVITO 中打开 `per_file/*_lc_labels.dump`。常用字段如下：

| 字段 | 含义 | 建议用途 |
|---|---|---|
| `lc_cluster` | 中心距 visual cluster id | 快速看空间 cluster |
| `lc_cluster_size` | visual cluster size | 看紧凑视觉团块 |
| `lc_aggregation_tier` | `0` none, `1` weak/transition, `2` aggregation, `3` core shoulder, `4` strict core | 推荐第一优先上色字段 |
| `lc_contact_degree` | mesogen 参与的 accepted contact 数 | 找 contact hub |
| `lc_min_pair_energy` | 该 mesogen 最强吸引 pair energy | 找强吸引区域 |
| `lc_mean_pair_energy` | 周围平均 pair energy | 看局部能量环境 |
| `lc_max_gb_strength` | 最大归一化吸引强度 | 找强 GB contact |
| `lc_mean_gb_strength` | 平均归一化吸引强度 | 区分孤立强边和整体强吸引 |
| `lc_core_contact_degree` | core-shoulder contact 数 | 看 `P2 > 0.71` 核心肩部 |
| `lc_core_nonlocal_degree` | 非局域 core-shoulder contact 数 | 看非局域支持 |
| `lc_is_core_particle` | 是否至少有一条 core-shoulder contact | 二值 core 上色 |
| `lc_strict_core_contact_degree` | strict-core contact 数 | 看最强核心区域 |
| `lc_is_strict_core_particle` | 是否至少有一条 strict-core contact | 二值 strict-core 上色 |
| `lc_domain` | weak 或 robust domain id | 局部束成员 |
| `lc_pearl` | 3D pearl id | 珠状 bead 成员 |
| `lc_state` | `0` none，`1` weak domain，`2` robust domain | 紧凑 domain 状态 |

如果要看 contact geometry，另行加载 `*_lc_contact_segments.vtk`。如果没有 OVITO Pro，可以先看 label dump，再单独打开代表帧的 VTK 文件。

### 如何分析结果

单条轨迹建议：

- 画 `robust_domain_count`、`weak_domain_count`、`pearl_count` 和 `largest_pearl_fraction`；
- 比较 `local_edge_fraction` 与 `nonlocal_edge_fraction`；
- 用 `domain_size_vs_s2.png` 审查 domain 质量；
- 用 `domain_size_vs_domain_count.png` 审查每帧中同 size domain 的数量分布；
- 在 OVITO 中核对代表帧。

外力序列建议：

- potential 不变时复用同一 potential cache；
- threshold prior 只有在代表性构型和温度范围足够时才复用；
- 比较稳态 `L_parallel`、`pearl_count`、`largest_pearl_fraction` 和 nonlocal support；
- 报告时间窗口或重复模拟的不确定性。

降温序列建议：

- 尽量用覆盖高温到低温构型变化的轨迹生成 threshold prior；
- 检查 2D lobe split 图；
- 用 robust domain count、pearl count 和 core-contact fraction 比较聚集起始。

### 关键参数

| 参数 | 默认或来源 | 目的 |
|---|---|---|
| `contact_mode` | `gayberne` | 使用重算 Gay-Berne pair energy |
| `mesogen_type` | `1` | 聚集成员 |
| `anchor_types` | `3` | 端点和拉伸轴支持 |
| `gb_off_strength` | threshold prior | gray/support edge 阈值 |
| `gb_on_strength` | threshold prior | strong edge 阈值 |
| `p2_cut` | threshold prior | pair 取向门槛 |
| `gb_core_strength` | threshold prior | diagnostic core-shoulder GB gate |
| `p2_core_cut` | `0.71` | diagnostic core-shoulder P2 gate |
| `gb_strict_core_strength` | threshold prior | diagnostic strict-core GB gate |
| `p2_strict_core_cut` | `0.80` | diagnostic strict-core P2 gate |
| `s_excl` | `1` | 局域链序支持范围 |
| `n_min` | `3` | robust-domain size evidence |
| `robust_min_s2` | `0.70` | robust-domain orientation evidence |
| `domain_min_lifetime` | `2` | persistence evidence |
| `pearl_gap_cut` | `auto` | pearl 合并最大 3D gap |
| `pearl_min_cross_contacts` | `2` | 跨域接触数下限 |
| `pearl_min_boundary_particles` | `2` | 边界支持粒子数下限 |
| `pearl_max_aspect_ratio` | `3.0` | bead-like 紧凑性限制 |
| `threshold_prior.global_frame_stride` | `1` | 阈值先验帧采样 stride |
| `workers` | `auto` | 并行核数，由 CPU、文件数和环境上限共同限制 |

### Debug 和可靠性检查

建议按以下顺序检查：

1. 确认 `lc_pearl_preflight/validation/verified_potential.json`。
2. 检查 `lc_pearl_preflight/thresholds/global_thresholds.json`。
3. 对比 `gb_strength_vs_p2_stream_hist.png`、`gb_strength_vs_p2_stream_dotgrid.png` 和 `gb_strength_vs_p2_stream_hexbin.png`。
4. 检查 `lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`。
5. 检查 `diagnostics/diagnostic_summary.json` 中实际应用的阈值。
6. 用 OVITO 检查代表帧 label。
7. 对比 `N_domain`、`N_pearl`、largest domain fraction、largest pearl fraction、`L_parallel` 和 `Rg_parallel/Rg_perp`。

如果阈值看起来不合理，只删除 `lc_pearl_preflight/thresholds/global_thresholds.json` 后重跑 threshold prior。若 potential 参数改变，则需要重新验证 potential 并重建 threshold prior。

### 文档

完整手册包括：

- [MANUAL.md](MANUAL.md)
- [docs/lc_pearl_algorithm_details.md](docs/lc_pearl_algorithm_details.md)
- [docs/lc_pearl_user_guide.md](docs/lc_pearl_user_guide.md)
- [docs/lc_pearl_operation_manual.md](docs/lc_pearl_operation_manual.md)
- [docs/LC_Pearl_v2.1_Academic_User_Manual.md](docs/LC_Pearl_v2.1_Academic_User_Manual.md)

历史过程文档和浏览器 HTML 副本保留在 `docs/` 中，但当前 v2.1 以 Markdown 文档为准。

### 发行版

本仓库当前发布版本为 `LC-Pearl v2.1.0`。旧版源码快照保存在 `releases/` 下，用于可复现回溯。
