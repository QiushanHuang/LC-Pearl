# LC-Pearl v2.1.0

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Release](https://img.shields.io/badge/release-v2.1.0-blue.svg)](RELEASE_NOTES_v2.1.0.md)

LC-Pearl is a reproducible analysis pipeline for quantifying liquid-crystal mesogen aggregation in LAMMPS dump trajectories. It is designed for single-chain or chain-like liquid-crystal simulations where Gay-Berne attraction, local orientational order, chain connectivity, temperature, and external stretching compete to form pearl-necklace-like morphologies.

LC-Pearl 将液晶聚集拆成三个层级：pair contact、domain 和 pearl。这样可以保留相邻 mesogen 的 weak local bundling，同时避免把链连接本身误判为稳定 robust aggregation。

## What v2.1.0 Provides

- Gay-Berne pair-energy reconstruction with validation/cache support.
- Streaming 2D `GB strength x P2` threshold prior before the main analysis.
- Weak and robust domain classification with size, orientation, persistence, nonlocal support, and parameter-stability evidence.
- 3D pearl merging for compact bead-like assemblies of robust domains.
- Diagnostic core tiers: `core shoulder` with `P2 > 0.71` and `strict core` with `P2 > 0.80`.
- OVITO label dumps, contact segments, visual cluster envelopes, and scalar fields for color coding.
- Mechanics observables: `L_parallel`, `Rg_parallel`, `Rg_perp`, and `S2_force`.
- New v2.1.0 diagnostic scatter: `domain_size_vs_domain_count.png`, where each point is one `(frame, domain size)` group and the y-axis is the number of domains of that size in that frame.

## Documentation

- Complete academic manual: [MANUAL.md](MANUAL.md)
- Browser user guide: [docs/lc_pearl_user_guide.html](docs/lc_pearl_user_guide.html)
- Browser algorithm reference: [docs/lc_pearl_algorithm_reference.html](docs/lc_pearl_algorithm_reference.html)
- Packaged manual pointer: [docs/LC_Pearl_v2.1_Academic_User_Manual.md](docs/LC_Pearl_v2.1_Academic_User_Manual.md)
- Release notes: [RELEASE_NOTES_v2.1.0.md](RELEASE_NOTES_v2.1.0.md)

## Installation

Use Python 3.11 or newer.

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 -m pip install -e .
```

Runtime dependencies are `numpy` and `matplotlib`.

## Recommended Current-Directory Workflow

Enter a dump output folder:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

On first use, LC-Pearl creates `lc_pearl_preflight/` and reports missing inputs. Put the original LAMMPS input, LAMMPS executable path, and optional topology data there, then run:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

`auto` performs:

1. potential validation or cache reuse;
2. threshold prior generation or cache reuse;
3. application of recommended thresholds before the main analysis;
4. domain, pearl, mechanics, OVITO, and diagnostics output.

## Preflight Inputs

| Path | Required | Purpose |
|---|---:|---|
| `lc_pearl_preflight/lammps/*.in` or `*.lmp` | yes | Source for `pair_style gayberne` and type 1-1 `pair_coeff` |
| `lc_pearl_preflight/lammps/lammps_executable.txt` | recommended | LAMMPS executable for run-0 / microstate validation |
| `lc_pearl_preflight/topology/*.data` | optional but recommended | Generates local/excluded pair tables |
| `lc_pearl_preflight/validation/verified_potential.json` | generated | Validated potential cache |
| `lc_pearl_preflight/thresholds/global_thresholds.json` | generated | Reusable threshold-prior cache |

## Config-Driven Workflow

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py run --config configs/quick_run.toml
```

Validation only:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py validate --config configs/quick_validate.toml
```

## Main Outputs

Output root defaults to `lc_domain_pearl_v2_output/`.

| Output | Meaning |
|---|---|
| `aggregation_timeseries.tsv` | Frame-level quantitative table |
| `per_file/*_lc_labels.dump` | OVITO-readable labels |
| `per_file/*_lc_contact_segments.vtk` | OVITO contact line segments |
| `diagnostics/domain_diagnostics.tsv` | Domain evidence table |
| `diagnostics/domain_size_vs_s2.png` | Domain size vs internal S2 scatter |
| `diagnostics/domain_size_frame_counts.tsv` | Per-frame, per-size domain counts |
| `diagnostics/domain_size_vs_domain_count.png` | Per-frame domain-count scatter by size |
| `diagnostics/diagnostic_summary.json` | Applied thresholds and provenance |
| `lc_pearl_preflight/thresholds/global_thresholds.json` | Threshold-prior artifact |

## OVITO Color Coding

Open `per_file/*_lc_labels.dump`.

Recommended fields:

- `lc_aggregation_tier`: `0` no contact, `1` weak/transition, `2` aggregation, `3` core shoulder, `4` strict core.
- `lc_domain`: weak or robust domain membership.
- `lc_pearl`: 3D pearl membership.
- `lc_state`: `0` none, `1` weak domain, `2` robust domain.
- `lc_contact_degree`: number of accepted contacts.
- `lc_core_contact_degree`: contacts passing `gb_core_strength` and `P2 > 0.71`.
- `lc_strict_core_contact_degree`: contacts passing `gb_strict_core_strength` and `P2 > 0.80`.

Load `*_lc_contact_segments.vtk` separately to inspect contact lines.

## Reliability Checklist

Before using results in a group meeting or publication, verify:

1. `verified_potential.json` matches the current LAMMPS potential.
2. `global_thresholds.json` was generated from representative data.
3. `gb_strength_vs_p2_stream_*` figures show an interpretable 2D structure.
4. `diagnostics/diagnostic_summary.json` contains the thresholds you expect.
5. OVITO labels agree with representative frames.
6. `domain_size_vs_s2.png` and `domain_size_vs_domain_count.png` are consistent with the observed morphology.

## 中文快速说明

LC-Pearl v2.1.0 的默认分析顺序是：

```text
dump -> 前置文件夹 -> GB potential 验证 -> 2D 阈值先验 -> domain/pearl 主分析 -> OVITO/diagnostics
```

核心概念：

- `domain`：局部 mesogen 束，分为 weak 和 robust。
- `pearl`：一个或多个 robust domains 在 3D 空间中组成的 bead-like 珠子。
- `gb_strength`：归一化 Gay-Berne 吸引强度。
- `P2`：pair 取向一致性。
- `core shoulder` / `strict core`：用于 OVITO 和审查的诊断层，不改变主 domain/pearl 定义。

完整中文/英文混合学术手册见 [MANUAL.md](MANUAL.md)。浏览器版见 [docs/lc_pearl_user_guide.html](docs/lc_pearl_user_guide.html) 和 [docs/lc_pearl_algorithm_reference.html](docs/lc_pearl_algorithm_reference.html)。
