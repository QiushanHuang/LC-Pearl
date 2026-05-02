# LC-Pearl v1.0.0

[![简体中文](https://img.shields.io/badge/语言-简体中文-1677ff)](README.zh-CN.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Release](https://img.shields.io/badge/release-v1.0.0-blue.svg)](RELEASE_NOTES_v1.0.0.md)

LC-Pearl is a reproducible analysis pipeline for quantifying liquid-crystal mesogen aggregation in LAMMPS dump trajectories. It was designed for single-chain or chain-like LC elastomer simulations where Gay-Berne attraction, local orientational order, chain connectivity, and external stretching compete to form pearl-necklace-like structures.

The v1.0.0 release defines a three-level hierarchy:

- `mesogen contact`: a pairwise attractive, orientationally relevant type-1 ellipsoid interaction.
- `domain`: a local mesogen bundle inferred from the contact graph and classified as weak or robust.
- `pearl`: a compact 3D bead-like assembly of one or more robust domains.

The pipeline also reports axial mechanics (`L_parallel`, `Rg_parallel`, `Rg_perp`, `S2_force`) so aggregation can be compared against force, stretch, temperature, and time.

## What Is Released

This repository contains the complete LC-Pearl v1 pipeline:

- Python analysis scripts in `scripts/`.
- Current-directory launcher `lc_pearl_here.py`.
- TOML-driven launcher `lc_pearl_cli.py`.
- Preflight template in `templates/lc_pearl_preflight/`.
- Quick configuration files in `configs/`.
- Academic user and algorithm manuals in `docs/`.
- Regression test for the 2D lobe threshold prior in `tests/`.

The software release is `v1.0.0`. The current threshold-prior artifact schema is `schema_version = 5`, named `LC Domain-Pearl V2 2D lobe streaming threshold prior`; this is an internal algorithm schema name, not the GitHub release number.

## Installation

Use Python 3.11 or newer. The minimal runtime dependencies are `numpy` and `matplotlib`.

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 -m pip install -e .
```

If you do not install the package, run scripts directly with the same Python environment:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

## Recommended Current-Directory Workflow

Enter a folder that contains LAMMPS dump files and run:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

On first use, LC-Pearl creates `lc_pearl_preflight/` and tells you what is missing. Put the required source files there:

- `lc_pearl_preflight/lammps/*.in` or `*.lmp`: original LAMMPS input containing `pair_style gayberne` and the type 1-1 `pair_coeff`.
- `lc_pearl_preflight/lammps/lammps_executable.txt`: one line with the LAMMPS executable path, for example `/Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp`.
- Optional `lc_pearl_preflight/topology/*.data` or `*.dat`: LAMMPS data file with `Atoms` and `Bonds` sections. LC-Pearl converts it into local-pair and excluded-pair tables.

Then run:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

`auto` performs the following sequence:

1. Reuse or create `lc_pearl_preflight/validation/verified_potential.json`.
2. Reuse or create `lc_pearl_preflight/thresholds/global_thresholds.json`.
3. Apply the 2D lobe threshold prior before the main analysis.
4. Run the domain, pearl, OVITO-label, and mechanics analysis.

## Config-Driven Workflow

Edit `configs/quick_run.toml`, then preview the command:

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

## Main Outputs

Output is written by default to `lc_domain_pearl_v2_output/`.

- `aggregation_timeseries.tsv`: frame-level quantitative table for aggregation, domains, pearls, and mechanics.
- `per_file/*_aggregation.tsv`: per-input-file time series.
- `per_file/*_summary.json`: compact per-file summary.
- `per_file/*_lc_labels.dump`: OVITO-readable particle labels, when OVITO label output is enabled.
- `per_file/*_lc_cluster_envelopes.dump`: visual cluster envelope particles.
- `per_file/*_lc_contact_edges.dump`: local-style edge records for attractive contacts.
- `per_file/*_lc_contact_segments.vtk`: OVITO-loadable contact line segments.
- `diagnostics/diagnostic_summary.json`: thresholds, counts, and diagnostic provenance.
- `diagnostics/gb_strength_vs_p2.png`: accepted-edge diagnostic under the current main-analysis gates.
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hist.png`: full streaming candidate-pair distribution used for threshold prior.
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`: dot-grid view of the same full 2D distribution with selected thresholds.
- `lc_pearl_preflight/thresholds/global_thresholds.json`: reusable threshold-prior artifact.

By default, LC-Pearl does not write the full `edge_diagnostics.tsv` or full per-frame JSONL debug records because these files can become multi-GB outputs. Enable them only when needed:

```toml
[analysis]
edge_diagnostics_table = "sample"  # off, sample, or full
edge_diagnostics_sample_size = 200000
write_frame_jsonl = true
```

## OVITO Visualization

Open `per_file/*_lc_labels.dump` in OVITO. Useful scalar fields include:

- `lc_cluster`: visual center-distance cluster id.
- `lc_cluster_size`: size of the visual cluster.
- `lc_contact_degree`: number of accepted contact edges incident on a mesogen.
- `lc_min_pair_energy`: most attractive pair energy involving that mesogen.
- `lc_mean_pair_energy`: mean attractive pair energy involving that mesogen.
- `lc_max_gb_strength`: strongest normalized attraction involving that mesogen.
- `lc_mean_gb_strength`: mean normalized attraction involving that mesogen.
- `lc_domain`: robust or weak domain id.
- `lc_pearl`: 3D pearl id.
- `lc_state`: compact state code for unassigned, weak-domain, robust-domain, and pearl-supported particles.

For contact geometry, load `*_lc_contact_segments.vtk` as an additional OVITO pipeline. If OVITO Pro is unavailable, inspect labels and envelopes first, then load line segments separately for representative frames.

## Algorithm Summary

LC-Pearl uses type-1 ellipsoids as mesogen members. Type-2 spheres and type-3 anchors are not aggregation members, but they can support chain ordering, endpoints, stretch direction, and topology reconstruction.

For each candidate E-E pair, the Gay-Berne pair energy is reconstructed from the dump coordinates, quaternions, shape axes, and the LAMMPS `pair_style/pair_coeff` source file. The normalized attraction strength is

```text
gb_strength = max(0, -U_GB / U_well)
```

The orientational score is

```text
P2 = (3 |u_i dot u_j|^2 - 1) / 2
```

The current threshold prior is selected from the full 2D distribution of `gb_strength x P2`, not from already accepted edges. It identifies the high-P2 weak-contact lobe, the high-P2 strong-attraction lobe, and the valley between them. This gives `gb_on`; `gb_off` is a conservative left-lobe shoulder for gray/support contacts; `p2_cut` is the orientation gate used consistently in the 2D split.

Domains are then built from the contact graph. Weak local domains are retained. Robust domains require size and evidence such as orientation, persistence, nonlocal support, and parameter stability. Pearls are defined at the next level: robust domains merge into a pearl only when their 3D gap, cross contacts, boundary support, and bead-like aspect ratio satisfy the pearl criteria.

## Key Parameters

- `contact_mode = "gayberne"`: use reconstructed Gay-Berne pair energy for contact strength.
- `gb_on_strength`: strong edge threshold from the threshold prior.
- `gb_off_strength`: gray/support edge threshold from the threshold prior.
- `p2_cut`: pair orientational threshold.
- `n_min`: minimum mesogens for robust-domain size evidence.
- `s_excl`: chain-sequence separation treated as local support.
- `domain_min_lifetime`: processed-frame age used as persistence evidence.
- `pearl_gap_cut`: maximum 3D domain-domain gap for pearl merging.
- `pearl_min_cross_contacts`: minimum cross-domain contacts for pearl merging.
- `pearl_min_boundary_particles`: minimum supported boundary particles on both domains.
- `pearl_max_aspect_ratio`: upper bound for bead-like compactness.
- `threshold_prior.global_frame_stride`: full or strided frame sampling for the threshold prior. `1` means all candidate frames.
- `workers = "auto"`: uses up to CPU count, input size, and the environment cap `LC_PEARL_MAX_WORKERS` or `LC_PEARL_MAX_AUTO_WORKERS`.

## Debugging and Reliability Checks

Use this order:

1. Confirm the potential cache: `lc_pearl_preflight/validation/verified_potential.json`.
2. Inspect `lc_pearl_preflight/thresholds/global_thresholds.json`.
3. Compare `gb_strength_vs_p2_stream_hist.png` and `gb_strength_vs_p2_stream_lobe_split_dotgrid.png`.
4. Check `diagnostics/diagnostic_summary.json` for actual applied thresholds.
5. Inspect OVITO labels for representative frames.
6. Compare `N_domain`, `N_pearl`, largest domain fraction, largest pearl fraction, `L_parallel`, and `Rg_parallel/Rg_perp` across force or temperature.

If threshold cuts look wrong, rebuild the threshold prior by deleting only `lc_pearl_preflight/thresholds/global_thresholds.json` and rerunning. If potential parameters change, rerun validation and rebuild the threshold prior.

## Documentation

The full academic manuals are:

- `docs/lc_pearl_user_guide.html`
- `docs/lc_pearl_algorithm_reference.html`

Historical process notes are kept separately in `docs/` and are not used as the current algorithm reference.

## Release

This repository is prepared as `LC-Pearl v1.0.0`, the first stable release of the domain-pearl aggregation analysis pipeline.
