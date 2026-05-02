# LC-Pearl v1.0.0 English README

[![简体中文](https://img.shields.io/badge/语言-简体中文-1677ff)](README.zh-CN.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

LC-Pearl is an open-source analysis pipeline for quantifying liquid-crystal mesogen aggregation in LAMMPS dump trajectories. It reconstructs or audits Gay-Berne pair attraction between type-1 ellipsoids and combines attraction strength, orientational order, chain topology, local/nonlocal support, and 3D bead geometry into reproducible domain-pearl metrics.

## Scientific Objective

LC-Pearl is designed for simulations where external stretching competes with mesogen attraction. The main goal is to convert OVITO-observed structures, such as several medium-size clusters packed together and partially pulled apart, into auditable quantitative observables.

The hierarchy has three levels:

- `mesogen contact`: whether two ellipsoidal mesogens have sufficient attraction and compatible orientation.
- `domain`: whether contact-connected mesogens form a local bundle. Weak local domains are retained, while robust domains require evidence from size, orientation, persistence, nonlocal support, and parameter stability.
- `pearl`: whether one or more robust domains form the same compact 3D bead-like aggregate.

## Installation

Use Python 3.11 or newer with `numpy` and `matplotlib`.

```bash
cd /path/to/LC-Pearl
python3 -m pip install -e .
```

The scripts can also be run directly without package installation.

## Minimal Workflow

Enter a directory containing LAMMPS dump files:

```bash
cd /path/to/dump_output_folder
python3 /path/to/LC-Pearl/lc_pearl_here.py
```

The first run creates `lc_pearl_preflight/` and reports missing inputs. Typical inputs are:

- `lc_pearl_preflight/lammps/*.in` or `*.lmp`: a LAMMPS input file containing `pair_style gayberne` and the type 1-1 `pair_coeff`.
- `lc_pearl_preflight/lammps/lammps_executable.txt`: one line containing the LAMMPS executable path.
- `lc_pearl_preflight/topology/*.data` or `*.dat`: optional LAMMPS data file used to generate local/exclude pair tables.

Then run:

```bash
python3 /path/to/LC-Pearl/lc_pearl_here.py auto
```

## Outputs

The default output directory is `lc_domain_pearl_v2_output/`.

- `aggregation_timeseries.tsv`: frame-level table for aggregation, domain, pearl, and mechanical observables.
- `per_file/*_lc_labels.dump`: OVITO-readable particle label dump.
- `per_file/*_lc_cluster_envelopes.dump`: auxiliary cluster-envelope visualization particles.
- `per_file/*_lc_contact_segments.vtk`: contact line segments.
- `diagnostics/diagnostic_summary.json`: thresholds and diagnostic provenance used by the main analysis.
- `lc_pearl_preflight/validation/verified_potential.json`: validated Gay-Berne potential cache.
- `lc_pearl_preflight/thresholds/global_thresholds.json`: reusable global threshold-prior artifact.
- `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`: visualization of the 2D lobe threshold prior.

## OVITO Fields

- `lc_cluster`: center-distance visual cluster id, not the final robust-domain definition.
- `lc_cluster_size`: number of mesogens in the visual cluster.
- `lc_contact_degree`: number of accepted contacts involving the mesogen.
- `lc_min_pair_energy`: strongest attractive pair energy involving the mesogen.
- `lc_mean_pair_energy`: mean contact energy involving the mesogen.
- `lc_max_gb_strength`: maximum normalized Gay-Berne attraction strength.
- `lc_mean_gb_strength`: mean normalized Gay-Berne attraction strength.
- `lc_domain`: domain id.
- `lc_pearl`: pearl id.
- `lc_state`: compact state code for unassigned, weak, robust, and pearl-supported particles.

## Documentation

- User guide: [docs/lc_pearl_user_guide.html](docs/lc_pearl_user_guide.html)
- Algorithm reference: [docs/lc_pearl_algorithm_reference.html](docs/lc_pearl_algorithm_reference.html)
- Release notes: [RELEASE_NOTES_v1.0.0.md](RELEASE_NOTES_v1.0.0.md)

## Contributor

Qiushan Huang.

## Open-Source Notice

LC-Pearl is released under the MIT License. You may use, modify, copy, and redistribute the software, provided that the copyright and license notice are retained. For academic use, cite `CITATION.cff` and report the Gay-Berne validation status, threshold prior, domain definition, and pearl definition used in the analysis.
