# LC-Pearl v1.0.0 Release Notes

LC-Pearl v1.0.0 is the first stable release of the liquid-crystal domain-pearl aggregation pipeline.

This release is open source under the MIT License. The listed contributor for the v1.0.0 release is Qiushan Huang / 黄秋山.

## Scientific Scope

- Quantifies mesogen aggregation in LAMMPS dump trajectories with type-1 ellipsoids as mesogen members.
- Separates aggregation into contact, weak/robust domain, and 3D pearl levels.
- Supports OVITO visualization through per-atom labels, cluster envelopes, contact-edge dumps, and VTK contact segments.
- Reports mechanics-compatible quantities including `L_parallel`, `Rg_parallel`, `Rg_perp`, and `S2_force`.

## Algorithm Highlights

- Gay-Berne pair energy reconstruction from dump coordinates, quaternions, shape axes, and LAMMPS input parameters.
- Full-streaming 2D threshold prior over `GB attraction strength x pair P2`.
- Threshold artifact schema v5: `LC Domain-Pearl V2 2D lobe streaming threshold prior`.
- Robust-domain evidence from size, orientation, persistence, nonlocal support, and parameter-stability checks.
- Pearl construction as 3D compact bead merging of robust domains, constrained by gap, cross contacts, boundary support, and aspect ratio.

## Engineering Defaults

- Post-run GB candidate-pair audit is not part of the default main pipeline.
- Full `edge_diagnostics.tsv` and full per-frame JSONL debug records are disabled by default to avoid multi-GB text outputs.
- `run_summary.json` is compact and records output locations rather than embedding all per-frame records.
- Preflight folders make potential validation and threshold-prior reuse portable across related simulations.

## Primary Documentation

- `docs/lc_pearl_user_guide.html`
- `docs/lc_pearl_algorithm_reference.html`

## Validation

The v1.0.0 source tree was checked with Python bytecode compilation and the 2D lobe threshold-prior regression test.
