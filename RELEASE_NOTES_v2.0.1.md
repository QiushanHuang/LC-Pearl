# LC-Pearl v2.0.1

Released: 2026-05-02

## Changes

- Adds a two-level diagnostic core tier on top of the main `gb_off / gb_on / p2_cut` contact graph.
- Redefines `gb_core_strength / p2_core_cut` as the high-P2 strong-lobe `core_shoulder` tier, estimated from the left shoulder of the strong-contact lobe.
- Adds `gb_strict_core_strength / p2_strict_core_cut` as a stricter `strict_core` tier, estimated from the high-P2 weighted-q25 region with default `p2_strict_core_cut = 0.80`.
- Propagates strict-core thresholds into the pipeline, main analysis, label dumps, contact-edge dumps, VTK contact segments, diagnostics, and OVITO Basic merge output.
- Keeps domain and pearl construction unchanged: `gb_core_*` and `gb_strict_core_*` remain visualization/audit tiers and do not alter robust-domain or pearl merging rules.
- Bumps the streaming threshold-prior schema to `7` so v2.0.1 priors are not confused with v2.0.0 priors.

## OVITO Interpretation

- `lc_state`: `0` non-domain, `1` weak domain, `2` robust domain.
- `lc_domain`: mesogen contact-domain id.
- `lc_pearl`: 3D bead-like pearl id after robust-domain merging.
- `lc_core_*`: core-shoulder contact tier.
- `lc_strict_core_*`: stricter strong-core contact tier.

