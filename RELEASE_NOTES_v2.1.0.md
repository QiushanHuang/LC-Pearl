# LC-Pearl v2.1.0

Released: 2026-05-03

## Release Summary

LC-Pearl v2.1.0 is the first fully documented public pipeline release for the current V2 aggregation algorithm. The core V2 method is a streaming two-dimensional `GB strength x P2` threshold split, followed by multi-tier contact classification, weak/robust domain construction, 3D pearl merging, mechanics, OVITO visualization, and diagnostics.

The central V2 algorithm is:

```text
validated Gay-Berne reconstruction
  -> streaming 2D GB-strength x P2 prior
  -> tiered contacts: none / weak-transition / aggregation / core shoulder / strict core
  -> weak and robust mesogen domains
  -> compact 3D pearls
  -> OVITO lc_aggregation_tier and quantitative diagnostics
```

The key visualization field is `lc_aggregation_tier`, a particle-level tier for OVITO:

| Tier | Meaning |
|---:|---|
| 0 | No accepted contact |
| 1 | Weak or transition support contact |
| 2 | Aggregation/strong contact |
| 3 | Core shoulder contact, using the high-P2 strong-lobe shoulder |
| 4 | Strict core contact, using the stricter high-P2 core subset |

`lc_aggregation_tier` is diagnostic and visual. Domain and pearl classification still use the graph/evidence hierarchy documented in the algorithm details.

## Changes

- Documents the V2 core algorithm explicitly as a 2D `GB strength x P2` split, not a one-dimensional cluster cutoff.
- Documents the multi-tier contact partition: no contact, weak/transition, aggregation, core shoulder, and strict core.
- Documents `lc_aggregation_tier` as the recommended first OVITO color-coding field.
- Adds the modular `lc_domain_size_counts.py` diagnostic tool.
- Adds `diagnostics/domain_size_frame_counts.tsv`, a per-frame and per-domain-size count table.
- Adds `diagnostics/domain_size_vs_domain_count.png`, a scatter plot where each point represents one `(frame, domain size)` group and the y-axis is the number of domains of that size in that frame.
- Integrates the new domain-size per-frame scatter diagnostic into the full pipeline.
- Adds Markdown documentation as the canonical documentation format:
  - `docs/lc_pearl_algorithm_details.md`
  - `docs/lc_pearl_user_guide.md`
  - `docs/lc_pearl_operation_manual.md`
- Adds the complete academic Markdown manual `MANUAL.md`.
- Adds `docs/LC_Pearl_v2.1_Academic_User_Manual.md` as a packaged documentation pointer.
- Updates README and citation metadata for the v2.1.0 release.

## Compatibility

- The core contact, domain, pearl, and threshold-prior logic remains compatible with v2.0.1.
- Existing v2.0.1 threshold-prior caches using the same schema can still be inspected, but v2.1.0 release documentation should be used for interpretation.
- The new diagnostic plot is derived from `domain_diagnostics.tsv` and does not require re-parsing dump trajectories when run standalone.
- Existing HTML documents are retained as browser-readable copies, but Markdown files are the canonical release documentation.

## Validation

The release was checked with Python bytecode compilation, the threshold-prior regression test, the new domain-size count regression test, and source distribution metadata checks.
