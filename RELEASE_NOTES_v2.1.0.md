# LC-Pearl v2.1.0

Released: 2026-05-03

## Release Summary

LC-Pearl v2.1.0 is the first fully documented public pipeline release for the current V2 aggregation algorithm. It packages the complete Gay-Berne validation, streaming threshold prior, weak/robust domain, 3D pearl, mechanics, OVITO visualization, and diagnostics workflow.

## Changes

- Adds the modular `lc_domain_size_counts.py` diagnostic tool.
- Adds `diagnostics/domain_size_frame_counts.tsv`, a per-frame and per-domain-size count table.
- Adds `diagnostics/domain_size_vs_domain_count.png`, a scatter plot where each point represents one `(frame, domain size)` group and the y-axis is the number of domains of that size in that frame.
- Integrates the new domain-size per-frame scatter diagnostic into the full pipeline.
- Updates the canonical user guide `docs/lc_pearl_user_guide.html`.
- Updates the canonical algorithm reference `docs/lc_pearl_algorithm_reference.html`.
- Adds the complete academic Markdown manual `MANUAL.md`.
- Adds `docs/LC_Pearl_v2.1_Academic_User_Manual.md` as a packaged documentation pointer.
- Updates README and citation metadata for the v2.1.0 release.

## Compatibility

- The core contact, domain, pearl, and threshold-prior logic remains compatible with v2.0.1.
- Existing v2.0.1 threshold-prior caches using the same schema can still be inspected, but v2.1.0 release documentation should be used for interpretation.
- The new diagnostic plot is derived from `domain_diagnostics.tsv` and does not require re-parsing dump trajectories when run standalone.

## Validation

The release was checked with Python bytecode compilation, the threshold-prior regression test, the new domain-size count regression test, and source distribution metadata checks.
