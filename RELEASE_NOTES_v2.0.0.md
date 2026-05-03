# LC-Pearl v2.0.0 Release Notes

Date: 2026-05-02

## Main Change

LC-Pearl v2.0.0 adds a diagnostic very-strong core-contact layer on top of the existing domain/pearl analysis.

- `gb_off_strength`, `gb_on_strength`, and `p2_cut` still define gray/support and strong contacts for domain and pearl construction.
- `gb_core_strength` and `p2_core_cut` mark stricter core contacts for audit and visualization only.
- The default strict orientation gate is `p2_core_cut = 0.71`, equivalent to `|u_i dot u_j| > 0.898` or within about 26 degrees of parallel/anti-parallel alignment.
- `gb_core_strength` is selected from the full streaming `GB strength x P2` histogram using the high-P2 strong-lobe shoulder, with a weighted-quantile fallback.

## Output Changes

- Threshold prior schema upgraded from `5` to `6`, so old threshold caches are rebuilt instead of silently reused.
- `gb_strength_vs_p2_stream_hist.png` now uses the dot-grid style and includes `gb_core` and `p2_core` lines.
- `lobe_split_preview/` now contains the main lobe/core split preview and `gb_core_slice_hist.png`.
- OVITO label dumps include core-contact fields such as `lc_core_contact_degree`, `lc_core_nonlocal_degree`, `lc_is_core_particle`, `lc_max_core_gb_strength`, and `lc_max_core_p2`.
- Contact edge dumps and VTK contact segments include `is_core_contact`, `is_core_seed`, `contact_tier_code`, `passes_gb_core`, `passes_p2_core`, and core margins.
- Timeseries and diagnostic summaries include core-contact counts and fractions.

## Compatibility

The v1.0.0 source snapshot is preserved in `releases/lc-pearl-v1.0.0-source-20260502.tar.gz`.

The core-contact layer is intentionally diagnostic-only in v2.0.0. It does not change robust-domain classification or pearl merging by default, so the main domain/pearl time series remains comparable with the previous V2 lobe-prior logic.
