# LC-Pearl v2.1 User Guide

This guide explains how to use the pipeline outputs to analyze liquid-crystal aggregation.

## 1. Fast Run

From a folder containing LAMMPS dump files:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

`auto` will:

1. create or reuse `lc_pearl_preflight/validation/verified_potential.json`;
2. create or reuse `lc_pearl_preflight/thresholds/global_thresholds.json`;
3. apply the 2D threshold prior before the main analysis;
4. generate labels, time series, diagnostics, and OVITO visualization files under `lc_domain_pearl_v2_output/`.

## 2. What To Check First

Start with these files:

| File | What to check |
|---|---|
| `lc_pearl_preflight/thresholds/global_thresholds.json` | Confirm applied `gb_off`, `gb_on`, `p2_cut`, `gb_core`, and `gb_strict_core`. |
| `lc_pearl_preflight/thresholds/gb_strength_vs_p2_stream_hexbin.png` | Check whether the 2D split matches the visible GB-strength x P2 distribution. |
| `lc_domain_pearl_v2_output/aggregation_timeseries.tsv` | Main time series: domains, pearls, contacts, mechanics. |
| `lc_domain_pearl_v2_output/diagnostics/domain_diagnostics.tsv` | Domain-level evidence: size, S2, persistence, nonlocal support, classification. |
| `lc_domain_pearl_v2_output/diagnostics/domain_size_vs_s2.png` | Whether robust domains are large and internally aligned. |
| `lc_domain_pearl_v2_output/diagnostics/domain_size_vs_domain_count.png` | How many domains of each size occur per frame. |
| `lc_domain_pearl_v2_output/per_file/*_lc_labels.dump` | Main OVITO label dump. |

## 3. OVITO Workflow

Open `per_file/*_lc_labels.dump`.

Recommended inspection order:

1. Color by `lc_aggregation_tier`.
2. Color by `lc_domain`.
3. Color by `lc_pearl`.
4. Inspect `lc_core_contact_degree` and `lc_strict_core_contact_degree` for strong core regions.
5. Compare with `lc_cluster` only as a simple spatial-cluster view.

`lc_aggregation_tier` values:

| Value | Meaning |
|---:|---|
| 0 | No accepted contact. |
| 1 | Weak/transition support contact. |
| 2 | Aggregation/strong contact. |
| 3 | Core shoulder: high-P2 strong-lobe shoulder. |
| 4 | Strict core: strongest high-P2 contact subset. |

This field is particle-level: a mesogen receives the maximum tier among its incident contacts.

## 4. Interpreting Domain And Pearl

Use `lc_domain` to ask:

```text
Which mesogens form a local LC bundle?
```

Use `lc_pearl` to ask:

```text
Which robust domains occupy the same compact 3D bead?
```

They are intentionally different. A pearl can contain multiple robust domains if the domains are spatially close, have enough cross-domain support, and remain bead-like after merging.

## 5. Interpreting The 2D Threshold Figures

Use these figures in the preflight threshold folder:

| Figure | Purpose |
|---|---|
| `gb_strength_vs_p2_stream_hist.png` | Dense heatmap of the full streaming distribution. |
| `gb_strength_vs_p2_stream_dotgrid.png` | Dot-grid view for nonzero histogram bins. |
| `gb_strength_vs_p2_stream_hexbin.png` | Hexbin-style view closest to the old audit plot style. |
| `lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png` | Actual threshold overlay. |
| `lobe_split_preview/gb_core_slice_hist.png` | Explains `gb_core` and `gb_strict_core`. |

The most important question is not whether a line "proves stable aggregation" by itself. The threshold lines define contact tiers used by the graph algorithm. Stability is then evaluated by domain size, S2, persistence, nonlocal support, and parameter stability.

## 6. Main Quantities For Analysis

| Quantity | Meaning |
|---|---|
| `robust_domain_count` | Number of robust local bundles in a frame. |
| `weak_domain_count` | Number of weak or transitional bundles. |
| `pearl_count` | Number of compact 3D beads. |
| `largest_domain_fraction` | Whether one local bundle dominates. |
| `largest_pearl_fraction` | Whether one pearl dominates. |
| `local_edge_fraction` | Fraction of accepted contacts that are local/adjacent. |
| `nonlocal_edge_fraction` | Fraction of contacts that support true nonlocal aggregation. |
| `core_contact_fraction` | Fraction of contacts in the core-shoulder tier. |
| `strict_core_contact_fraction` | Fraction of contacts in the strict-core tier. |
| `L_parallel` | Extension along stretch axis. |
| `Rg_parallel`, `Rg_perp` | Chain shape along and perpendicular to the stretch axis. |
| `S2_force` | Mesogen alignment relative to stretch direction. |

For a force series, compare time-window averages of `L_parallel`, `pearl_count`, `largest_pearl_fraction`, `nonlocal_edge_fraction`, and core fractions.

## 7. How To Judge Reliability

A result is ready for discussion only if:

- the potential cache is valid for the current LAMMPS input;
- the threshold prior was generated from representative data or deliberately reused;
- the 2D threshold figures show separable or at least interpretable regions;
- `diagnostic_summary.json` reports the expected applied thresholds;
- OVITO `lc_aggregation_tier` matches visual aggregation in representative frames;
- `lc_domain` and `lc_pearl` do not contradict the physical picture;
- domain-size plots are consistent with the expected fragmentation/coarsening behavior.

## 中文使用要点

最先看 `lc_aggregation_tier`，它把粒子分成无接触、weak/transition、aggregation、core shoulder、strict core 五类。然后看 `lc_domain` 判断局部液晶束，看 `lc_pearl` 判断多个 robust domain 是否形成同一个 3D 珠子。阈值图的作用是定义 contact tier，不是单独证明稳定聚集；稳定性由后续 domain/pearl 证据共同判断。
