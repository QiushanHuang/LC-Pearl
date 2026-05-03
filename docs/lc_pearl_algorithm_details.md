# LC-Pearl v2.1 Algorithm Details

This document is the canonical Markdown algorithm-detail document for LC-Pearl v2.1. The central v2 algorithm is not a one-dimensional distance cluster. It is a two-dimensional contact classification and a hierarchical aggregation pipeline:

```text
LAMMPS dump + preflight inputs
  -> validated Gay-Berne reconstruction
  -> streaming 2D GB-strength x P2 prior
  -> multi-tier mesogen contacts
  -> weak/robust domain graph
  -> 3D pearl merge
  -> mechanics, OVITO labels, diagnostics
```

## 1. Physical Objects

LC-Pearl separates four concepts that should not be mixed:

| Object | Meaning | Why it exists |
|---|---|---|
| Mesogen | Type-1 ellipsoid particle. This is the aggregation member. | Only mesogen-mesogen attraction defines liquid-crystal aggregation. |
| Contact | Pair-level attraction plus orientational compatibility. | Quantifies whether two mesogens support aggregation. |
| Domain | Local bundle of mesogens connected by accepted contacts. | Measures local LC bundling, including weak adjacent-chain bundles. |
| Pearl | Compact 3D bead-like object made from one or more robust domains. | Measures pearl-necklace bead morphology, not just graph connectivity. |

Type-2 spheres and type-3 anchors are not domain or pearl members by default. They can still be used to infer chain order, endpoints, stretch direction, and topology-derived local/excluded pairs.

## 2. Pair Quantities

For each candidate mesogen pair `i,j`, LC-Pearl computes:

```text
r_i = unwrapped mesogen position
u_i = mesogen long-axis direction from quaternion
q_ij = |u_i dot u_j|
P2_ij = (3 q_ij^2 - 1) / 2
U_GB,ij = reconstructed Gay-Berne pair energy
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
delta_s_ij = mesogen sequence separation, if topology/order is available
```

`P2` is a liquid-crystal-style pair alignment score. `P2 = 1` means parallel or antiparallel long axes, and `P2 = -0.5` means perpendicular axes. The absolute dot product is used because ellipsoid axes are head-tail symmetric.

`gb_strength` is not a new potential. It is a normalized attraction score. If the pair is repulsive or non-attractive, the score is 0. If the pair approaches the orientation-dependent attractive well, the score approaches 1. Values slightly above 1 can appear from numerical/binning details or if the reconstructed well normalization differs from a sampled local minimum; the plotting code bins values at the upper edge instead of dropping them.

## 3. Validated Gay-Berne Energy

The production contact mode is `contact_mode = gayberne`.

LC-Pearl reads the LAMMPS input file from `lc_pearl_preflight/lammps/gb_param_source.in` and extracts the `pair_style gayberne` and type 1-1 `pair_coeff` parameters. Python then reconstructs the GB pair energy from:

- positions `xu,yu,zu`;
- quaternions;
- ellipsoid shape axes;
- GB epsilon/sigma/mu/nu/cutoff parameters;
- type mapping for mesogens.

The validation layer creates or reuses:

```text
lc_pearl_preflight/validation/verified_potential.json
```

Use the cache when the potential settings, units, shape axes, type mapping, and quaternion convention are unchanged. Revalidate when any of those change. Temperature, stretching force, and trajectory length do not by themselves require a new validation cache if the potential is unchanged.

## 4. Core V2 Threshold Algorithm: 2D GB-Strength x P2 Split

The v2 threshold prior is the core algorithmic upgrade. Thresholds are selected before the main analysis from the global or strided candidate-pair distribution.

The prior module `scripts/lc_threshold_prior.py` streams over frames and pairs. It does not need to write a huge all-pair table. It accumulates a 2D histogram:

```text
x-axis: gb_strength
y-axis: P2
weight: candidate-pair count
```

This 2D distribution is used because aggregation is not defined by energy alone or orientation alone. A pair should support LC aggregation only when it is both attractive and sufficiently aligned.

## 5. Multi-Tier Contact Partition

The v2 algorithm partitions contacts into ordered tiers:

| Tier | Edge meaning | Pair rule | Main use |
|---:|---|---|---|
| 0 | No accepted contact | Fails gray/support gate | Background particles in OVITO. |
| 1 | Weak/transition support | `gb_strength >= gb_off_strength` and `P2 >= p2_cut`, but not strong | Boundary support, weak local domain evidence. |
| 2 | Aggregation contact | `gb_strength >= gb_on_strength` and `P2 >= p2_cut` | Strong edge used to seed domain graph. |
| 3 | Core shoulder | `gb_strength >= gb_core_strength` and `P2 > 0.71` | Diagnostic high-alignment strong-lobe shoulder. |
| 4 | Strict core | `gb_strength >= gb_strict_core_strength` and `P2 > 0.80` | Diagnostic strongest high-alignment core subset. |

Important: tiers 3 and 4 are diagnostic. They do not replace the main strong/gray graph rules. They make OVITO and diagnostics clearer by separating normal aggregation contacts from the strongest high-P2 core contacts.

## 6. How Thresholds Are Chosen

The threshold-prior artifact is:

```text
lc_pearl_preflight/thresholds/global_thresholds.json
```

It stores the applied values, confidence/status fields, and provenance. The key thresholds are:

| Parameter | How it is selected | Why |
|---|---|---|
| `p2_cut` | Valley between low-orientation background and high-orientation lobe, clamped to a physically useful range when needed. | Prevents random or perpendicular contacts from joining domains. |
| `gb_on_strength` | High-P2 GB-strength valley between weak-contact lobe and strong-attraction lobe; fallback is high-P2 weighted quantile. | Defines strong aggregation edges. |
| `gb_off_strength` | Conservative left-lobe shoulder below `gb_on_strength`. | Keeps weak/boundary support without making it a strong seed. |
| `p2_core_cut` | Fixed at `0.71`. | Requires high pair alignment for core-shoulder diagnostics. |
| `gb_core_strength` | Left shoulder of the high-P2 strong lobe above `gb_on_strength`. | Marks the beginning of the visibly strong core-contact region. |
| `p2_strict_core_cut` | Fixed at `0.80`. | Stricter high-alignment diagnostic gate. |
| `gb_strict_core_strength` | Weighted high-strength subset of the strict high-P2 strong lobe. | Marks the strongest core contacts. |

`gb_off`, `gb_on`, `gb_core`, and `gb_strict_core` are ordered thresholds. Their intended interpretation is:

```text
below gb_off: no accepted aggregation support
gb_off to gb_on: transition / weak support
gb_on to gb_core: aggregation contact
gb_core to gb_strict_core with high P2: core shoulder
above gb_strict_core with stricter P2: strict core
```

The exact numeric values are data-driven, not guessed constants. A cache built from a representative cooling trajectory can be reused for force or temperature folders if the potential and morphology range are compatible. If a new system changes the GB parameters, particle shapes, density regime, or mesogen type mapping, rebuild the threshold prior.

## 7. Diagnostic Figures For The 2D Split

The threshold-prior stage writes:

| File | Purpose |
|---|---|
| `stream_hist_2d_gb_p2.tsv` | Numeric 2D histogram backing the figures. |
| `gb_strength_vs_p2_stream_hist.png` | Dense heatmap view of the full streaming 2D histogram. |
| `gb_strength_vs_p2_stream_dotgrid.png` | Bin-center dot-grid view, useful for nonzero-bin inspection. |
| `gb_strength_vs_p2_stream_hexbin.png` | Hexbin-style visual view, similar to the older audit plot style. |
| `lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png` | Threshold overlay for the actual lobe split. |
| `lobe_split_preview/gb_core_slice_hist.png` | GB slices for `P2 > 0.71` and `P2 > 0.80`, explaining `gb_core` and `gb_strict_core`. |

These figures should show why the threshold lines are placed. If the distribution has no visible separation, the thresholds are lower-confidence and should be reviewed with OVITO.

## 8. Contact Graph Construction

After thresholds are known, the main analysis builds a contact graph:

```text
strong edge = tier 2 or higher
gray edge   = tier 1 support edge
```

Strong edges create seed components. Gray/support edges can attach particles to existing seeds but should not by themselves create large robust aggregation. Local contacts are retained, not deleted. Their role is tracked separately through `delta_s`, local-pair tables, and nonlocal-support counts.

This is why adjacent mesogens are not ignored. Adjacent local bundles remain visible as weak local domains. They are upgraded to robust domains only when there is additional evidence such as size, internal alignment, persistence, nonlocal support, or parameter stability.

## 9. Domain Rules

A domain is a local mesogen bundle. LC-Pearl first finds connected components from strong edges, then attaches supported gray/local particles.

Robust-domain evidence includes:

| Evidence | Meaning |
|---|---|
| Size | Default robust size floor is `n_min = 3`. Dimers are kept but are not automatically robust. |
| Internal alignment | Domain-level `S2` must pass the robust orientation rule when enabled. |
| Persistence | Domain track must survive beyond transient single-frame noise when persistence is used. |
| Nonlocal support | At least one non-adjacent support contact can be required or audited. |
| Parameter stability | Domain overlap under threshold perturbation is tracked as stability evidence. |

`lc_state` in OVITO summarizes domain state:

```text
0 = no domain
1 = weak domain
2 = robust domain
```

## 10. Pearl Rules

A pearl is a 3D bead-like object made from one or more robust domains. It is not the same as a larger domain.

Domain merging into a pearl requires:

- small 3D envelope gap;
- enough cross-domain contacts;
- boundary support on both domains, not just one bridge edge;
- merged aspect ratio below `pearl_max_aspect_ratio`.

This distinction is necessary because a large domain means one connected local contact bundle, while a pearl means a compact bead-like spatial object. Several robust domains can occupy the same bead and should be counted as one pearl if the 3D criteria are satisfied.

## 11. OVITO Scalar Logic

The most important v2 visualization field is:

```text
lc_aggregation_tier
```

It is a particle-level maximum over incident contact tiers:

| `lc_aggregation_tier` | Meaning in OVITO | Suggested color meaning |
|---:|---|---|
| 0 | No accepted incident contact | background / unaggregated |
| 1 | At least one weak/transition support contact | weak boundary or possible local support |
| 2 | At least one aggregation/strong contact | aggregated mesogen |
| 3 | At least one core-shoulder contact | strong high-P2 core region |
| 4 | At least one strict-core contact | strongest core subset |

Use `lc_aggregation_tier` first for quick inspection, then switch to:

- `lc_domain` to see local bundle membership;
- `lc_pearl` to see bead membership;
- `lc_cluster` to see simple spatial cluster labels;
- `lc_core_contact_degree` and `lc_strict_core_contact_degree` to see how many core contacts each mesogen has;
- `lc_max_gb_strength` and `lc_mean_gb_strength` to distinguish isolated strong contacts from broad attraction.

## 12. New v2.1 Domain-Size Count Statistic

`scripts/lc_domain_size_counts.py` adds:

```text
diagnostics/domain_size_frame_counts.tsv
diagnostics/domain_size_vs_domain_count.png
```

This statistic is not a global histogram. Each point in the plot is one `(frame, domain size)` group:

```text
x = domain size
y = number of domains of that size in that frame
```

This makes it possible to see whether a trajectory has many repeated small domains in the same frame, a few large domains, or mixed fragmentation states.

## 13. Main Architecture

| Module | Role |
|---|---|
| `lc_pearl_here.py` | Current-directory launcher. Creates and uses `lc_pearl_preflight/`. |
| `lc_pearl_cli.py` | TOML-driven command builder for reproducible runs. |
| `scripts/lc_domain_pearl_pipeline.py` | Orchestrates validation, threshold prior, main analysis, and diagnostics. |
| `scripts/lc_threshold_prior.py` | Streaming 2D threshold-prior module. |
| `scripts/liquid_crystal_aggregation.py` | Main contact/domain/pearl/mechanics/OVITO analysis engine. |
| `scripts/lc_domain_size_counts.py` | Standalone and pipeline-integrated domain-size count diagnostic. |
| `scripts/gb_pair_audit.py` | Optional legacy/targeted pair audit; no longer the primary threshold engine. |
| `scripts/lammps_pairlocal_compare.py` | Optional LAMMPS comparison/validation helper. |
| `scripts/lc_topology_from_data.py` | Converts topology data into local/excluded pair tables when available. |

## 14. Reliability Rules

Treat a run as reliable only when:

- `verified_potential.json` matches the current potential setup;
- `global_thresholds.json` was generated before the main run or deliberately reused;
- the main output reports the same applied thresholds as the prior;
- 2D split figures are physically interpretable;
- OVITO `lc_aggregation_tier`, `lc_domain`, and `lc_pearl` agree with representative frames;
- topology-dependent local/nonlocal claims use a reliable topology file or are labeled as inferred.

## 中文摘要

LC-Pearl v2 的核心不是简单距离聚类，也不是单独取向聚类，而是：

```text
GB strength x P2 二维分布切割
  -> 0/1/2/3/4 多层 contact tier
  -> weak/robust domain
  -> 3D pearl
  -> OVITO lc_aggregation_tier 和统计图
```

其中 `gb_off` 到 `gb_on` 是过渡/weak support，`gb_on` 以上是主聚集边，`gb_core + P2>0.71` 是 core shoulder 可视化层，`gb_strict_core + P2>0.80` 是更强的 strict core 可视化层。`lc_aggregation_tier` 是 OVITO 中最推荐优先看的上色字段。
