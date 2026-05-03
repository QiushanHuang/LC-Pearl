# LC-Pearl v2.1.0 Academic User Manual

LC-Pearl is a reproducible pipeline for quantifying liquid-crystal mesogen aggregation in LAMMPS dump trajectories. It was built for single-chain and chain-like liquid-crystal simulations where Gay-Berne attraction, local orientational order, chain connectivity, temperature, and external stretching compete to form pearl-necklace-like morphologies.

This manual is the canonical Markdown reference for the v2.1.0 release. The focused Markdown manuals are:

- `docs/lc_pearl_algorithm_details.md`
- `docs/lc_pearl_user_guide.md`
- `docs/lc_pearl_operation_manual.md`

Browser-ready HTML copies are retained for convenience, but the Markdown documents are the canonical release documentation.

## 1. Scope And Scientific Model

LC-Pearl separates three physical questions that are often mixed together in simple cluster analysis:

1. Pair contact: are two mesogens attractively and orientationally coupled?
2. Domain: do contacts form a local mesogen bundle?
3. Pearl: do one or more robust domains form the same compact 3D bead-like object?

The analysis members are type-1 ellipsoids by default. Type-2 spheres and type-3 anchors are not domain or pearl members, but they can support chain ordering, endpoint detection, stretch-axis inference, and topology reconstruction.

This distinction matters because adjacent chain segments can form real local bundles, but chain connectivity alone should not automatically create a large robust domain. LC-Pearl therefore keeps weak local domains while requiring additional evidence before a domain is classified as robust.

## 2. Pipeline Overview

The recommended v2 pipeline order is:

```text
LAMMPS dump files
  -> preflight input discovery
  -> Gay-Berne potential validation or cache reuse
  -> streaming 2D GB-strength x P2 threshold prior
  -> main mesogen contact graph
  -> weak/robust domain classification
  -> 3D pearl merging
  -> mechanics and OVITO outputs
  -> diagnostics and post-run review
```

The key design rule is that threshold selection happens before the main analysis. LC-Pearl does not first run the full domain/pearl analysis using fixed default thresholds and then retroactively recommend new thresholds. The main analysis receives the recommended prior values before contact graph construction.

The core V2 algorithm is the 2D `GB strength x P2` split followed by a multi-tier contact partition:

| Tier | Meaning | Rule |
|---:|---|---|
| 0 | no accepted contact | below gray/support gate |
| 1 | weak or transition support | `gb_strength >= gb_off_strength` and `P2 >= p2_cut`, but not strong |
| 2 | aggregation/strong contact | `gb_strength >= gb_on_strength` and `P2 >= p2_cut` |
| 3 | core shoulder | `gb_strength >= gb_core_strength` and `P2 > 0.71` |
| 4 | strict core | `gb_strength >= gb_strict_core_strength` and `P2 > 0.80` |

This tier model is also exposed in OVITO as `lc_aggregation_tier`. The field is particle-level and equals the maximum tier among incident contacts for each mesogen. It is the recommended first color-coding field for visual review.

## 3. Installation

Use Python 3.11 or newer. The runtime dependencies are `numpy` and `matplotlib`.

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 -m pip install -e .
```

Direct source-tree execution also works:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

## 4. Recommended Current-Directory Workflow

Enter a folder containing LAMMPS dump files and run:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

On first use, the wrapper creates `lc_pearl_preflight/` and reports missing inputs. Put the required source files there:

| Preflight path | What you provide | What LC-Pearl extracts |
|---|---|---|
| `lc_pearl_preflight/lammps/*.in` or `*.lmp` | Original LAMMPS input file containing `pair_style gayberne` and type 1-1 `pair_coeff` | Gay-Berne parameters used by Python energy reconstruction |
| `lc_pearl_preflight/lammps/lammps_executable.txt` | One line with the LAMMPS executable path | `run 0` / microstate validation |
| `lc_pearl_preflight/topology/*.data` | Optional LAMMPS data file with `Atoms` and `Bonds` sections | Local-pair and excluded-pair tables |
| `lc_pearl_preflight/validation/verified_potential.json` | Normally generated automatically | Reusable potential-validation cache |
| `lc_pearl_preflight/thresholds/global_thresholds.json` | Normally generated automatically | Reusable threshold-prior cache |

Run the complete workflow:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

Run validation only:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py validate
```

Run analysis using existing caches:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py run
```

## 5. Config-Driven Workflow

Edit `configs/quick_run.toml`, then preview the generated command:

```bash
cd /Users/joshua/Desktop/MD/LC-Pearl
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py print-run --config configs/quick_run.toml
```

Run:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py run --config configs/quick_run.toml
```

Validate only:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 lc_pearl_cli.py validate --config configs/quick_validate.toml
```

Use TOML configs when you want reproducible batch runs across different force, temperature, density, or chain-length folders.

## 6. Potential Validation

In `contact_mode = "gayberne"`, LC-Pearl reconstructs type-1 ellipsoid pair energies from dump coordinates, quaternions, shape axes, and LAMMPS Gay-Berne parameters. The normalized attraction strength is

```text
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
```

where `U_GB,ij` is the reconstructed pair energy and `U_well,ij` is the orientation-dependent attractive well depth. This normalization makes a contact score comparable across different relative orientations.

Potential validation is used to check the Python reconstruction against LAMMPS. The validation layer can use a representative dump frame and sampled two-particle microstates. After successful validation, LC-Pearl writes:

```text
lc_pearl_preflight/validation/verified_potential.json
```

Reuse this cache when the Gay-Berne potential, shape parameters, unit convention, and type mapping have not changed.

Revalidate when:

- `pair_style`, `pair_coeff`, units, shape axes, or quaternion interpretation changes.
- A new LAMMPS build changes aspherical pair behavior.
- You analyze a different mesogen type mapping.
- You suspect the current `verified_potential.json` came from a different input deck.

You usually do not need to revalidate just because force, stretch length, or temperature changes, provided the potential settings are identical.

## 7. Threshold Prior

The threshold prior is computed from the full or strided candidate-pair distribution before the main analysis. For each candidate type-1 ellipsoid pair, LC-Pearl computes:

```text
q_ij = |u_i dot u_j|
P2_ij = (3 q_ij^2 - 1) / 2
gb_strength_ij = max(0, -U_GB,ij / U_well,ij)
```

`P2 = 1` means parallel or antiparallel long axes, and `P2 = -0.5` means perpendicular axes. The absolute dot product is used because ellipsoid long axes are head-tail symmetric.

The streaming prior accumulates a 2D histogram of `gb_strength x P2`. It then estimates:

| Parameter | Source | Interpretation |
|---|---|---|
| `p2_cut` | Valley or stable split separating low-orientation background and high-orientation lobe | Pair-orientation gate for contact edges |
| `gb_off_strength` | Conservative shoulder on the weak-contact side | Gray/support edge threshold |
| `gb_on_strength` | Valley between weak-contact and strong-attraction high-P2 lobes | Strong edge threshold |
| `p2_core_cut = 0.71` | Fixed diagnostic high-orientation gate | Core-shoulder visualization gate |
| `gb_core_strength` | Strong-lobe left shoulder at `P2 > 0.71` | Core-shoulder contact tier |
| `p2_strict_core_cut = 0.80` | Fixed stricter high-orientation gate | Strict-core visualization gate |
| `gb_strict_core_strength` | Strong-lobe high-strength quantile at `P2 > 0.80` | Strict-core contact tier |

The core tiers are diagnostic. They do not replace `gb_on_strength` and do not change robust-domain or pearl definitions. They make OVITO coloring and review more interpretable by distinguishing no/transition/aggregation/core/strict-core contact classes.

The main threshold-prior artifact is:

```text
lc_pearl_preflight/thresholds/global_thresholds.json
```

The most useful review figures are:

- `gb_strength_vs_p2_stream_hist.png`
- `gb_strength_vs_p2_stream_dotgrid.png`
- `gb_strength_vs_p2_stream_hexbin.png`
- `lobe_split_preview/gb_strength_vs_p2_stream_lobe_split_dotgrid.png`
- `lobe_split_preview/gb_core_slice_hist.png`

Use `threshold_prior.global_frame_stride = 1` for full frame sampling. Use `10` or `100` for quick checks on directories with many one-frame dump files. For final parameters, prefer a representative trajectory covering the relevant temperature and morphology range.

## 8. Contact Graph

The main contact graph uses two edge classes:

```text
strong edge: gb_strength >= gb_on_strength  and P2 >= p2_cut
gray edge:   gb_strength >= gb_off_strength and P2 >= p2_cut
```

Strong edges form candidate seeds. Gray edges provide support and attachment. Local chain contacts are retained but separately marked using `delta_s` and topology information.

`delta_s` is the mesogen sequence separation. If reliable topology exists, LC-Pearl uses preflight-generated local and excluded pair tables. If topology is missing, local/nonlocal status is inferred from chain order and IDs, which is less exact.

## 9. Weak And Robust Domains

A domain is a local mesogen bundle inferred from the contact graph. LC-Pearl distinguishes weak local domains from robust domains.

Weak local domains:

- retain adjacent or near-adjacent liquid-crystal bundling;
- are important for OVITO interpretation and transitional states;
- are not discarded just because they are local.

Robust domains require evidence. The default evidence categories are:

| Evidence | Typical parameter | Meaning |
|---|---|---|
| Size | `n_min = 3` | More than a pure dimer |
| Orientation | `robust_min_s2 = 0.70` | Internally aligned mesogen bundle |
| Persistence | `domain_min_lifetime = 2` | Survives beyond a single transient frame |
| Nonlocal support | topology/local-pair logic | Not only adjacent-chain support |
| Parameter stability | perturbation scan | Stable under stricter thresholds |

The classification rule is deliberately conservative: adjacent liquid-crystal aggregation is preserved as weak evidence, but a robust domain requires more than local adjacency alone.

## 10. Pearl Definition

A pearl is not just a larger domain. A domain answers whether a local mesogen bundle exists. A pearl answers whether one or more robust domains occupy the same compact 3D bead-like region.

Two robust domains merge into one pearl only if they satisfy all relevant criteria:

- 3D envelope gap is below `pearl_gap_cut`;
- cross-domain contact count is at least `pearl_min_cross_contacts`;
- both domain boundaries have support, avoiding a single bridge edge;
- merged object remains bead-like with `aspect_ratio <= pearl_max_aspect_ratio`.

Default pearl parameters:

| Parameter | Default | Purpose |
|---|---:|---|
| `pearl_gap_cut` | `auto` | Spatial merge gap |
| `pearl_min_cross_contacts` | `2` | Avoid one-edge bridge artifacts |
| `pearl_min_boundary_particles` | `2` | Require support on both sides |
| `pearl_max_aspect_ratio` | `3.0` | Reject elongated connector-like objects |

Use `pearl_count`, `domains per pearl`, `largest_pearl_fraction`, and connector lengths together. Do not interpret `pearl_count` alone.

## 11. Mechanics Layer

LC-Pearl computes axial and shape observables:

| Quantity | Meaning | Use |
|---|---|---|
| `L_parallel` | projected end-to-end length along stretch axis | Main force-extension observable |
| `Rg_parallel` | parallel radius-of-gyration component | Shape response |
| `Rg_perp` | perpendicular radius-of-gyration component | Lateral compaction |
| `S2_force` | mesogen alignment relative to stretch axis | Stretch-induced orientational order |

For multiple constant-force simulations, compare plateau or steady-state averages across force. If appropriate, estimate compliance by

```text
C_parallel = Delta <L_parallel> / Delta F
```

`Rg_parallel` and `Rg_perp` are supporting shape descriptors. They should not replace `L_parallel` as the primary extension measure.

## 12. Main Output Files

The default output root is:

```text
lc_domain_pearl_v2_output/
```

Important outputs:

| Output | Purpose |
|---|---|
| `aggregation_timeseries.tsv` | Frame-level domain, pearl, cluster, contact, and mechanics table |
| `per_file/*_aggregation.tsv` | Per-input-file time series |
| `per_file/*_summary.json` | Compact per-file summary |
| `per_file/*_lc_labels.dump` | OVITO particle labels |
| `per_file/*_lc_cluster_envelopes.dump` | Visual cluster envelope particles |
| `per_file/*_lc_contact_edges.dump` | Contact-edge records for visualization/debugging |
| `per_file/*_lc_contact_segments.vtk` | OVITO-loadable contact line segments |
| `diagnostics/domain_diagnostics.tsv` | Domain-level evidence and classification table |
| `diagnostics/domain_size_vs_s2.png` | Scatter of domain size vs internal S2 |
| `diagnostics/domain_size_frame_counts.tsv` | Per-frame, per-size domain counts |
| `diagnostics/domain_size_vs_domain_count.png` | Scatter of domain size vs domain count per frame |
| `diagnostics/diagnostic_summary.json` | Applied thresholds, counts, and provenance |
| `diagnostics/pearl_candidate_diagnostics.tsv` | Domain-domain pearl merge candidates |

Large debug outputs are disabled by default:

```toml
[analysis]
edge_diagnostics_table = "off"  # off, sample, full
write_frame_jsonl = false
accepted_edge_audit = false
```

Enable them only for targeted debugging.

## 13. OVITO Visualization

Open `per_file/*_lc_labels.dump` in OVITO. Useful scalar properties:

| Property | Meaning | Recommended use |
|---|---|---|
| `lc_cluster` | center-distance visual cluster id | Fast spatial cluster view |
| `lc_cluster_size` | visual cluster size | Identify compact groups |
| `lc_aggregation_tier` | 0 no contact, 1 weak/transition, 2 aggregation, 3 core shoulder, 4 strict core | Best quick color-coding field |
| `lc_contact_degree` | accepted contacts incident on a mesogen | Identify central contact hubs |
| `lc_min_pair_energy` | most attractive pair energy | Locate strong local attraction |
| `lc_mean_pair_energy` | mean pair energy around a mesogen | Review local energetic environment |
| `lc_max_gb_strength` | strongest normalized attraction | Review strongest contact |
| `lc_mean_gb_strength` | mean normalized attraction | Distinguish isolated strong edges from broad attraction |
| `lc_core_contact_degree` | contacts passing `gb_core_strength` and `P2 > 0.71` | Core-shoulder regions |
| `lc_core_nonlocal_degree` | nonlocal core-shoulder contacts | Nonlocal aggregation support |
| `lc_is_core_particle` | has at least one core-shoulder contact | Binary core-shoulder coloring |
| `lc_strict_core_contact_degree` | contacts passing `gb_strict_core_strength` and `P2 > 0.80` | Strongest core contacts |
| `lc_is_strict_core_particle` | has at least one strict-core contact | Binary strict-core coloring |
| `lc_domain` | weak or robust domain id | Local bundle membership |
| `lc_pearl` | 3D bead id | Pearl-necklace bead membership |
| `lc_state` | 0 none, 1 weak domain, 2 robust domain | Compact state view |

For line contacts, load `*_lc_contact_segments.vtk` as an additional pipeline. If OVITO Pro is unavailable, inspect `*_lc_labels.dump` first, then open representative VTK files separately.

## 14. Analysis Strategies

### 14.1 Single trajectory

Use:

- `N_domain(t)` and `N_pearl(t)` to detect fragmentation or coarsening;
- `largest_domain_fraction` and `largest_pearl_fraction` to distinguish one dominant aggregate from many small beads;
- `local_edge_fraction` and `nonlocal_edge_fraction` to separate chain-adjacent bundling from nonlocal aggregation;
- `domain_size_vs_s2.png` to review domain quality;
- `domain_size_vs_domain_count.png` to see how many same-sized domains occur in each frame.

### 14.2 Force series

For multiple constant-force simulations:

- use the same validated potential cache if the potential is unchanged;
- use the same threshold prior if it was built from a representative morphology and temperature range;
- compare steady-state windows of `L_parallel`, `pearl_count`, `largest_pearl_fraction`, and `nonlocal_edge_fraction`;
- report uncertainty across time windows or replicate simulations.

### 14.3 Temperature or cooling series

For cooling trajectories:

- build the threshold prior from a trajectory covering the full high-temperature to low-temperature transition when possible;
- inspect `gb_strength_vs_p2_stream_*` plots for lobe separability;
- compare aggregation onset temperature using `robust_domain_count`, `pearl_count`, and `core_contact_fraction`.

### 14.4 Debugging a suspicious frame

1. Find the timestep from `aggregation_timeseries.tsv`.
2. Open the corresponding label dump in OVITO.
3. Color by `lc_aggregation_tier`, then by `lc_domain`, then by `lc_pearl`.
4. Load contact segments for the same frame if needed.
5. Inspect `domain_diagnostics.tsv` rows for that timestep.
6. If edge-level evidence is needed, rerun with `edge_diagnostics_table = "sample"` or `"full"` for a smaller frame subset.

## 15. Parameter Reference

| Parameter | Default or source | Role | Notes |
|---|---|---|---|
| `contact_mode` | `gayberne` | Contact score source | Recommended production mode |
| `mesogen_type` | `1` | Aggregation members | Type-1 ellipsoids |
| `anchor_types` | `3` | Stretch-axis/end support | Not aggregation members |
| `gb_off_strength` | threshold prior | Gray/support edge threshold | Conservative support contacts |
| `gb_on_strength` | threshold prior | Strong edge threshold | Main energy gate |
| `p2_cut` | threshold prior | Pair orientation gate | Main orientation gate |
| `gb_core_strength` | threshold prior | Diagnostic core-shoulder GB gate | Visualization/audit only |
| `p2_core_cut` | `0.71` | Diagnostic core-shoulder P2 gate | Roughly high alignment |
| `gb_strict_core_strength` | threshold prior | Diagnostic strict-core GB gate | Visualization/audit only |
| `p2_strict_core_cut` | `0.80` | Diagnostic strict-core P2 gate | Stricter alignment |
| `s_excl` | `1` | Local sequence exclusion/support range | Use topology when available |
| `n_min` | `3` | Robust size evidence | Avoid pure dimers as robust cores |
| `robust_min_s2` | `0.70` | Robust orientation evidence | Internal domain order |
| `domain_min_lifetime` | `2` | Persistence evidence | Processed-frame lifetime |
| `stable_overlap_fraction` | implementation default | Parameter stability evidence | Perturbation overlap criterion |
| `pearl_gap_cut` | `auto` | 3D pearl merge distance | Derived from contact length scale if not set |
| `pearl_min_cross_contacts` | `2` | Domain-domain support count | Avoid single bridge artifacts |
| `pearl_min_boundary_particles` | `2` | Boundary support on both domains | More robust pearl merging |
| `pearl_max_aspect_ratio` | `3.0` | Compact bead constraint | Reject long connector-like merges |
| `threshold_prior.global_frame_stride` | `1` | Prior frame sampling | `1` = all candidate frames |
| `threshold_prior.workers` | `auto` | Prior parallelism | Capped by CPU, chunk count, environment |
| `analysis.workers` | `auto` | Main analysis parallelism | Cross-file tracking can require sequential processing |
| `edge_diagnostics_table` | `off` | Edge table output | Use `sample` or `full` only for debug |
| `write_frame_jsonl` | `false` | Full JSON debug output | Can be very large |

## 16. Reliability Checklist

Before using results in a group meeting or paper, verify:

1. `verified_potential.json` exists and matches the current LAMMPS potential.
2. `global_thresholds.json` exists, has the expected schema, and was generated from representative data.
3. `gb_strength_vs_p2_stream_*` figures show a physically interpretable distribution.
4. `diagnostics/diagnostic_summary.json` reports the same thresholds you expect.
5. OVITO label colors agree with representative snapshots.
6. `domain_size_vs_s2.png` does not show unexpected low-S2 robust domains.
7. `domain_size_vs_domain_count.png` is consistent with the expected number of same-size domains per frame.
8. `aggregation_timeseries.tsv` trends are consistent with visual inspection.
9. Any topology-dependent local/nonlocal claims are supported by a reliable data file or clearly labeled as inferred.

## 17. Known Limits

- Exact LAMMPS `special_bonds` semantics require reliable topology and, for final validation, LAMMPS comparison.
- Threshold priors are statistical. If the 2D distribution lacks a clear lobe structure, recommended thresholds should be treated as lower-confidence.
- Pearl merging is a 3D bead definition. It intentionally differs from domain construction.
- Diagnostic core tiers improve visualization but do not change the main robust-domain/pearl classification.
- Full edge tables and JSONL records can be extremely large; defaults keep outputs compact.

## 18. Citation

Use `CITATION.cff` for software citation. In publications or group reports, also state:

- LC-Pearl version;
- LAMMPS potential validation status;
- threshold-prior source trajectory and sampling settings;
- applied `gb_off_strength`, `gb_on_strength`, `p2_cut`, `gb_core_strength`, and `gb_strict_core_strength`;
- robust-domain evidence settings;
- pearl merge settings;
- whether topology data were used.
