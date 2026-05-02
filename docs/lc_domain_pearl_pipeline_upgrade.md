# LC Domain/Pearl Pipeline Upgrade Notes

This note records the current project-level plan for making the LC aggregation analysis auditable rather than manually thresholded.

## Current Pipeline

1. `Coding/liquid_crystal_aggregation.py`
   - Runs domain/pearl/mechanics analysis.
   - Writes `aggregation_timeseries.tsv`, OVITO label dumps, contact edge tables, and diagnostics.
   - Supports file-level multiprocessing only when cross-file tracking is disabled.

2. `Coding/gb_pair_audit.py`
   - New calibration layer.
   - Recomputes unscreened type-1/type-1 Gay-Berne candidate pairs inside the GB cutoff.
   - Writes `gb_candidate_pairs.tsv` and `threshold_auto/threshold_recommendations.json`.

3. `Coding/lc_threshold_recommend.py`
   - Recommends `gb_off_strength`, `gb_on_strength`, `p2_cut`, and `robust_min_s2` from diagnostics or candidate pair data.
   - Uses deterministic 1D k-means boundaries and writes readable density plots instead of unreadable full-screen scatter plots.

4. `Coding/lc_topology_exclusions.py`
   - Converts a LAMMPS data file `Bonds` section into explicit `exclude_pairs.tsv` and `local_special_pairs.tsv`.
   - Intended to replace heuristic `adjacent_id_gap` when a real topology is available.

5. `Coding/lammps_pairlocal_compare.py`
   - Compares Python GB audit energy tables against a trusted LAMMPS-generated pair table.
   - Also writes a template for a LAMMPS validation input, but the final validation route must be checked against the exact LAMMPS build.

6. `Coding/lc_ovito_basic_merge.py`
   - Creates a single OVITO Basic-compatible dump containing original labeled atoms, halo marker atoms, and dotted contact markers.
   - Avoids `Add to scene`, which requires OVITO Pro.

7. `Coding/lc_domain_pearl_pipeline.py`
   - Thin wrapper that runs the legacy analysis, threshold recommendation, and optional GB candidate audit.

## Important Precision Corrections

- `shapex/shapey/shapez` in the dump are diameters. The GB reconstruction now converts them to semiaxes before evaluating the pair formula.
- LAMMPS `pair_gayberne` uses an effective `upsilon/2` exponent internally. Parsed LAMMPS input now stores that effective exponent.
- The Python GB energy is still a reconstruction and should be validated against LAMMPS `run 0` or an equivalent trusted pair list before being described as exact.

## Threshold Policy

The analysis should not require visual judgment of `0.12 / 0.30`.

Recommended workflow:

```bash
python /Users/joshua/Desktop/MD/Coding/gb_pair_audit.py \
  traj.force_clamp_aligned.0.dump traj.force_clamp_aligned.1013000.dump \
  --gb-param-file /path/to/in.lmp \
  --output-root gb_pair_audit
```

Then inspect:

- `gb_pair_audit/gb_candidate_pairs.tsv`
- `gb_pair_audit/threshold_auto/threshold_recommendations.json`
- `gb_pair_audit/threshold_auto/gb_strength_vs_p2_hexbin_recommended.png`

The recommended thresholds are data-derived analysis parameters, not universal physical constants.

## Topology / special_bonds Policy

Do not treat sequence adjacency as LAMMPS `special_bonds`.

Preferred order:

1. Exact LAMMPS data file or restart converted to data.
2. `lc_topology_exclusions.py from-data topology.data`.
3. Pass `--exclude-pair-file exclude_pairs.tsv` and, if needed, `--local-pair-file local_special_pairs.tsv` to the analysis.
4. Only if no topology exists, fall back to `s_excl` and inferred chain order, and mark the result as heuristic.

For the current a-E-a-7S single-chain system, type-1 ellipsoids are not directly in harmonic bonds, so `special_bonds` is expected to have little direct effect on type-1/type-1 GB contacts. The sequence local map is still important for classifying local-support vs nonlocal contacts.

## Recommended Future LAMMPS Dump Additions

These are not required for the current scripts but will improve future reliability:

- Add `mol` and explicit mesogen sequence index, e.g. `mesogen_s`.
- Save or export topology as a LAMMPS data file alongside each production folder.
- Keep `quatw quati quatj quatk` and `shapex shapey shapez` in every analysis dump.
- Record the exact LAMMPS input and resolved variables used for each output folder.

## Parallel Strategy

- Final robust-domain statistics with cross-file persistence: sequential over files.
- Quick visual scans: use `--no-track-across-files --workers auto`.
- GB pair audit: file-level multiprocessing is safe because it does not need persistence.
- Reporting/OVITO Basic merge: safe to parallelize per file in a later version.
