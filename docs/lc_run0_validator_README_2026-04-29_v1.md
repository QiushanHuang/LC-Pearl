# LC LAMMPS Run0 Validator v1

This is an immutable process note for the standalone run0 validator added on 2026-04-29.

## Purpose

`lc_lammps_run0_validate.py` validates the Python type1-type1 Gay-Berne energy reconstruction against LAMMPS `run 0`.

It intentionally validates only the type-1 mesogen pair energy used by the aggregation algorithm. It does not validate sphere/anchor interactions, bond energy, rigid-body integration, or domain/pearl thresholds.

## Required Inputs

The validator is standalone. It can be run in the simulation directory, in a scratch directory, or on a copied representative dump frame. It needs:

- A LAMMPS dump frame containing `id`, `type`, coordinates (`xu/yu/zu` preferred, `x/y/z` accepted), `quatw/quati/quatj/quatk`, and `shapex/shapey/shapez`.
- The LAMMPS input file or parameter file containing the relevant `pair_style gayberne`, `pair_coeff 1 1`, `pair_modify`, variables used by those lines, and `special_bonds`.
- A working LAMMPS executable. Your usual command can be expressed as `--lammps-executable /Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp --mpi-prefix "mpiexec -np 4"`.

Optional but recommended:

- `--frame-index N` to validate a representative nontrivial frame, not necessarily the first frame.
- `--local-pair-file` and `--exclude-pair-file` from `lc_topology_prepare.py` so the verified-potential cache fingerprint records the same topology context used by analysis.
- `--min-pairs`, `--min-attractive-pairs`, and `--min-abs-python-total` to prevent an empty or near-zero frame from being accepted as proof.
- `--microstate-sample-count 12` controls how many two-particle pair microstates are sampled by default.
- `--microstate-sample-percent P` samples at least `P%` of candidate pairs. Because the default count is 12, use `--microstate-sample-count 0 --microstate-sample-percent P` if you want a pure percentage below 12 pairs. Use `--microstate-sample-percent 100` for full pair sampling.

## Why It Does Not Use `pair/local eng`

Common LAMMPS `pair_gayberne` builds do not expose per-pair energy through `compute pair/local eng`. The validator therefore builds a temporary type1-only system from one dump frame and compares:

```text
Python sum over type1-type1 GB pair energies
vs.
LAMMPS pair energy after run 0, with `thermo_modify norm no`
```

It then samples candidate pairs and runs isolated two-particle LAMMPS microstates. Each microstate has exactly one type1-type1 pair, so it checks individual pair energies without relying on `pair/local eng`.

Default behavior is total-energy validation plus 12 sampled pair microstates. To fully sample every candidate pair, use:

```bash
--microstate-sample-percent 100
```

For pure fractional sampling, for example 10% without the default 12-pair floor, use:

```bash
--microstate-sample-count 0 --microstate-sample-percent 10
```

Full-sampling validation command:

```bash
/Users/joshua/Desktop/MD/venv/bin/python /Users/joshua/Desktop/MD/Coding/lc_lammps_run0_validate.py \
  --dump-file traj.force_clamp_aligned.0.dump \
  --gb-param-file ../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp \
  --output-root lc_run0_validation_full \
  --local-pair-file lc_topology_pairs/local_special_pairs.tsv \
  --exclude-pair-file lc_topology_pairs/exclude_pairs.tsv \
  --lammps-executable /Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp \
  --mpi-prefix "mpiexec -np 4" \
  --microstate-sample-percent 100 \
  --max-abs-delta 1e-5
```

Internally, each sampled pair is rebuilt from the Python minimum-image pair vector `r12` in a padded nonperiodic two-particle box. This avoids confusing pair-level validation with the original frame's periodic images or unrelated neighbors.

## Dry Run

Use dry-run first to generate the LAMMPS input and Python pair table without claiming validation:

```bash
/Users/joshua/Desktop/MD/venv/bin/python /Users/joshua/Desktop/MD/Coding/lc_lammps_run0_validate.py \
  --dump-file traj.force_clamp_aligned.0.dump \
  --gb-param-file ../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp \
  --output-root lc_run0_validation \
  --local-pair-file lc_topology_pairs/local_special_pairs.tsv \
  --exclude-pair-file lc_topology_pairs/exclude_pairs.tsv \
  --lammps-executable /Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp \
  --mpi-prefix "mpiexec -np 4" \
  --dry-run
```

Dry-run writes:

```text
lc_run0_validation/
  validate_type1_gb_run0.in
  validation_atom_id_map.tsv
  python_type1_pair_energies.tsv
  run0_validation_manifest.json
  verified_potential.json   # status = validation_required
```

`microstate_checks/` is written only during actual validation, not during dry-run.

The generated LAMMPS input deliberately uses two LAMMPS-specific safeguards:

- Dump `quatw/quati/quatj/quatk` is converted to the axis-angle form required by `set atom quat`.
- `thermo_modify norm no` is placed after `thermo_style` so the reported energy is total energy, not per-atom normalized energy.

## Actual Validation

```bash
/Users/joshua/Desktop/MD/venv/bin/python /Users/joshua/Desktop/MD/Coding/lc_lammps_run0_validate.py \
  --dump-file traj.force_clamp_aligned.0.dump \
  --gb-param-file ../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp \
  --output-root lc_run0_validation \
  --local-pair-file lc_topology_pairs/local_special_pairs.tsv \
  --exclude-pair-file lc_topology_pairs/exclude_pairs.tsv \
  --lammps-executable /Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp \
  --mpi-prefix "mpiexec -np 4" \
  --microstate-sample-count 12 \
  --max-abs-delta 1e-5
```

Only if both the whole-frame total-energy comparison and all sampled microstate checks pass does `verified_potential.json` get:

```json
{
  "status": "validated"
}
```

The main pipeline can then reuse it:

```bash
--verified-potential-file lc_run0_validation/verified_potential.json \
--gb-param-file ../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp \
-- \
traj.force_clamp_aligned.0.dump \
--contact-mode gayberne \
--local-pair-file lc_topology_pairs/local_special_pairs.tsv \
--exclude-pair-file lc_topology_pairs/exclude_pairs.tsv
```

The topology pair-file flags must match between validation and analysis. If you validate without `--local-pair-file` / `--exclude-pair-file` and later analyze with them, the cache fingerprint will intentionally miss and request revalidation.

The cache acceptance fingerprint is based on potential settings, mesogen type, dump schema, topology pair-file hashes, analysis code hash, run0 validator code hash, and validation contract. The LAMMPS executable/version is recorded as provenance, but it is not required to match byte-for-byte when reusing a standalone artifact.

## Status Meaning

- `validated`: the selected frame was nontrivial, LAMMPS ran, total type1-type1 GB pair energy matched, and all sampled two-particle microstates matched within tolerance. The main pipeline may reuse this file.
- `validation_required`: dry-run or cache miss. This is not proof and should not be used as a trusted potential.
- `failed`: LAMMPS launch failed, the marker was missing, the frame was uninformative, or Python and LAMMPS energies disagreed. Do not use this artifact as validation.

`microstate_summary.json` records `passed`, `candidate_pair_count`, `selected_pair_count`, `passed_pair_count`, `failed_pair_count`, `sample_count`, `sample_percent`, `seed`, `max_abs_tolerance`, `max_abs_delta`, `mean_abs_delta`, `rmse`, `table`, `results_table`, `output_root`, `selection`, `worst_pair`, and per-pair `rows`.

## When To Revalidate

Revalidate when any of these changes:

- `pair_style gayberne`
- `pair_coeff 1 1`
- `pair_modify`
- LAMMPS executable/build/suffix. This is recorded in the artifact as provenance; if you intentionally switch LAMMPS builds, rerun validation manually.
- Python GB energy code
- dump quaternion or shape columns
- mesogen type or shape convention
- topology/exclusion files used in the analysis fingerprint

You do not need raw-energy revalidation only because `gb_on/gb_off/p2_cut` changes. Those are threshold calibration parameters, not potential parameters.

## Current Important Finding

The first actual run0 check on the F0.38 representative frame did not validate the Python GB reconstruction. That failure was useful: it exposed three validator/reconstruction issues that have now been fixed in code.

- Dump `shapex/shapey/shapez` are diameters, while `pair_gayberne` uses internal semiaxes; the Python energy path now converts dump shape values to semiaxes.
- `set atom quat` requires axis-angle input, not raw quaternion components; the validator now converts dump quaternions before writing the run0 input.
- LAMMPS thermo output can be per-atom normalized; the validator now writes `thermo_modify norm no` after `thermo_style`.

After these fixes, the F0.38 representative frame validated with:

```text
Python type1-type1 GB total = -95.14282713795289
LAMMPS run0 pair total      = -95.1428271379528
absolute delta              = 8.53e-14
sampled microstates         = 12/12 passed by default
```

## If Validation Fails

Use `python_type1_pair_energies.tsv`, `validate_type1_gb_run0.in`, `lammps_stdout.txt`, `lammps_stderr.txt`, `run0_validation_summary.json`, and `microstate_checks/microstate_pair_checks.tsv` to locate the source. A total-energy failure usually points to global reconstruction or LAMMPS input setup; a microstate failure points to a specific pair, quaternion, shape, cutoff, or shift mismatch.
