# LC Topology Preparation v1

This is an immutable process note for the topology preparation module added on 2026-04-29.

## Purpose

The aggregation members are still only type-1 mesogens. The topology module does not add sphere or anchor atoms to domains or pearls. It produces reusable helper files so the analysis can distinguish:

- true LAMMPS `special_bonds` exclusions,
- local chain support contacts,
- mesogen sequence order along the chain,
- rigid a-E-a molecule information used only for analysis mapping.

## Input

Preferred input is a LAMMPS `write_data` file containing `Atoms` and `Bonds`, for example:

```text
TOPOLOGY.rho030.data
```

Your heating-chain input already writes this file:

```lammps
write_data ${Tstar_tag}/${rho_tag}/TOPOLOGY.${rho_tag}.data
```

If only a restart exists, first write a conversion template:

```bash
/Users/joshua/Desktop/MD/venv/bin/python /Users/joshua/Desktop/MD/Coding/lc_topology_prepare.py write-restart-template \
  --restart Restart.GB.rho030.1000000 \
  --output write_topology_data.in
```

Then run LAMMPS:

```bash
mpiexec -np 4 /Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp -in write_topology_data.in
```

## Command

```bash
/Users/joshua/Desktop/MD/venv/bin/python /Users/joshua/Desktop/MD/Coding/lc_topology_prepare.py from-data \
  TOPOLOGY.rho030.data \
  --special-lj 0,1,1 \
  --mesogen-type 1 \
  --anchor-type 3 \
  --output-root lc_topology_pairs
```

## Output

```text
lc_topology_pairs/
  exclude_pairs.tsv
  local_special_pairs.tsv
  chain_index_map.tsv
  topology_manifest.json
```

`exclude_pairs.tsv` should be passed to:

```bash
--exclude-pair-file lc_topology_pairs/exclude_pairs.tsv
```

`local_special_pairs.tsv` should be passed to:

```bash
--local-pair-file lc_topology_pairs/local_special_pairs.tsv
```

`chain_index_map.tsv` is currently an analysis manifest. The main analysis still infers chain order internally, but this file is the stable target for future multi-chain support.

## Pipeline Placement

The topology files are forwarded after the wrapper's `--` separator, because they are arguments for `liquid_crystal_aggregation.py`:

```bash
/Users/joshua/Desktop/MD/venv/bin/python /Users/joshua/Desktop/MD/Coding/lc_domain_pearl_pipeline.py run \
  --output-root lc_domain_pearl_v1_output \
  --gb-param-file ../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp \
  --verified-potential-file lc_run0_validation/verified_potential.json \
  -- \
  . \
  --pattern "traj.force_clamp_aligned.*.dump" \
  --contact-mode gayberne \
  --local-pair-file lc_topology_pairs/local_special_pairs.tsv \
  --exclude-pair-file lc_topology_pairs/exclude_pairs.tsv
```

Use the same pair-file paths when running `lc_lammps_run0_validate.py`; otherwise the verified-potential fingerprint will not match the analysis run.

## Important Semantics

In your current model, a-E-a is rigid by shared molecule ID and `fix rigid/nvt/small`, not by harmonic bonds. Therefore rigid a-E-a links are useful for chain sequence mapping but must not be treated as LAMMPS `special_bonds`.

For your typical input:

```lammps
special_bonds lj 0.0 1.0 1.0
```

Only true 1-2 bonded pairs are excluded from LJ/GB nonbonded interaction. 1-3 and 1-4 pairs remain full nonbonded interactions and should not be put into `exclude_pairs.tsv`.
