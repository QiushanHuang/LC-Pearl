# Validation Cache

LC-Pearl writes these files automatically:

- `verified_potential.json`: portable validation cache.
- `run0_validation/`: detailed LAMMPS run 0 and microstate validation outputs.

You can copy this folder with the rest of `lc_pearl_preflight/` to a similar dump directory. LC-Pearl will reuse it only if the fingerprint still matches.
