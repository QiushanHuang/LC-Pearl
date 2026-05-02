# LAMMPS Inputs

Put these files here:

- Any original `.in` or `.lmp` file that contains the Gay-Berne potential settings. `lc_pearl_here.py` auto-detects the file with `pair_style gayberne` and `pair_coeff`, then copies it to the canonical internal name `gb_param_source.in`.
- `lammps_executable.txt`: required unless `[paths].lammps_executable` is set in `../lc_pearl_config.toml`.

`gb_param_source.in` is the standardized internal replacement for the old command-line `--gb-param-file`; you normally do not need to create that exact filename by hand.
