# Optional Topology Input

Recommended input:

- Put one raw LAMMPS data file here, for example `topology_for_lc_analysis.data`.
- The data file should contain `Atoms` and `Bonds` sections.
- `lc_pearl_here.py` reads the raw data file and automatically generates `local_pairs.tsv` and `exclude_pairs.tsv`.

Generated internal files:

- `local_pairs.tsv`: two columns, `atom_i atom_j`, for local chain-neighbor support pairs.
- `exclude_pairs.tsv`: two columns, `atom_i atom_j`, for pairs excluded from contact analysis.

You normally should not write those TSV files by hand. If they already exist, LC-Pearl uses them and includes their SHA256 hashes in the validation fingerprint.
