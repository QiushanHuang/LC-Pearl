# LC-Pearl v2.1.0

LC-Pearl is a reproducible pipeline for liquid-crystal mesogen aggregation analysis in LAMMPS dump trajectories.

See the main bilingual README and full Markdown manuals:

- `README.md`
- `MANUAL.md`
- `docs/lc_pearl_algorithm_details.md`
- `docs/lc_pearl_user_guide.md`
- `docs/lc_pearl_operation_manual.md`

Core workflow:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

LC-Pearl validates Gay-Berne potential reconstruction, builds a streaming 2D `GB strength x P2` threshold prior, partitions contacts into none / weak-transition / aggregation / core shoulder / strict core tiers, exposes `lc_aggregation_tier` for OVITO, then computes weak/robust domains, 3D pearls, mechanics, labels, and diagnostics.
