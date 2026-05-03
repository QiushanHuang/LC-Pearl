# LC-Pearl v2.1.0

LC-Pearl is a reproducible pipeline for liquid-crystal mesogen aggregation analysis in LAMMPS dump trajectories.

See the main bilingual README and full manual:

- `README.md`
- `MANUAL.md`
- `docs/lc_pearl_user_guide.html`
- `docs/lc_pearl_algorithm_reference.html`

Core workflow:

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

LC-Pearl validates Gay-Berne potential reconstruction, builds a streaming 2D `GB strength x P2` threshold prior, then computes weak/robust domains, 3D pearls, mechanics, OVITO labels, and diagnostics.
