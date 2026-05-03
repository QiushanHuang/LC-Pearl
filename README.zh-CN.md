# LC-Pearl v2.1.0

LC-Pearl 是用于 LAMMPS dump 轨迹中液晶 mesogen 聚集分析的可复现 pipeline。

主要文档：

- `README.md`
- `MANUAL.md`
- `docs/lc_pearl_algorithm_details.md`
- `docs/lc_pearl_user_guide.md`
- `docs/lc_pearl_operation_manual.md`

推荐运行：

```bash
cd /path/to/dump_output_folder
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

LC-Pearl 会验证 Gay-Berne potential 重算，生成 streaming 2D `GB strength x P2` 阈值先验，将 contact 分为 none / weak-transition / aggregation / core shoulder / strict core 多层，并在 OVITO 中输出 `lc_aggregation_tier`，然后计算 weak/robust domain、3D pearl、力学量、label 和 diagnostics。
