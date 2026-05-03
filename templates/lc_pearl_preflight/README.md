# LC-Pearl Preflight Folder

把整个 `lc_pearl_preflight/` 放进一个 dump 输出目录，然后在该 dump 目录运行：

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

默认行为是：如果 `validation/verified_potential.json` 存在且和当前参数匹配，就复用验证；否则先跑 LAMMPS run 0 / microstate 验证。然后检查 `thresholds/global_thresholds.json`：只有在 GB potential 已验证且 prior fingerprint 匹配时，才把其中的 `gb_off/gb_on/p2_cut` 以及 diagnostic `gb_core/p2_core`、`gb_strict_core/p2_strict_core` 应用到主分析；如果 prior 不存在，就先用 V2.1.0 streaming prior 生成它。

## 必须放进去的文件

| 文件 | 必须性 | LC-Pearl 从里面提取什么 | 例子 |
|---|---:|---|---|
| `lammps/*.in` 或 `lammps/*.lmp` | 必须 | 自动识别含 `pair_style gayberne` 和 `pair_coeff` 的原始 LAMMPS input，并复制/标准化为 `lammps/gb_param_source.in`。随后从标准文件读取 `units`、`atom_style`、`pair_style gayberne`、`pair_modify`、type 1-1 的 `pair_coeff`、`special_bonds`。 | 你以前命令里的 `--gb-param-file "../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp"`，现在直接复制到 `lammps/` 里即可，不必手动改名。 |
| `lammps/lammps_executable.txt` | 必须，除非已在 config 写绝对路径 | 第一条非注释行作为 LAMMPS 可执行文件路径，用来跑 run 0 验证。 | `/Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp` |
| 当前 dump 目录中的代表帧 dump | 必须 | `id/type/xu yu zu/quaternion/shape/box`，用于 run 0 验证。默认自动选第一个匹配 `*.dump` 的文件；也可在 config 指定。 | `traj.force_clamp_aligned.6055000.dump` |

## 可选但推荐的文件

| 文件 | 用途 |
|---|---|
| `topology/*.data` 或 `topology/*.dat` | 原始 LAMMPS data 文件，包含 `Atoms` 和 `Bonds` section。preflight 会自动根据 Bonds 和 `special_bonds lj` 生成 `topology/local_pairs.tsv` 与 `topology/exclude_pairs.tsv`。 |
| `topology/local_pairs.tsv` | 自动生成的内部标准文件；两列 atom id，表示链局部/相邻支持 pair。你通常不需要手写。 |
| `topology/exclude_pairs.tsv` | 自动生成的内部标准文件；两列 atom id，表示不应参与 contact 判断的 pair。你通常不需要手写。 |
| `thresholds/global_thresholds.json` | V2.1.0 阈值先验。若存在且 fingerprint 匹配、GB potential 已验证，wrapper 会在主分析前自动应用其中的 `gb_off_strength`、`gb_on_strength`、`p2_cut`、`gb_core_strength`、`p2_core_cut`、`gb_strict_core_strength`、`p2_strict_core_cut`。若不存在会自动生成。 |
| `thresholds/threshold_recommendations.json` | 与 `global_thresholds.json` 同步保存的推荐记录，方便人工复查。 |
| `validation/verified_potential.json` | 自动生成。迁移到相似项目时，如果 GB 参数、拓扑排除、dump schema、代码 fingerprint 匹配，就可跳过验证。 |

## 对应旧命令

旧命令：

```bash
python lc_lammps_run0_validate.py \
  --dump-file traj.force_clamp_aligned.6055000.dump \
  --gb-param-file "../in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp"
```

现在改成：

```text
dump_output_folder/
  traj.force_clamp_aligned.6055000.dump
  lc_pearl_preflight/
    lc_pearl_config.toml
    lammps/
      in.single_chain_constant_force_aligned_box-F0.38-tail_v0_no_thermal-xuyu-zu.lmp
                                  # 复制旧 --gb-param-file 指向的 in/lmp 文件；会自动标准化为 gb_param_source.in
      lammps_executable.txt       # 写 LAMMPS 可执行文件路径
    topology/
      topology_for_lc_analysis.data # 可选；放原始 data 文件后自动生成 local/exclude TSV
```

如果你希望固定代表帧，在 `lc_pearl_config.toml` 写：

```toml
[validation]
representative_dump = "traj.force_clamp_aligned.6055000.dump"
```

如果保持 `representative_dump = "auto"`，LC-Pearl 会自动选当前 dump 文件夹中第一个匹配 `[paths].pattern` 的 dump。
