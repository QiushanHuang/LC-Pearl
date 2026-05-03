# LC-Pearl v2.1 Operation Manual

This manual focuses on running, rerunning, and debugging the pipeline.

## 1. Required Inputs

In a dump output folder, LC-Pearl expects or creates:

```text
lc_pearl_preflight/
  lammps/
  topology/
  validation/
  thresholds/
```

You provide:

| Input | Put it here | Purpose |
|---|---|---|
| LAMMPS input deck | `lc_pearl_preflight/lammps/*.in` or `*.lmp` | Extract Gay-Berne parameters. |
| LAMMPS executable path | `lc_pearl_preflight/lammps/lammps_executable.txt` | Run validation when needed. |
| LAMMPS data/topology file | `lc_pearl_preflight/topology/*.data` or `*.dat` | Convert bonds into local/excluded pair tables. |

LC-Pearl generates:

| Generated file | Purpose |
|---|---|
| `lc_pearl_preflight/lammps/gb_param_source.in` | Normalized copied input used by the pipeline. |
| `lc_pearl_preflight/topology/local_pairs.tsv` | Local/adjacent mesogen pair table. |
| `lc_pearl_preflight/topology/exclude_pairs.tsv` | Excluded/special pair table when topology is available. |
| `lc_pearl_preflight/validation/verified_potential.json` | Potential-validation cache. |
| `lc_pearl_preflight/thresholds/global_thresholds.json` | Reusable 2D threshold-prior cache. |

## 2. Commands

Create/check preflight:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py
```

Full run:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py auto
```

Validation only:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py validate
```

Run using existing caches:

```bash
python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_here.py run
```

Preview a TOML-driven run:

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 /Users/joshua/Desktop/MD/LC-Pearl/lc_pearl_cli.py print-run --config /Users/joshua/Desktop/MD/LC-Pearl/configs/quick_run.toml
```

## 3. Sampling And Parallelism

Important configuration fields:

```toml
[threshold_prior]
global_frame_stride = 1        # 1 = all frames; 10 or 100 = faster preview
global_frame_budget = 0        # 0 = no budget cap
workers = "auto"               # capped by CPU, chunk count, and LC_PEARL_MAX_WORKERS

[analysis]
workers = "auto"
edge_diagnostics_table = "off" # off, sample, full
write_frame_jsonl = false
accepted_edge_audit = false
```

Use full frame sampling for final threshold priors. Use stride 10 or 100 for fast checks. If memory pressure or SIGKILL occurs, increase stride, reduce chunk size, or reuse an already validated representative `global_thresholds.json`.

## 4. Output Root

Default v2 output:

```text
lc_domain_pearl_v2_output/
```

Important subfolders:

| Folder | Meaning |
|---|---|
| `per_file/` | Per-dump labels, summaries, optional edges/segments. |
| `diagnostics/` | Domain/pearl/threshold/mechanics diagnostics. |
| `gb_pair_audit/` | Optional legacy/targeted pair audit, not the primary v2 threshold engine. |

The primary threshold-prior evidence is in:

```text
lc_pearl_preflight/thresholds/
```

not in `gb_pair_audit/`.

## 5. Standalone Domain-Size Count Diagnostic

This diagnostic does not reparse dump files. It reads `diagnostics/domain_diagnostics.tsv`.

```bash
/Users/joshua/Desktop/MD/venv/bin/python3 \
  /Users/joshua/Desktop/MD/LC-Pearl/scripts/lc_domain_size_counts.py \
  /path/to/lc_domain_pearl_v2_output
```

Outputs:

```text
diagnostics/domain_size_frame_counts.tsv
diagnostics/domain_size_vs_domain_count.png
```

## 6. Debugging Checklist

If labels look too broad:

1. Check `global_thresholds.json` and confirm the applied values are not old defaults.
2. Inspect `gb_strength_vs_p2_stream_hexbin.png` and lobe split preview.
3. Color by `lc_aggregation_tier` in OVITO, not only `lc_domain`.
4. Compare `lc_cluster` with `lc_domain`; `lc_cluster` is simple spatial clustering and is intentionally not the final robust-domain definition.
5. If local contacts dominate, verify topology files and `local_pairs.tsv`.
6. If necessary, rerun a small frame subset with `edge_diagnostics_table = "sample"` or `"full"`.

If the run is slow:

1. Confirm whether `threshold_prior.global_frame_stride = 1` is intentionally full sampling.
2. Reuse a representative `global_thresholds.json` instead of rebuilding it for every similar folder.
3. Keep full edge tables and JSONL debug output disabled.
4. Use `LC_PEARL_MAX_WORKERS=10` or a lower value appropriate for the machine.

## 中文操作要点

正式分析建议先用一个覆盖高温到低温、或覆盖主要构型变化的代表轨迹生成 `global_thresholds.json`，然后把整个 `lc_pearl_preflight/` 迁移到相同 potential 的其他力或温度目录。这样后续目录可以直接复用验证和阈值先验，避免重复做全量 2D prior。
