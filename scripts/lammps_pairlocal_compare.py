#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def normalized_pair(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a <= b else (b, a)


def read_table(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        first = handle.readline()
        delimiter = "\t" if "\t" in first else None
        handle.seek(0)
        return list(csv.DictReader(handle, delimiter=delimiter))


def find_column(row: Dict[str, str], candidates: Tuple[str, ...]) -> Optional[str]:
    lowered = {key.lower(): key for key in row}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def pair_key_from_row(row: Dict[str, str]) -> Tuple[int, int, int]:
    timestep_col = find_column(row, ("timestep", "step", "TimeStep"))
    atom_i_col = find_column(row, ("atom_i", "id1", "i", "atom1", "c_pid[1]"))
    atom_j_col = find_column(row, ("atom_j", "id2", "j", "atom2", "c_pid[2]"))
    if atom_i_col is None or atom_j_col is None:
        raise ValueError("pair table must contain atom_i/atom_j or id1/id2 columns")
    timestep = int(float(row[timestep_col])) if timestep_col is not None and row.get(timestep_col, "") else 0
    atom_i = int(float(row[atom_i_col]))
    atom_j = int(float(row[atom_j_col]))
    left, right = normalized_pair(atom_i, atom_j)
    return timestep, left, right


def energy_from_row(row: Dict[str, str]) -> float:
    energy_col = find_column(row, ("pair_energy", "eng", "energy", "c_pair[1]", "c_plocal[1]"))
    if energy_col is None:
        raise ValueError("pair table must contain pair_energy/eng/energy column")
    return float(row[energy_col])


def compare_pair_tables(python_table: Path, lammps_table: Path, output_root: Path) -> Dict[str, object]:
    output_root.mkdir(parents=True, exist_ok=True)
    python_rows = read_table(python_table)
    lammps_rows = read_table(lammps_table)
    python_map: Dict[Tuple[int, int, int], float] = {
        pair_key_from_row(row): energy_from_row(row)
        for row in python_rows
    }
    lammps_map: Dict[Tuple[int, int, int], float] = {
        pair_key_from_row(row): energy_from_row(row)
        for row in lammps_rows
    }
    shared = sorted(set(python_map) & set(lammps_map))
    diffs: List[float] = []
    with (output_root / "pairlocal_comparison.tsv").open("w", encoding="utf-8") as handle:
        handle.write("timestep\tatom_i\tatom_j\tpython_pair_energy\tlammps_pair_energy\tdelta_python_minus_lammps\tabs_delta\n")
        for timestep, atom_i, atom_j in shared:
            py_energy = python_map[(timestep, atom_i, atom_j)]
            lm_energy = lammps_map[(timestep, atom_i, atom_j)]
            delta = py_energy - lm_energy
            diffs.append(delta)
            handle.write(
                f"{timestep}\t{atom_i}\t{atom_j}\t{py_energy:.12g}\t{lm_energy:.12g}\t"
                f"{delta:.12g}\t{abs(delta):.12g}\n"
            )
    if diffs:
        rmse = math.sqrt(sum(delta * delta for delta in diffs) / len(diffs))
        max_abs = max(abs(delta) for delta in diffs)
        mean_abs = sum(abs(delta) for delta in diffs) / len(diffs)
    else:
        rmse = max_abs = mean_abs = math.nan
    summary = {
        "python_rows": len(python_rows),
        "lammps_rows": len(lammps_rows),
        "matched_pairs": len(shared),
        "missing_in_lammps": len(set(python_map) - set(lammps_map)),
        "missing_in_python": len(set(lammps_map) - set(python_map)),
        "rmse": rmse,
        "mean_abs_delta": mean_abs,
        "max_abs_delta": max_abs,
        "interpretation": (
            "Small deltas support the Python GB reconstruction. Large deltas usually mean pair_modify shift, "
            "special_bonds exclusions, quaternion conventions, dump precision, or pair/local table formatting differ."
        ),
    }
    (output_root / "pairlocal_comparison_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def write_pairlocal_template(output: Path, dump_file: str, topology_file: str = "topology_for_lc_analysis.data") -> None:
    output.write_text(
        "\n".join(
            [
                "# LAMMPS pair/local validation template for LC GB audit.",
                "# This is intentionally a template: verify atom_style, pair_style, pair_coeff, and rerun fields",
                "# against the simulation input before using the resulting table as ground truth.",
                "clear",
                f"read_data {topology_file}",
                "# Re-add the exact pair_style, pair_modify, and pair_coeff lines from the production input here.",
                "# Example:",
                "# pair_style hybrid gayberne 1.0 3.0 1.0 5.0 lj/cut 1.122462048309373",
                "# pair_modify shift yes",
                "# pair_coeff 1 1 gayberne 1.0 1.0 1.0 1.0 0.2 1.0 1.0 0.2 5.0",
                "# pair_gayberne does not support compute pair/local eng in common LAMMPS builds.",
                "# Use pair/local only for the neighbor pair list/geometry, and compare Python energy sums",
                "# against thermo pe or compute pe/atom pair reduced over atoms.",
                "compute pid all property/local patom1 patom2",
                "compute geom all pair/local dist dx dy dz",
                "compute peatom all pe/atom pair",
                "compute pesum all reduce sum c_peatom",
                "thermo_style custom step pe c_pesum",
                "dump plocal all local 1 lammps_pairlocal_geometry.tsv c_pid[1] c_pid[2] c_geom[1] c_geom[2] c_geom[3] c_geom[4]",
                f"rerun {dump_file} dump x y z box yes replace yes",
                "undump plocal",
                "",
                "# If this LAMMPS build cannot rerun quaternion/aspherical orientation from the dump, do not use",
                "# rerun for final GB energy validation. Generate per-frame atom_style ellipsoid data from the dump",
                "# and run 0 instead, then compare Python totals against LAMMPS pe/c_pesum.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare Python GB audit pair energies against a LAMMPS pair/local table.")
    sub = parser.add_subparsers(dest="command", required=True)
    cmp_parser = sub.add_parser("compare")
    cmp_parser.add_argument("--python-table", type=Path, required=True)
    cmp_parser.add_argument("--lammps-table", type=Path, required=True)
    cmp_parser.add_argument("--output-root", type=Path, default=Path("pairlocal_compare_output"))
    tmpl = sub.add_parser("write-template")
    tmpl.add_argument("--dump-file", required=True)
    tmpl.add_argument("--topology-file", default="topology_for_lc_analysis.data")
    tmpl.add_argument("--output", type=Path, default=Path("rerun_pairlocal_template.in"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "write-template":
        write_pairlocal_template(args.output, dump_file=args.dump_file, topology_file=args.topology_file)
        print(args.output)
        return
    summary = compare_pair_tables(args.python_table, args.lammps_table, args.output_root)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
