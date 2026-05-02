#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


def normalized_pair(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a <= b else (b, a)


def parse_lammps_data_bonds(data_path: Path) -> List[Tuple[int, int, int, int]]:
    bonds: List[Tuple[int, int, int, int]] = []
    in_bonds = False
    with data_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            head = line.split()[0]
            if head == "Bonds":
                in_bonds = True
                continue
            if in_bonds and head.isalpha():
                break
            if not in_bonds:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                bond_id = int(parts[0])
                bond_type = int(parts[1])
                atom_i = int(parts[2])
                atom_j = int(parts[3])
            except ValueError:
                continue
            bonds.append((bond_id, bond_type, atom_i, atom_j))
    return bonds


def special_pairs_from_bonds(
    bonds: Sequence[Tuple[int, int, int, int]],
    lj_factors: Tuple[float, float, float],
) -> Dict[str, set[Tuple[int, int]]]:
    adjacency: Dict[int, set[int]] = {}
    for _bond_id, _bond_type, atom_i, atom_j in bonds:
        adjacency.setdefault(atom_i, set()).add(atom_j)
        adjacency.setdefault(atom_j, set()).add(atom_i)

    by_distance: Dict[int, set[Tuple[int, int]]] = {1: set(), 2: set(), 3: set()}
    for start in adjacency:
        queue: deque[Tuple[int, int]] = deque([(start, 0)])
        visited = {start}
        while queue:
            atom, distance = queue.popleft()
            if distance >= 3:
                continue
            for neighbor in adjacency.get(atom, set()):
                if neighbor in visited:
                    continue
                visited.add(neighbor)
                next_distance = distance + 1
                if 1 <= next_distance <= 3 and start != neighbor:
                    by_distance[next_distance].add(normalized_pair(start, neighbor))
                queue.append((neighbor, next_distance))

    excluded: set[Tuple[int, int]] = set()
    local: set[Tuple[int, int]] = set()
    for distance, factor in enumerate(lj_factors, start=1):
        pairs = by_distance[distance]
        if factor == 0.0:
            excluded.update(pairs)
        else:
            local.update(pairs)
    return {
        "one_two": by_distance[1],
        "one_three": by_distance[2],
        "one_four": by_distance[3],
        "excluded": excluded,
        "local_nonzero_special": local,
    }


def write_pair_file(path: Path, pairs: Iterable[Tuple[int, int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("atom_i\tatom_j\n")
        for atom_i, atom_j in sorted(set(pairs)):
            handle.write(f"{atom_i}\t{atom_j}\n")


def write_restart_to_data_template(path: Path, restart_path: str = "Restart.cooldown.52000000") -> None:
    path.write_text(
        "\n".join(
            [
                "# Generated helper for topology extraction.",
                "# Run with LAMMPS in the directory that contains the restart file.",
                "# This does not run dynamics; it converts restart topology to a readable data file.",
                "clear",
                f"read_restart {restart_path}",
                "write_data topology_for_lc_analysis.data",
                "",
            ]
        ),
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate explicit pair exclusion/local maps from a LAMMPS data Bonds section.")
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("from-data", help="Read LAMMPS data file and write pair maps.")
    build.add_argument("data_file", type=Path)
    build.add_argument("--special-lj", default="0,1,1", help="special_bonds lj factors for 1-2,1-3,1-4.")
    build.add_argument("--output-root", type=Path, default=Path("lc_topology_pairs"))
    tmpl = sub.add_parser("write-restart-template", help="Write a LAMMPS input that converts a restart to a data file.")
    tmpl.add_argument("--restart", default="Restart.cooldown.52000000")
    tmpl.add_argument("--output", type=Path, default=Path("write_topology_data.in"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "write-restart-template":
        write_restart_to_data_template(args.output, restart_path=args.restart)
        print(args.output)
        return

    factors = tuple(float(item.strip()) for item in args.special_lj.split(",") if item.strip())
    if len(factors) != 3:
        raise SystemExit("--special-lj must contain three comma-separated values")
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    bonds = parse_lammps_data_bonds(args.data_file)
    pair_sets = special_pairs_from_bonds(bonds, factors)  # type: ignore[arg-type]
    write_pair_file(output / "exclude_pairs.tsv", pair_sets["excluded"])
    write_pair_file(output / "local_special_pairs.tsv", pair_sets["local_nonzero_special"])
    summary = {
        "data_file": str(args.data_file),
        "bond_count": len(bonds),
        "special_lj": list(factors),
        "one_two_pairs": len(pair_sets["one_two"]),
        "one_three_pairs": len(pair_sets["one_three"]),
        "one_four_pairs": len(pair_sets["one_four"]),
        "excluded_pairs": len(pair_sets["excluded"]),
        "local_nonzero_special_pairs": len(pair_sets["local_nonzero_special"]),
    }
    (output / "topology_pair_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
