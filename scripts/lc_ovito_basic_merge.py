#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple


Frame = Tuple[int, List[str], List[Tuple[float, float]], List[str], List[List[str]]]


def read_atom_dump(path: Path) -> Iterator[Frame]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    idx = 0
    while idx < len(lines):
        if not lines[idx].startswith("ITEM: TIMESTEP"):
            idx += 1
            continue
        timestep = int(lines[idx + 1].strip())
        n_atoms = int(lines[idx + 3].strip())
        box_header = lines[idx + 4].split()[3:]
        bounds = [
            tuple(float(v) for v in lines[idx + 5].split()[:2]),
            tuple(float(v) for v in lines[idx + 6].split()[:2]),
            tuple(float(v) for v in lines[idx + 7].split()[:2]),
        ]
        columns = lines[idx + 8].split()[2:]
        rows = [lines[idx + 9 + row_idx].split() for row_idx in range(n_atoms)]
        yield timestep, box_header, bounds, columns, rows
        idx += 9 + n_atoms


def read_contact_edges(path: Path) -> Dict[int, List[Dict[str, str]]]:
    frames: Dict[int, List[Dict[str, str]]] = {}
    if not path.exists():
        return frames
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    idx = 0
    while idx < len(lines):
        if not lines[idx].startswith("ITEM: TIMESTEP"):
            idx += 1
            continue
        timestep = int(lines[idx + 1].strip())
        n_entries = int(lines[idx + 3].strip())
        columns = lines[idx + 8].split()[2:]
        rows = []
        for row_idx in range(n_entries):
            parts = lines[idx + 9 + row_idx].split()
            rows.append({column: parts[col_idx] for col_idx, column in enumerate(columns) if col_idx < len(parts)})
        frames[timestep] = rows
        idx += 9 + n_entries
    return frames


def set_value(row: List[str], col_index: Dict[str, int], names: Iterable[str], value: object) -> None:
    for name in names:
        if name in col_index:
            row[col_index[name]] = str(value)


def get_float(row: List[str], col_index: Dict[str, int], names: Iterable[str], default: float = 0.0) -> float:
    for name in names:
        if name in col_index:
            try:
                return float(row[col_index[name]])
            except ValueError:
                return default
    return default


def get_int(row: List[str], col_index: Dict[str, int], name: str, default: int = 0) -> int:
    if name not in col_index:
        return default
    try:
        return int(float(row[col_index[name]]))
    except ValueError:
        return default


def merged_columns(label_columns: List[str]) -> List[str]:
    columns = list(label_columns)
    for name in ("vis_kind", "vis_edge_type", "source_atom_id"):
        if name not in columns:
            columns.append(name)
    return columns


def pad_row(row: List[str], old_columns: List[str], new_columns: List[str]) -> List[str]:
    mapped = {column: row[idx] for idx, column in enumerate(old_columns) if idx < len(row)}
    return [mapped.get(column, "0") for column in new_columns]


def make_marker_from_atom(
    source_row: List[str],
    columns: List[str],
    *,
    new_id: int,
    marker_type: int,
    x: float,
    y: float,
    z: float,
    radius: float,
    vis_kind: int,
    edge_type: int = 0,
    source_atom_id: int = 0,
) -> List[str]:
    col_index = {name: idx for idx, name in enumerate(columns)}
    row = list(source_row)
    set_value(row, col_index, ("id",), new_id)
    set_value(row, col_index, ("type",), marker_type)
    set_value(row, col_index, ("x", "xu"), f"{x:.10g}")
    set_value(row, col_index, ("y", "yu"), f"{y:.10g}")
    set_value(row, col_index, ("z", "zu"), f"{z:.10g}")
    set_value(row, col_index, ("quatw",), "1")
    set_value(row, col_index, ("quati", "quatj", "quatk"), "0")
    set_value(row, col_index, ("shapex", "shapey", "shapez"), f"{radius:.10g}")
    set_value(row, col_index, ("mass",), "1")
    set_value(row, col_index, ("vis_kind",), vis_kind)
    set_value(row, col_index, ("vis_edge_type",), edge_type)
    set_value(row, col_index, ("source_atom_id",), source_atom_id)
    return row


def write_basic_vis_dump(
    label_dump: Path,
    output_dump: Path,
    *,
    contact_edges: Path | None = None,
    halo_radius: float = 1.9,
    contact_radius: float = 0.25,
    contact_samples: int = 5,
) -> None:
    contact_by_timestep = read_contact_edges(contact_edges) if contact_edges else {}
    with output_dump.open("w", encoding="utf-8") as handle:
        for timestep, box_header, bounds, label_columns, label_rows in read_atom_dump(label_dump):
            out_columns = merged_columns(label_columns)
            col_index = {name: idx for idx, name in enumerate(out_columns)}
            padded_rows = [pad_row(row, label_columns, out_columns) for row in label_rows]
            atom_by_id = {get_int(row, col_index, "id"): row for row in padded_rows}
            output_rows: List[List[str]] = []
            max_id = max(atom_by_id, default=0)
            for row in padded_rows:
                set_value(row, col_index, ("vis_kind",), 0)
                set_value(row, col_index, ("vis_edge_type",), 0)
                set_value(row, col_index, ("source_atom_id",), get_int(row, col_index, "id"))
                output_rows.append(row)

            for row in padded_rows:
                if get_int(row, col_index, "lc_cluster") <= 0:
                    continue
                max_id += 1
                atom_id = get_int(row, col_index, "id")
                output_rows.append(
                    make_marker_from_atom(
                        row,
                        out_columns,
                        new_id=max_id,
                        marker_type=99,
                        x=get_float(row, col_index, ("x", "xu")),
                        y=get_float(row, col_index, ("y", "yu")),
                        z=get_float(row, col_index, ("z", "zu")),
                        radius=halo_radius,
                        vis_kind=1,
                        source_atom_id=atom_id,
                    )
                )

            for edge in contact_by_timestep.get(timestep, []):
                try:
                    atom_i = int(float(edge.get("atom_i", "0")))
                    atom_j = int(float(edge.get("atom_j", "0")))
                    edge_type = int(float(edge.get("edge_type_code", "0")))
                except ValueError:
                    continue
                left = atom_by_id.get(atom_i)
                right = atom_by_id.get(atom_j)
                if left is None or right is None:
                    continue
                lx = get_float(left, col_index, ("x", "xu"))
                ly = get_float(left, col_index, ("y", "yu"))
                lz = get_float(left, col_index, ("z", "zu"))
                rx = get_float(right, col_index, ("x", "xu"))
                ry = get_float(right, col_index, ("y", "yu"))
                rz = get_float(right, col_index, ("z", "zu"))
                for sample in range(1, max(2, contact_samples + 1)):
                    frac = sample / float(contact_samples + 1)
                    max_id += 1
                    output_rows.append(
                        make_marker_from_atom(
                            left,
                            out_columns,
                            new_id=max_id,
                            marker_type=98,
                            x=lx + frac * (rx - lx),
                            y=ly + frac * (ry - ly),
                            z=lz + frac * (rz - lz),
                            radius=contact_radius,
                            vis_kind=2,
                            edge_type=edge_type,
                            source_atom_id=atom_i,
                        )
                    )

            handle.write("ITEM: TIMESTEP\n")
            handle.write(f"{timestep}\n")
            handle.write("ITEM: NUMBER OF ATOMS\n")
            handle.write(f"{len(output_rows)}\n")
            handle.write("ITEM: BOX BOUNDS " + " ".join(box_header) + "\n")
            for lo, hi in bounds:
                handle.write(f"{lo:.12g} {hi:.12g}\n")
            handle.write("ITEM: ATOMS " + " ".join(out_columns) + "\n")
            for row in output_rows:
                handle.write(" ".join(row) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge LC labels, halo markers, and contact markers into one OVITO Basic-compatible dump.")
    parser.add_argument("--label-dump", type=Path, required=True)
    parser.add_argument("--contact-edges", type=Path, default=None)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--halo-radius", type=float, default=1.9)
    parser.add_argument("--contact-radius", type=float, default=0.25)
    parser.add_argument("--contact-samples", type=int, default=5)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    write_basic_vis_dump(
        args.label_dump,
        args.output,
        contact_edges=args.contact_edges,
        halo_radius=float(args.halo_radius),
        contact_radius=float(args.contact_radius),
        contact_samples=int(args.contact_samples),
    )
    print(args.output)


if __name__ == "__main__":
    main()
