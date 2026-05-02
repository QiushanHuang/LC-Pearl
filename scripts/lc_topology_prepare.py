#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from lc_topology_exclusions import (
    parse_lammps_data_bonds,
    special_pairs_from_bonds,
    write_pair_file,
    write_restart_to_data_template,
)


SCHEMA_VERSION = 1


@dataclass(frozen=True)
class AtomRecord:
    atom_id: int
    atom_type: int
    molecule_id: Optional[int]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def section_name(line: str) -> Optional[str]:
    stripped = line.strip()
    if not stripped:
        return None
    head = stripped.split("#", 1)[0].strip().split()
    if not head:
        return None
    if head[0].isalpha():
        return head[0]
    return None


def parse_lammps_data_atoms(data_path: Path) -> List[AtomRecord]:
    atoms: List[AtomRecord] = []
    in_atoms = False
    with data_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            name = section_name(line)
            if name == "Atoms":
                in_atoms = True
                continue
            if in_atoms and name is not None:
                break
            if not in_atoms or line.startswith("#"):
                continue
            parts = line.split("#", 1)[0].split()
            if len(parts) < 2:
                continue
            try:
                atom_id = int(parts[0])
                atom_type = int(parts[1])
            except ValueError:
                continue
            molecule_id: Optional[int] = None
            # LAMMPS write_data for atom_style hybrid ellipsoid bond commonly writes:
            # id type x y z ellipsoidflag density molecule-ID ix iy iz
            if len(parts) >= 8:
                try:
                    molecule_id = int(float(parts[7]))
                except ValueError:
                    molecule_id = None
            atoms.append(AtomRecord(atom_id=atom_id, atom_type=atom_type, molecule_id=molecule_id))
    return atoms


def parse_special_lj(value: str) -> Tuple[float, float, float]:
    factors = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    if len(factors) != 3:
        raise ValueError("--special-lj must contain exactly three comma-separated values")
    return factors  # type: ignore[return-value]


def infer_mesogen_chain_map(
    atoms: Sequence[AtomRecord],
    mesogen_type: int,
    anchor_type: int,
) -> List[Dict[str, object]]:
    mesogens = [atom for atom in atoms if atom.atom_type == int(mesogen_type)]
    if not mesogens:
        return []
    anchors_by_mol: Dict[int, List[int]] = {}
    for atom in atoms:
        if atom.atom_type == int(anchor_type) and atom.molecule_id is not None:
            anchors_by_mol.setdefault(int(atom.molecule_id), []).append(atom.atom_id)
    with_mol = [atom for atom in mesogens if atom.molecule_id is not None and atom.molecule_id > 0]
    if len(with_mol) == len(mesogens) and len({atom.molecule_id for atom in mesogens}) == len(mesogens):
        ordered = sorted(mesogens, key=lambda atom: (int(atom.molecule_id or 0), atom.atom_id))
        method = "rigid_molecule_id_sorted"
        confidence = "high"
    else:
        ordered = sorted(mesogens, key=lambda atom: atom.atom_id)
        method = "atom_id_sorted"
        confidence = "medium"
    rows: List[Dict[str, object]] = []
    for chain_index, atom in enumerate(ordered):
        anchors = sorted(anchors_by_mol.get(int(atom.molecule_id or -1), []))
        rows.append(
            {
                "mesogen_atom_id": atom.atom_id,
                "chain_id": 1,
                "chain_s": chain_index,
                "atom_type": atom.atom_type,
                "rigid_mol_id": atom.molecule_id if atom.molecule_id is not None else "",
                "left_anchor_id": anchors[0] if anchors else "",
                "right_anchor_id": anchors[-1] if len(anchors) >= 2 else "",
                "method": method,
                "confidence": confidence,
            }
        )
    return rows


def write_chain_index_map(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("mesogen_atom_id\tchain_id\tchain_s\tatom_type\trigid_mol_id\tleft_anchor_id\tright_anchor_id\tmethod\tconfidence\n")
        for row in rows:
            handle.write(
                f"{row['mesogen_atom_id']}\t{row['chain_id']}\t{row['chain_s']}\t"
                f"{row['atom_type']}\t{row['rigid_mol_id']}\t{row['left_anchor_id']}\t"
                f"{row['right_anchor_id']}\t{row['method']}\t{row['confidence']}\n"
            )


def prepare_topology(
    *,
    data_file: Path,
    output_root: Path,
    special_lj: Tuple[float, float, float],
    mesogen_type: int = 1,
    anchor_type: int = 3,
) -> Dict[str, object]:
    output = output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    bonds = parse_lammps_data_bonds(data_file)
    atoms = parse_lammps_data_atoms(data_file)
    pair_sets = special_pairs_from_bonds(bonds, special_lj)
    chain_rows = infer_mesogen_chain_map(atoms, mesogen_type=mesogen_type, anchor_type=anchor_type)
    exclude_path = output / "exclude_pairs.tsv"
    local_path = output / "local_special_pairs.tsv"
    chain_path = output / "chain_index_map.tsv"
    write_pair_file(exclude_path, pair_sets["excluded"])
    write_pair_file(local_path, pair_sets["local_nonzero_special"])
    write_chain_index_map(chain_path, chain_rows)

    warnings: List[str] = []
    if not atoms:
        warnings.append("No Atoms section was parsed; chain_index_map.tsv is empty.")
    if not bonds:
        warnings.append("No Bonds section was parsed; pair exclusion maps are empty.")
    if not chain_rows:
        warnings.append("No mesogen atoms were found for the requested mesogen type.")
    if chain_rows and any(row.get("method") == "rigid_molecule_id_sorted" for row in chain_rows):
        warnings.append("Rigid a-E-a molecule IDs were used only for sequence mapping; they are not LAMMPS special_bonds.")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_file": str(data_file),
        "data_file_sha256": sha256_file(data_file),
        "mesogen_type": int(mesogen_type),
        "anchor_type": int(anchor_type),
        "special_lj": list(special_lj),
        "outputs": {
            "exclude_pair_file": str(exclude_path),
            "local_pair_file": str(local_path),
            "chain_index_map": str(chain_path),
        },
        "counts": {
            "atom_count": len(atoms),
            "bond_count": len(bonds),
            "mesogen_count": len(chain_rows),
            "one_two_pairs": len(pair_sets["one_two"]),
            "one_three_pairs": len(pair_sets["one_three"]),
            "one_four_pairs": len(pair_sets["one_four"]),
            "excluded_pairs": len(pair_sets["excluded"]),
            "local_nonzero_special_pairs": len(pair_sets["local_nonzero_special"]),
        },
        "warnings": warnings,
    }
    (output / "topology_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare reusable LC topology files for aggregation analysis.")
    sub = parser.add_subparsers(dest="command", required=True)
    from_data = sub.add_parser("from-data", help="Read a LAMMPS data file and write topology helper files.")
    from_data.add_argument("data_file", type=Path)
    from_data.add_argument("--special-lj", default="0,1,1", help="special_bonds lj factors for 1-2,1-3,1-4.")
    from_data.add_argument("--mesogen-type", type=int, default=1)
    from_data.add_argument("--anchor-type", type=int, default=3)
    from_data.add_argument("--output-root", type=Path, default=Path("lc_topology_pairs"))
    tmpl = sub.add_parser("write-restart-template", help="Write a restart-to-data conversion template.")
    tmpl.add_argument("--restart", default="Restart.GB.rho030.0")
    tmpl.add_argument("--output", type=Path, default=Path("write_topology_data.in"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "write-restart-template":
        write_restart_to_data_template(args.output, restart_path=args.restart)
        print(args.output)
        return
    manifest = prepare_topology(
        data_file=args.data_file,
        output_root=args.output_root,
        special_lj=parse_special_lj(args.special_lj),
        mesogen_type=int(args.mesogen_type),
        anchor_type=int(args.anchor_type),
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
