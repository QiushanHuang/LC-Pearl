#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Sequence


SCHEMA_VERSION = 1


def sha256_file(path: Optional[Path]) -> Optional[str]:
    if path is None or not Path(path).exists():
        return None
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> Dict[str, object]:
    resolved = Path(path).resolve()
    return {
        "path": str(resolved),
        "exists": resolved.exists(),
        "sha256": sha256_file(resolved),
        "bytes": int(resolved.stat().st_size) if resolved.exists() else None,
    }


def write_sampling_manifest(
    output_path: Path,
    *,
    inputs: Sequence[Path],
    gb_param_file: Optional[Path],
    output_root: Path,
    frame_strategy: str,
    frame_limit: int,
    pair_strategy: str,
    sample_pairs_per_frame: int,
    seed: int,
    thresholds: Dict[str, object],
    topology_files: Dict[str, Optional[Path | str]],
    every: int = 1,
    min_pairs_per_stratum: Optional[int] = None,
    notes: Optional[Sequence[str]] = None,
) -> Dict[str, object]:
    topology_records: Dict[str, object] = {}
    for key, value in topology_files.items():
        path = Path(str(value)).resolve() if value else None
        topology_records[key] = file_record(path) if path is not None else None
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "output_root": str(Path(output_root).resolve()),
        "inputs": [file_record(path) for path in inputs],
        "gb_param_file": file_record(Path(gb_param_file)) if gb_param_file is not None else None,
        "topology_files": topology_records,
        "sampling": {
            "frame_strategy": str(frame_strategy),
            "frame_limit": int(frame_limit),
            "every": int(every),
            "pair_strategy": str(pair_strategy),
            "sample_pairs_per_frame": int(sample_pairs_per_frame),
            "min_pairs_per_stratum": int(min_pairs_per_stratum) if min_pairs_per_stratum is not None else None,
        },
        "thresholds": dict(thresholds),
        "random_seed": int(seed),
        "notes": list(notes or []),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Write an auditable sampling manifest for LC aggregation calibration runs.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--gb-param-file", type=Path, default=None)
    parser.add_argument("--frame-strategy", default="stratified")
    parser.add_argument("--frame-limit", type=int, default=0)
    parser.add_argument("--every", type=int, default=1)
    parser.add_argument("--pair-strategy", default="stratified")
    parser.add_argument("--sample-pairs-per-frame", type=int, default=0)
    parser.add_argument("--min-pairs-per-stratum", type=int, default=None)
    parser.add_argument("--seed", type=int, default=20260429)
    parser.add_argument("--gb-off-strength", type=float, default=0.12)
    parser.add_argument("--gb-on-strength", type=float, default=0.30)
    parser.add_argument("--p2-cut", type=float, default=0.70)
    parser.add_argument("--local-pair-file", type=Path, default=None)
    parser.add_argument("--exclude-pair-file", type=Path, default=None)
    parser.add_argument("inputs", nargs="*", type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    manifest = write_sampling_manifest(
        args.output,
        inputs=args.inputs,
        gb_param_file=args.gb_param_file,
        output_root=args.output_root,
        frame_strategy=args.frame_strategy,
        frame_limit=args.frame_limit,
        every=args.every,
        pair_strategy=args.pair_strategy,
        sample_pairs_per_frame=args.sample_pairs_per_frame,
        min_pairs_per_stratum=args.min_pairs_per_stratum,
        seed=args.seed,
        thresholds={
            "gb_off_strength": float(args.gb_off_strength),
            "gb_on_strength": float(args.gb_on_strength),
            "p2_cut": float(args.p2_cut),
        },
        topology_files={
            "local_pair_file": args.local_pair_file,
            "exclude_pair_file": args.exclude_pair_file,
        },
    )
    print(json.dumps({"manifest": str(args.output), "schema_version": manifest["schema_version"]}, indent=2))


if __name__ == "__main__":
    main()
