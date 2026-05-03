#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence


DEFAULT_MAX_WORKERS = 10


@dataclass(frozen=True)
class DomainSizeFrameCount:
    source_file: str
    timestep: int
    size: int
    total_count: int
    weak_count: int
    robust_count: int
    other_count: int


def parse_workers(raw_value: str | int, task_count: int) -> int:
    if isinstance(raw_value, int):
        requested = raw_value
    else:
        raw = str(raw_value).strip().lower()
        if raw == "auto":
            max_workers = max(1, int(os.environ.get("LC_PEARL_MAX_WORKERS", DEFAULT_MAX_WORKERS)))
            return max(1, min(os.cpu_count() or 1, max_workers, max(1, task_count)))
        requested = int(raw)
    if requested <= 0:
        raise ValueError("--workers must be a positive integer or 'auto'")
    return max(1, min(requested, max(1, task_count)))


def resolve_domain_table(input_path: Path) -> Path:
    """Resolve a diagnostics dir, output root, or direct TSV to domain_diagnostics.tsv."""
    path = Path(input_path).expanduser().resolve()
    if path.is_file():
        return path
    if (path / "domain_diagnostics.tsv").exists():
        return path / "domain_diagnostics.tsv"
    if (path / "diagnostics" / "domain_diagnostics.tsv").exists():
        return path / "diagnostics" / "domain_diagnostics.tsv"
    raise FileNotFoundError(
        f"cannot find domain_diagnostics.tsv from {path}; pass the TSV, diagnostics dir, or LC-Pearl output root"
    )


def default_output_dir(input_path: Path, domain_table: Path) -> Path:
    path = Path(input_path).expanduser().resolve()
    if path.is_file():
        return path.parent
    if (path / "domain_diagnostics.tsv").exists():
        return path
    if (path / "diagnostics" / "domain_diagnostics.tsv").exists():
        return path / "diagnostics"
    return domain_table.parent


def _safe_int(raw_value: object, default: int = 0) -> int:
    try:
        return int(float(str(raw_value).strip()))
    except (TypeError, ValueError):
        return default


def _record_to_key(record: dict[str, object]) -> Optional[tuple[str, int, int]]:
    size = _safe_int(record.get("size"), default=0)
    if size <= 0:
        return None
    source_file = str(record.get("source_file", "")).strip()
    timestep = _safe_int(record.get("timestep"), default=0)
    return source_file, timestep, size


def _counts_to_rows(counts: dict[tuple[str, int, int], dict[str, int]]) -> List[DomainSizeFrameCount]:
    rows: List[DomainSizeFrameCount] = []
    for (source_file, timestep, size), values in sorted(counts.items(), key=lambda item: (item[0][1], item[0][0], item[0][2])):
        rows.append(
            DomainSizeFrameCount(
                source_file=source_file,
                timestep=int(timestep),
                size=int(size),
                total_count=int(values["total"]),
                weak_count=int(values["weak"]),
                robust_count=int(values["robust"]),
                other_count=int(values["other"]),
            )
        )
    return rows


def collect_domain_size_frame_counts(domain_table: Path) -> List[DomainSizeFrameCount]:
    counts: dict[tuple[str, int, int], dict[str, int]] = {}
    with Path(domain_table).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        required = {"source_file", "timestep", "size", "classification"}
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"{domain_table} is missing required column(s): {', '.join(missing)}")

        for row in reader:
            key = _record_to_key(row)
            if key is None:
                continue
            bucket = counts.setdefault(key, {"total": 0, "weak": 0, "robust": 0, "other": 0})
            classification = str(row.get("classification", "")).strip().lower()
            bucket["total"] += 1
            if classification == "weak":
                bucket["weak"] += 1
            elif classification == "robust":
                bucket["robust"] += 1
            else:
                bucket["other"] += 1
    return _counts_to_rows(counts)


def collect_domain_size_frame_counts_from_records(records: Iterable[dict[str, object]]) -> List[DomainSizeFrameCount]:
    counts: dict[tuple[str, int, int], dict[str, int]] = {}
    for record in records:
        key = _record_to_key(record)
        if key is None:
            continue
        bucket = counts.setdefault(key, {"total": 0, "weak": 0, "robust": 0, "other": 0})
        classification = str(record.get("classification", "")).strip().lower()
        bucket["total"] += 1
        if classification == "weak":
            bucket["weak"] += 1
        elif classification == "robust":
            bucket["robust"] += 1
        else:
            bucket["other"] += 1
    return _counts_to_rows(counts)


def write_domain_size_frame_count_table(rows: Sequence[DomainSizeFrameCount], output_path: Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["source_file", "timestep", "size", "total_count", "weak_count", "robust_count", "other_count"])
        for row in rows:
            writer.writerow(
                [
                    row.source_file,
                    row.timestep,
                    row.size,
                    row.total_count,
                    row.weak_count,
                    row.robust_count,
                    row.other_count,
                ]
            )
    return output_path


def plot_domain_size_frame_counts(
    rows: Sequence[DomainSizeFrameCount],
    output_path: Path,
    *,
    size_floor: Optional[int] = None,
    title: str = "Domain size vs per-frame domain count",
    log_y: bool = False,
) -> Path:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "domain_size_vs_domain_count.png requires numpy and matplotlib. "
            "Use the LC-Pearl environment, for example: "
            "/Users/joshua/Desktop/MD/venv/bin/python3 scripts/lc_domain_size_counts.py <diagnostics_dir>"
        ) from exc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.4, 5.0), dpi=160)
    if rows:
        sizes = np.array([row.size for row in rows], dtype=float)
        total = np.array([row.total_count for row in rows], dtype=float)
        weak = np.array([row.weak_count for row in rows], dtype=float)
        robust = np.array([row.robust_count for row in rows], dtype=float)
        other = np.array([row.other_count for row in rows], dtype=float)

        ax.scatter(sizes, total, s=22, alpha=0.34, c="#315f9f", label="all domains/frame")
        robust_mask = robust > 0
        if np.any(robust_mask):
            ax.scatter(sizes[robust_mask], robust[robust_mask], s=24, alpha=0.58, c="#9b2c2c", label="robust domains/frame")
        weak_mask = weak > 0
        if np.any(weak_mask):
            ax.scatter(sizes[weak_mask], weak[weak_mask], s=20, alpha=0.38, c="#66716b", label="weak domains/frame")
        other_mask = other > 0
        if np.any(other_mask):
            ax.scatter(sizes[other_mask], other[other_mask], s=18, alpha=0.40, c="#7c3aed", label="other domains/frame")
        if size_floor is not None:
            ax.axvline(float(size_floor), color="#315f9f", lw=1.5, label=f"size_floor={size_floor}")
        ax.set_xlim(float(np.min(sizes)) - 0.8, float(np.max(sizes)) + 0.8)
        if log_y:
            ax.set_yscale("log")
    else:
        ax.text(0.5, 0.5, "No domain rows", ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("domain size")
    ax.set_ylabel("domain count per frame")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def write_domain_size_frame_count_outputs(
    domain_table: Path,
    output_dir: Path,
    *,
    size_floor: Optional[int] = None,
    log_y: bool = False,
) -> dict[str, object]:
    rows = collect_domain_size_frame_counts(Path(domain_table))
    output_dir = Path(output_dir)
    table_path = write_domain_size_frame_count_table(rows, output_dir / "domain_size_frame_counts.tsv")
    plot_path = plot_domain_size_frame_counts(
        rows,
        output_dir / "domain_size_vs_domain_count.png",
        size_floor=size_floor,
        log_y=log_y,
    )
    return {
        "domain_table": str(Path(domain_table)),
        "frame_count_table": str(table_path),
        "plot": str(plot_path),
        "n_points": len(rows),
        "total_domains": int(sum(row.total_count for row in rows)),
    }


def process_one(input_path: Path, output_dir: Optional[Path], size_floor: Optional[int], log_y: bool) -> dict[str, object]:
    domain_table = resolve_domain_table(input_path)
    resolved_output_dir = Path(output_dir).expanduser().resolve() if output_dir is not None else default_output_dir(input_path, domain_table)
    return write_domain_size_frame_count_outputs(domain_table, resolved_output_dir, size_floor=size_floor, log_y=log_y)


def process_many(
    input_paths: Sequence[Path],
    *,
    output_dir: Optional[Path] = None,
    size_floor: Optional[int] = None,
    log_y: bool = False,
    workers: str | int = "auto",
) -> List[dict[str, object]]:
    if output_dir is not None and len(input_paths) != 1:
        raise ValueError("--output-dir can only be used with one input path")
    worker_count = parse_workers(workers, len(input_paths))
    jobs = [(Path(path), output_dir, size_floor, log_y) for path in input_paths]
    if worker_count == 1:
        return [process_one(*job) for job in jobs]
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        return list(executor.map(_process_one_star, jobs))


def _process_one_star(args: tuple[Path, Optional[Path], Optional[int], bool]) -> dict[str, object]:
    return process_one(*args)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate domain_size_vs_domain_count scatter diagnostics from existing LC-Pearl domain_diagnostics.tsv files."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="domain_diagnostics.tsv, diagnostics directory, or LC-Pearl output root.",
    )
    parser.add_argument("--output-dir", type=Path, help="Output directory for a single input. Defaults to the diagnostics directory.")
    parser.add_argument("--size-floor", type=int, help="Optional vertical size threshold line, usually max(n_min, robust_min_size).")
    parser.add_argument("--log-y", action="store_true", help="Use a log-scaled per-frame count axis.")
    parser.add_argument("--workers", default="auto", help="Parallel workers for multiple input paths. Default: auto.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    outputs = process_many(
        args.inputs,
        output_dir=args.output_dir,
        size_floor=args.size_floor,
        log_y=bool(args.log_y),
        workers=args.workers,
    )
    for item in outputs:
        print(
            f"[OK] domain_size_vs_domain_count: domains={item['total_domains']} points={item['n_points']} "
            f"plot={item['plot']}"
        )


if __name__ == "__main__":
    main()
