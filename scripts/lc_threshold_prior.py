#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import gb_pair_audit
import lc_calibrate_thresholds as cal
import liquid_crystal_aggregation as lca
from lc_threshold_recommend import annotate_hline, annotate_vline, threshold_label


SCHEMA_VERSION = 5
METHOD_NAME = "LC Domain-Pearl V2 2D lobe streaming threshold prior"
DEFAULT_GB_BINS = 120
DEFAULT_P2_BINS = 120
DEFAULT_MAX_AUTO_WORKERS = 10
DEFAULT_FILE_CHUNK_SIZE = 500
DEFAULT_MAX_BLOCK_HISTOGRAMS = 256
FRAME_SUMMARY_COLUMNS = [
    "source_file",
    "candidate_frames",
    "selected_frames",
    "histogram_support_frames",
    "skipped_frames",
    "candidate_pairs",
    "included_attractive_pairs_weighted",
    "out_of_hist_range_pairs_weighted",
]
SELECTED_FRAME_COLUMNS = [
    "source_file",
    "timestep",
    "local_frame_index",
    "candidate_frame_index_in_file",
    "global_candidate_frame_index",
    "frame_stride",
    "frame_offset",
    "frame_sampling_weight",
    "frame_sampling_reason",
    "used",
]


def sha256_file(path: Optional[Path]) -> Optional[str]:
    if path is None or not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_dump_frames(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                count += 1
    return count


def make_frame_sampling_plan(
    files: Sequence[Path],
    *,
    every: int,
    global_frame_budget: int,
    global_frame_stride: int = 1,
    sample_seed: int = 20260429,
) -> Dict[str, object]:
    frame_counts = [count_dump_frames(path) for path in files]
    candidate_frame_counts = [int(math.ceil(count / int(every))) for count in frame_counts]
    candidate_frame_count_total = int(sum(candidate_frame_counts))
    requested_stride = max(1, int(global_frame_stride))
    frame_stride_source = "all_frames"
    if requested_stride > 1:
        frame_stride = requested_stride
        frame_offset = 0
        frame_stride_source = "global_frame_stride"
    elif int(global_frame_budget) > 0:
        frame_stride = max(1, int(math.ceil(max(candidate_frame_count_total, 1) / max(int(global_frame_budget), 1))))
        frame_offset = int(sample_seed) % frame_stride if frame_stride > 1 else 0
        frame_stride_source = "global_frame_budget"
    else:
        frame_stride = 1
        frame_offset = 0
    if frame_stride > 1:
        selected_frame_count_estimate = sum(
            1 for idx in range(candidate_frame_count_total)
            if (idx - frame_offset) % frame_stride == 0
        )
    else:
        selected_frame_count_estimate = candidate_frame_count_total
    frame_weight = (
        float(candidate_frame_count_total) / float(selected_frame_count_estimate)
        if selected_frame_count_estimate > 0
        else 1.0
    )
    offsets: List[int] = []
    offset = 0
    for count in candidate_frame_counts:
        offsets.append(offset)
        offset += int(count)
    return {
        "raw_frame_counts": frame_counts,
        "candidate_frame_counts": candidate_frame_counts,
        "raw_frame_count_total": int(sum(frame_counts)),
        "candidate_frame_count_estimate": int(candidate_frame_count_total),
        "selected_frame_count_estimate": int(selected_frame_count_estimate),
        "frame_stride": int(frame_stride),
        "frame_offset": int(frame_offset),
        "frame_stride_source": str(frame_stride_source),
        "requested_global_frame_stride": int(requested_stride),
        "frame_weight": float(frame_weight),
        "global_offsets": offsets,
    }


def resolve_worker_count(value: object, task_count: int) -> int:
    if str(value).lower() == "auto":
        max_auto = max(1, int(os.environ.get("LC_PEARL_MAX_AUTO_WORKERS", DEFAULT_MAX_AUTO_WORKERS)))
        return max(1, min(os.cpu_count() or 1, max(1, int(task_count)), max_auto))
    parsed = int(value)
    if parsed < 1:
        raise ValueError("--workers must be >= 1 or auto")
    return max(1, min(parsed, max(1, int(task_count))))


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def append_tsv_body(source: Path, target_handle) -> None:
    with source.open("r", encoding="utf-8", newline="") as handle:
        next(handle, None)
        for line in handle:
            target_handle.write(line)


def chunked(sequence: Sequence[object], chunk_size: int) -> Iterable[List[object]]:
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    for start in range(0, len(sequence), chunk_size):
        yield list(sequence[start : start + chunk_size])


def append_block_histogram(block_hist2d: List[np.ndarray], block: np.ndarray, max_saved: int) -> None:
    if max_saved <= 0:
        return
    if len(block_hist2d) < max_saved:
        block_hist2d.append(np.asarray(block, dtype=float).copy())


def histogram_rows_1d(edges: np.ndarray, counts: np.ndarray, value_name: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for idx, count in enumerate(np.asarray(counts, dtype=float)):
        if count <= 0:
            continue
        rows.append(
            {
                f"{value_name}_low": float(edges[idx]),
                f"{value_name}_high": float(edges[idx + 1]),
                f"{value_name}_center": float(0.5 * (edges[idx] + edges[idx + 1])),
                "weighted_count": float(count),
            }
        )
    return rows


def histogram_rows_2d(gb_edges: np.ndarray, p2_edges: np.ndarray, counts: np.ndarray) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for gi in range(counts.shape[0]):
        for pi in range(counts.shape[1]):
            count = float(counts[gi, pi])
            if count <= 0:
                continue
            rows.append(
                {
                    "gb_low": float(gb_edges[gi]),
                    "gb_high": float(gb_edges[gi + 1]),
                    "gb_center": float(0.5 * (gb_edges[gi] + gb_edges[gi + 1])),
                    "p2_low": float(p2_edges[pi]),
                    "p2_high": float(p2_edges[pi + 1]),
                    "p2_center": float(0.5 * (p2_edges[pi] + p2_edges[pi + 1])),
                    "weighted_count": count,
                }
            )
    return rows


def nonzero_histogram_samples(gb_edges: np.ndarray, p2_edges: np.ndarray, hist2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    gb_centers = 0.5 * (gb_edges[:-1] + gb_edges[1:])
    p2_centers = 0.5 * (p2_edges[:-1] + p2_edges[1:])
    gi, pi = np.nonzero(hist2d > 0)
    strengths = gb_centers[gi].astype(float)
    p2_values = p2_centers[pi].astype(float)
    weights = hist2d[gi, pi].astype(float)
    return strengths, p2_values, weights


def smooth_1d(values: np.ndarray, window: int = 5) -> np.ndarray:
    data = np.asarray(values, dtype=float)
    if data.size == 0:
        return data
    width = max(1, int(window))
    if width <= 1:
        return data.copy()
    if width % 2 == 0:
        width += 1
    if data.size < width:
        return data.copy()
    kernel = np.ones(width, dtype=float) / float(width)
    return np.convolve(data, kernel, mode="same")


def weighted_quantile_from_hist(centers: np.ndarray, counts: np.ndarray, quantile: float) -> float:
    x = np.asarray(centers, dtype=float)
    w = np.asarray(counts, dtype=float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0.0)
    if not np.any(valid):
        return float("nan")
    x = x[valid]
    w = w[valid]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cumulative = np.cumsum(w)
    target = min(max(float(quantile), 0.0), 1.0) * float(cumulative[-1])
    idx = int(np.searchsorted(cumulative, target, side="left"))
    return float(x[min(idx, x.size - 1)])


def strongest_two_lobe_valley(
    centers: np.ndarray,
    counts: np.ndarray,
    *,
    left_range: Tuple[float, float],
    right_range: Tuple[float, float],
    smooth_window: int = 5,
) -> Dict[str, object]:
    x = np.asarray(centers, dtype=float)
    y = np.asarray(counts, dtype=float)
    if x.size == 0 or y.size != x.size:
        return {"available": False, "reason": "empty_histogram"}
    left = np.where((x >= float(left_range[0])) & (x <= float(left_range[1])) & np.isfinite(y) & (y > 0.0))[0]
    right = np.where((x >= float(right_range[0])) & (x <= float(right_range[1])) & np.isfinite(y) & (y > 0.0))[0]
    if left.size == 0 or right.size == 0:
        return {"available": False, "reason": "missing_left_or_right_lobe"}
    left_peak = int(left[np.argmax(y[left])])
    right_peak = int(right[np.argmax(y[right])])
    if right_peak <= left_peak + 2:
        return {"available": False, "reason": "lobes_not_separated"}
    score = smooth_1d(np.log10(np.maximum(y, 0.0) + 1.0), window=smooth_window)
    valley_candidates = np.arange(left_peak + 1, right_peak, dtype=int)
    valley = int(valley_candidates[np.argmin(score[valley_candidates])])
    valley_count = float(max(y[valley], 0.0))
    left_peak_count = float(max(y[left_peak], 0.0))
    right_peak_count = float(max(y[right_peak], 0.0))
    weaker_peak = max(min(left_peak_count, right_peak_count), 1.0)
    contrast = float(weaker_peak / max(valley_count, 1.0))
    return {
        "available": True,
        "left_peak": float(x[left_peak]),
        "left_peak_count": left_peak_count,
        "right_peak": float(x[right_peak]),
        "right_peak_count": right_peak_count,
        "valley": float(x[valley]),
        "valley_count": valley_count,
        "valley_contrast": contrast,
    }


def estimate_lobe_thresholds_from_histograms(
    *,
    gb_edges: np.ndarray,
    p2_edges: np.ndarray,
    hist2d: np.ndarray,
    current: Dict[str, float],
    min_pairs: int = 100,
    min_oriented_pairs: int = 60,
) -> Dict[str, object]:
    """Select GB/P2 cuts from the full 2D histogram, not a screened edge table.

    The intended LC-Pearl V2 threshold target is the visible high-P2 lobe split:
    weak/high-orientation contacts sit near low GB attraction strength, while
    compact pearl-like contacts form a separate high-strength high-P2 lobe.
    """
    hist = np.asarray(hist2d, dtype=float)
    gb_centers = 0.5 * (np.asarray(gb_edges[:-1], dtype=float) + np.asarray(gb_edges[1:], dtype=float))
    p2_centers = 0.5 * (np.asarray(p2_edges[:-1], dtype=float) + np.asarray(p2_edges[1:], dtype=float))
    total_weight = float(np.sum(hist))
    current_clean = {
        "gb_off_strength": float(current["gb_off_strength"]),
        "gb_on_strength": float(current["gb_on_strength"]),
        "p2_cut": float(current["p2_cut"]),
    }
    if total_weight <= 0.0:
        return {
            "calibration_status": "low",
            "apply_allowed": False,
            "current": current_clean,
            "recommended": current_clean,
            "parameters": {},
            "sample_sizes": {"weighted_pair_rows": 0.0},
            "warnings": ["No attractive candidate-pair histogram entries were available."],
            "notes": ["2D lobe split did not run because the histogram was empty."],
        }

    p2_counts = hist.sum(axis=0)
    p2_valley_detail = strongest_two_lobe_valley(
        p2_centers,
        p2_counts,
        left_range=(-0.50, 0.35),
        right_range=(0.70, 1.00),
        smooth_window=5,
    )
    if bool(p2_valley_detail.get("available")):
        p2_cut_raw = float(p2_valley_detail["valley"])
        p2_cut = min(0.75, max(0.50, p2_cut_raw))
        p2_decision = "2d_hist_p2_marginal_valley_clamped_to_domain_gate"
    else:
        p2_cut = 0.50
        p2_decision = "fallback_default_pair_p2_gate"

    high_p2_gate = float(p2_cut)
    high_mask = p2_centers >= high_p2_gate
    gb_counts = hist[:, high_mask].sum(axis=1)
    if float(np.sum(gb_counts)) < float(min_oriented_pairs):
        high_p2_gate = max(0.70, float(p2_cut))
        high_mask = p2_centers >= high_p2_gate
        gb_counts = hist[:, high_mask].sum(axis=1)
    if float(np.sum(gb_counts)) < float(min_oriented_pairs):
        gb_counts = hist.sum(axis=1)
        high_p2_gate = float("nan")

    gb_valley_detail = strongest_two_lobe_valley(
        gb_centers,
        gb_counts,
        left_range=(0.0, 0.25),
        right_range=(0.25, min(float(gb_edges[-1]), 1.0)),
        smooth_window=5,
    )
    warnings: List[str] = []
    if bool(gb_valley_detail.get("available")):
        gb_on = float(gb_valley_detail["valley"])
        gb_on_decision = "2d_high_p2_lobe_valley"
    else:
        gb_on = weighted_quantile_from_hist(gb_centers, gb_counts, 0.90)
        if not math.isfinite(gb_on):
            gb_on = float(current_clean["gb_on_strength"])
        gb_on_decision = "fallback_high_p2_weighted_q90"
        warnings.append("Could not identify two separated high-P2 GB lobes; gb_on used a high-P2 weighted quantile fallback.")

    below_on = gb_centers < gb_on
    shoulder = weighted_quantile_from_hist(gb_centers[below_on], gb_counts[below_on], 0.99)
    if not math.isfinite(shoulder):
        shoulder = max(gb_on - 0.02, 0.02)
    gb_off = float(shoulder)
    gb_off = min(max(gb_off, 0.02), max(gb_on - 0.02, 0.02))
    if gb_on <= gb_off:
        gb_on = min(float(gb_edges[-1]), gb_off + 0.02)

    contrast = float(gb_valley_detail.get("valley_contrast", 0.0)) if isinstance(gb_valley_detail, dict) else 0.0
    status = "high" if total_weight >= float(min_pairs) and contrast >= 50.0 else "medium" if total_weight >= float(min_pairs) and contrast >= 5.0 else "low"
    apply_allowed = bool(gb_on > gb_off)
    if not apply_allowed:
        warnings.append("2D lobe split produced non-monotonic thresholds; inspect the lobe split preview before applying.")

    recommended = {
        "gb_off_strength": cal.rounded(float(gb_off)),
        "gb_on_strength": cal.rounded(float(gb_on)),
        "p2_cut": cal.rounded(float(p2_cut)),
    }
    return {
        "calibration_status": status,
        "apply_allowed": bool(apply_allowed),
        "current": current_clean,
        "recommended": recommended,
        "delta_recommended_minus_current": {
            key: cal.rounded(float(recommended[key]) - float(current_clean[key]))
            for key in recommended
        },
        "parameters": {
            "method": "2d_lobe_split",
            "p2_cut": {
                **p2_valley_detail,
                "estimate": float(p2_cut),
                "decision": p2_decision,
            },
            "gb_on_strength": {
                **gb_valley_detail,
                "estimate": float(gb_on),
                "decision": gb_on_decision,
                "orientation_gate_for_gb_histogram": high_p2_gate,
            },
            "gb_off_strength": {
                "estimate": float(gb_off),
                "decision": "left_lobe_q99_shoulder_bounded_below_gb_on",
                "left_lobe_q99_shoulder": float(shoulder),
            },
        },
        "parameter_confidence": {
            "gb_off_strength": status,
            "gb_on_strength": status,
            "p2_cut": "medium" if bool(p2_valley_detail.get("available")) else "low",
            "joint2d": status,
        },
        "joint_2d": {
            "status": status,
            "method": "direct_2d_lobe_split_from_stream_histogram",
            "selected": recommended,
        },
        "sample_sizes": {
            "nonzero_histogram_cells": int(np.count_nonzero(hist > 0.0)),
            "weighted_pair_rows": float(round(total_weight, 6)),
            "high_p2_weighted_pair_rows": float(round(float(np.sum(gb_counts)), 6)),
        },
        "warnings": warnings,
        "notes": [
            "Thresholds were selected from the full streaming GB-strength x P2 histogram.",
            "gb_on is the high-P2 GB-strength valley between the weak-contact lobe and the strong-attraction lobe.",
            "gb_off is a conservative left-lobe shoulder for gray/support contacts, not the main robust-domain boundary.",
            "apply_allowed now means only that recommended thresholds are numerically usable; low/medium/high status is descriptive, not a time/block gate.",
        ],
    }


def estimate_thresholds_from_histograms(
    *,
    gb_edges: np.ndarray,
    p2_edges: np.ndarray,
    hist2d: np.ndarray,
    current: Dict[str, float],
    min_pairs: int = 100,
    min_oriented_pairs: int = 60,
    independent_pair_blocks: int = 2,
    min_pair_blocks: int = 1,
) -> Dict[str, object]:
    _ = independent_pair_blocks, min_pair_blocks
    return estimate_lobe_thresholds_from_histograms(
        gb_edges=gb_edges,
        p2_edges=p2_edges,
        hist2d=hist2d,
        current=current,
        min_pairs=min_pairs,
        min_oriented_pairs=min_oriented_pairs,
    )


def draw_threshold_lines(ax, recommended: Dict[str, object]) -> None:
    if not recommended:
        return
    annotate_vline(ax, float(recommended["gb_off_strength"]), threshold_label("gb_off", recommended["gb_off_strength"]), color="#f59e0b")
    annotate_vline(ax, float(recommended["gb_on_strength"]), threshold_label("gb_on", recommended["gb_on_strength"]), color="#ef4444", ymax=0.82)
    annotate_hline(ax, float(recommended["p2_cut"]), threshold_label("p2_cut", recommended["p2_cut"]), color="#38bdf8")


def draw_stream_dotgrid(ax, gb_edges: np.ndarray, p2_edges: np.ndarray, hist2d: np.ndarray, *, title: str):
    gb_centers = 0.5 * (gb_edges[:-1] + gb_edges[1:])
    p2_centers = 0.5 * (p2_edges[:-1] + p2_edges[1:])
    gi, pi = np.nonzero(hist2d > 0.0)
    values = np.log10(hist2d[gi, pi] + 1.0)
    scatter = ax.scatter(
        gb_centers[gi],
        p2_centers[pi],
        c=values,
        s=18,
        marker="s",
        linewidths=0.0,
        cmap="viridis",
    )
    ax.set_xlabel("GB attraction strength")
    ax.set_ylabel("pair P2")
    ax.set_title(title)
    return scatter


def write_lobe_split_plots(output_dir: Path, gb_edges: np.ndarray, p2_edges: np.ndarray, hist2d: np.ndarray, recommendation: Dict[str, object]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        (output_dir / "plot_error.txt").write_text(str(exc), encoding="utf-8")
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    recommended = recommendation.get("recommended", {})
    fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=180)
    scatter = draw_stream_dotgrid(ax, gb_edges, p2_edges, hist2d, title="LC Domain-Pearl V2 2D lobe split")
    draw_threshold_lines(ax, recommended)
    fig.colorbar(scatter, ax=ax, label="log10(weighted count + 1)")
    fig.tight_layout()
    fig.savefig(output_dir / "gb_strength_vs_p2_stream_lobe_split_dotgrid.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.3), dpi=180)
    mesh = axes[0].pcolormesh(gb_edges, p2_edges, np.log10(hist2d.T + 1.0), cmap="viridis", shading="auto")
    draw_threshold_lines(axes[0], recommended)
    axes[0].set_xlabel("GB attraction strength")
    axes[0].set_ylabel("pair P2")
    axes[0].set_title("Full streaming histogram")
    scatter = draw_stream_dotgrid(axes[1], gb_edges, p2_edges, hist2d, title="Dot-grid view of same full histogram")
    draw_threshold_lines(axes[1], recommended)
    fig.colorbar(mesh, ax=axes[0], label="log10(weighted count + 1)")
    fig.colorbar(scatter, ax=axes[1], label="log10(weighted count + 1)")
    fig.tight_layout()
    fig.savefig(output_dir / "gb_strength_vs_p2_stream_lobe_split_comparison.png")
    plt.close(fig)


def write_stream_plots(output_dir: Path, gb_edges: np.ndarray, p2_edges: np.ndarray, hist2d: np.ndarray, recommendation: Dict[str, object]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        (output_dir / "plot_error.txt").write_text(str(exc), encoding="utf-8")
        return
    recommended = recommendation.get("recommended", {})
    fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=180)
    mesh = ax.pcolormesh(gb_edges, p2_edges, np.log10(hist2d.T + 1.0), cmap="viridis", shading="auto")
    draw_threshold_lines(ax, recommended)
    ax.set_xlabel("GB attraction strength")
    ax.set_ylabel("pair P2")
    ax.set_title("LC Domain-Pearl V2 streaming GB strength x P2 histogram")
    fig.colorbar(mesh, ax=ax, label="log10(weighted count + 1)")
    fig.tight_layout()
    fig.savefig(output_dir / "gb_strength_vs_p2_stream_hist.png")
    plt.close(fig)
    write_lobe_split_plots(output_dir, gb_edges, p2_edges, hist2d, recommendation)


def reservoir_update(
    sample: List[Dict[str, object]],
    rows: Sequence[Dict[str, object]],
    *,
    seen_count: int,
    max_items: int,
    rng: np.random.Generator,
) -> int:
    if max_items <= 0:
        return seen_count + len(rows)
    for row in rows:
        seen_count += 1
        if len(sample) < max_items:
            sample.append(dict(row))
            continue
        idx = int(rng.integers(0, seen_count))
        if idx < max_items:
            sample[idx] = dict(row)
    return seen_count


def weighted_reservoir_update(
    sample: List[Dict[str, object]],
    rows: Sequence[Dict[str, object]],
    weights: Sequence[float],
    *,
    seen_weight: float,
    max_items: int,
    rng: np.random.Generator,
) -> float:
    """Merge chunk-level reservoirs into an approximately global weighted reservoir."""
    total_weight = float(seen_weight)
    if max_items <= 0:
        return total_weight + float(np.sum(np.asarray(weights, dtype=float))) if weights else total_weight
    for row, raw_weight in zip(rows, weights):
        weight = float(raw_weight)
        if not math.isfinite(weight) or weight <= 0.0:
            continue
        total_weight += weight
        key = math.log(max(float(rng.random()), 1e-300)) / weight
        candidate = dict(row)
        candidate["_weighted_reservoir_key"] = key
        if len(sample) < max_items:
            sample.append(candidate)
            continue
        min_idx, min_item = min(enumerate(sample), key=lambda item: float(item[1].get("_weighted_reservoir_key", -math.inf)))
        if key > float(min_item.get("_weighted_reservoir_key", -math.inf)):
            sample[min_idx] = candidate
    return total_weight


def build_config(args: argparse.Namespace) -> lca.AggregationConfig:
    config = gb_pair_audit.build_config(args)
    return config


def run_chunk_jobs_bounded(
    chunk_jobs: Sequence[Tuple[Sequence[Tuple[Path, int]], lca.AggregationConfig, Dict[str, object]]],
    *,
    worker_count: int,
    merge_result,
) -> None:
    """Submit only a small window of chunk jobs so completed results cannot pile up."""
    if worker_count <= 1 or len(chunk_jobs) <= 1:
        for job in chunk_jobs:
            merge_result(process_file_chunk_for_streaming(job))
        return
    max_pending = max(1, int(worker_count) * 2)
    iterator = iter(chunk_jobs)
    pending = set()
    with ProcessPoolExecutor(max_workers=int(worker_count)) as executor:
        for _ in range(min(max_pending, len(chunk_jobs))):
            try:
                pending.add(executor.submit(process_file_chunk_for_streaming, next(iterator)))
            except StopIteration:
                break
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                merge_result(future.result())
            for _ in range(min(len(done), max_pending - len(pending))):
                try:
                    pending.add(executor.submit(process_file_chunk_for_streaming, next(iterator)))
                except StopIteration:
                    break


def histogram_bin(edges: np.ndarray, value: float) -> Optional[int]:
    if not math.isfinite(value) or value < float(edges[0]) or value > float(edges[-1]):
        return None
    if value == float(edges[-1]):
        return int(edges.size - 2)
    idx = int(np.searchsorted(edges, value, side="right") - 1)
    if idx < 0 or idx >= edges.size - 1:
        return None
    return idx


def collect_frame_streaming_histograms(
    *,
    source_file: Path,
    timestep: int,
    box: lca.BoxSpec,
    data: np.ndarray,
    col_index: Dict[str, int],
    config: lca.AggregationConfig,
    gb_edges: np.ndarray,
    p2_edges: np.ndarray,
    local_pairs: set[Tuple[int, int]],
    excluded_pairs: set[Tuple[int, int]],
    params: lca.GayBerneParams,
    broad_cut: float,
    frame_index: int,
    frame_weight: float,
    frame_sampling_reason: str,
    sample: List[Dict[str, object]],
    sample_seen: int,
    audit_examples: int,
    rng: np.random.Generator,
) -> Dict[str, object]:
    all_pos = lca.extract_positions(data, col_index)
    all_ids = lca.extract_particle_ids(data, col_index)
    all_types = lca.extract_particle_types(data, col_index)
    all_shapes = lca.extract_shape_axes(data, col_index)
    if all_shapes is None:
        raise RuntimeError(f"{source_file}: GB threshold prior requires shapex/shapey/shapez columns")
    all_quaternions = lca.extract_quaternions(data, col_index)
    mesogen_indices = lca.select_mesogen_indices(all_types, config.mesogen_type)
    hist_gb = np.zeros(gb_edges.size - 1, dtype=float)
    hist_p2 = np.zeros(p2_edges.size - 1, dtype=float)
    hist2d = np.zeros((gb_edges.size - 1, p2_edges.size - 1), dtype=float)
    local_hist2d = np.zeros_like(hist2d)
    nonlocal_hist2d = np.zeros_like(hist2d)
    if mesogen_indices.size == 0:
        return {
            "candidate_pairs": 0,
            "included_attractive_pairs_weighted": 0.0,
            "out_of_hist_range_pairs_weighted": 0.0,
            "histogram_support": False,
            "hist_gb": hist_gb,
            "hist_p2": hist_p2,
            "hist2d": hist2d,
            "local_hist2d": local_hist2d,
            "nonlocal_hist2d": nonlocal_hist2d,
            "sample_seen": int(sample_seen),
        }
    pos = all_pos[mesogen_indices, :]
    ids = all_ids[mesogen_indices].astype(int)
    shapes = all_shapes[mesogen_indices, :]
    quats = all_quaternions[mesogen_indices, :]
    u, _axis_used = lca.extract_orientations(data[mesogen_indices, :], col_index, config.axis)
    chain_indices = lca.infer_chain_indices(
        data,
        col_index,
        mesogen_indices,
        all_ids,
        positions=all_pos,
        types=all_types,
        anchor_types=config.anchor_types,
    )
    weight_semantics, inclusion_probability_known = gb_pair_audit.sampling_semantics("global_stream", frame_sampling_reason)
    candidate_pairs = 0
    included_attractive_pairs = 0.0
    out_of_hist_range_pairs = 0.0
    for i in range(pos.shape[0] - 1):
        delta = pos[i + 1 :, :] - pos[i, :]
        delta[:, 0] = lca.minimum_image(delta[:, 0], box[0])
        delta[:, 1] = lca.minimum_image(delta[:, 1], box[1])
        delta[:, 2] = lca.minimum_image(delta[:, 2], box[2])
        distances = np.linalg.norm(delta, axis=1)
        close = np.nonzero((distances > 1e-12) & (distances <= broad_cut))[0] + i + 1
        for neighbor in close:
            j = int(neighbor)
            delta_vec = lca.minimum_image_vector(pos[j, :] - pos[i, :], box)
            metrics = lca.gayberne_pair_metrics(
                r12=delta_vec,
                quat_i=quats[i, :],
                quat_j=quats[j, :],
                shape_i=shapes[i, :],
                shape_j=shapes[j, :],
                params=params,
            )
            if metrics is None:
                continue
            candidate_pairs += 1
            pair_energy, well_depth = metrics
            atom_i = int(ids[i])
            atom_j = int(ids[j])
            pair_key = lca.normalized_pair(atom_i, atom_j)
            q_score = float(abs(np.clip(np.dot(u[i, :], u[j, :]), -1.0, 1.0)))
            p2_score = float(lca.p2(q_score))
            delta_s = int(abs(int(chain_indices[i]) - int(chain_indices[j])))
            is_local = int(delta_s <= config.s_excl or pair_key in local_pairs)
            is_excluded = int(pair_key in excluded_pairs)
            attraction_strength = max(0.0, -float(pair_energy) / max(float(well_depth), 1e-24))
            if (
                is_excluded
                or not math.isfinite(attraction_strength)
                or not math.isfinite(p2_score)
                or attraction_strength <= 0.0
            ):
                continue
            weight = float(frame_weight)
            gb_idx = histogram_bin(gb_edges, attraction_strength)
            p2_idx = histogram_bin(p2_edges, p2_score)
            if gb_idx is None or p2_idx is None:
                out_of_hist_range_pairs += weight
            else:
                hist_gb[gb_idx] += weight
                hist_p2[p2_idx] += weight
                hist2d[gb_idx, p2_idx] += weight
                if is_local:
                    local_hist2d[gb_idx, p2_idx] += weight
                else:
                    nonlocal_hist2d[gb_idx, p2_idx] += weight
                included_attractive_pairs += weight
            row = {
                "source_file": str(source_file),
                "timestep": int(timestep),
                "atom_i": atom_i,
                "atom_j": atom_j,
                "distance": float(np.linalg.norm(delta_vec)),
                "pair_energy": float(pair_energy),
                "well_depth": float(well_depth),
                "attraction_strength": float(attraction_strength),
                "q_score": q_score,
                "p2_score": p2_score,
                "delta_s": delta_s,
                "is_local": is_local,
                "is_excluded": is_excluded,
                "frame_index": int(frame_index),
                "frame_weight": float(frame_weight),
                "frame_sampling_reason": str(frame_sampling_reason),
                "frame_event_score": 0.0,
                "weight_semantics": weight_semantics,
                "inclusion_probability_known": inclusion_probability_known,
                "sampling_strata": "global_stream_all_pairs",
                "sampling_weight": weight,
                "sampling_probability": 1.0 if inclusion_probability_known else "",
                "pair_sampling_rate": 1.0,
                "candidate_pairs_in_stratum": 1,
                "selected_pairs_in_stratum": 1,
                "sampling_selected_reason": "stream_reservoir",
            }
            sample_seen = reservoir_update(sample, [row], seen_count=sample_seen, max_items=audit_examples, rng=rng)
    return {
        "candidate_pairs": int(candidate_pairs),
        "included_attractive_pairs_weighted": float(included_attractive_pairs),
        "out_of_hist_range_pairs_weighted": float(out_of_hist_range_pairs),
        "histogram_support": bool(np.sum(hist2d) > 0.0),
        "hist_gb": hist_gb,
        "hist_p2": hist_p2,
        "hist2d": hist2d,
        "local_hist2d": local_hist2d,
        "nonlocal_hist2d": nonlocal_hist2d,
        "sample_seen": int(sample_seen),
    }


def process_file_for_streaming(args_tuple: Tuple[Path, lca.AggregationConfig, Dict[str, object]]) -> Dict[str, object]:
    dump_path, config, options = args_tuple
    if "chunk_manifest_dir" not in options:
        options = dict(options)
        options["chunk_manifest_dir"] = str(Path.cwd() / "_threshold_prior_single_file_manifest")
        options.setdefault("chunk_id", 0)
    return process_file_chunk_for_streaming(([(dump_path, int(options.get("global_offset", 0)))], config, options))


def process_file_chunk_for_streaming(
    args_tuple: Tuple[Sequence[Tuple[Path, int]], lca.AggregationConfig, Dict[str, object]]
) -> Dict[str, object]:
    file_jobs, config, options = args_tuple
    gb_edges = np.asarray(options["gb_edges"], dtype=float)
    p2_edges = np.asarray(options["p2_edges"], dtype=float)
    every = int(options["every"])
    frame_stride = int(options["frame_stride"])
    frame_offset = int(options["frame_offset"])
    frame_stride_source = str(options.get("frame_stride_source", "all_frames"))
    frame_weight = float(options["frame_weight"])
    block_size = int(options["block_size"])
    audit_examples = int(options["audit_examples"])
    max_block_histograms = int(options.get("max_block_histograms", DEFAULT_MAX_BLOCK_HISTOGRAMS))
    sample_seed = int(options["sample_seed"])
    hist_gb = np.zeros(gb_edges.size - 1, dtype=float)
    hist_p2 = np.zeros(p2_edges.size - 1, dtype=float)
    hist2d = np.zeros((gb_edges.size - 1, p2_edges.size - 1), dtype=float)
    local_hist2d = np.zeros_like(hist2d)
    nonlocal_hist2d = np.zeros_like(hist2d)
    block_hist2d: List[np.ndarray] = []
    block_count = 0
    current_block = np.zeros_like(hist2d)
    current_block_frames = 0
    sample: List[Dict[str, object]] = []
    sample_seen = 0
    seed_basis = [str(path.resolve()) for path, _offset in file_jobs[:3]]
    if file_jobs:
        seed_basis.append(str(file_jobs[-1][0].resolve()))
    rng = np.random.default_rng(gb_pair_audit.stable_sample_seed("threshold_prior_chunk", tuple(seed_basis), len(file_jobs), sample_seed))
    chunk_id = int(options.get("chunk_id", 0))
    manifest_dir = Path(str(options["chunk_manifest_dir"]))
    manifest_dir.mkdir(parents=True, exist_ok=True)
    frame_summary_path = manifest_dir / f"frame_stream_summary_chunk_{chunk_id:06d}.tsv"
    selected_frame_path = manifest_dir / f"selected_frame_manifest_chunk_{chunk_id:06d}.tsv"
    local_pairs = lca.read_pair_list(config.local_pair_file)
    excluded_pairs = lca.read_pair_list(config.exclude_pair_file)
    params = config.gayberne_params
    if params is None:
        raise RuntimeError("threshold prior requires gayberne_params")
    broad_cut = float(config.r_energy_cap or params.cutoff)
    frame_summary_handle = frame_summary_path.open("w", encoding="utf-8", newline="")
    try:
        selected_frame_handle = selected_frame_path.open("w", encoding="utf-8", newline="")
    except Exception:
        frame_summary_handle.close()
        raise
    frame_summary_writer = csv.DictWriter(frame_summary_handle, fieldnames=FRAME_SUMMARY_COLUMNS, delimiter="\t")
    selected_frame_writer = csv.DictWriter(selected_frame_handle, fieldnames=SELECTED_FRAME_COLUMNS, delimiter="\t")
    try:
        frame_summary_writer.writeheader()
        selected_frame_writer.writeheader()
        for dump_path, global_offset in file_jobs:
            candidate_frames = 0
            selected_frames = 0
            histogram_support_frames = 0
            candidate_pairs = 0
            included_attractive_pairs = 0.0
            out_of_hist_range_pairs = 0.0
            skipped_frames = 0
            local_frame_index = 0
            for timestep, box, _columns, col_index, data in lca.parse_dump_frames(dump_path):
                local_frame_index += 1
                if (local_frame_index - 1) % every != 0:
                    continue
                candidate_frames += 1
                global_frame_index = int(global_offset) + candidate_frames - 1
                if frame_stride > 1 and (global_frame_index - frame_offset) % frame_stride != 0:
                    continue
                frame_sampling_reason = (
                    f"global_stream_stride_{frame_stride_source}"
                    if frame_stride > 1
                    else f"global_stream_{frame_stride_source}"
                )
                frame_stats = collect_frame_streaming_histograms(
                    source_file=dump_path,
                    timestep=int(timestep),
                    box=box,
                    data=data,
                    col_index=col_index,
                    config=config,
                    gb_edges=gb_edges,
                    p2_edges=p2_edges,
                    local_pairs=local_pairs,
                    excluded_pairs=excluded_pairs,
                    params=params,
                    broad_cut=broad_cut,
                    frame_index=local_frame_index,
                    frame_weight=float(frame_weight),
                    frame_sampling_reason=frame_sampling_reason,
                    sample=sample,
                    sample_seen=sample_seen,
                    audit_examples=audit_examples,
                    rng=rng,
                )
                sample_seen = int(frame_stats["sample_seen"])
                selected_frame_writer.writerow(
                    {
                        "source_file": str(dump_path),
                        "timestep": int(timestep),
                        "local_frame_index": int(local_frame_index),
                        "candidate_frame_index_in_file": int(candidate_frames - 1),
                        "global_candidate_frame_index": int(global_frame_index),
                        "frame_stride": int(frame_stride),
                        "frame_offset": int(frame_offset),
                        "frame_sampling_weight": float(frame_weight),
                        "frame_sampling_reason": frame_sampling_reason,
                        "used": 1,
                    }
                )
                selected_frames += 1
                candidate_pairs += int(frame_stats["candidate_pairs"])
                included_attractive_pairs += float(frame_stats["included_attractive_pairs_weighted"])
                out_of_hist_range_pairs += float(frame_stats["out_of_hist_range_pairs_weighted"])
                if bool(frame_stats["histogram_support"]):
                    histogram_support_frames += 1
                hist_gb += frame_stats["hist_gb"]
                hist_p2 += frame_stats["hist_p2"]
                frame_hist2d = frame_stats["hist2d"]
                hist2d += frame_hist2d
                local_hist2d += frame_stats["local_hist2d"]
                nonlocal_hist2d += frame_stats["nonlocal_hist2d"]
                if bool(frame_stats["histogram_support"]):
                    current_block += frame_hist2d
                    current_block_frames += 1
                    if current_block_frames >= block_size:
                        block_count += 1
                        append_block_histogram(block_hist2d, current_block, max_block_histograms)
                        current_block = np.zeros_like(hist2d)
                        current_block_frames = 0
            frame_summary_writer.writerow(
                {
                    "source_file": str(dump_path),
                    "candidate_frames": int(candidate_frames),
                    "selected_frames": int(selected_frames),
                    "histogram_support_frames": int(histogram_support_frames),
                    "skipped_frames": int(skipped_frames),
                    "candidate_pairs": int(candidate_pairs),
                    "included_attractive_pairs_weighted": float(included_attractive_pairs),
                    "out_of_hist_range_pairs_weighted": float(out_of_hist_range_pairs),
                }
            )
    finally:
        frame_summary_handle.close()
        selected_frame_handle.close()
    if current_block_frames:
        block_count += 1
        append_block_histogram(block_hist2d, current_block, max_block_histograms)
    return {
        "chunk_file_count": int(len(file_jobs)),
        "frame_summary_path": str(frame_summary_path),
        "selected_frame_manifest_path": str(selected_frame_path),
        "hist_gb": hist_gb,
        "hist_p2": hist_p2,
        "hist2d": hist2d,
        "local_hist2d": local_hist2d,
        "nonlocal_hist2d": nonlocal_hist2d,
        "block_hist2d": block_hist2d,
        "block_count": int(block_count),
        "audit_sample": sample,
        "audit_sample_seen": int(sample_seen),
    }


def build_streaming_prior(args: argparse.Namespace) -> Dict[str, object]:
    inputs = args.inputs if args.inputs else [Path.cwd()]
    files = lca.iter_input_files(inputs, pattern=args.pattern, recursive=bool(args.recursive))
    if not files:
        raise SystemExit("No dump files matched the given inputs.")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    gb_edges = np.linspace(0.0, float(args.gb_strength_max), int(args.gb_bins) + 1)
    p2_edges = np.linspace(-0.5, 1.0, int(args.p2_bins) + 1)
    frame_plan = make_frame_sampling_plan(
        files,
        every=int(args.every),
        global_frame_budget=int(args.global_frame_budget),
        global_frame_stride=int(args.global_frame_stride),
        sample_seed=int(args.sample_seed),
    )
    global_offsets = list(frame_plan["global_offsets"])
    frame_stride = int(frame_plan["frame_stride"])
    frame_offset = int(frame_plan["frame_offset"])
    frame_stride_source = str(frame_plan["frame_stride_source"])
    frame_weight = float(frame_plan["frame_weight"])
    config = build_config(args)
    options = {
        "gb_edges": gb_edges,
        "p2_edges": p2_edges,
        "every": int(args.every),
        "frame_stride": int(frame_stride),
        "frame_offset": int(frame_offset),
        "frame_stride_source": str(frame_stride_source),
        "frame_weight": float(frame_weight),
        "block_size": int(args.block_size_frames),
        "audit_examples": int(args.audit_example_pairs),
        "sample_seed": int(args.sample_seed),
    }
    file_jobs = list(zip(files, global_offsets))
    chunk_size = int(args.file_chunk_size)
    file_chunks = list(chunked(file_jobs, chunk_size))
    per_chunk_block_cap = (
        0
        if int(args.max_block_histograms) <= 0
        else max(1, int(math.ceil(int(args.max_block_histograms) / max(len(file_chunks), 1))))
    )
    chunk_manifest_dir = output_dir / "_threshold_prior_chunks"
    chunk_manifest_dir.mkdir(parents=True, exist_ok=True)
    for stale in chunk_manifest_dir.glob("*.tsv"):
        stale.unlink()
    chunk_jobs = []
    for chunk_id, chunk in enumerate(file_chunks):
        per_chunk_options = dict(options)
        per_chunk_options["chunk_id"] = int(chunk_id)
        per_chunk_options["chunk_manifest_dir"] = str(chunk_manifest_dir)
        per_chunk_options["max_block_histograms"] = int(per_chunk_block_cap)
        chunk_jobs.append((chunk, config, per_chunk_options))
    worker_count = resolve_worker_count(args.workers, len(chunk_jobs))
    print(
        "[LC-Pearl threshold-prior] "
        f"files={len(files)} chunks={len(chunk_jobs)} chunk_size={chunk_size} "
        f"workers={worker_count} candidate_frames={frame_plan['candidate_frame_count_estimate']} "
        f"selected_frames_estimate={frame_plan['selected_frame_count_estimate']} "
        f"frame_stride={frame_stride} frame_stride_source={frame_stride_source}"
    )
    hist_gb = np.zeros(gb_edges.size - 1, dtype=float)
    hist_p2 = np.zeros(p2_edges.size - 1, dtype=float)
    hist2d = np.zeros((gb_edges.size - 1, p2_edges.size - 1), dtype=float)
    local_hist2d = np.zeros_like(hist2d)
    nonlocal_hist2d = np.zeros_like(hist2d)
    block_hist2d: List[np.ndarray] = []
    chunk_block_count = 0
    candidate_frame_count_total = 0
    selected_frame_count_total = 0
    histogram_support_frame_count_total = 0
    candidate_pairs_total = 0
    audit_sample: List[Dict[str, object]] = []
    audit_reservoir_seen_weight = 0.0
    audit_seen_total = 0
    out_of_hist_range_pairs = 0.0
    rng = np.random.default_rng(int(args.sample_seed))
    frame_summary_handle = (output_dir / "frame_stream_summary.tsv").open("w", encoding="utf-8", newline="")
    try:
        selected_frame_handle = (output_dir / "selected_frame_manifest.tsv").open("w", encoding="utf-8", newline="")
    except Exception:
        frame_summary_handle.close()
        raise
    frame_summary_writer = csv.DictWriter(frame_summary_handle, fieldnames=FRAME_SUMMARY_COLUMNS, delimiter="\t")
    selected_frame_writer = csv.DictWriter(selected_frame_handle, fieldnames=SELECTED_FRAME_COLUMNS, delimiter="\t")
    try:
        frame_summary_writer.writeheader()
        selected_frame_writer.writeheader()
    except Exception:
        frame_summary_handle.close()
        selected_frame_handle.close()
        raise

    def merge_result(result: Dict[str, object]) -> None:
        nonlocal hist_gb, hist_p2, hist2d, local_hist2d, nonlocal_hist2d
        nonlocal chunk_block_count, audit_reservoir_seen_weight, audit_seen_total, out_of_hist_range_pairs
        nonlocal candidate_frame_count_total, selected_frame_count_total, histogram_support_frame_count_total, candidate_pairs_total
        hist_gb += result["hist_gb"]
        hist_p2 += result["hist_p2"]
        hist2d += result["hist2d"]
        local_hist2d += result["local_hist2d"]
        nonlocal_hist2d += result["nonlocal_hist2d"]
        chunk_block_count += int(result.get("block_count", 0))
        remaining_blocks = max(0, int(args.max_block_histograms) - len(block_hist2d))
        if remaining_blocks:
            block_hist2d.extend(result["block_hist2d"][:remaining_blocks])
        append_tsv_body(Path(str(result["selected_frame_manifest_path"])), selected_frame_handle)
        with Path(str(result["frame_summary_path"])).open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                frame_summary_writer.writerow({column: row.get(column, "") for column in FRAME_SUMMARY_COLUMNS})
                candidate_frame_count_total += int(row["candidate_frames"])
                selected_frame_count_total += int(row["selected_frames"])
                histogram_support_frame_count_total += int(row.get("histogram_support_frames", 0))
                candidate_pairs_total += int(row["candidate_pairs"])
                out_of_hist_range_pairs += float(row["out_of_hist_range_pairs_weighted"])
        for key in ("frame_summary_path", "selected_frame_manifest_path"):
            try:
                Path(str(result[key])).unlink()
            except OSError:
                pass
        chunk_seen = int(result.get("audit_sample_seen", 0))
        audit_seen_total += chunk_seen
        chunk_sample = list(result["audit_sample"])
        chunk_weight = float(chunk_seen) / max(len(chunk_sample), 1) if chunk_sample else 0.0
        audit_reservoir_seen_weight = weighted_reservoir_update(
            audit_sample,
            chunk_sample,
            [chunk_weight] * len(chunk_sample),
            seen_weight=audit_reservoir_seen_weight,
            max_items=int(args.audit_example_pairs),
            rng=rng,
        )

    try:
        run_chunk_jobs_bounded(chunk_jobs, worker_count=worker_count, merge_result=merge_result)
    finally:
        frame_summary_handle.close()
        selected_frame_handle.close()
    try:
        chunk_manifest_dir.rmdir()
    except OSError:
        pass
    confidence_block_count = max(1, int(histogram_support_frame_count_total) // max(1, int(args.block_size_frames)))
    recommendation = estimate_thresholds_from_histograms(
        gb_edges=gb_edges,
        p2_edges=p2_edges,
        hist2d=hist2d,
        current={
            "gb_off_strength": float(args.current_gb_off),
            "gb_on_strength": float(args.current_gb_on),
            "p2_cut": float(args.current_p2_cut),
        },
        min_pairs=int(args.min_pairs),
        min_oriented_pairs=int(args.min_oriented_pairs),
        independent_pair_blocks=int(confidence_block_count),
        min_pair_blocks=2,
    )
    hist_range_fraction = (
        float(out_of_hist_range_pairs) / float(out_of_hist_range_pairs + np.sum(hist2d))
        if float(out_of_hist_range_pairs + np.sum(hist2d)) > 0.0
        else 0.0
    )
    if hist_range_fraction > 0.01:
        recommendation.setdefault("warnings", []).append(
            f"{hist_range_fraction:.3%} of weighted attractive pairs fell outside histogram ranges; increase --gb-strength-max or inspect stream_histograms.npz. "
            "LC-Pearl V2 records this as coverage metadata but does not use it as a threshold-application gate."
        )
    artifact = {
        "schema_version": SCHEMA_VERSION,
        "method_name": METHOD_NAME,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        **recommendation,
        "algorithm_fingerprint": {
            "schema_version": SCHEMA_VERSION,
            "method_name": METHOD_NAME,
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "gb_bins": int(args.gb_bins),
            "p2_bins": int(args.p2_bins),
            "gb_strength_max": float(args.gb_strength_max),
            "min_pairs": int(args.min_pairs),
            "min_oriented_pairs": int(args.min_oriented_pairs),
            "calibration_mode": "joint2d_streaming_histogram",
        },
        "build_parameters": {
            "every": int(args.every),
            "global_frame_budget": int(args.global_frame_budget),
            "global_frame_stride": int(args.global_frame_stride),
            "block_size_frames": int(args.block_size_frames),
            "file_chunk_size": int(args.file_chunk_size),
            "max_block_histograms": int(args.max_block_histograms),
            "audit_example_pairs": int(args.audit_example_pairs),
            "sample_seed": int(args.sample_seed),
        },
        "input_provenance": {
            "inputs": [str(path.resolve()) for path in inputs],
            "input_file_count": int(len(files)),
            "raw_frame_count_total": int(frame_plan["raw_frame_count_total"]),
            "candidate_frame_count_estimate": int(frame_plan["candidate_frame_count_estimate"]),
            "pattern": str(args.pattern),
            "recursive": bool(args.recursive),
            "gb_param_file": str(args.gb_param_file.resolve()),
            "gb_param_file_sha256": sha256_file(args.gb_param_file),
            "mesogen_type": int(args.mesogen_type),
            "anchor_types": str(args.anchor_types),
            "axis": str(args.axis),
            "s_excl": int(args.s_excl),
            "r_energy_cap": str(args.r_energy_cap),
            "local_pair_file": str(args.local_pair_file.resolve()) if args.local_pair_file else None,
            "local_pair_file_sha256": sha256_file(args.local_pair_file),
            "exclude_pair_file": str(args.exclude_pair_file.resolve()) if args.exclude_pair_file else None,
            "exclude_pair_file_sha256": sha256_file(args.exclude_pair_file),
        },
        "streaming": {
            "candidate_frame_count_total": int(candidate_frame_count_total),
            "selected_frame_count_total": int(selected_frame_count_total),
            "histogram_support_frame_count_total": int(histogram_support_frame_count_total),
            "selected_frame_count_estimate": int(frame_plan["selected_frame_count_estimate"]),
            "frame_stride": int(frame_stride),
            "frame_offset": int(frame_offset),
            "frame_stride_source": str(frame_stride_source),
            "frame_sampling_weight": float(frame_weight),
            "global_frame_budget": int(args.global_frame_budget),
            "global_frame_stride": int(args.global_frame_stride),
            "frame_budget_semantics": "global_frame_stride=1 means no explicit global stride; 10 means every 10th global candidate frame; 100 means every 100th. If global_frame_stride > 1 it takes precedence over global_frame_budget. global_frame_budget=0 means all remaining frames; positive budget derives a seeded deterministic stride over ordered input files/frames.",
            "candidate_pairs_total": int(candidate_pairs_total),
            "included_attractive_pairs_weighted": float(round(float(np.sum(hist2d)), 6)),
            "out_of_hist_range_pairs_weighted": float(round(float(out_of_hist_range_pairs), 6)),
            "out_of_hist_range_pair_fraction": float(round(hist_range_fraction, 6)),
            "gb_bins": int(args.gb_bins),
            "p2_bins": int(args.p2_bins),
            "gb_strength_max": float(args.gb_strength_max),
            "block_count": int(confidence_block_count),
            "block_count_semantics": "Full histogram-supported frame blocks recorded for provenance only. LC-Pearl V2 threshold application is not gated by time/block confidence.",
            "chunk_local_block_count": int(chunk_block_count),
            "saved_block_histogram_count": int(len(block_hist2d)),
            "max_saved_block_histograms": int(args.max_block_histograms),
            "block_size_frames": int(args.block_size_frames),
            "file_chunk_size": int(args.file_chunk_size),
            "chunk_count": int(len(chunk_jobs)),
            "audit_example_pairs": int(len(audit_sample)),
            "audit_example_seen_pairs": int(audit_seen_total),
            "audit_example_reservoir_seen_weight": float(round(float(audit_reservoir_seen_weight), 6)),
            "workers": int(worker_count),
        },
    }
    threshold_file = args.global_threshold_file or (output_dir / "global_thresholds.json")
    threshold_file.parent.mkdir(parents=True, exist_ok=True)
    threshold_file.write_text(json.dumps(artifact, indent=2, ensure_ascii=False), encoding="utf-8")
    (output_dir / "threshold_recommendations.json").write_text(json.dumps(artifact, indent=2, ensure_ascii=False), encoding="utf-8")
    write_tsv(output_dir / "stream_hist_1d_gb.tsv", histogram_rows_1d(gb_edges, hist_gb, "gb_strength"), ["gb_strength_low", "gb_strength_high", "gb_strength_center", "weighted_count"])
    write_tsv(output_dir / "stream_hist_1d_p2.tsv", histogram_rows_1d(p2_edges, hist_p2, "p2"), ["p2_low", "p2_high", "p2_center", "weighted_count"])
    write_tsv(output_dir / "stream_hist_2d_gb_p2.tsv", histogram_rows_2d(gb_edges, p2_edges, hist2d), ["gb_low", "gb_high", "gb_center", "p2_low", "p2_high", "p2_center", "weighted_count"])
    if audit_sample:
        write_tsv(output_dir / "audit_pair_examples.tsv", audit_sample, gb_pair_audit.PAIR_COLUMNS)
    np.savez_compressed(
        output_dir / "stream_histograms.npz",
        gb_edges=gb_edges,
        p2_edges=p2_edges,
        hist_gb=hist_gb,
        hist_p2=hist_p2,
        hist2d=hist2d,
        local_hist2d=local_hist2d,
        nonlocal_hist2d=nonlocal_hist2d,
        block_hist2d=np.asarray(block_hist2d, dtype=float) if block_hist2d else np.zeros((0, hist2d.shape[0], hist2d.shape[1]), dtype=float),
    )
    (output_dir / "streaming_manifest.json").write_text(json.dumps(artifact["streaming"] | artifact["input_provenance"], indent=2, ensure_ascii=False), encoding="utf-8")
    write_stream_plots(output_dir, gb_edges, p2_edges, hist2d, artifact)
    print(f"Global threshold prior: {threshold_file}")
    return artifact


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="LC Domain-Pearl V2 streaming global threshold prior builder.")
    parser.add_argument("inputs", nargs="*", type=Path, help="Dump files or directories. Defaults to current directory.")
    parser.add_argument("--pattern", default="*.dump")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--gb-param-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("lc_threshold_prior_output"))
    parser.add_argument("--global-threshold-file", type=Path, default=None)
    parser.add_argument("--axis", default="auto", choices=["auto", "x", "y", "z"])
    parser.add_argument("--mesogen-type", type=int, default=1)
    parser.add_argument("--anchor-types", default="2,3")
    parser.add_argument("--s-excl", type=int, default=1)
    parser.add_argument("--local-pair-file", type=Path, default=None)
    parser.add_argument("--exclude-pair-file", type=Path, default=None)
    parser.add_argument("--r-energy-cap", default="auto")
    parser.add_argument("--every", type=int, default=1)
    parser.add_argument("--global-frame-stride", type=int, default=1, help="Explicit global frame stride after --every. 1 means all candidate frames; 10 means every 10th global candidate frame; 100 means every 100th. Takes precedence over --global-frame-budget when > 1.")
    parser.add_argument("--global-frame-budget", type=int, default=0, help="0 means stream all frames; positive values use deterministic global stride.")
    parser.add_argument("--block-size-frames", type=int, default=100)
    parser.add_argument("--file-chunk-size", type=int, default=DEFAULT_FILE_CHUNK_SIZE, help="Number of dump files processed per worker task; prevents per-file task explosions for single-frame dump directories.")
    parser.add_argument("--max-block-histograms", type=int, default=DEFAULT_MAX_BLOCK_HISTOGRAMS, help="Maximum block histograms saved to stream_histograms.npz; block_count is still computed from the full stream.")
    parser.add_argument("--audit-example-pairs", type=int, default=5000)
    parser.add_argument("--gb-bins", type=int, default=DEFAULT_GB_BINS)
    parser.add_argument("--p2-bins", type=int, default=DEFAULT_P2_BINS)
    parser.add_argument("--gb-strength-max", type=float, default=1.5)
    parser.add_argument("--workers", default="auto", help="Parallel file-chunk workers. 'auto' uses min(CPU cores, chunk count, LC_PEARL_MAX_AUTO_WORKERS); default cap is 10.")
    parser.add_argument("--sample-seed", type=int, default=20260429)
    parser.add_argument("--current-gb-off", type=float, default=0.12)
    parser.add_argument("--current-gb-on", type=float, default=0.30)
    parser.add_argument("--current-p2-cut", type=float, default=0.70)
    parser.add_argument("--current-s2-cut", type=float, default=0.70)
    parser.add_argument("--n-min", type=int, default=3)
    parser.add_argument("--min-pairs", type=int, default=100)
    parser.add_argument("--min-oriented-pairs", type=int, default=60)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.every <= 0:
        raise SystemExit("--every must be positive")
    if args.global_frame_stride <= 0:
        raise SystemExit("--global-frame-stride must be positive")
    if args.global_frame_budget < 0:
        raise SystemExit("--global-frame-budget must be non-negative")
    if args.block_size_frames <= 0:
        raise SystemExit("--block-size-frames must be positive")
    if args.file_chunk_size <= 0:
        raise SystemExit("--file-chunk-size must be positive")
    if args.max_block_histograms < 0:
        raise SystemExit("--max-block-histograms must be non-negative")
    if args.gb_bins < 10 or args.p2_bins < 10:
        raise SystemExit("--gb-bins and --p2-bins must be at least 10")
    if args.gb_strength_max <= 0.0:
        raise SystemExit("--gb-strength-max must be positive")
    build_streaming_prior(args)


if __name__ == "__main__":
    main()
