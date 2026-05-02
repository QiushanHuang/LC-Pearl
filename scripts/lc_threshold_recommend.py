#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(row: Dict[str, str], key: str, default: float = math.nan) -> float:
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return default


def row_is_excluded(row: Dict[str, str]) -> bool:
    value = as_float(row, "is_excluded", 0.0)
    return math.isfinite(value) and int(value) == 1


def clamp(value: float, low: float, high: float) -> float:
    return max(low, min(high, float(value)))


def one_dimensional_kmeans(values: Sequence[float], k: int, max_iter: int = 80) -> Tuple[np.ndarray, np.ndarray]:
    array = np.asarray([float(v) for v in values if math.isfinite(float(v))], dtype=float)
    if array.size < k:
        raise ValueError(f"need at least {k} finite values")
    quantiles = np.linspace(0.0, 1.0, k + 2)[1:-1]
    centers = np.quantile(array, quantiles)
    centers = np.asarray(centers, dtype=float)
    labels = np.zeros(array.shape[0], dtype=int)
    for _ in range(max_iter):
        distances = np.abs(array[:, None] - centers[None, :])
        new_labels = np.argmin(distances, axis=1)
        new_centers = centers.copy()
        for idx in range(k):
            group = array[new_labels == idx]
            if group.size:
                new_centers[idx] = float(np.mean(group))
        order = np.argsort(new_centers)
        remap = {int(old): int(new) for new, old in enumerate(order)}
        sorted_centers = new_centers[order]
        sorted_labels = np.array([remap[int(label)] for label in new_labels], dtype=int)
        if np.array_equal(sorted_labels, labels) and np.allclose(sorted_centers, centers):
            centers = sorted_centers
            labels = sorted_labels
            break
        centers = sorted_centers
        labels = sorted_labels
    return centers, labels


def cluster_boundary(centers: Sequence[float], left_idx: int, right_idx: int) -> float:
    left = float(centers[left_idx])
    right = float(centers[right_idx])
    return 0.5 * (left + right)


def separation_confidence(centers: Sequence[float], values: Sequence[float]) -> str:
    array = np.asarray([float(v) for v in values if math.isfinite(float(v))], dtype=float)
    if array.size < 20:
        return "low"
    spread = float(np.std(array))
    if spread <= 1e-12:
        return "low"
    sep = float(np.max(centers) - np.min(centers)) / spread
    if sep >= 1.25:
        return "high"
    if sep >= 0.75:
        return "medium"
    return "low"


def recommend_from_rows(
    edge_rows: Sequence[Dict[str, str]],
    domain_rows: Sequence[Dict[str, str]],
    *,
    current_gb_off: float = 0.12,
    current_gb_on: float = 0.30,
    current_p2_cut: float = 0.70,
    current_s2_cut: float = 0.70,
    n_min: int = 3,
) -> Dict[str, object]:
    calibration_edge_rows = [row for row in edge_rows if not row_is_excluded(row)]
    strengths = np.array(
        [as_float(row, "attraction_strength") for row in calibration_edge_rows],
        dtype=float,
    )
    p2_values = np.array(
        [as_float(row, "p2_score") for row in calibration_edge_rows],
        dtype=float,
    )
    valid_pair_mask = np.isfinite(strengths) & np.isfinite(p2_values) & (strengths > 0.0)
    strengths = strengths[valid_pair_mask]
    p2_values = p2_values[valid_pair_mask]

    notes: List[str] = []
    if strengths.size < 50:
        notes.append("Pair sample is small; recommendations should be treated as preliminary.")

    try:
        p2_centers, _ = one_dimensional_kmeans(p2_values, k=2)
        recommended_p2 = clamp(cluster_boundary(p2_centers, 0, 1), 0.50, 0.90)
        p2_confidence = separation_confidence(p2_centers, p2_values)
    except ValueError:
        p2_centers = np.array([current_p2_cut, current_p2_cut], dtype=float)
        recommended_p2 = current_p2_cut
        p2_confidence = "low"
        notes.append("Could not infer P2 threshold from data; kept current p2_cut.")

    oriented_strengths = strengths[p2_values >= recommended_p2]
    if oriented_strengths.size < 30 and strengths.size:
        oriented_strengths = strengths[p2_values >= current_p2_cut]
    if oriented_strengths.size < 30 and strengths.size:
        oriented_strengths = strengths
        notes.append("Too few oriented pairs; strength thresholds use all attractive candidate pairs.")

    try:
        strength_centers, _ = one_dimensional_kmeans(oriented_strengths, k=3)
        recommended_gb_off = clamp(cluster_boundary(strength_centers, 0, 1), 0.02, 0.80)
        recommended_gb_on = clamp(cluster_boundary(strength_centers, 1, 2), recommended_gb_off + 0.02, 0.95)
        strength_confidence = separation_confidence(strength_centers, oriented_strengths)
    except ValueError:
        strength_centers = np.array([current_gb_off, current_gb_on], dtype=float)
        recommended_gb_off = current_gb_off
        recommended_gb_on = current_gb_on
        strength_confidence = "low"
        notes.append("Could not infer GB strength thresholds from data; kept current gb_off/gb_on.")

    domain_s2 = []
    domain_sizes = []
    for row in domain_rows:
        size = as_float(row, "size")
        s2 = as_float(row, "s2")
        if math.isfinite(size) and math.isfinite(s2):
            domain_sizes.append(size)
            domain_s2.append(s2)
    domain_s2_array = np.asarray(domain_s2, dtype=float)
    domain_sizes_array = np.asarray(domain_sizes, dtype=float)
    eligible_s2 = domain_s2_array[domain_sizes_array >= float(n_min)]
    try:
        s2_centers, _ = one_dimensional_kmeans(eligible_s2, k=2)
        recommended_s2 = clamp(cluster_boundary(s2_centers, 0, 1), 0.45, 0.90)
        s2_confidence = separation_confidence(s2_centers, eligible_s2)
    except ValueError:
        s2_centers = np.array([current_s2_cut, current_s2_cut], dtype=float)
        recommended_s2 = current_s2_cut
        s2_confidence = "low"
        notes.append("Could not infer robust_min_s2 from domain data; kept current robust_min_s2.")

    def rounded(value: float) -> float:
        return float(round(float(value), 4))

    current = {
        "gb_off_strength": current_gb_off,
        "gb_on_strength": current_gb_on,
        "p2_cut": current_p2_cut,
        "robust_min_s2": current_s2_cut,
    }
    recommended = {
        "gb_off_strength": rounded(recommended_gb_off),
        "gb_on_strength": rounded(recommended_gb_on),
        "p2_cut": rounded(recommended_p2),
        "robust_min_s2": rounded(recommended_s2),
    }
    deltas = {
        key: rounded(float(recommended[key]) - float(current[key]))
        for key in current
    }

    return {
        "method": {
            "pair_p2": "1D k-means with k=2 on pair P2 scores; boundary between low/high orientation clusters.",
            "gb_strength": "1D k-means with k=3 on oriented attractive-pair gb_strength; boundaries define gray and strong thresholds.",
            "domain_s2": "1D k-means with k=2 on domain S2 for domains with size >= n_min.",
        },
        "current": current,
        "recommended": recommended,
        "delta_recommended_minus_current": deltas,
        "confidence": {
            "p2_cut": p2_confidence,
            "gb_strength": strength_confidence,
            "robust_min_s2": s2_confidence,
        },
        "sample_sizes": {
            "input_pair_rows": int(len(edge_rows)),
            "excluded_pair_rows": int(len(edge_rows) - len(calibration_edge_rows)),
            "attractive_pair_rows": int(strengths.size),
            "oriented_pair_rows_used_for_strength": int(oriented_strengths.size),
            "domain_rows": int(domain_s2_array.size),
            "eligible_domain_rows_for_s2": int(eligible_s2.size),
        },
        "cluster_centers": {
            "p2": [rounded(v) for v in p2_centers.tolist()],
            "gb_strength": [rounded(v) for v in strength_centers.tolist()],
            "domain_s2": [rounded(v) for v in s2_centers.tolist()],
        },
        "notes": notes,
    }


def write_pair_bin_summary(edge_rows: Sequence[Dict[str, str]], output_path: Path, bins: int = 20) -> None:
    strengths = np.array([as_float(row, "attraction_strength") for row in edge_rows], dtype=float)
    p2_values = np.array([as_float(row, "p2_score") for row in edge_rows], dtype=float)
    is_local = np.array([int(as_float(row, "is_local", 0.0)) for row in edge_rows], dtype=int)
    mask = np.isfinite(strengths) & np.isfinite(p2_values)
    strengths = strengths[mask]
    p2_values = p2_values[mask]
    is_local = is_local[mask]
    strength_edges = np.linspace(0.0, max(1.0, float(np.max(strengths)) if strengths.size else 1.0), bins + 1)
    p2_edges = np.linspace(-0.5, 1.0, bins + 1)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("strength_low\tstrength_high\tp2_low\tp2_high\tcount\tlocal_count\tnonlocal_count\n")
        for si in range(bins):
            s_mask = (strengths >= strength_edges[si]) & (strengths < strength_edges[si + 1])
            if si == bins - 1:
                s_mask = (strengths >= strength_edges[si]) & (strengths <= strength_edges[si + 1])
            for pi in range(bins):
                p_mask = (p2_values >= p2_edges[pi]) & (p2_values < p2_edges[pi + 1])
                if pi == bins - 1:
                    p_mask = (p2_values >= p2_edges[pi]) & (p2_values <= p2_edges[pi + 1])
                selected = s_mask & p_mask
                count = int(np.sum(selected))
                if count == 0:
                    continue
                local_count = int(np.sum(is_local[selected] == 1))
                handle.write(
                    f"{strength_edges[si]:.8g}\t{strength_edges[si+1]:.8g}\t"
                    f"{p2_edges[pi]:.8g}\t{p2_edges[pi+1]:.8g}\t"
                    f"{count}\t{local_count}\t{count - local_count}\n"
                )


def threshold_label(name: str, value: object, digits: int = 3) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return f"{name}={value}"
    return f"{name}={numeric:.{digits}f}"


def annotate_vline(ax: object, value: float, label: str, *, color: str, ymax: float = 0.96, lw: float = 1.8) -> None:
    ax.axvline(float(value), color=color, lw=lw, label=label)
    ax.annotate(
        label,
        xy=(float(value), ymax),
        xycoords=("data", "axes fraction"),
        xytext=(4, -4),
        textcoords="offset points",
        rotation=90,
        va="top",
        ha="left",
        color=color,
        fontsize=8,
    )


def annotate_hline(ax: object, value: float, label: str, *, color: str, xmax: float = 0.98, lw: float = 1.8) -> None:
    ax.axhline(float(value), color=color, lw=lw, label=label)
    ax.annotate(
        label,
        xy=(xmax, float(value)),
        xycoords=("axes fraction", "data"),
        xytext=(-4, 4),
        textcoords="offset points",
        va="bottom",
        ha="right",
        color=color,
        fontsize=8,
    )


def write_plots(edge_rows: Sequence[Dict[str, str]], domain_rows: Sequence[Dict[str, str]], output_dir: Path, recommendation: Dict[str, object]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        (output_dir / "plot_error.txt").write_text(str(exc), encoding="utf-8")
        return

    recommended = recommendation["recommended"]  # type: ignore[index]
    strengths = np.array([as_float(row, "attraction_strength") for row in edge_rows], dtype=float)
    p2_values = np.array([as_float(row, "p2_score") for row in edge_rows], dtype=float)
    mask = np.isfinite(strengths) & np.isfinite(p2_values)
    strengths = strengths[mask]
    p2_values = p2_values[mask]
    if strengths.size:
        fig, ax = plt.subplots(figsize=(7.0, 5.2), dpi=180)
        hb = ax.hexbin(strengths, p2_values, gridsize=42, bins="log", cmap="viridis", mincnt=1)
        annotate_vline(ax, float(recommended["gb_off_strength"]), threshold_label("gb_off", recommended["gb_off_strength"]), color="#f59e0b")
        annotate_vline(ax, float(recommended["gb_on_strength"]), threshold_label("gb_on", recommended["gb_on_strength"]), color="#ef4444", ymax=0.82)
        annotate_hline(ax, float(recommended["p2_cut"]), threshold_label("p2_cut", recommended["p2_cut"]), color="#38bdf8")
        ax.set_xlabel("GB attraction strength")
        ax.set_ylabel("pair P2")
        ax.legend(frameon=False, loc="lower right")
        fig.colorbar(hb, ax=ax, label="log10(count)")
        fig.tight_layout()
        fig.savefig(output_dir / "gb_strength_vs_p2_hexbin_recommended.png")
        plt.close(fig)

    sizes = np.array([as_float(row, "size") for row in domain_rows], dtype=float)
    s2 = np.array([as_float(row, "s2") for row in domain_rows], dtype=float)
    mask = np.isfinite(sizes) & np.isfinite(s2)
    sizes = sizes[mask]
    s2 = s2[mask]
    if sizes.size:
        fig, ax = plt.subplots(figsize=(7.0, 5.2), dpi=180)
        jitter = np.linspace(-0.08, 0.08, sizes.size) if sizes.size > 1 else np.array([0.0])
        ax.scatter(sizes + jitter, s2, s=22, alpha=0.65, color="#475569")
        annotate_hline(ax, float(recommended["robust_min_s2"]), threshold_label("robust_min_s2", recommended["robust_min_s2"]), color="#ef4444")
        annotate_vline(ax, 3.0, "n_min=3", color="#38bdf8", lw=1.4)
        ax.set_xlabel("domain size")
        ax.set_ylabel("domain S2")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "domain_size_vs_s2_recommended.png")
        plt.close(fig)


def write_recommendation_report(
    diagnostics_dir: Path,
    output_dir: Optional[Path] = None,
    *,
    current_gb_off: float = 0.12,
    current_gb_on: float = 0.30,
    current_p2_cut: float = 0.70,
    current_s2_cut: float = 0.70,
    n_min: int = 3,
) -> Dict[str, object]:
    output = output_dir or diagnostics_dir / "threshold_auto"
    output.mkdir(parents=True, exist_ok=True)
    edge_path = diagnostics_dir / "edge_diagnostics.tsv"
    domain_path = diagnostics_dir / "domain_diagnostics.tsv"
    if not edge_path.exists():
        raise FileNotFoundError(edge_path)
    if not domain_path.exists():
        raise FileNotFoundError(domain_path)
    edge_rows = read_tsv(edge_path)
    domain_rows = read_tsv(domain_path)
    recommendation = recommend_from_rows(
        edge_rows,
        domain_rows,
        current_gb_off=current_gb_off,
        current_gb_on=current_gb_on,
        current_p2_cut=current_p2_cut,
        current_s2_cut=current_s2_cut,
        n_min=n_min,
    )
    (output / "threshold_recommendations.json").write_text(
        json.dumps(recommendation, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_pair_bin_summary(edge_rows, output / "pair_bin_summary.tsv")
    write_plots(edge_rows, domain_rows, output, recommendation)
    return recommendation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Recommend LC/GB analysis thresholds from diagnostics tables.")
    parser.add_argument("diagnostics_dir", type=Path, help="Directory containing edge_diagnostics.tsv and domain_diagnostics.tsv.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory. Defaults to diagnostics_dir/threshold_auto.")
    parser.add_argument("--current-gb-off", type=float, default=0.12)
    parser.add_argument("--current-gb-on", type=float, default=0.30)
    parser.add_argument("--current-p2-cut", type=float, default=0.70)
    parser.add_argument("--current-s2-cut", type=float, default=0.70)
    parser.add_argument("--n-min", type=int, default=3)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    recommendation = write_recommendation_report(
        args.diagnostics_dir,
        args.output_dir,
        current_gb_off=float(args.current_gb_off),
        current_gb_on=float(args.current_gb_on),
        current_p2_cut=float(args.current_p2_cut),
        current_s2_cut=float(args.current_s2_cut),
        n_min=int(args.n_min),
    )
    print(json.dumps(recommendation["recommended"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
