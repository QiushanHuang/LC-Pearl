#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from lc_threshold_recommend import as_float, write_pair_bin_summary, write_plots


SCHEMA_VERSION = 1
PARAMETER_KEYS = ("gb_off_strength", "gb_on_strength", "p2_cut", "robust_min_s2")
CALIBRATION_MODES = ("1d", "joint2d")
DEFAULT_MAX_AUTO_WORKERS = 10


def resolve_worker_count(value: object, task_count: int) -> int:
    if str(value).lower() == "auto":
        max_auto = int(os.environ.get("LC_PEARL_MAX_AUTO_CALIBRATION_WORKERS", DEFAULT_MAX_AUTO_WORKERS))
        return max(1, min(int(task_count), os.cpu_count() or 1, max_auto))
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise ValueError("--workers must be an integer or 'auto'") from None
    if parsed < 1:
        raise ValueError("--workers must be >= 1 or 'auto'")
    return max(1, min(parsed, max(int(task_count), 1)))


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Optional[Path]) -> Optional[str]:
    if path is None or not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_float(value: object, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def stable_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def clamp(value: float, low: float, high: float) -> float:
    return max(float(low), min(float(high), float(value)))


def rounded(value: float, digits: int = 4) -> float:
    if not math.isfinite(float(value)):
        return float("nan")
    return float(round(float(value), digits))


def finite_array(values: Iterable[object]) -> np.ndarray:
    arr = np.asarray([stable_float(value) for value in values], dtype=float)
    return arr[np.isfinite(arr)]


def filter_pair_rows(rows: Sequence[Dict[str, str]], *, include_excluded: bool = False) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    filtered: List[Dict[str, str]] = []
    excluded = 0
    invalid = 0
    non_attractive = 0
    for row in rows:
        if not include_excluded and stable_int(row.get("is_excluded", 0)) == 1:
            excluded += 1
            continue
        strength = stable_float(row.get("attraction_strength"))
        p2_score = stable_float(row.get("p2_score"))
        if not (math.isfinite(strength) and math.isfinite(p2_score)):
            invalid += 1
            continue
        if strength <= 0.0:
            non_attractive += 1
            continue
        filtered.append(row)
    return filtered, {
        "input_pair_rows": int(len(rows)),
        "excluded_pair_rows": int(excluded),
        "invalid_pair_rows": int(invalid),
        "non_attractive_pair_rows": int(non_attractive),
        "calibration_pair_rows": int(len(filtered)),
    }


def pair_arrays(rows: Sequence[Dict[str, str]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    strengths = np.asarray([stable_float(row.get("attraction_strength")) for row in rows], dtype=float)
    p2_values = np.asarray([stable_float(row.get("p2_score")) for row in rows], dtype=float)
    weights = np.asarray([sampling_weight(row) for row in rows], dtype=float)
    mask = np.isfinite(strengths) & np.isfinite(p2_values) & np.isfinite(weights) & (weights > 0.0) & (strengths > 0.0)
    return strengths[mask], p2_values[mask], weights[mask]


def sampling_weight(row: Dict[str, str]) -> float:
    value = stable_float(row.get("sampling_weight"), 1.0)
    if not math.isfinite(value) or value <= 0.0:
        return 1.0
    return float(value)


def kish_effective_size(weights: Sequence[float]) -> float:
    arr = np.asarray([float(value) for value in weights if math.isfinite(float(value)) and float(value) > 0.0], dtype=float)
    if arr.size == 0:
        return 0.0
    return float((np.sum(arr) ** 2) / max(float(np.sum(arr * arr)), 1e-12))


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: Sequence[float]) -> np.ndarray:
    if values.size == 0:
        return np.asarray([], dtype=float)
    order = np.argsort(values)
    x = np.asarray(values[order], dtype=float)
    w = np.asarray(weights[order], dtype=float)
    cumulative = np.cumsum(w)
    total = float(cumulative[-1])
    if total <= 0.0:
        return np.quantile(x, quantiles)
    return np.interp(np.asarray(quantiles, dtype=float) * total, cumulative, x)


def candidate_values(
    values: np.ndarray,
    weights: np.ndarray,
    *,
    anchors: Sequence[float],
    quantiles: Sequence[float],
    clamp_range: Tuple[float, float],
) -> List[float]:
    candidates: List[float] = []
    if values.size and weights.size == values.size:
        candidates.extend(float(value) for value in weighted_quantile(values, weights, quantiles))
    candidates.extend(float(value) for value in anchors)
    low, high = float(clamp_range[0]), float(clamp_range[1])
    cleaned = sorted({
        rounded(clamp(float(value), low, high), 4)
        for value in candidates
        if math.isfinite(float(value))
    })
    return [float(value) for value in cleaned if low <= float(value) <= high]


def joint2d_calibration(
    strengths: np.ndarray,
    p2_values: np.ndarray,
    weights: np.ndarray,
    *,
    current: Dict[str, float],
    recommended: Dict[str, float],
    min_pairs: int,
    min_oriented_pairs: int,
    independent_pair_blocks: int = 0,
    min_pair_blocks: int = 2,
) -> Dict[str, object]:
    mask = np.isfinite(strengths) & np.isfinite(p2_values) & np.isfinite(weights) & (weights > 0.0)
    x = np.asarray(strengths[mask], dtype=float)
    y = np.asarray(p2_values[mask], dtype=float)
    w = np.asarray(weights[mask], dtype=float)
    if x.size == 0:
        return {
            "status": "low",
            "reason": "No valid attractive pair rows for joint 2D calibration.",
            "grid_evaluated": 0,
            "independent_pair_blocks": int(independent_pair_blocks),
            "minimum_pair_blocks": int(min_pair_blocks),
            "selected": {},
            "best_unconstrained": {},
            "grid": [],
        }

    off_base = float(recommended.get("gb_off_strength", current["gb_off_strength"]))
    on_base = float(recommended.get("gb_on_strength", current["gb_on_strength"]))
    p2_base = float(recommended.get("p2_cut", current["p2_cut"]))
    off_candidates = candidate_values(
        x,
        w,
        anchors=[current["gb_off_strength"], off_base, off_base - 0.04, off_base + 0.04],
        quantiles=[0.30, 0.40, 0.50, 0.60],
        clamp_range=(0.02, 0.80),
    )
    on_candidates = candidate_values(
        x,
        w,
        anchors=[current["gb_on_strength"], on_base, on_base - 0.05, on_base + 0.05],
        quantiles=[0.60, 0.70, 0.80, 0.90],
        clamp_range=(0.04, 0.95),
    )
    p2_candidates = candidate_values(
        y,
        w,
        anchors=[current["p2_cut"], p2_base, p2_base - 0.04, p2_base + 0.04],
        quantiles=[0.45, 0.55, 0.65, 0.75],
        clamp_range=(0.50, 0.90),
    )

    rows: List[Dict[str, object]] = []
    total_eff = kish_effective_size(w)
    best_unconstrained: Optional[Dict[str, object]] = None
    best_feasible: Optional[Dict[str, object]] = None
    min_background = max(5.0, 0.5 * float(min_oriented_pairs))
    min_strong = max(5.0, 0.5 * float(min_oriented_pairs))
    for p2_cut in p2_candidates:
        oriented = y >= p2_cut
        oriented_eff = kish_effective_size(w[oriented])
        for gb_off in off_candidates:
            if gb_off >= 0.93:
                continue
            gray = oriented & (x >= gb_off)
            gray_eff = kish_effective_size(w[gray])
            for gb_on in on_candidates:
                if gb_on <= gb_off + 0.02:
                    continue
                strong = oriented & (x >= gb_on)
                strong_eff = kish_effective_size(w[strong])
                background = ~(gray)
                background_eff = kish_effective_size(w[background])
                if np.any(strong):
                    strong_strength = float(np.average(x[strong], weights=w[strong]))
                    strong_p2 = float(np.average(y[strong], weights=w[strong]))
                else:
                    strong_strength = 0.0
                    strong_p2 = 0.0
                if np.any(background):
                    background_strength = float(np.average(x[background], weights=w[background]))
                    background_p2 = float(np.average(y[background], weights=w[background]))
                else:
                    background_strength = 0.0
                    background_p2 = 0.0
                contrast = max(0.0, strong_strength - background_strength) + 0.5 * max(0.0, strong_p2 - background_p2)
                support_penalty = 0.0
                if gray_eff < min_oriented_pairs:
                    support_penalty += (float(min_oriented_pairs) - gray_eff) / max(float(min_oriented_pairs), 1.0)
                if strong_eff < max(5.0, 0.5 * float(min_oriented_pairs)):
                    support_penalty += (max(5.0, 0.5 * float(min_oriented_pairs)) - strong_eff) / max(float(min_oriented_pairs), 1.0)
                distance_penalty = (
                    abs(gb_off - off_base) / 0.20
                    + abs(gb_on - on_base) / 0.25
                    + abs(p2_cut - p2_base) / 0.20
                )
                score = (
                    math.sqrt(max(gray_eff, 0.0))
                    + math.sqrt(max(strong_eff, 0.0))
                    + 10.0 * contrast
                    - 5.0 * support_penalty
                    - 0.15 * distance_penalty
                )
                row = {
                    "gb_off_strength": rounded(gb_off),
                    "gb_on_strength": rounded(gb_on),
                    "p2_cut": rounded(p2_cut),
                    "oriented_effective_n": rounded(oriented_eff, 3),
                    "gray_effective_n": rounded(gray_eff, 3),
                    "strong_effective_n": rounded(strong_eff, 3),
                    "background_effective_n": rounded(background_eff, 3),
                    "contrast": rounded(contrast, 6),
                    "score": rounded(score, 6),
                    "feasible": bool(
                        total_eff >= min_pairs
                        and oriented_eff >= min_oriented_pairs
                        and gray_eff >= min_oriented_pairs
                        and strong_eff >= min_strong
                        and background_eff >= min_background
                        and independent_pair_blocks >= min_pair_blocks
                    ),
                }
                rows.append(row)
                if best_unconstrained is None or float(row["score"]) > float(best_unconstrained["score"]):
                    best_unconstrained = row
                if bool(row["feasible"]) and (best_feasible is None or float(row["score"]) > float(best_feasible["score"])):
                    best_feasible = row

    if best_unconstrained is None:
        return {
            "status": "low",
            "reason": "No valid joint 2D grid candidate satisfied gb_on > gb_off.",
            "grid_evaluated": 0,
            "independent_pair_blocks": int(independent_pair_blocks),
            "minimum_pair_blocks": int(min_pair_blocks),
            "selected": {},
            "best_unconstrained": {},
            "grid": [],
        }
    selected = best_feasible
    reference = best_feasible or best_unconstrained
    strong_eff = float(reference["strong_effective_n"])
    gray_eff = float(reference["gray_effective_n"])
    status = "high"
    reason = "Weighted 2D rectangular gate has joint gray/strong support."
    if best_feasible is None:
        status = "low"
        reason = "Joint 2D gate did not meet effective support requirements."
    elif gray_eff < 2.0 * float(min_oriented_pairs):
        status = "medium"
    return {
        "status": status,
        "reason": reason,
        "total_effective_n": rounded(total_eff, 3),
        "independent_pair_blocks": int(independent_pair_blocks),
        "minimum_pair_blocks": int(min_pair_blocks),
        "minimum_background_effective_n": rounded(min_background, 3),
        "minimum_strong_effective_n": rounded(min_strong, 3),
        "grid_evaluated": int(len(rows)),
        "grid_rows_stored": int(min(len(rows), 250)),
        "selected": selected or {},
        "best_unconstrained": best_unconstrained,
        "grid": sorted(rows, key=lambda row: float(row["score"]), reverse=True)[:250],
    }


def write_joint2d_grid(path: Path, joint: Dict[str, object]) -> None:
    grid = joint.get("grid", [])
    columns = [
        "gb_off_strength",
        "gb_on_strength",
        "p2_cut",
        "oriented_effective_n",
        "gray_effective_n",
        "strong_effective_n",
        "background_effective_n",
        "contrast",
        "score",
        "feasible",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        if not isinstance(grid, list):
            return
        for row in grid:
            if not isinstance(row, dict):
                continue
            handle.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")


def summarize_sampling(rows: Sequence[Dict[str, str]]) -> Dict[str, object]:
    strata: Dict[str, int] = {}
    reasons: Dict[str, int] = {}
    weighted_total = 0.0
    weighted_rows = 0
    non_unbiased_rows = 0
    for row in rows:
        stratum = str(row.get("sampling_strata", ""))
        reason = str(row.get("sampling_selected_reason", ""))
        if stratum:
            strata[stratum] = strata.get(stratum, 0) + 1
        if reason:
            reasons[reason] = reasons.get(reason, 0) + 1
        if "sampling_weight" in row:
            weighted_rows += 1
        weighted_total += sampling_weight(row)
        semantics = str(row.get("weight_semantics", ""))
        if "coverage_expansion_not_unbiased" in semantics or "global_candidate_cap_weight" in semantics:
            non_unbiased_rows += 1
    return {
        "selected_pair_rows": int(len(rows)),
        "weighted_pair_rows": float(round(weighted_total, 6)),
        "rows_with_sampling_weight": int(weighted_rows),
        "rows_with_non_unbiased_event_weight": int(non_unbiased_rows),
        "unique_sampling_strata": int(len(strata)),
        "top_sampling_strata": sorted(strata.items(), key=lambda item: (-item[1], item[0]))[:20],
        "selection_reasons": reasons,
    }


def domain_arrays(rows: Sequence[Dict[str, str]], n_min: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sizes = np.asarray([stable_float(row.get("size")) for row in rows], dtype=float)
    s2_values = np.asarray([stable_float(row.get("s2")) for row in rows], dtype=float)
    mask = np.isfinite(sizes) & np.isfinite(s2_values)
    sizes = sizes[mask]
    s2_values = s2_values[mask]
    eligible = s2_values[sizes >= float(n_min)]
    return sizes, s2_values, eligible


def logsumexp_2d(values: np.ndarray, axis: int = 1) -> np.ndarray:
    max_values = np.max(values, axis=axis, keepdims=True)
    return np.squeeze(max_values, axis=axis) + np.log(np.sum(np.exp(values - max_values), axis=axis))


def fit_gmm_1d(
    values: Sequence[float],
    k: int,
    *,
    sample_weights: Optional[Sequence[float]] = None,
    max_iter: int = 120,
    tol: float = 1e-7,
) -> Dict[str, object]:
    raw_x = np.asarray([float(v) for v in values], dtype=float)
    raw_w = np.ones(raw_x.size, dtype=float) if sample_weights is None else np.asarray([float(v) for v in sample_weights], dtype=float)
    if raw_w.size != raw_x.size:
        raise ValueError("sample_weights length must match values length")
    mask = np.isfinite(raw_x) & np.isfinite(raw_w) & (raw_w > 0.0)
    x = raw_x[mask]
    data_weights = raw_w[mask]
    if x.size < k:
        raise ValueError(f"need at least {k} finite values")
    if float(np.max(x) - np.min(x)) <= 1e-12:
        raise ValueError("values are effectively constant")
    data_weight_sum = max(float(np.sum(data_weights)), 1e-12)
    means = np.asarray(weighted_quantile(x, data_weights, np.linspace(0.0, 1.0, k + 2)[1:-1]), dtype=float)
    weighted_mean = float(np.sum(data_weights * x) / data_weight_sum)
    weighted_var = float(np.sum(data_weights * (x - weighted_mean) ** 2) / data_weight_sum)
    variance_floor = max(weighted_var * 1e-6, 1e-10)
    variances = np.full(k, max(weighted_var, variance_floor), dtype=float)
    component_weights = np.full(k, 1.0 / float(k), dtype=float)
    previous_ll = -float("inf")
    labels = np.zeros(x.size, dtype=int)

    for _iter in range(max_iter):
        log_prob = np.empty((x.size, k), dtype=float)
        for component in range(k):
            var = max(float(variances[component]), variance_floor)
            log_prob[:, component] = (
                math.log(max(float(component_weights[component]), 1e-15))
                - 0.5 * math.log(2.0 * math.pi * var)
                - ((x - float(means[component])) ** 2) / (2.0 * var)
            )
        log_norm = logsumexp_2d(log_prob, axis=1)
        ll = float(np.sum(data_weights * log_norm))
        resp = np.exp(log_prob - log_norm[:, None])
        weighted_resp = resp * data_weights[:, None]
        component_weights = np.sum(weighted_resp, axis=0) / data_weight_sum
        for component in range(k):
            total = float(np.sum(weighted_resp[:, component]))
            if total <= 1e-12:
                means[component] = float(weighted_quantile(x, data_weights, [(component + 1.0) / (k + 1.0)])[0])
                variances[component] = max(weighted_var, variance_floor)
                component_weights[component] = 1.0 / float(k)
                continue
            means[component] = float(np.sum(weighted_resp[:, component] * x) / total)
            variances[component] = max(float(np.sum(weighted_resp[:, component] * (x - means[component]) ** 2) / total), variance_floor)
        if abs(ll - previous_ll) <= tol * max(1.0, abs(previous_ll)):
            break
        previous_ll = ll

    order = np.argsort(means)
    means = means[order]
    variances = variances[order]
    component_weights = component_weights[order]

    log_prob = np.empty((x.size, k), dtype=float)
    for component in range(k):
        var = max(float(variances[component]), variance_floor)
        log_prob[:, component] = (
            math.log(max(float(component_weights[component]), 1e-15))
            - 0.5 * math.log(2.0 * math.pi * var)
            - ((x - float(means[component])) ** 2) / (2.0 * var)
        )
    log_norm = logsumexp_2d(log_prob, axis=1)
    ll = float(np.sum(data_weights * log_norm))
    labels = np.argmax(log_prob - log_norm[:, None], axis=1)
    n_params = (k - 1) + k + k
    n_eff = kish_effective_size(data_weights)
    bic = float(n_params * math.log(max(float(n_eff), 1.0)) - 2.0 * ll)
    return {
        "k": int(k),
        "weights": component_weights,
        "means": means,
        "variances": variances,
        "log_likelihood": ll,
        "bic": bic,
        "labels": labels,
        "n": int(x.size),
        "n_eff": float(n_eff),
    }


def mixture_log_prob(grid: np.ndarray, model: Dict[str, object]) -> np.ndarray:
    means = np.asarray(model["means"], dtype=float)
    variances = np.asarray(model["variances"], dtype=float)
    weights = np.asarray(model["weights"], dtype=float)
    log_prob = np.empty((grid.size, means.size), dtype=float)
    for component in range(means.size):
        var = max(float(variances[component]), 1e-12)
        log_prob[:, component] = (
            math.log(max(float(weights[component]), 1e-15))
            - 0.5 * math.log(2.0 * math.pi * var)
            - ((grid - float(means[component])) ** 2) / (2.0 * var)
        )
    return log_prob


def posterior_boundary(model: Dict[str, object], left_idx: int, right_idx: int) -> float:
    means = np.asarray(model["means"], dtype=float)
    left = float(means[left_idx])
    right = float(means[right_idx])
    if right <= left:
        return 0.5 * (left + right)
    grid = np.linspace(left, right, 512)
    log_prob = mixture_log_prob(grid, model)
    log_norm = logsumexp_2d(log_prob, axis=1)
    post = np.exp(log_prob - log_norm[:, None])
    diff = np.abs(post[:, left_idx] - post[:, right_idx])
    return float(grid[int(np.argmin(diff))])


def kde_bandwidth(values: np.ndarray, weights: Optional[np.ndarray] = None) -> float:
    x = np.asarray(values, dtype=float)
    if x.size < 2:
        return 1.0
    w = np.ones(x.size, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w) & (w > 0.0)
    x = x[mask]
    w = w[mask]
    if x.size < 2:
        return 1.0
    w_sum = max(float(np.sum(w)), 1e-12)
    mean = float(np.sum(w * x) / w_sum)
    n_eff = max(kish_effective_size(w), 1.0)
    std = math.sqrt(max(float(np.sum(w * (x - mean) ** 2) / w_sum), 0.0))
    q75, q25 = weighted_quantile(x, w, [0.75, 0.25])
    robust = min(std, float((q75 - q25) / 1.349)) if std > 0.0 else float((q75 - q25) / 1.349)
    if robust <= 1e-12 or not math.isfinite(robust):
        robust = max(std, float(np.max(x) - np.min(x)) / 6.0, 1e-3)
    return max(0.9 * robust * (float(n_eff) ** (-0.2)), 1e-4)


def kde_valley(values: Sequence[float], left: float, right: float, *, sample_weights: Optional[Sequence[float]] = None) -> Optional[float]:
    raw_x = np.asarray([float(v) for v in values], dtype=float)
    raw_w = np.ones(raw_x.size, dtype=float) if sample_weights is None else np.asarray([float(v) for v in sample_weights], dtype=float)
    if raw_w.size != raw_x.size:
        return None
    mask = np.isfinite(raw_x) & np.isfinite(raw_w) & (raw_w > 0.0)
    x = raw_x[mask]
    w = raw_w[mask]
    if x.size < 5 or right <= left:
        return None
    bandwidth = kde_bandwidth(x, w)
    grid = np.linspace(left, right, 512)
    scaled = (grid[:, None] - x[None, :]) / bandwidth
    density = np.sum(w[None, :] * np.exp(-0.5 * scaled * scaled), axis=1)
    density = density / (max(float(np.sum(w)), 1e-12) * bandwidth * math.sqrt(2.0 * math.pi))
    return float(grid[int(np.argmin(density))])


def estimate_boundary(
    values: Sequence[float],
    *,
    k: int,
    left_idx: int,
    right_idx: int,
    clamp_range: Tuple[float, float],
    current: float,
    min_samples: int,
    sample_weights: Optional[Sequence[float]] = None,
) -> Dict[str, object]:
    raw_x = np.asarray([float(v) for v in values], dtype=float)
    raw_w = np.ones(raw_x.size, dtype=float) if sample_weights is None else np.asarray([float(v) for v in sample_weights], dtype=float)
    if raw_w.size != raw_x.size:
        raw_w = np.ones(raw_x.size, dtype=float)
    mask = np.isfinite(raw_x) & np.isfinite(raw_w) & (raw_w > 0.0)
    x = raw_x[mask]
    weights = raw_w[mask]
    n_eff = kish_effective_size(weights)
    notes: List[str] = []
    if n_eff < min_samples:
        return {
            "estimate": float(current),
            "decision": "keep_current",
            "confidence": "low",
            "fallback_used": True,
            "reason": f"effective sample size {rounded(n_eff, 3)} < {min_samples}",
            "model": {},
            "kde_valley": None,
            "gmm_boundary": None,
            "notes": [f"Kept current value because effective sample size was below {min_samples}."],
        }
    try:
        model = fit_gmm_1d(x, k=k, sample_weights=weights)
        gmm_boundary = posterior_boundary(model, left_idx, right_idx)
    except ValueError as exc:
        return {
            "estimate": float(current),
            "decision": "keep_current",
            "confidence": "low",
            "fallback_used": True,
            "reason": str(exc),
            "model": {},
            "kde_valley": None,
            "gmm_boundary": None,
            "notes": ["Kept current value because GMM fitting failed."],
        }

    means = np.asarray(model["means"], dtype=float)
    valley = kde_valley(x, float(means[left_idx]), float(means[right_idx]), sample_weights=weights)
    candidates = [float(gmm_boundary)]
    if valley is not None and math.isfinite(float(valley)):
        candidates.append(float(valley))
    raw_estimate = float(np.median(np.asarray(candidates, dtype=float)))
    estimate = clamp(raw_estimate, clamp_range[0], clamp_range[1])
    value_range = max(float(np.max(x) - np.min(x)), 1e-9)
    gmm_kde_gap = abs(float(gmm_boundary) - float(valley)) if valley is not None else 0.0
    if valley is None:
        notes.append("KDE valley unavailable; used GMM posterior boundary only.")
    elif gmm_kde_gap > 0.18 * value_range:
        notes.append("GMM boundary and KDE valley differ noticeably; confidence is downgraded.")

    component_gap = abs(float(means[right_idx]) - float(means[left_idx]))
    pooled_std = math.sqrt(0.5 * (float(model["variances"][left_idx]) + float(model["variances"][right_idx])))
    separation = component_gap / max(pooled_std, 1e-9)
    confidence = "high"
    if n_eff < max(200, min_samples * 2) or separation < 1.8 or gmm_kde_gap > 0.18 * value_range:
        confidence = "medium"
    if n_eff < min_samples or separation < 0.9 or estimate in {float(clamp_range[0]), float(clamp_range[1])}:
        confidence = "low"

    model_summary = {
        "k": int(model["k"]),
        "n": int(model["n"]),
        "n_eff": rounded(float(model.get("n_eff", n_eff)), 3),
        "bic": rounded(float(model["bic"]), 6),
        "weights": [rounded(v, 6) for v in np.asarray(model["weights"], dtype=float).tolist()],
        "means": [rounded(v, 6) for v in means.tolist()],
        "variances": [rounded(v, 8) for v in np.asarray(model["variances"], dtype=float).tolist()],
        "component_separation": rounded(separation, 6),
    }
    return {
        "estimate": float(estimate),
        "decision": "apply" if confidence != "low" else "keep_current",
        "confidence": confidence,
        "fallback_used": False,
        "reason": "GMM posterior boundary cross-checked by KDE valley.",
        "model": model_summary,
        "kde_valley": rounded(valley, 6) if valley is not None else None,
        "gmm_boundary": rounded(gmm_boundary, 6),
        "notes": notes,
    }


def make_blocks(rows: Sequence[Dict[str, str]]) -> List[List[Dict[str, str]]]:
    blocks: Dict[str, List[Dict[str, str]]] = {}
    for idx, row in enumerate(rows):
        source = row.get("source_file", "")
        timestep = row.get("timestep", "")
        key = f"{source}|{timestep}" if source or timestep else f"rowblock:{idx // 250}"
        blocks.setdefault(key, []).append(row)
    return list(blocks.values())


def resample_blocks(blocks: Sequence[Sequence[Dict[str, str]]], rng: np.random.Generator) -> List[Dict[str, str]]:
    if not blocks:
        return []
    sampled: List[Dict[str, str]] = []
    for block_idx in rng.integers(0, len(blocks), size=len(blocks)):
        sampled.extend(blocks[int(block_idx)])
    return sampled


BootstrapDraw = Tuple[Tuple[int, ...], Tuple[int, ...]]


def resample_blocks_by_indices(blocks: Sequence[Sequence[Dict[str, str]]], indices: Sequence[int]) -> List[Dict[str, str]]:
    sampled: List[Dict[str, str]] = []
    for block_idx in indices:
        sampled.extend(blocks[int(block_idx)])
    return sampled


def bootstrap_draws(pair_block_count: int, domain_block_count: int, *, seed: int, samples: int) -> List[BootstrapDraw]:
    rng = np.random.default_rng(int(seed))
    draws: List[BootstrapDraw] = []
    for _ in range(int(samples)):
        pair_indices = tuple(int(value) for value in rng.integers(0, pair_block_count, size=pair_block_count).tolist()) if pair_block_count else tuple()
        domain_indices = tuple(int(value) for value in rng.integers(0, domain_block_count, size=domain_block_count).tolist()) if domain_block_count else tuple()
        draws.append((pair_indices, domain_indices))
    return draws


def bootstrap_draw_chunk(
    pair_blocks: Sequence[Sequence[Dict[str, str]]],
    domain_blocks: Sequence[Sequence[Dict[str, str]]],
    estimate_kwargs: Dict[str, object],
    draws: Sequence[BootstrapDraw],
) -> Tuple[Dict[str, List[float]], int]:
    values: Dict[str, List[float]] = {key: [] for key in PARAMETER_KEYS}
    failed = 0
    for pair_indices, domain_indices in draws:
        sampled_pairs = resample_blocks_by_indices(pair_blocks, pair_indices)
        sampled_domains = resample_blocks_by_indices(domain_blocks, domain_indices) if domain_blocks else []
        try:
            result = estimate_once(sampled_pairs, sampled_domains, **estimate_kwargs)
        except Exception:
            failed += 1
            continue
        recommended = result["recommended"]
        for key in PARAMETER_KEYS:
            value = stable_float(recommended.get(key))
            if math.isfinite(value):
                values[key].append(value)
    return values, failed


def _bootstrap_draw_chunk_star(args: Tuple[Sequence[Sequence[Dict[str, str]]], Sequence[Sequence[Dict[str, str]]], Dict[str, object], Sequence[BootstrapDraw]]) -> Tuple[Dict[str, List[float]], int]:
    return bootstrap_draw_chunk(*args)


def merge_bootstrap_chunks(chunks: Sequence[Tuple[Dict[str, List[float]], int]]) -> Tuple[Dict[str, List[float]], int]:
    values: Dict[str, List[float]] = {key: [] for key in PARAMETER_KEYS}
    failed = 0
    for chunk_values, chunk_failed in chunks:
        failed += int(chunk_failed)
        for key in PARAMETER_KEYS:
            values[key].extend(chunk_values.get(key, []))
    return values, failed


def estimate_once(
    pair_rows: Sequence[Dict[str, str]],
    domain_rows: Sequence[Dict[str, str]],
    *,
    current_gb_off: float,
    current_gb_on: float,
    current_p2_cut: float,
    current_s2_cut: float,
    n_min: int,
    min_pairs: int,
    min_oriented_pairs: int,
    min_domains: int,
    calibration_mode: str = "1d",
) -> Dict[str, object]:
    strengths, p2_values, pair_weights = pair_arrays(pair_rows)
    sizes, domain_s2, eligible_s2 = domain_arrays(domain_rows, n_min)
    notes: List[str] = []
    current = {
        "gb_off_strength": float(current_gb_off),
        "gb_on_strength": float(current_gb_on),
        "p2_cut": float(current_p2_cut),
        "robust_min_s2": float(current_s2_cut),
    }
    p2_detail = estimate_boundary(
        p2_values,
        k=2,
        left_idx=0,
        right_idx=1,
        clamp_range=(0.50, 0.90),
        current=current_p2_cut,
        min_samples=min_pairs,
        sample_weights=pair_weights,
    )
    p2_cut = float(p2_detail["estimate"])
    oriented_mask = p2_values >= p2_cut
    oriented_strengths = strengths[oriented_mask]
    oriented_weights = pair_weights[oriented_mask]
    if kish_effective_size(oriented_weights) < min_oriented_pairs and strengths.size:
        oriented_mask = p2_values >= float(current_p2_cut)
        oriented_strengths = strengths[oriented_mask]
        oriented_weights = pair_weights[oriented_mask]
        notes.append("Too few pairs above recommended p2_cut; strength thresholds used current p2_cut as orientation gate.")
    if kish_effective_size(oriented_weights) < min_oriented_pairs and strengths.size:
        oriented_strengths = strengths
        oriented_weights = pair_weights
        notes.append("Too few oriented pairs; strength thresholds used all attractive candidate pairs.")

    gb_off_detail = estimate_boundary(
        oriented_strengths,
        k=3,
        left_idx=0,
        right_idx=1,
        clamp_range=(0.02, 0.80),
        current=current_gb_off,
        min_samples=min_oriented_pairs,
        sample_weights=oriented_weights,
    )
    off_value = float(gb_off_detail["estimate"])
    gb_on_detail = estimate_boundary(
        oriented_strengths,
        k=3,
        left_idx=1,
        right_idx=2,
        clamp_range=(off_value + 0.02, 0.95),
        current=max(current_gb_on, off_value + 0.02),
        min_samples=min_oriented_pairs,
        sample_weights=oriented_weights,
    )

    s2_detail = estimate_boundary(
        eligible_s2,
        k=2,
        left_idx=0,
        right_idx=1,
        clamp_range=(0.45, 0.90),
        current=current_s2_cut,
        min_samples=min_domains,
    )
    if not domain_rows:
        s2_detail["notes"] = list(s2_detail.get("notes", [])) + [
            "No domain diagnostics were provided; robust_min_s2 is not changed during pre-analysis calibration."
        ]

    parameter_details = {
        "p2_cut": p2_detail,
        "gb_off_strength": gb_off_detail,
        "gb_on_strength": gb_on_detail,
        "robust_min_s2": s2_detail,
    }
    recommended = {
        "gb_off_strength": rounded(float(gb_off_detail["estimate"])),
        "gb_on_strength": rounded(float(gb_on_detail["estimate"])),
        "p2_cut": rounded(p2_cut),
        "robust_min_s2": rounded(float(s2_detail["estimate"])),
    }
    joint_2d: Optional[Dict[str, object]] = None
    if calibration_mode == "joint2d":
        joint_2d = joint2d_calibration(
            strengths,
            p2_values,
            pair_weights,
            current=current,
            recommended={key: float(value) for key, value in recommended.items()},
            min_pairs=int(min_pairs),
            min_oriented_pairs=int(min_oriented_pairs),
            independent_pair_blocks=len(make_blocks(pair_rows)),
        )
        if joint_2d.get("status") in {"high", "medium"} and isinstance(joint_2d.get("selected"), dict):
            selected = joint_2d["selected"]  # type: ignore[index]
            for key in ("gb_off_strength", "gb_on_strength", "p2_cut"):
                selected_value = rounded(stable_float(selected.get(key), float(recommended[key])))  # type: ignore[union-attr]
                recommended[key] = selected_value
                detail = parameter_details.get(key, {})
                if isinstance(detail, dict):
                    detail["estimate"] = float(selected_value)
                    detail["joint2d_estimate"] = float(selected_value)
                    detail["confidence"] = str(joint_2d.get("status"))
                    detail["decision"] = "apply"
                    detail["reason"] = "Weighted joint 2D gate cross-check selected this scalar threshold."
                    detail["notes"] = list(detail.get("notes", [])) + [
                        f"Joint 2D calibration status={joint_2d.get('status')}; selected weighted rectangular gate."
                    ]
        else:
            for key in ("gb_off_strength", "gb_on_strength", "p2_cut"):
                recommended[key] = rounded(float(current[key]))
                detail = parameter_details.get(key, {})
                if isinstance(detail, dict):
                    detail["confidence"] = "low"
                    detail["decision"] = "keep_current"
                    detail["notes"] = list(detail.get("notes", [])) + [
                        f"Joint 2D calibration failed support gate: {joint_2d.get('reason') if joint_2d else 'unavailable'}"
                    ]
    if recommended["gb_on_strength"] <= recommended["gb_off_strength"]:
        recommended["gb_on_strength"] = rounded(float(recommended["gb_off_strength"]) + 0.02)
        parameter_details["gb_on_strength"]["confidence"] = "low"
        parameter_details["gb_on_strength"]["decision"] = "keep_current"
        parameter_details["gb_on_strength"]["notes"] = list(parameter_details["gb_on_strength"].get("notes", [])) + [
            "Adjusted gb_on to preserve gb_on > gb_off; do not auto-apply this estimate."
        ]

    return {
        "current": current,
        "recommended": recommended,
        "parameters": parameter_details,
        "joint_2d": joint_2d,
        "sample_sizes": {
            "attractive_pair_rows": int(strengths.size),
            "attractive_pair_effective_n": rounded(kish_effective_size(pair_weights), 3),
            "oriented_pair_rows_used_for_strength": int(oriented_strengths.size),
            "oriented_pair_effective_n_used_for_strength": rounded(kish_effective_size(oriented_weights), 3),
            "domain_rows": int(domain_s2.size),
            "eligible_domain_rows_for_s2": int(eligible_s2.size),
        },
        "notes": notes,
    }


def bootstrap_estimates(
    pair_rows: Sequence[Dict[str, str]],
    domain_rows: Sequence[Dict[str, str]],
    *,
    seed: int,
    samples: int,
    estimate_kwargs: Dict[str, object],
    workers: object = 1,
) -> Dict[str, object]:
    if samples <= 0:
        return {"samples": 0, "workers": 0, "parameter_ci": {}, "failed_fit_count": 0}
    pair_blocks = make_blocks(pair_rows)
    domain_blocks = make_blocks(domain_rows)
    resolved_workers = resolve_worker_count(workers, int(samples))
    draws = bootstrap_draws(len(pair_blocks), len(domain_blocks), seed=int(seed), samples=int(samples))
    if resolved_workers <= 1 or int(samples) <= 1:
        values, failed = bootstrap_draw_chunk(pair_blocks, domain_blocks, estimate_kwargs, draws)
    else:
        chunks: List[List[BootstrapDraw]] = [[] for _ in range(resolved_workers)]
        for idx, draw in enumerate(draws):
            chunks[idx % resolved_workers].append(draw)
        tasks = [
            (pair_blocks, domain_blocks, estimate_kwargs, chunk)
            for chunk in chunks
            if chunk
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=resolved_workers) as executor:
            chunk_results = list(executor.map(_bootstrap_draw_chunk_star, tasks))
        values, failed = merge_bootstrap_chunks(chunk_results)

    parameter_ci: Dict[str, object] = {}
    for key, key_values in values.items():
        arr = np.asarray(key_values, dtype=float)
        if arr.size == 0:
            parameter_ci[key] = {"n": 0, "median": None, "ci95": None, "iqr": None}
            continue
        q025, q25, q50, q75, q975 = np.quantile(arr, [0.025, 0.25, 0.50, 0.75, 0.975])
        parameter_ci[key] = {
            "n": int(arr.size),
            "median": rounded(q50),
            "ci95": [rounded(q025), rounded(q975)],
            "iqr": rounded(float(q75 - q25)),
        }
    return {
        "samples": int(samples),
        "workers": int(resolved_workers),
        "seed": int(seed),
        "block_count": {
            "pair_blocks": int(len(pair_blocks)),
            "domain_blocks": int(len(domain_blocks)),
        },
        "parameter_ci": parameter_ci,
        "failed_fit_count": int(failed),
    }


def downgrade_confidence(value: str) -> str:
    if value == "high":
        return "medium"
    if value == "medium":
        return "low"
    return "low"


def apply_bootstrap_confidence(result: Dict[str, object], bootstrap: Dict[str, object]) -> None:
    parameter_ci = bootstrap.get("parameter_ci", {})
    parameters = result.get("parameters", {})
    for key, detail in parameters.items():
        if not isinstance(detail, dict):
            continue
        ci_info = parameter_ci.get(key, {}) if isinstance(parameter_ci, dict) else {}
        ci = ci_info.get("ci95") if isinstance(ci_info, dict) else None
        if not ci or len(ci) != 2:
            detail["confidence"] = downgrade_confidence(str(detail.get("confidence", "low")))
            detail["decision"] = "keep_current"
            detail["notes"] = list(detail.get("notes", [])) + ["Bootstrap CI unavailable; automatic application disabled for this parameter."]
            continue
        width = float(ci[1]) - float(ci[0])
        scale = 1.0 if key != "gb_off_strength" and key != "gb_on_strength" else 0.95
        if width > 0.28 * scale:
            detail["confidence"] = downgrade_confidence(str(detail.get("confidence", "low")))
            detail["notes"] = list(detail.get("notes", [])) + ["Bootstrap CI is wide; confidence downgraded."]
        if str(detail.get("confidence")) == "low":
            detail["decision"] = "keep_current"


def overall_status(parameters: Dict[str, object]) -> Tuple[str, bool]:
    required = ["p2_cut", "gb_off_strength", "gb_on_strength"]
    confidences = [str(parameters.get(key, {}).get("confidence", "low")) for key in required if isinstance(parameters.get(key), dict)]
    decisions = [str(parameters.get(key, {}).get("decision", "keep_current")) for key in required if isinstance(parameters.get(key), dict)]
    if not confidences or any(item == "low" for item in confidences):
        return "low", False
    if len(decisions) != len(required) or any(item != "apply" for item in decisions):
        return "low", False
    if any(item == "medium" for item in confidences):
        return "medium", True
    return "high", True


def write_bootstrap_table(path: Path, bootstrap: Dict[str, object]) -> None:
    parameter_ci = bootstrap.get("parameter_ci", {})
    with path.open("w", encoding="utf-8") as handle:
        handle.write("parameter\tn\tmedian\tci95_low\tci95_high\tiqr\n")
        if not isinstance(parameter_ci, dict):
            return
        for key in PARAMETER_KEYS:
            info = parameter_ci.get(key, {})
            if not isinstance(info, dict) or not info.get("ci95"):
                handle.write(f"{key}\t0\t\t\t\t\n")
                continue
            ci = info["ci95"]
            handle.write(f"{key}\t{info.get('n', 0)}\t{info.get('median')}\t{ci[0]}\t{ci[1]}\t{info.get('iqr')}\n")


def calibrate_thresholds(
    *,
    candidate_pairs: Path,
    domain_diagnostics: Optional[Path],
    output_dir: Path,
    current_gb_off: float,
    current_gb_on: float,
    current_p2_cut: float,
    current_s2_cut: float,
    n_min: int,
    min_pairs: int,
    min_oriented_pairs: int,
    min_domains: int,
    bootstrap_samples: int,
    seed: int,
    workers: object = 1,
    include_excluded: bool = False,
    input_kind: str = "candidate_pairs",
    potential_validation_status: str = "unknown",
    sampling_manifest: Optional[Path] = None,
    use_sampling_weights: bool = True,
    calibration_mode: str = "joint2d",
) -> Dict[str, object]:
    if calibration_mode not in CALIBRATION_MODES:
        raise ValueError(f"calibration_mode must be one of {', '.join(CALIBRATION_MODES)}")
    raw_pair_rows = read_tsv(candidate_pairs)
    domain_rows = read_tsv(domain_diagnostics)
    pair_rows, filter_summary = filter_pair_rows(raw_pair_rows, include_excluded=include_excluded)
    sampling_summary = summarize_sampling(pair_rows)
    fit_pair_rows = list(pair_rows)
    if not use_sampling_weights:
        fit_pair_rows = [dict(row, sampling_weight="1.0") for row in fit_pair_rows]
    estimate_kwargs: Dict[str, object] = {
        "current_gb_off": float(current_gb_off),
        "current_gb_on": float(current_gb_on),
        "current_p2_cut": float(current_p2_cut),
        "current_s2_cut": float(current_s2_cut),
        "n_min": int(n_min),
        "min_pairs": int(min_pairs),
        "min_oriented_pairs": int(min_oriented_pairs),
        "min_domains": int(min_domains),
        "calibration_mode": str(calibration_mode),
    }
    result = estimate_once(fit_pair_rows, domain_rows, **estimate_kwargs)
    raw_recommended = dict(result["recommended"])  # type: ignore[arg-type]
    unweighted_check: Optional[Dict[str, object]] = None
    if use_sampling_weights and sampling_summary.get("rows_with_sampling_weight", 0):
        unweighted_rows = [dict(row, sampling_weight="1.0") for row in fit_pair_rows]
        try:
            unweighted_result = estimate_once(unweighted_rows, domain_rows, **estimate_kwargs)
            weighted_recommended = result.get("recommended", {})
            unweighted_recommended = unweighted_result.get("recommended", {})
            deltas = {
                key: rounded(
                    stable_float(weighted_recommended.get(key), 0.0)
                    - stable_float(unweighted_recommended.get(key), 0.0)
                )
                for key in PARAMETER_KEYS
            }
            unweighted_check = {
                "unweighted_recommended": unweighted_recommended,
                "weighted_minus_unweighted": deltas,
            }
        except Exception as exc:
            unweighted_check = {"error": str(exc)}
    bootstrap = bootstrap_estimates(
        fit_pair_rows,
        domain_rows,
        seed=int(seed),
        samples=int(bootstrap_samples),
        estimate_kwargs=estimate_kwargs,
        workers=workers,
    )
    apply_bootstrap_confidence(result, bootstrap)
    status, apply_allowed = overall_status(result["parameters"])  # type: ignore[arg-type]
    warnings: List[str] = []
    if input_kind != "candidate_pairs":
        warnings.append("Input table is not an unscreened candidate-pair table; use this as a diagnostic report, not as a primary threshold source.")
    if potential_validation_status not in {"validated", "off", "not_required"}:
        warnings.append("GB potential is not marked validated in this artifact; automatic thresholds are statistical recommendations, not proof of potential reconstruction accuracy.")
        if status == "high":
            status = "medium"
        apply_allowed = False
        for detail in result["parameters"].values():  # type: ignore[union-attr]
            if isinstance(detail, dict):
                detail["decision"] = "keep_current"
                detail["notes"] = list(detail.get("notes", [])) + [
                    "Potential validation status is not validated/off/not_required; automatic application is disabled."
                ]
    if not apply_allowed:
        warnings.append("At least one required parameter has low confidence; pipeline should keep current values unless explicitly forced.")
    if int(sampling_summary.get("rows_with_non_unbiased_event_weight", 0)) > 0:
        warnings.append("Event-aware rows are event-enriched coverage samples, not an unbiased random sample; sampling weights are audit weights, not strict inclusion-probability weights.")
        if status == "high":
            status = "medium"
        apply_allowed = False
        for detail in result["parameters"].values():  # type: ignore[union-attr]
            if isinstance(detail, dict):
                detail["decision"] = "keep_current"
                detail["notes"] = list(detail.get("notes", [])) + [
                    "Event-aware calibration rows are coverage-enriched rather than unbiased; automatic application is disabled."
                ]

    output_dir.mkdir(parents=True, exist_ok=True)
    current = result["current"]
    recommended = dict(result["recommended"])  # type: ignore[arg-type]
    for key in PARAMETER_KEYS:
        detail = result["parameters"].get(key, {}) if isinstance(result.get("parameters"), dict) else {}
        if isinstance(detail, dict) and detail.get("decision") != "apply":
            recommended[key] = rounded(float(current[key]))  # type: ignore[index]
    if float(recommended["gb_on_strength"]) <= float(recommended["gb_off_strength"]):
        off_value = min(float(recommended["gb_off_strength"]), 0.93)
        recommended["gb_off_strength"] = rounded(off_value)
        recommended["gb_on_strength"] = rounded(off_value + 0.02)
    deltas = {key: rounded(float(recommended[key]) - float(current[key])) for key in PARAMETER_KEYS}  # type: ignore[index]
    artifact = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "calibration_status": status,
        "apply_allowed": bool(apply_allowed),
        "input_kind": input_kind,
        "random_seed": int(seed),
        "input_provenance": {
            "candidate_pairs": str(candidate_pairs),
            "candidate_pairs_sha256": sha256_file(candidate_pairs),
            "domain_diagnostics": str(domain_diagnostics) if domain_diagnostics else None,
            "domain_diagnostics_sha256": sha256_file(domain_diagnostics),
            "sampling_manifest": str(sampling_manifest) if sampling_manifest else None,
            "sampling_manifest_sha256": sha256_file(sampling_manifest),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "potential_validation_status": potential_validation_status,
        },
        "filters": {
            **filter_summary,
            "include_excluded": bool(include_excluded),
            "min_pairs": int(min_pairs),
            "min_oriented_pairs": int(min_oriented_pairs),
            "min_domains": int(min_domains),
            "use_sampling_weights": bool(use_sampling_weights),
            "fit_pair_rows": int(len(fit_pair_rows)),
        },
        "method": {
            "calibration_mode": str(calibration_mode),
            "p2_cut": "Weighted 1D Gaussian mixture (k=2) posterior boundary cross-checked by weighted 1D Gaussian KDE valley.",
            "gb_strength": "Weighted 1D Gaussian mixture (k=3) on oriented attractive pairs; low/mid boundary gives gb_off and mid/high boundary gives gb_on.",
            "joint2d": "Optional weighted 2D rectangular gate grid in (GB strength, P2) space; records joint gray/strong support and only adjusts scalar thresholds when support is sufficient.",
            "robust_min_s2": "1D Gaussian mixture (k=2) on domain S2 for size >= n_min when domain diagnostics are available.",
            "bootstrap": "Block bootstrap by source_file+timestep; fixed RNG seed; confidence downgraded when CI is unavailable or wide.",
        },
        "current": current,
        "raw_recommended": raw_recommended,
        "recommended": recommended,
        "delta_recommended_minus_current": deltas,
        "parameters": result["parameters"],
        "sample_sizes": {
            **filter_summary,
            **result["sample_sizes"],  # type: ignore[arg-type]
        },
        "joint_2d": result.get("joint_2d"),
        "sampling_summary": sampling_summary,
        "weighted_vs_unweighted_check": unweighted_check,
        "bootstrap": bootstrap,
        "warnings": warnings,
        "notes": result["notes"],
    }
    artifact_path = output_dir / "calibration_artifact.json"
    artifact_path.write_text(json.dumps(artifact, indent=2, ensure_ascii=False), encoding="utf-8")
    if isinstance(artifact.get("joint_2d"), dict):
        write_joint2d_grid(output_dir / "joint2d_grid.tsv", artifact["joint_2d"])  # type: ignore[arg-type]

    recommendation = {
        "method": artifact["method"],
        "current": current,
        "raw_recommended": raw_recommended,
        "recommended": recommended,
        "delta_recommended_minus_current": deltas,
        "confidence": {key: artifact["parameters"][key]["confidence"] for key in PARAMETER_KEYS},  # type: ignore[index]
        "decisions": {key: artifact["parameters"][key]["decision"] for key in PARAMETER_KEYS},  # type: ignore[index]
        "sample_sizes": artifact["sample_sizes"],
        "sampling_summary": sampling_summary,
        "bootstrap": bootstrap.get("parameter_ci", {}),
        "calibration_status": status,
        "apply_allowed": bool(apply_allowed),
        "calibration_artifact": str(artifact_path),
        "joint_2d": artifact.get("joint_2d"),
        "joint2d_grid": str(output_dir / "joint2d_grid.tsv") if isinstance(artifact.get("joint_2d"), dict) else None,
        "warnings": warnings,
        "notes": result["notes"],
    }
    (output_dir / "threshold_recommendations.json").write_text(
        json.dumps(recommendation, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_bootstrap_table(output_dir / "threshold_bootstrap_summary.tsv", bootstrap)
    write_pair_bin_summary(pair_rows, output_dir / "candidate_pair_bin_summary.tsv")
    write_plots(pair_rows, domain_rows, output_dir, recommendation)
    return artifact


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Precision calibration for LC domain/pearl threshold parameters.")
    parser.add_argument("--candidate-pairs", type=Path, required=True, help="Unscreened GB candidate-pair TSV, or screened edge table for diagnostic-only reports.")
    parser.add_argument("--domain-diagnostics", type=Path, default=None, help="Optional domain_diagnostics.tsv for robust_min_s2 calibration.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--current-gb-off", type=float, default=0.12)
    parser.add_argument("--current-gb-on", type=float, default=0.30)
    parser.add_argument("--current-p2-cut", type=float, default=0.70)
    parser.add_argument("--current-s2-cut", type=float, default=0.70)
    parser.add_argument("--n-min", type=int, default=3)
    parser.add_argument("--min-pairs", type=int, default=100)
    parser.add_argument("--min-oriented-pairs", type=int, default=60)
    parser.add_argument("--min-domains", type=int, default=25)
    parser.add_argument("--bootstrap-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260429)
    parser.add_argument("--workers", default="1", help="Parallel bootstrap workers. Use an integer or 'auto'. Default: 1.")
    parser.add_argument("--include-excluded", action="store_true", help="Include is_excluded=1 rows. Default excludes them from calibration.")
    parser.add_argument("--input-kind", choices=["candidate_pairs", "screened_edges"], default="candidate_pairs")
    parser.add_argument("--potential-validation-status", default="unknown")
    parser.add_argument("--sampling-manifest", type=Path, default=None)
    parser.add_argument("--no-sampling-weights", dest="use_sampling_weights", action="store_false", default=True)
    parser.add_argument("--calibration-mode", choices=list(CALIBRATION_MODES), default="joint2d")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.bootstrap_samples < 0:
        raise SystemExit("--bootstrap-samples must be non-negative")
    try:
        resolve_worker_count(args.workers, max(int(args.bootstrap_samples), 1))
    except ValueError as exc:
        raise SystemExit(str(exc)) from None
    artifact = calibrate_thresholds(
        candidate_pairs=args.candidate_pairs,
        domain_diagnostics=args.domain_diagnostics,
        output_dir=args.output_dir,
        current_gb_off=float(args.current_gb_off),
        current_gb_on=float(args.current_gb_on),
        current_p2_cut=float(args.current_p2_cut),
        current_s2_cut=float(args.current_s2_cut),
        n_min=int(args.n_min),
        min_pairs=int(args.min_pairs),
        min_oriented_pairs=int(args.min_oriented_pairs),
        min_domains=int(args.min_domains),
        bootstrap_samples=int(args.bootstrap_samples),
        seed=int(args.seed),
        workers=args.workers,
        include_excluded=bool(args.include_excluded),
        input_kind=args.input_kind,
        potential_validation_status=args.potential_validation_status,
        sampling_manifest=args.sampling_manifest,
        use_sampling_weights=bool(args.use_sampling_weights),
        calibration_mode=str(args.calibration_mode),
    )
    print(json.dumps({
        "calibration_status": artifact["calibration_status"],
        "apply_allowed": artifact["apply_allowed"],
        "recommended": artifact["recommended"],
        "warnings": artifact["warnings"],
        "joint_2d": {
            "status": artifact.get("joint_2d", {}).get("status") if isinstance(artifact.get("joint_2d"), dict) else None,
            "selected": artifact.get("joint_2d", {}).get("selected") if isinstance(artifact.get("joint_2d"), dict) else None,
        },
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
