#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

import liquid_crystal_aggregation as lca
import lc_sampling_manifest as sampling_manifest
from lc_threshold_recommend import recommend_from_rows, write_pair_bin_summary, write_plots


PAIR_COLUMNS = [
    "source_file",
    "timestep",
    "atom_i",
    "atom_j",
    "distance",
    "pair_energy",
    "well_depth",
    "attraction_strength",
    "q_score",
    "p2_score",
    "delta_s",
    "is_local",
    "is_excluded",
    "frame_index",
    "frame_weight",
    "frame_sampling_reason",
    "frame_event_score",
    "weight_semantics",
    "inclusion_probability_known",
    "sampling_strata",
    "sampling_weight",
    "sampling_probability",
    "pair_sampling_rate",
    "candidate_pairs_in_stratum",
    "selected_pairs_in_stratum",
    "sampling_selected_reason",
]


def iter_input_files(inputs: Sequence[Path], pattern: str, recursive: bool) -> List[Path]:
    return lca.iter_input_files(inputs, pattern=pattern, recursive=recursive)


def finite_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def safe_int(value: object, default: int = 0) -> int:
    parsed = finite_float(value)
    if not math.isfinite(parsed):
        return default
    return int(parsed)


def deterministic_sample_indices(n_items: int, max_items: int, seed: int) -> np.ndarray:
    if max_items <= 0 or n_items <= max_items:
        return np.arange(n_items, dtype=int)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_items, size=max_items, replace=False))


def stable_sample_seed(*parts: object) -> int:
    encoded = "\u241f".join(str(part) for part in parts).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:16], 16) % (2**32)


def pair_sampling_stratum(
    row: Dict[str, object],
    *,
    gb_off: float,
    gb_on: float,
    p2_cut: float,
    distance_cap: float,
) -> str:
    strength = finite_float(row.get("attraction_strength"))
    p2_score = finite_float(row.get("p2_score"))
    distance = finite_float(row.get("distance"))
    if not math.isfinite(strength):
        strength_bin = "strength_invalid"
    elif strength >= gb_on:
        strength_bin = "strength_strong"
    elif strength >= gb_off:
        strength_bin = "strength_gray"
    elif strength > 0.0:
        strength_bin = "strength_weak"
    else:
        strength_bin = "strength_zero"
    if math.isfinite(strength) and (
        abs(strength - gb_off) <= max(0.015, 0.08 * max(gb_off, 1e-6))
        or abs(strength - gb_on) <= max(0.02, 0.08 * max(gb_on, 1e-6))
    ):
        strength_bin += "_near_threshold"

    if not math.isfinite(p2_score):
        p2_bin = "p2_invalid"
    elif abs(p2_score - p2_cut) <= 0.05:
        p2_bin = "p2_near_cut"
    elif p2_score >= p2_cut:
        p2_bin = "p2_high"
    elif p2_score >= 0.30:
        p2_bin = "p2_mid"
    else:
        p2_bin = "p2_low"

    if not math.isfinite(distance) or distance_cap <= 0.0:
        distance_bin = "distance_invalid"
    else:
        frac = distance / max(distance_cap, 1e-12)
        if frac >= 0.90:
            distance_bin = "distance_near_cutoff"
        elif frac >= 0.60:
            distance_bin = "distance_mid"
        else:
            distance_bin = "distance_close"

    local_bin = "local" if safe_int(row.get("is_local")) == 1 else "nonlocal"
    excluded_bin = "excluded" if safe_int(row.get("is_excluded")) == 1 else "included"
    return "|".join((strength_bin, p2_bin, distance_bin, local_bin, excluded_bin))


def quantile_bin_labels(rows: Sequence[Dict[str, object]]) -> List[str]:
    fields = [
        ("attraction_strength", "strength_q"),
        ("p2_score", "p2_q"),
        ("distance", "distance_q"),
    ]
    values_by_field: Dict[str, np.ndarray] = {}
    for field, _prefix in fields:
        arr = np.asarray([finite_float(row.get(field)) for row in rows], dtype=float)
        finite = arr[np.isfinite(arr)]
        values_by_field[field] = finite
    edges_by_field: Dict[str, np.ndarray] = {}
    for field, _prefix in fields:
        finite = values_by_field[field]
        if finite.size >= 5 and float(np.max(finite) - np.min(finite)) > 1e-12:
            edges_by_field[field] = np.unique(np.quantile(finite, [0.2, 0.4, 0.6, 0.8]))
        else:
            edges_by_field[field] = np.asarray([], dtype=float)
    labels: List[str] = []
    for row in rows:
        parts: List[str] = []
        for field, prefix in fields:
            value = finite_float(row.get(field))
            edges = edges_by_field[field]
            if not math.isfinite(value):
                parts.append(f"{prefix}invalid")
            elif edges.size == 0:
                parts.append(f"{prefix}all")
            else:
                parts.append(f"{prefix}{int(np.searchsorted(edges, value, side='right'))}")
        labels.append("|".join(parts))
    return labels


def stratified_sample_pair_rows(
    rows: Sequence[Dict[str, object]],
    *,
    max_items: int,
    seed: int,
    gb_off: float,
    gb_on: float,
    p2_cut: float,
    distance_cap: float,
    min_per_stratum: int = 2,
    frame_weight: float = 1.0,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    prepared: List[Dict[str, object]] = []
    strata: Dict[str, List[int]] = {}
    near_threshold: List[int] = []
    quantile_labels = quantile_bin_labels(rows)
    for idx, row in enumerate(rows):
        enriched = dict(row)
        threshold_stratum = pair_sampling_stratum(
            enriched,
            gb_off=gb_off,
            gb_on=gb_on,
            p2_cut=p2_cut,
            distance_cap=distance_cap,
        )
        stratum = f"{threshold_stratum}|{quantile_labels[idx]}"
        enriched["sampling_strata"] = stratum
        enriched["frame_weight"] = float(frame_weight)
        enriched["sampling_weight"] = float(frame_weight)
        enriched["sampling_probability"] = "" if str(enriched.get("inclusion_probability_known", "")) == "0" else 1.0
        enriched["pair_sampling_rate"] = 1.0
        enriched["candidate_pairs_in_stratum"] = 1
        enriched["selected_pairs_in_stratum"] = 1
        enriched["sampling_selected_reason"] = "all"
        prepared.append(enriched)
        strata.setdefault(stratum, []).append(idx)
        strength = finite_float(row.get("attraction_strength"))
        p2_score = finite_float(row.get("p2_score"))
        distance = finite_float(row.get("distance"))
        if (
            math.isfinite(strength)
            and (
                abs(strength - gb_off) <= max(0.015, 0.08 * max(gb_off, 1e-6))
                or abs(strength - gb_on) <= max(0.02, 0.08 * max(gb_on, 1e-6))
            )
        ) or (math.isfinite(p2_score) and abs(p2_score - p2_cut) <= 0.05) or (
            math.isfinite(distance) and distance_cap > 0.0 and distance >= 0.90 * distance_cap
        ):
            near_threshold.append(idx)

    if max_items <= 0 or len(prepared) <= max_items:
        for stratum, indices in strata.items():
            candidate_count = len(indices)
            for idx in indices:
                prepared[idx]["candidate_pairs_in_stratum"] = int(candidate_count)
                prepared[idx]["selected_pairs_in_stratum"] = int(candidate_count)
                prepared[idx]["pair_sampling_rate"] = 1.0
        return prepared, {
            "strategy": "all",
            "candidate_pairs": int(len(prepared)),
            "selected_pairs": int(len(prepared)),
            "strata_total": int(len(strata)),
            "strata_selected": int(len(strata)),
        }

    rng = np.random.default_rng(int(seed))
    selected: set[int] = set()
    reasons: Dict[int, str] = {}

    strong_indices = [
        idx for idx, row in enumerate(prepared)
        if finite_float(row.get("attraction_strength")) >= gb_on
    ]
    for pool, reason in (
        (near_threshold, "near_threshold"),
        (strong_indices, "strong_extreme"),
    ):
        if not pool:
            continue
        take = min(len(pool), max(1, max_items // 5), max_items - len(selected))
        if take <= 0:
            break
        chosen = rng.choice(np.asarray(sorted(set(pool)), dtype=int), size=take, replace=False)
        for idx in chosen.tolist():
            selected.add(int(idx))
            reasons[int(idx)] = reason

    critical_strata = [
        (stratum, indices)
        for stratum, indices in strata.items()
        if (
            "near_threshold" in stratum
            or "p2_near_cut" in stratum
            or "distance_near_cutoff" in stratum
            or "strength_strong" in stratum
            or "strength_q0" in stratum
            or "strength_q4" in stratum
            or "p2_q0" in stratum
            or "p2_q4" in stratum
        )
    ]
    for stratum, indices in sorted(critical_strata, key=lambda item: (len(item[1]), item[0])):
        remaining = [idx for idx in indices if idx not in selected]
        if not remaining or len(selected) >= max_items:
            continue
        chosen = rng.choice(np.asarray(remaining, dtype=int), size=1, replace=False)
        for idx in chosen.tolist():
            selected.add(int(idx))
            reasons.setdefault(int(idx), "critical_stratum")

    for stratum, indices in sorted(strata.items(), key=lambda item: (-len(item[1]), item[0])):
        remaining = [idx for idx in indices if idx not in selected]
        if not remaining or len(selected) >= max_items:
            continue
        take = min(len(remaining), max(1, int(min_per_stratum)), max_items - len(selected))
        chosen = rng.choice(np.asarray(remaining, dtype=int), size=take, replace=False)
        for idx in chosen.tolist():
            selected.add(int(idx))
            reasons.setdefault(int(idx), "stratum_floor")

    if len(selected) < max_items:
        remaining = [idx for idx in range(len(prepared)) if idx not in selected]
        take = min(len(remaining), max_items - len(selected))
        if take > 0:
            chosen = rng.choice(np.asarray(remaining, dtype=int), size=take, replace=False)
            for idx in chosen.tolist():
                selected.add(int(idx))
                reasons.setdefault(int(idx), "random_fill")

    selected_indices = sorted(selected)
    selected_by_stratum: Dict[str, int] = {}
    for idx in selected_indices:
        selected_by_stratum[str(prepared[idx]["sampling_strata"])] = selected_by_stratum.get(str(prepared[idx]["sampling_strata"]), 0) + 1
    sampled: List[Dict[str, object]] = []
    for idx in selected_indices:
        row = dict(prepared[idx])
        stratum = str(row["sampling_strata"])
        candidate_count = len(strata[stratum])
        selected_count = max(selected_by_stratum.get(stratum, 1), 1)
        pair_weight = float(candidate_count / selected_count)
        row["sampling_weight"] = float(frame_weight) * pair_weight
        rate = float(selected_count / max(candidate_count, 1))
        row["sampling_probability"] = "" if str(row.get("inclusion_probability_known", "")) == "0" else rate
        row["pair_sampling_rate"] = rate
        row["candidate_pairs_in_stratum"] = int(candidate_count)
        row["selected_pairs_in_stratum"] = int(selected_count)
        row["sampling_selected_reason"] = reasons.get(idx, "selected")
        sampled.append(row)
    return sampled, {
        "strategy": "stratified",
        "candidate_pairs": int(len(prepared)),
        "selected_pairs": int(len(sampled)),
        "strata_total": int(len(strata)),
        "strata_selected": int(len(selected_by_stratum)),
        "near_threshold_candidates": int(len(set(near_threshold))),
        "strong_candidates": int(len(strong_indices)),
        "min_per_stratum": int(min_per_stratum),
    }


FrameRecord = Tuple[int, int, lca.BoxSpec, Dict[str, int], np.ndarray]
SelectedFrameRecord = Tuple[int, int, lca.BoxSpec, Dict[str, int], np.ndarray, float, str, float]


def sampling_semantics(frame_sample_strategy: str, frame_sampling_reason: str) -> Tuple[str, int]:
    if str(frame_sample_strategy) == "event-aware" and not str(frame_sampling_reason).startswith("event_static_fallback_"):
        return "coverage_expansion_not_unbiased", 0
    return "inverse_frame_pair_stratum_weight", 1


def frame_feature_value(record: FrameRecord) -> float:
    _frame_index, _timestep, box, col_index, data = record
    try:
        positions = lca.extract_positions(data, col_index)
        finite = positions[np.isfinite(positions).all(axis=1)]
        if finite.size:
            span = np.ptp(finite, axis=0)
            rg = np.linalg.norm(finite - np.mean(finite, axis=0), axis=1)
            return float(np.linalg.norm(span) + np.mean(rg))
    except Exception:
        pass
    arr = np.asarray(data, dtype=float)
    finite_values = arr[np.isfinite(arr)]
    if finite_values.size == 0:
        return 0.0
    return float(np.mean(finite_values))


def frame_event_scores(records: Sequence[FrameRecord]) -> List[float]:
    if not records:
        return []
    values = np.asarray([frame_feature_value(record) for record in records], dtype=float)
    if values.size == 1:
        return [0.0]
    scores: List[float] = []
    for idx, value in enumerate(values):
        previous_delta = abs(float(value - values[idx - 1])) if idx > 0 else 0.0
        next_delta = abs(float(values[idx + 1] - value)) if idx + 1 < values.size else 0.0
        scores.append(float(max(previous_delta, next_delta)))
    return scores


def event_center_rank(records: Sequence[FrameRecord], scores: Sequence[float]) -> List[int]:
    values = np.asarray([frame_feature_value(record) for record in records], dtype=float)
    ranked: List[Tuple[float, float, int]] = []
    for idx, score in enumerate(scores):
        left = values[idx - 1] if idx > 0 else values[idx]
        right = values[idx + 1] if idx + 1 < values.size else values[idx]
        local_extremum = abs(float(values[idx] - 0.5 * (left + right)))
        ranked.append((-float(score), -local_extremum, idx))
    return [idx for _neg_score, _neg_extremum, idx in sorted(ranked)]


def select_frame_records(records: Sequence[FrameRecord], frame_limit: int, strategy: str) -> List[SelectedFrameRecord]:
    if frame_limit <= 0 or len(records) <= frame_limit:
        return [(idx, timestep, box, col_index, data, 1.0, "all", 0.0) for idx, timestep, box, col_index, data in records]
    if strategy == "first":
        selected_indices = list(range(frame_limit))
        frame_weight = float(len(records) / max(len(selected_indices), 1))
        return [(*records[idx], frame_weight, "first", 0.0) for idx in selected_indices]
    if strategy == "event-aware":
        if frame_limit < 3:
            raise ValueError("event-aware frame sampling requires frame_limit >= 3 to keep event center and neighbors")
        scores = frame_event_scores(records)
        max_score = max(scores) if scores else 0.0
        if max_score <= 1e-12:
            return [
                (
                    frame_index,
                    timestep,
                    box,
                    col_index,
                    data,
                    frame_weight,
                    f"event_static_fallback_{reason}",
                    event_score,
                )
                for frame_index, timestep, box, col_index, data, frame_weight, reason, event_score
                in select_frame_records(records, frame_limit, "stratified")
            ]
        selected: set[int] = set()
        ranked = event_center_rank(records, scores)
        reason_by_idx: Dict[int, str] = {}
        for idx in ranked:
            if len(selected) >= frame_limit:
                break
            event_group = [idx]
            for neighbor in (idx - 1, idx + 1):
                if 0 <= neighbor < len(records) and len(event_group) < frame_limit:
                    event_group.append(neighbor)
            event_group = sorted(set(event_group))
            if len(selected | set(event_group)) > frame_limit and selected:
                continue
            for neighbor in event_group:
                selected.add(neighbor)
                if neighbor == idx:
                    reason_by_idx[neighbor] = "event_center"
                else:
                    reason_by_idx.setdefault(neighbor, "event_neighbor")
            if len(selected) >= frame_limit:
                break
        stratified_indices = np.linspace(0, len(records) - 1, num=frame_limit)
        for value in stratified_indices:
            if len(selected) >= frame_limit:
                break
            idx = int(round(value))
            selected.add(idx)
            reason_by_idx.setdefault(idx, "time_stratified")
        cursor = 0
        while len(selected) < frame_limit and cursor < len(records):
            selected.add(cursor)
            reason_by_idx.setdefault(cursor, "fill")
            cursor += 1
        if len(selected) > frame_limit:
            centers = [idx for idx, reason in reason_by_idx.items() if reason == "event_center"]
            keep = set(sorted(centers, key=lambda idx: (-scores[idx], idx))[:1])
            for idx in sorted(selected, key=lambda idx: (idx not in keep, abs(idx - next(iter(keep))) if keep else idx, idx)):
                if len(keep) >= frame_limit:
                    break
                keep.add(idx)
            selected = keep
        selected_indices = sorted(selected)
        return [
            (
                *records[idx],
                1.0,
                reason_by_idx.get(idx, "fill"),
                float(scores[idx]) if scores else 0.0,
            )
            for idx in selected_indices
        ]
    if strategy != "stratified":
        raise ValueError(f"unknown frame sample strategy: {strategy}")
    indices = np.linspace(0, len(records) - 1, num=frame_limit)
    chosen = sorted({int(round(value)) for value in indices})
    cursor = 0
    while len(chosen) < frame_limit and cursor < len(records):
        if cursor not in chosen:
            chosen.append(cursor)
        cursor += 1
    selected_indices = sorted(chosen[:frame_limit])
    frame_weight = float(len(records) / max(len(selected_indices), 1))
    return [(*records[idx], frame_weight, "time_stratified", 0.0) for idx in selected_indices]


def collect_candidate_pairs_for_frame(
    *,
    source_file: Path,
    timestep: int,
    box: lca.BoxSpec,
    data: np.ndarray,
    col_index: Dict[str, int],
    config: lca.AggregationConfig,
    sample_pairs_per_frame: int = 0,
    pair_sample_strategy: str = "stratified",
    sample_seed: int = 20260429,
    min_per_stratum: int = 2,
    frame_index: int = 0,
    frame_weight: float = 1.0,
    frame_sampling_reason: str = "all",
    frame_event_score: float = 0.0,
    frame_sample_strategy: str = "all",
) -> List[Dict[str, object]]:
    all_pos = lca.extract_positions(data, col_index)
    all_ids = lca.extract_particle_ids(data, col_index)
    all_types = lca.extract_particle_types(data, col_index)
    all_shapes = lca.extract_shape_axes(data, col_index)
    if all_shapes is None:
        raise RuntimeError(f"{source_file}: GB audit requires shapex/shapey/shapez columns")
    all_quaternions = lca.extract_quaternions(data, col_index)
    mesogen_indices = lca.select_mesogen_indices(all_types, config.mesogen_type)
    if mesogen_indices.size == 0:
        return []
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
    local_pairs = lca.read_pair_list(config.local_pair_file)
    excluded_pairs = lca.read_pair_list(config.exclude_pair_file)
    params = config.gayberne_params
    if params is None:
        raise RuntimeError("GB audit requires gayberne_params")
    broad_cut = float(config.r_energy_cap or params.cutoff)

    pair_indices: List[Tuple[int, int]] = []
    for i in range(pos.shape[0] - 1):
        delta = pos[i + 1 :, :] - pos[i, :]
        delta[:, 0] = lca.minimum_image(delta[:, 0], box[0])
        delta[:, 1] = lca.minimum_image(delta[:, 1], box[1])
        delta[:, 2] = lca.minimum_image(delta[:, 2], box[2])
        distances = np.linalg.norm(delta, axis=1)
        close = np.nonzero((distances > 1e-12) & (distances <= broad_cut))[0] + i + 1
        pair_indices.extend((i, int(j)) for j in close)

    rows: List[Dict[str, object]] = []
    for i, j in pair_indices:
        atom_i = int(ids[i])
        atom_j = int(ids[j])
        pair_key = lca.normalized_pair(atom_i, atom_j)
        delta_vec = lca.minimum_image_vector(pos[j, :] - pos[i, :], box)
        distance = float(np.linalg.norm(delta_vec))
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
        pair_energy, well_depth = metrics
        q_score = float(abs(np.clip(np.dot(u[i, :], u[j, :]), -1.0, 1.0)))
        p2_score = float(lca.p2(q_score))
        delta_s = int(abs(int(chain_indices[i]) - int(chain_indices[j])))
        is_local = int(delta_s <= config.s_excl or pair_key in local_pairs)
        is_excluded = int(pair_key in excluded_pairs)
        weight_semantics, inclusion_probability_known = sampling_semantics(frame_sample_strategy, frame_sampling_reason)
        rows.append(
            {
                "source_file": str(source_file),
                "timestep": int(timestep),
                "atom_i": atom_i,
                "atom_j": atom_j,
                "distance": distance,
                "pair_energy": float(pair_energy),
                "well_depth": float(well_depth),
                "attraction_strength": max(0.0, -float(pair_energy) / max(float(well_depth), 1e-24)),
                "q_score": q_score,
                "p2_score": p2_score,
                "delta_s": delta_s,
                "is_local": is_local,
                "is_excluded": is_excluded,
                "frame_index": int(frame_index),
                "frame_weight": float(frame_weight),
                "frame_sampling_reason": str(frame_sampling_reason),
                "frame_event_score": float(frame_event_score),
                "weight_semantics": weight_semantics,
                "inclusion_probability_known": inclusion_probability_known,
            }
        )
    if pair_sample_strategy == "uniform":
        selected_indices = deterministic_sample_indices(
            len(rows),
            sample_pairs_per_frame,
            seed=stable_sample_seed("gb_pair_audit_uniform", source_file.resolve(), int(timestep), len(rows), sample_seed),
        )
        selected_rows = []
        for idx in selected_indices:
            row = dict(rows[int(idx)])
            row["sampling_strata"] = "uniform"
            pair_weight = float(len(rows) / max(len(selected_indices), 1)) if sample_pairs_per_frame > 0 else 1.0
            row["sampling_weight"] = float(frame_weight) * pair_weight
            rate = float(len(selected_indices) / max(len(rows), 1))
            row["sampling_probability"] = "" if str(row.get("inclusion_probability_known", "")) == "0" else rate
            row["pair_sampling_rate"] = rate
            row["candidate_pairs_in_stratum"] = int(len(rows))
            row["selected_pairs_in_stratum"] = int(len(selected_indices))
            row["sampling_selected_reason"] = "uniform"
            selected_rows.append(row)
        return selected_rows
    if pair_sample_strategy != "stratified":
        raise ValueError(f"unknown pair sample strategy: {pair_sample_strategy}")
    selected_rows, _summary = stratified_sample_pair_rows(
        rows,
        max_items=int(sample_pairs_per_frame),
        seed=stable_sample_seed("gb_pair_audit_stratified", source_file.resolve(), int(timestep), len(rows), sample_seed),
        gb_off=float(config.gb_off_strength),
        gb_on=float(config.gb_on_strength),
        p2_cut=float(config.p2_cut),
        distance_cap=float(broad_cut),
        min_per_stratum=int(min_per_stratum),
        frame_weight=float(frame_weight),
    )
    return selected_rows


def frame_selection_records_for_output(dump_path: Path, selected: Sequence[SelectedFrameRecord], candidate_frame_count: int, strategy: str, every: int) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []
    selected_count = len(selected)
    for frame_index, timestep, _box, _col_index, _data, frame_weight, reason, event_score in selected:
        weight_semantics, inclusion_probability_known = sampling_semantics(strategy, reason)
        records.append(
            {
                "source_file": str(dump_path),
                "frame_index": int(frame_index),
                "timestep": int(timestep),
                "frame_sample_strategy": str(strategy),
                "audit_every": int(every),
                "candidate_frame_count": int(candidate_frame_count),
                "selected_frame_count": int(selected_count),
                "frame_sampling_rate_among_processed_frames": float(selected_count / max(candidate_frame_count, 1)),
                "frame_weight": float(frame_weight),
                "frame_sampling_reason": str(reason),
                "frame_event_score": float(event_score),
                "weight_semantics": weight_semantics,
                "inclusion_probability_known": inclusion_probability_known,
            }
        )
    return records


def audit_dump_file(
    dump_path: Path,
    config: lca.AggregationConfig,
    every: int,
    frame_limit: int,
    sample_pairs_per_frame: int,
    frame_sample_strategy: str = "stratified",
    pair_sample_strategy: str = "stratified",
    sample_seed: int = 20260429,
    min_per_stratum: int = 2,
) -> Tuple[Path, List[Dict[str, object]], List[Dict[str, object]]]:
    rows: List[Dict[str, object]] = []
    frame_selection_rows: List[Dict[str, object]] = []
    frame_index = 0
    frame_records: List[FrameRecord] = []
    for timestep, box, _columns, col_index, data in lca.parse_dump_frames(dump_path):
        frame_index += 1
        if (frame_index - 1) % every != 0:
            continue
        if frame_limit > 0:
            frame_records.append((int(frame_index), int(timestep), box, dict(col_index), np.array(data, copy=True)))
            continue
        rows.extend(
            collect_candidate_pairs_for_frame(
                source_file=dump_path,
                timestep=timestep,
                box=box,
                data=data,
                col_index=col_index,
                config=config,
                sample_pairs_per_frame=sample_pairs_per_frame,
                pair_sample_strategy=pair_sample_strategy,
                sample_seed=sample_seed,
                min_per_stratum=min_per_stratum,
                frame_index=frame_index,
                frame_weight=1.0,
                frame_sampling_reason="all",
                frame_event_score=0.0,
                frame_sample_strategy=frame_sample_strategy,
            )
        )
    if frame_limit > 0:
        selected_frames = select_frame_records(frame_records, frame_limit, frame_sample_strategy)
        frame_selection_rows.extend(
            frame_selection_records_for_output(dump_path, selected_frames, len(frame_records), frame_sample_strategy, every)
        )
        for selected_frame_index, timestep, box, col_index, data, frame_weight, frame_reason, event_score in selected_frames:
            rows.extend(
                collect_candidate_pairs_for_frame(
                    source_file=dump_path,
                    timestep=timestep,
                    box=box,
                    data=data,
                    col_index=col_index,
                    config=config,
                    sample_pairs_per_frame=sample_pairs_per_frame,
                    pair_sample_strategy=pair_sample_strategy,
                    sample_seed=sample_seed,
                    min_per_stratum=min_per_stratum,
                    frame_index=selected_frame_index,
                    frame_weight=frame_weight,
                    frame_sampling_reason=frame_reason,
                    frame_event_score=event_score,
                    frame_sample_strategy=frame_sample_strategy,
                )
            )
    return dump_path, rows, frame_selection_rows


def audit_dump_file_star(args: Tuple[Path, lca.AggregationConfig, int, int, int, str, str, int, int]) -> Tuple[Path, List[Dict[str, object]], List[Dict[str, object]]]:
    return audit_dump_file(*args)


def write_frame_selection_manifest(path: Path, records: Sequence[Dict[str, object]]) -> None:
    columns = [
        "source_file",
        "frame_index",
        "timestep",
        "frame_sample_strategy",
        "audit_every",
        "candidate_frame_count",
        "selected_frame_count",
        "frame_sampling_rate_among_processed_frames",
        "frame_weight",
        "frame_sampling_reason",
        "frame_event_score",
        "weight_semantics",
        "inclusion_probability_known",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow({column: row.get(column, "") for column in columns})


def summarize_frame_selection(records: Sequence[Dict[str, object]]) -> Dict[str, object]:
    reasons: Dict[str, int] = {}
    timesteps: List[int] = []
    scores: List[float] = []
    rates: List[float] = []
    for row in records:
        reason = str(row.get("frame_sampling_reason", ""))
        if reason:
            reasons[reason] = reasons.get(reason, 0) + 1
        timestep = safe_int(row.get("timestep"), default=-1)
        if timestep >= 0:
            timesteps.append(timestep)
        score = finite_float(row.get("frame_event_score"))
        if math.isfinite(score):
            scores.append(score)
        rate = finite_float(row.get("frame_sampling_rate_among_processed_frames"))
        if math.isfinite(rate):
            rates.append(rate)
    return {
        "selected_frame_rows": int(len(records)),
        "selected_timesteps": timesteps[:200],
        "frame_sampling_reasons": reasons,
        "event_score_range": [float(min(scores)), float(max(scores))] if scores else [],
        "frame_sampling_rate_range": [float(min(rates)), float(max(rates))] if rates else [],
        "weight_semantics": "coverage_expansion_not_unbiased for true event-aware selections; inverse_frame_pair_stratum_weight for stratified/first/all and event_static_fallback_* rows.",
    }


def write_rows(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in PAIR_COLUMNS})


def cap_candidate_rows_for_calibration(
    rows: Sequence[Dict[str, object]],
    max_total: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    input_count = int(len(rows))
    if max_total <= 0 or input_count <= int(max_total):
        return list(rows), {
            "cap_applied": False,
            "input_candidate_pairs": input_count,
            "written_candidate_pairs": input_count,
            "dropped_candidate_pairs": 0,
            "max_candidate_pairs_total": int(max_total),
            "selection_method": "none",
        }
    target = max(1, int(max_total))
    indices = np.linspace(0, input_count - 1, num=target, dtype=int)
    selected_indices: List[int] = []
    seen = set()
    for raw_idx in indices:
        idx = int(raw_idx)
        if idx not in seen:
            selected_indices.append(idx)
            seen.add(idx)
    cursor = 0
    while len(selected_indices) < target and cursor < input_count:
        if cursor not in seen:
            selected_indices.append(cursor)
            seen.add(cursor)
        cursor += 1
    selected_indices.sort()
    capped = [dict(rows[idx]) for idx in selected_indices[:target]]
    return capped, {
        "cap_applied": True,
        "input_candidate_pairs": input_count,
        "written_candidate_pairs": int(len(capped)),
        "dropped_candidate_pairs": int(input_count - len(capped)),
        "max_candidate_pairs_total": int(max_total),
        "selection_method": "deterministic_strided_across_sorted_candidate_rows",
        "selected_index_preview": selected_indices[:20],
    }


def append_rows_with_rolling_cap(
    current: List[Dict[str, object]],
    new_rows: Sequence[Dict[str, object]],
    max_total: int,
) -> List[Dict[str, object]]:
    if not new_rows:
        return current
    current.extend(dict(row) for row in new_rows)
    if max_total > 0 and len(current) > max_total:
        current, _ = cap_candidate_rows_for_calibration(current, max_total)
    return current


def apply_global_candidate_cap_weights(
    rows: Sequence[Dict[str, object]],
    cap_summary: Dict[str, object],
) -> List[Dict[str, object]]:
    if not cap_summary.get("cap_applied") or not rows:
        return [dict(row) for row in rows]
    input_count = max(1.0, float(cap_summary.get("input_candidate_pairs", len(rows))))
    written_count = max(1.0, float(cap_summary.get("written_candidate_pairs", len(rows))))
    weight_factor = input_count / written_count
    cap_summary["global_cap_weight_factor"] = float(weight_factor)
    weighted_rows: List[Dict[str, object]] = []
    for row in rows:
        updated = dict(row)
        weight = finite_float(updated.get("sampling_weight"))
        if not math.isfinite(weight) or weight <= 0.0:
            weight = 1.0
        updated["sampling_weight"] = float(weight * weight_factor)
        semantics = str(updated.get("weight_semantics", "")).strip()
        if semantics:
            updated["weight_semantics"] = f"{semantics}+global_candidate_cap_weight"
        else:
            updated["weight_semantics"] = "global_candidate_cap_weight"
        updated["inclusion_probability_known"] = 0
        weighted_rows.append(updated)
    return weighted_rows


def pair_sampling_audit_records(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    strata: Dict[str, Dict[str, object]] = {}
    for row in rows:
        key = str(row.get("sampling_strata", "unspecified"))
        record = strata.setdefault(
            key,
            {
                "sampling_strata": key,
                "selected_pairs": 0,
                "estimated_candidate_pairs": 0.0,
                "sampling_rate": 0.0,
                "sampling_weight_min": math.inf,
                "sampling_weight_max": 0.0,
                "inclusion_probability_known_rows": 0,
                "inclusion_probability_unknown_rows": 0,
                "weight_semantics": {},
                "selection_reasons": {},
            },
        )
        weight = finite_float(row.get("sampling_weight"))
        if not math.isfinite(weight) or weight <= 0.0:
            weight = 1.0
        record["selected_pairs"] = int(record["selected_pairs"]) + 1
        record["estimated_candidate_pairs"] = float(record["estimated_candidate_pairs"]) + weight
        record["sampling_weight_min"] = min(float(record["sampling_weight_min"]), weight)
        record["sampling_weight_max"] = max(float(record["sampling_weight_max"]), weight)
        if safe_int(row.get("inclusion_probability_known"), default=1) == 1:
            record["inclusion_probability_known_rows"] = int(record["inclusion_probability_known_rows"]) + 1
        else:
            record["inclusion_probability_unknown_rows"] = int(record["inclusion_probability_unknown_rows"]) + 1
        semantics = record["weight_semantics"]
        if isinstance(semantics, dict):
            label = str(row.get("weight_semantics", "unspecified"))
            semantics[label] = semantics.get(label, 0) + 1
        reason = str(row.get("sampling_selected_reason", ""))
        if reason:
            reasons = record["selection_reasons"]
            if isinstance(reasons, dict):
                reasons[reason] = reasons.get(reason, 0) + 1
    output: List[Dict[str, object]] = []
    for record in strata.values():
        selected = int(record["selected_pairs"])
        estimated = float(record["estimated_candidate_pairs"])
        record["sampling_rate"] = float(selected / max(estimated, 1e-12))
        if not math.isfinite(float(record["sampling_weight_min"])):
            record["sampling_weight_min"] = 0.0
        record["weight_semantics"] = json.dumps(record["weight_semantics"], sort_keys=True)
        record["selection_reasons"] = json.dumps(record["selection_reasons"], sort_keys=True)
        output.append(record)
    return sorted(output, key=lambda row: str(row["sampling_strata"]))


def pair_sampling_manifest_records(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    output: List[Dict[str, object]] = []
    for row in rows:
        weight = finite_float(row.get("sampling_weight"))
        if not math.isfinite(weight) or weight <= 0.0:
            weight = 1.0
        output.append(
            {
                "source_file": row.get("source_file", ""),
                "timestep": safe_int(row.get("timestep")),
                "frame_index": safe_int(row.get("frame_index")),
                "atom_i": safe_int(row.get("atom_i")),
                "atom_j": safe_int(row.get("atom_j")),
                "sampling_strata": row.get("sampling_strata", ""),
                "sampling_selected_reason": row.get("sampling_selected_reason", ""),
                "sampling_weight": float(weight),
                "sampling_probability": "" if safe_int(row.get("inclusion_probability_known"), default=0) == 0 else finite_float(row.get("sampling_probability")),
                "candidate_pairs_in_stratum": safe_int(row.get("candidate_pairs_in_stratum"), default=max(1, int(round(weight)))),
                "selected_pairs_in_stratum": safe_int(row.get("selected_pairs_in_stratum"), default=1),
                "estimated_stratum_candidate_pairs": float(safe_int(row.get("candidate_pairs_in_stratum"), default=max(1, int(round(weight))))),
                "pair_sampling_rate": finite_float(row.get("pair_sampling_rate")) if math.isfinite(finite_float(row.get("pair_sampling_rate"))) else float(1.0 / max(weight, 1e-12)),
                "frame_sampling_reason": row.get("frame_sampling_reason", ""),
                "frame_event_score": finite_float(row.get("frame_event_score")),
                "weight_semantics": row.get("weight_semantics", ""),
                "inclusion_probability_known": safe_int(row.get("inclusion_probability_known"), default=0),
            }
        )
    return output


def write_pair_sampling_manifest(path: Path, records: Sequence[Dict[str, object]]) -> None:
    columns = [
        "source_file",
        "timestep",
        "frame_index",
        "atom_i",
        "atom_j",
        "sampling_strata",
        "sampling_selected_reason",
        "sampling_weight",
        "sampling_probability",
        "candidate_pairs_in_stratum",
        "selected_pairs_in_stratum",
        "estimated_stratum_candidate_pairs",
        "pair_sampling_rate",
        "frame_sampling_reason",
        "frame_event_score",
        "weight_semantics",
        "inclusion_probability_known",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_pair_sampling_audit(path: Path, records: Sequence[Dict[str, object]]) -> None:
    columns = [
        "sampling_strata",
        "selected_pairs",
        "estimated_candidate_pairs",
        "sampling_rate",
        "sampling_weight_min",
        "sampling_weight_max",
        "inclusion_probability_known_rows",
        "inclusion_probability_unknown_rows",
        "weight_semantics",
        "selection_reasons",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow({column: row.get(column, "") for column in columns})


def summarize_rows(rows: Sequence[Dict[str, object]], recommendation: Dict[str, object]) -> Dict[str, object]:
    strengths = np.array([finite_float(row.get("attraction_strength")) for row in rows], dtype=float)
    p2_values = np.array([finite_float(row.get("p2_score")) for row in rows], dtype=float)
    energies = np.array([finite_float(row.get("pair_energy")) for row in rows], dtype=float)
    mask = np.isfinite(strengths) & np.isfinite(p2_values) & np.isfinite(energies)
    strengths = strengths[mask]
    p2_values = p2_values[mask]
    energies = energies[mask]
    recommended = recommendation.get("recommended", {})
    gb_off = float(recommended.get("gb_off_strength", 0.12))
    gb_on = float(recommended.get("gb_on_strength", 0.30))
    p2_cut = float(recommended.get("p2_cut", 0.70))
    return {
        "pair_rows": int(strengths.size),
        "attractive_pair_rows": int(np.sum(energies < 0.0)),
        "excluded_pair_rows": int(sum(1 for row in rows if safe_int(row.get("is_excluded")) == 1)),
        "calibration_pair_rows": int(sum(1 for row in rows if safe_int(row.get("is_excluded")) != 1)),
        "recommended_counts": {
            "gray_or_stronger": int(np.sum((strengths >= gb_off) & (p2_values >= p2_cut))),
            "strong": int(np.sum((strengths >= gb_on) & (p2_values >= p2_cut))),
        },
        "quantiles": {
            "gb_strength": [float(v) for v in np.quantile(strengths, [0.1, 0.25, 0.5, 0.75, 0.9]).tolist()] if strengths.size else [],
            "p2_score": [float(v) for v in np.quantile(p2_values, [0.1, 0.25, 0.5, 0.75, 0.9]).tolist()] if p2_values.size else [],
            "pair_energy": [float(v) for v in np.quantile(energies, [0.1, 0.25, 0.5, 0.75, 0.9]).tolist()] if energies.size else [],
        },
        "sampling": summarize_sampling(rows),
    }


def summarize_sampling(rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    strata: Dict[str, int] = {}
    reasons: Dict[str, int] = {}
    timesteps: Dict[str, int] = {}
    frame_weights: List[float] = []
    frame_reasons: Dict[str, int] = {}
    weight_sum = 0.0
    weighted = 0
    for row in rows:
        stratum = str(row.get("sampling_strata", ""))
        reason = str(row.get("sampling_selected_reason", ""))
        if stratum:
            strata[stratum] = strata.get(stratum, 0) + 1
        if reason:
            reasons[reason] = reasons.get(reason, 0) + 1
        timestep = str(row.get("timestep", ""))
        if timestep:
            timesteps[timestep] = timesteps.get(timestep, 0) + 1
        frame_weight = finite_float(row.get("frame_weight"))
        if math.isfinite(frame_weight):
            frame_weights.append(frame_weight)
        frame_reason = str(row.get("frame_sampling_reason", ""))
        if frame_reason:
            frame_reasons[frame_reason] = frame_reasons.get(frame_reason, 0) + 1
        weight = finite_float(row.get("sampling_weight"))
        if math.isfinite(weight):
            weight_sum += weight
            weighted += 1
    return {
        "selected_rows": int(len(rows)),
        "weighted_pair_rows": float(weight_sum),
        "rows_with_sampling_weight": int(weighted),
        "unique_sampling_strata": int(len(strata)),
        "selected_timesteps": int(len(timesteps)),
        "frame_weight_range": [float(min(frame_weights)), float(max(frame_weights))] if frame_weights else [],
        "frame_sampling_reasons": frame_reasons,
        "top_sampling_strata": sorted(strata.items(), key=lambda item: (-item[1], item[0]))[:20],
        "selection_reasons": reasons,
    }


def build_config(args: argparse.Namespace) -> lca.AggregationConfig:
    params = lca.parse_gayberne_params_from_lmp(args.gb_param_file)
    r_energy_cap = lca.parse_optional_positive_float(args.r_energy_cap)
    return lca.AggregationConfig(
        axis=args.axis,
        every=1,
        p2_cut=float(args.current_p2_cut),
        min_core_neighbors=1,
        cutoff_bins=60,
        cutoff_frames=1,
        r_cut_mode="manual",
        manual_r_cut=float(params.cutoff),
        mesogen_type=int(args.mesogen_type),
        anchor_types=tuple(int(item.strip()) for item in args.anchor_types.split(",") if item.strip()),
        contact_mode="gayberne",
        gayberne_params=params,
        gb_on_strength=float(args.current_gb_on),
        gb_off_strength=float(args.current_gb_off),
        robust_min_s2=float(args.current_s2_cut),
        r_energy_cap=r_energy_cap,
        s_excl=int(args.s_excl),
        local_pair_file=str(args.local_pair_file.resolve()) if args.local_pair_file else None,
        exclude_pair_file=str(args.exclude_pair_file.resolve()) if args.exclude_pair_file else None,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Collect unscreened Gay-Berne candidate pair diagnostics from LAMMPS dumps.")
    parser.add_argument("inputs", nargs="*", type=Path, help="Dump files or directories. Defaults to current directory.")
    parser.add_argument("--pattern", default="*.dump")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--gb-param-file", type=Path, required=True)
    parser.add_argument("--axis", default="auto", choices=["auto", "x", "y", "z"])
    parser.add_argument("--mesogen-type", type=int, default=1)
    parser.add_argument("--anchor-types", default="2,3")
    parser.add_argument("--s-excl", type=int, default=1)
    parser.add_argument("--local-pair-file", type=Path, default=None)
    parser.add_argument("--exclude-pair-file", type=Path, default=None)
    parser.add_argument("--r-energy-cap", default="auto")
    parser.add_argument("--every", type=int, default=1, help="Process every k-th frame in each dump.")
    parser.add_argument("--frame-limit", type=int, default=0, help="Maximum processed frames per dump; 0 means no limit.")
    parser.add_argument("--frame-sample-strategy", choices=["stratified", "first", "event-aware"], default="stratified", help="How --frame-limit selects frames. event-aware enriches high-change frames using cheap frame features.")
    parser.add_argument("--sample-pairs-per-frame", type=int, default=0, help="Deterministically sample at most N candidate pairs per frame; 0 means all.")
    parser.add_argument("--pair-sample-strategy", choices=["stratified", "uniform"], default="stratified", help="stratified preserves threshold, distance, strength, orientation, local/nonlocal, and excluded strata.")
    parser.add_argument("--sample-seed", type=int, default=20260429)
    parser.add_argument("--min-pairs-per-stratum", type=int, default=2)
    parser.add_argument(
        "--max-candidate-pairs-total",
        type=int,
        default=1_000_000,
        help="Maximum rows written to the combined gb_candidate_pairs.tsv used by calibration; 0 disables the cap.",
    )
    parser.add_argument("--workers", type=lca.parse_worker_count, default=1)
    parser.add_argument("--output-root", type=Path, default=Path("gb_pair_audit_output"))
    parser.add_argument("--current-gb-off", type=float, default=0.12)
    parser.add_argument("--current-gb-on", type=float, default=0.30)
    parser.add_argument("--current-p2-cut", type=float, default=0.70)
    parser.add_argument("--current-s2-cut", type=float, default=0.70)
    parser.add_argument("--n-min", type=int, default=3)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.every <= 0:
        raise SystemExit("--every must be positive")
    if args.frame_limit < 0:
        raise SystemExit("--frame-limit must be non-negative")
    if args.sample_pairs_per_frame < 0:
        raise SystemExit("--sample-pairs-per-frame must be non-negative")
    if args.min_pairs_per_stratum <= 0:
        raise SystemExit("--min-pairs-per-stratum must be positive")
    if args.max_candidate_pairs_total < 0:
        raise SystemExit("--max-candidate-pairs-total must be non-negative")

    input_paths = args.inputs if args.inputs else [Path.cwd()]
    files = iter_input_files(input_paths, pattern=args.pattern, recursive=bool(args.recursive))
    if not files:
        raise SystemExit("No dump files matched the given inputs.")
    config = build_config(args)
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    per_file_dir = output_root / "per_file"
    per_file_dir.mkdir(parents=True, exist_ok=True)

    worker_count = lca.resolve_worker_count(int(args.workers), len(files))
    all_rows: List[Dict[str, object]] = []
    all_frame_selection_rows: List[Dict[str, object]] = []
    raw_candidate_pairs = 0
    max_candidate_pairs_total = int(args.max_candidate_pairs_total)
    per_file_row_cap = 0
    if max_candidate_pairs_total > 0:
        per_file_row_cap = max(1, int(math.ceil(max_candidate_pairs_total / max(len(files), 1))))

    def process_job(job: Tuple[Path, List[Dict[str, object]], List[Dict[str, object]]]) -> None:
        nonlocal all_rows, raw_candidate_pairs
        dump_path, rows, frame_selection_rows = job
        raw_candidate_pairs += int(len(rows))
        all_rows = append_rows_with_rolling_cap(all_rows, rows, max_candidate_pairs_total)
        all_frame_selection_rows.extend(frame_selection_rows)
        safe = lca.build_output_stems([dump_path])[dump_path]
        per_file_rows = rows
        if per_file_row_cap > 0:
            per_file_rows, _ = cap_candidate_rows_for_calibration(rows, per_file_row_cap)
        write_rows(per_file_dir / f"{safe}_gb_candidate_pairs.tsv", per_file_rows)
        print(f"[OK] {dump_path.name}: candidate_pairs={len(rows)}, written_per_file={len(per_file_rows)}")

    job_args = [
        (
            path,
            config,
            args.every,
            args.frame_limit,
            args.sample_pairs_per_frame,
            args.frame_sample_strategy,
            args.pair_sample_strategy,
            args.sample_seed,
            args.min_pairs_per_stratum,
        )
        for path in files
    ]
    if worker_count <= 1 or len(files) <= 1:
        for job_arg in job_args:
            process_job(audit_dump_file_star(job_arg))
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            for job in executor.map(audit_dump_file_star, job_args):
                process_job(job)

    combined_path = output_root / "gb_candidate_pairs.tsv"
    frame_selection_path = output_root / "frame_selection_manifest.tsv"
    pair_sampling_path = output_root / "pair_sampling_audit.tsv"
    pair_sampling_manifest_path = output_root / "pair_sampling_manifest.tsv"
    capped_rows = list(all_rows)
    cap_summary = {
        "cap_applied": bool(max_candidate_pairs_total > 0 and raw_candidate_pairs > len(capped_rows)),
        "input_candidate_pairs": int(raw_candidate_pairs),
        "written_candidate_pairs": int(len(capped_rows)),
        "dropped_candidate_pairs": int(max(raw_candidate_pairs - len(capped_rows), 0)),
        "max_candidate_pairs_total": int(max_candidate_pairs_total),
        "selection_method": "rolling_deterministic_strided_across_ordered_files" if max_candidate_pairs_total > 0 else "none",
        "per_file_candidate_pair_cap": int(per_file_row_cap),
    }
    capped_rows = apply_global_candidate_cap_weights(capped_rows, cap_summary)
    if cap_summary.get("cap_applied"):
        print(
            "[WARN] candidate pair table capped: "
            f"{cap_summary['written_candidate_pairs']} / {cap_summary['input_candidate_pairs']} rows written to {combined_path.name}"
        )
    write_rows(combined_path, capped_rows)
    write_frame_selection_manifest(frame_selection_path, all_frame_selection_rows)
    write_pair_sampling_audit(pair_sampling_path, pair_sampling_audit_records(capped_rows))
    write_pair_sampling_manifest(pair_sampling_manifest_path, pair_sampling_manifest_records(capped_rows))
    calibration_rows = [row for row in capped_rows if safe_int(row.get("is_excluded")) != 1]
    if not calibration_rows and capped_rows:
        calibration_rows = list(capped_rows)
    recommendation = recommend_from_rows(
        calibration_rows,  # type: ignore[arg-type]
        [],
        current_gb_off=float(args.current_gb_off),
        current_gb_on=float(args.current_gb_on),
        current_p2_cut=float(args.current_p2_cut),
        current_s2_cut=float(args.current_s2_cut),
        n_min=int(args.n_min),
    )
    recommendation.setdefault("notes", [])
    recommendation["notes"].append(
        "GB pair audit calibrates pair thresholds only. robust_min_s2 is kept at the current value unless domain diagnostics are provided to lc_calibrate_thresholds.py."
    )
    recommendation["notes"].append(
        "This threshold_auto report uses the legacy unweighted recommender and is diagnostic when pair sampling is enabled; pipeline precision calibration uses sampling weights."
    )
    recommendation["input_kind"] = "legacy_unweighted_gb_pair_audit"
    recommendation["recommended"]["robust_min_s2"] = float(args.current_s2_cut)
    threshold_dir = output_root / "threshold_auto"
    threshold_dir.mkdir(exist_ok=True)
    (threshold_dir / "threshold_recommendations.json").write_text(
        json.dumps(recommendation, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_pair_bin_summary(capped_rows, threshold_dir / "candidate_pair_bin_summary.tsv")  # type: ignore[arg-type]
    write_plots(capped_rows, [], threshold_dir, recommendation)  # type: ignore[arg-type]
    summary = summarize_rows(capped_rows, recommendation)
    summary["candidate_pair_cap"] = cap_summary
    summary["raw_candidate_pairs_before_cap"] = int(len(all_rows))
    summary["threshold_recommendation"] = recommendation
    summary["frame_selection_manifest"] = str(frame_selection_path)
    summary["frame_selection_summary"] = summarize_frame_selection(all_frame_selection_rows)
    summary["pair_sampling_audit"] = str(pair_sampling_path)
    summary["pair_sampling_manifest"] = str(pair_sampling_manifest_path)
    manifest_path = output_root / "sampling_manifest.json"
    manifest = sampling_manifest.write_sampling_manifest(
        manifest_path,
        inputs=files,
        gb_param_file=args.gb_param_file,
        output_root=output_root,
        frame_strategy=args.frame_sample_strategy,
        frame_limit=int(args.frame_limit),
        every=int(args.every),
        pair_strategy=args.pair_sample_strategy,
        sample_pairs_per_frame=int(args.sample_pairs_per_frame),
        min_pairs_per_stratum=int(args.min_pairs_per_stratum),
        seed=int(args.sample_seed),
        thresholds={
            "gb_off_strength": float(args.current_gb_off),
            "gb_on_strength": float(args.current_gb_on),
            "p2_cut": float(args.current_p2_cut),
            "robust_min_s2": float(args.current_s2_cut),
            "n_min": int(args.n_min),
            "r_energy_cap": str(args.r_energy_cap),
            "max_candidate_pairs_total": int(args.max_candidate_pairs_total),
        },
        topology_files={
            "local_pair_file": args.local_pair_file,
            "exclude_pair_file": args.exclude_pair_file,
        },
        notes=[
            "gb_pair_audit collects unscreened GB candidate pairs for threshold calibration.",
            "Stratified pair sampling preserves near-threshold, strong, orientation, distance, local/nonlocal, and excluded strata when sampling is enabled.",
            "event-aware frame sampling is event-enriched coverage, not an unbiased random sample; see frame_selection_manifest.tsv.",
            "gb_candidate_pairs.tsv may be globally capped for calibration safety; see gb_pair_audit_summary.json candidate_pair_cap.",
        ],
    )
    summary["sampling_manifest"] = str(manifest_path)
    summary["sampling_manifest_sha256"] = sampling_manifest.sha256_file(manifest_path)
    summary["sampling_manifest_record"] = manifest
    (output_root / "gb_pair_audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"Candidate pair table: {combined_path}")
    print(f"Recommended thresholds: {threshold_dir / 'threshold_recommendations.json'}")


if __name__ == "__main__":
    main()
