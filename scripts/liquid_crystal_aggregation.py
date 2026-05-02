#!/usr/bin/env python3
"""Liquid-crystal domain and pearl-necklace analysis for LAMMPS dump files.

Algorithm v1 separates three questions that were previously mixed together:
support contacts define local LC bundles, robust domains require additional
evidence, and pearls are compact 3D bead-like groups of robust domains.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from lc_threshold_recommend import annotate_hline, annotate_vline, threshold_label


DEFAULT_MAX_WORKERS = 10


class AxisLength(float):
    """Box-axis length carrying whether minimum-image wrapping is valid."""

    def __new__(cls, value: float, periodic: bool = True) -> "AxisLength":
        obj = float.__new__(cls, value)
        obj.periodic = bool(periodic)
        return obj


@dataclass(frozen=True)
class BoxSpec:
    lengths: Tuple[float, float, float]
    bounds: Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]
    boundary_tokens: Tuple[str, str, str] = ("pp", "pp", "pp")

    @property
    def periodic(self) -> Tuple[bool, bool, bool]:
        return tuple(token.startswith("p") for token in self.boundary_tokens)  # type: ignore[return-value]

    def __getitem__(self, index: int) -> AxisLength:
        return AxisLength(self.lengths[index], self.periodic[index])

    def __iter__(self) -> Iterator[AxisLength]:
        for index in range(3):
            yield self[index]

    def __len__(self) -> int:
        return 3


@dataclass(frozen=True)
class GayBerneParams:
    gamma: float = 1.0
    upsilon: float = 3.0
    mu: float = 1.0
    cutoff: float = 5.0
    epsilon: float = 1.0
    sigma: float = 1.0
    eps_i: Tuple[float, float, float] = (1.0, 1.0, 0.2)
    eps_j: Tuple[float, float, float] = (1.0, 1.0, 0.2)


@dataclass(frozen=True)
class AggregationConfig:
    axis: str
    every: int
    p2_cut: float
    min_core_neighbors: int
    cutoff_bins: int
    cutoff_frames: int
    r_cut_mode: str
    manual_r_cut: Optional[float]
    max_auto_r_cut: Optional[float] = None
    auto_r_cut_shape_factor: float = 1.8
    mesogen_type: int = 1
    anchor_types: Tuple[int, ...] = (2, 3)
    contact_mode: str = "center"
    gayberne_params: Optional[GayBerneParams] = None
    gb_threshold_mode: str = "relative"
    gb_require_pair_orientation: bool = True
    u_on: float = -0.30
    u_off: float = -0.12
    gb_on_strength: float = 0.30
    gb_off_strength: float = 0.12
    r_energy_cap: Optional[float] = None
    g_on: float = 1.0
    g_off: float = 1.25
    q_on: Optional[float] = None
    q_off: Optional[float] = None
    s_excl: int = 1
    n_min: int = 3
    robust_min_size: int = 3
    robust_min_s2: float = 0.70
    robust_min_evidence: int = 4
    robust_require_orientation: bool = True
    robust_require_nonlocal: bool = False
    domain_min_lifetime: int = 2
    adjacent_id_gap: int = 1
    local_pair_file: Optional[str] = None
    exclude_pair_file: Optional[str] = None
    perturbation_r_cut_scale: float = 0.95
    perturbation_p2_margin: float = 0.05
    stable_overlap_fraction: float = 0.75
    pearl_gap_cut: Optional[float] = None
    pearl_min_cross_contacts: int = 2
    pearl_min_boundary_particles: int = 2
    pearl_max_aspect_ratio: float = 3.0
    track_jaccard: float = 0.50
    consensus_threshold: float = 0.70
    enable_robustness_scan: bool = False
    write_ovito_labels: bool = True
    write_cluster_envelopes: bool = True
    write_contact_edges: bool = True
    write_contact_segments: bool = True
    write_frame_jsonl: bool = False
    write_diagnostics: bool = True
    edge_diagnostics_table: str = "off"
    edge_diagnostics_sample_size: int = 200_000
    shared_r_cut: bool = True
    track_across_files: bool = True
    workers: int = 1
    cluster_cut: Optional[float] = None
    cluster_cut_shape_factor: float = 1.35
    cluster_min_size: int = 2
    cluster_envelope_padding: float = 0.4


@dataclass
class FrameAggregationResult:
    source_file: str
    timestep: int
    axis_used: str
    n_particles: int
    r_cut: float
    aggregation_degree: float
    largest_cluster_fraction: float
    clustered_fraction: float
    core_fraction: float
    mean_cluster_s2: float
    max_cluster_s2: float
    global_s2: float
    n_clusters: int
    qualified_pairs: int
    qualified_pair_fraction: float
    energy_edge_count: int
    min_pair_energy: float
    mean_pair_energy: float
    cluster_sizes: List[int]
    cluster_s2: List[float]
    visual_cluster_count: int
    largest_visual_cluster_fraction: float
    visual_cluster_sizes: List[int]
    weak_domain_count: int
    robust_domain_count: int
    weak_domain_fraction: float
    robust_domain_fraction: float
    weak_domain_sizes: List[int]
    robust_domain_sizes: List[int]
    robust_domain_s2: List[float]
    pearl_count: int
    largest_pearl_fraction: float
    pearl_sizes: List[int]
    pearl_domain_counts: List[int]
    pearl_compactness: List[float]
    pearl_aspect_ratios: List[float]
    pearl_axis_centers: List[float]
    pearl_axis_widths: List[float]
    connector_lengths: List[float]
    axial_axis: List[float]
    local_edge_fraction: float
    nonlocal_edge_fraction: float
    ambiguous_mesogen_fraction: float
    stable_core_pair_fraction: float
    l_parallel: float
    rg_parallel: float
    rg_perp: float
    s2_force: float
    stretch_axis: List[float]
    particle_contact_records: List[Dict[str, object]]
    contact_edge_records: List[Dict[str, object]]
    visual_cluster_records: List[Dict[str, object]]
    domain_records: List[Dict[str, object]]
    pearl_records: List[Dict[str, object]]
    pearl_candidate_records: List[Dict[str, object]]


@dataclass(frozen=True)
class SupportEdge:
    i: int
    j: int
    atom_i: int
    atom_j: int
    distance: float
    g_value: float
    q_score: float
    p2_score: float
    delta_s: int
    edge_type: str
    can_seed: bool
    is_local: bool
    pair_energy: Optional[float] = None
    well_depth: Optional[float] = None
    attraction_strength: Optional[float] = None
    contact_mode: str = "center"


@dataclass
class DomainCandidate:
    domain_id: int
    members: List[int]
    atom_ids: List[int]
    s2: float
    director: np.ndarray
    edge_count: int
    adjacent_edge_count: int
    nonlocal_edge_count: int
    stable_under_perturbation: bool
    track_id: int = -1
    age: int = 1
    evidence: Optional[Dict[str, bool]] = None
    evidence_count: int = 0
    classification: str = "weak"


@dataclass
class PearlCandidate:
    pearl_id: int
    domain_ids: List[int]
    members: List[int]
    atom_ids: List[int]
    center: np.ndarray
    radius_of_gyration: float
    max_radius: float
    compactness: float
    aspect_ratio: float


class DomainTracker:
    """Track domains across processed frames by atom-id overlap."""

    def __init__(self, min_jaccard: float) -> None:
        self.min_jaccard = float(min_jaccard)
        self._next_track_id = 1
        self._tracks: Dict[int, Tuple[set[int], int]] = {}

    def update(self, domains: Sequence[DomainCandidate]) -> None:
        used_tracks: set[int] = set()

        for domain in sorted(domains, key=lambda item: (-len(item.atom_ids), item.atom_ids)):
            current = set(domain.atom_ids)
            best_track: Optional[int] = None
            best_score = 0.0
            best_age = 0
            for track_id, (previous, age) in self._tracks.items():
                if track_id in used_tracks:
                    continue
                union = current | previous
                if not union:
                    continue
                score = len(current & previous) / len(union)
                if score > best_score:
                    best_score = score
                    best_track = track_id
                    best_age = age

            if best_track is not None and best_score >= self.min_jaccard:
                domain.track_id = best_track
                domain.age = best_age + 1
                used_tracks.add(best_track)
            else:
                domain.track_id = self._next_track_id
                domain.age = 1
                self._next_track_id += 1

        self._tracks = {
            domain.track_id: (set(domain.atom_ids), int(domain.age))
            for domain in domains
            if domain.track_id > 0
        }


def axis_from_quat(w: float, x: float, y: float, z: float, axis: str) -> Tuple[float, float, float]:
    """Rotate a body-frame basis axis into the lab frame using quaternion q=(w,x,y,z)."""
    if axis == "x":
        nx = 1.0 - 2.0 * (y * y + z * z)
        ny = 2.0 * (x * y + w * z)
        nz = 2.0 * (x * z - w * y)
    elif axis == "y":
        nx = 2.0 * (x * y - w * z)
        ny = 1.0 - 2.0 * (x * x + z * z)
        nz = 2.0 * (y * z + w * x)
    elif axis == "z":
        nx = 2.0 * (x * z + w * y)
        ny = 2.0 * (y * z - w * x)
        nz = 1.0 - 2.0 * (x * x + y * y)
    else:
        raise ValueError("axis must be one of {'x', 'y', 'z'}")

    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm == 0.0:
        return 0.0, 0.0, 0.0
    return nx / norm, ny / norm, nz / norm


def quat_to_mat_trans_lammps(quat: np.ndarray) -> np.ndarray:
    q = np.asarray(quat, dtype=float)
    norm = float(np.linalg.norm(q))
    if norm <= 1e-12:
        q = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    else:
        q = q / norm
    w, i, j, k = (float(value) for value in q)
    w2 = w * w
    i2 = i * i
    j2 = j * j
    k2 = k * k
    twoij = 2.0 * i * j
    twoik = 2.0 * i * k
    twojk = 2.0 * j * k
    twoiw = 2.0 * i * w
    twojw = 2.0 * j * w
    twokw = 2.0 * k * w
    return np.array(
        [
            [w2 + i2 - j2 - k2, twoij + twokw, twoik - twojw],
            [twoij - twokw, w2 - i2 + j2 - k2, twojk + twoiw],
            [twojw + twoik, twojk - twoiw, w2 - i2 - j2 + k2],
        ],
        dtype=float,
    )


def p2(values: np.ndarray) -> np.ndarray:
    return 0.5 * (3.0 * values * values - 1.0)


def compute_Q_and_S(u: np.ndarray) -> Tuple[float, np.ndarray]:
    """Return the largest eigenvalue S2 and its eigenvector for unit vectors."""
    if u.size == 0:
        return 0.0, np.array([0.0, 0.0, 1.0], dtype=float)
    uu = (u[:, :, None] * u[:, None, :]).mean(axis=0)
    q_tensor = 1.5 * uu - 0.5 * np.eye(3)
    evals, evecs = np.linalg.eigh(q_tensor)
    idx = int(np.argmax(evals))
    s2 = float(evals[idx])
    director = evecs[:, idx]
    if director[2] < 0.0:
        director = -director
    return s2, director


def minimum_image(dx: np.ndarray, box_length: float) -> np.ndarray:
    if not getattr(box_length, "periodic", True):
        return dx
    return dx - box_length * np.round(dx / box_length)


def minimum_image_vector(delta: np.ndarray, box: Tuple[float, float, float]) -> np.ndarray:
    adjusted = np.array(delta, dtype=float, copy=True)
    adjusted[0] = minimum_image(adjusted[0], box[0])
    adjusted[1] = minimum_image(adjusted[1], box[1])
    adjusted[2] = minimum_image(adjusted[2], box[2])
    return adjusted


def pbc_distance(a: np.ndarray, b: np.ndarray, box: Tuple[float, float, float]) -> float:
    return float(np.linalg.norm(minimum_image_vector(np.asarray(b, dtype=float) - np.asarray(a, dtype=float), box)))


def box_axis_bounds(box: Tuple[float, float, float] | BoxSpec, index: int) -> Tuple[float, float]:
    if isinstance(box, BoxSpec):
        return box.bounds[index]
    return 0.0, float(box[index])


def write_box_bounds(handle, box: Tuple[float, float, float] | BoxSpec) -> None:
    if isinstance(box, BoxSpec):
        handle.write(f"ITEM: BOX BOUNDS {' '.join(box.boundary_tokens)}\n")
        for lower, upper in box.bounds:
            handle.write(f"{lower:.12g} {upper:.12g}\n")
        return
    handle.write("ITEM: BOX BOUNDS pp pp pp\n")
    handle.write(f"0.0 {box[0]:.12g}\n")
    handle.write(f"0.0 {box[1]:.12g}\n")
    handle.write(f"0.0 {box[2]:.12g}\n")


def parse_dump_frames(dump_path: Path) -> Iterator[Tuple[int, BoxSpec, List[str], Dict[str, int], np.ndarray]]:
    """Yield frames from a text LAMMPS dump with an orthorhombic box."""
    with dump_path.open("r", encoding="utf-8", errors="ignore") as handle:
        while True:
            line = handle.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(handle.readline().strip())

            line = handle.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError(f"{dump_path}: expected 'ITEM: NUMBER OF ATOMS' near timestep {timestep}")
            n_atoms = int(handle.readline().strip())

            line = handle.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError(f"{dump_path}: expected 'ITEM: BOX BOUNDS' near timestep {timestep}")
            boundary_tokens = tuple(line.strip().split()[3:6])
            if len(boundary_tokens) != 3:
                boundary_tokens = ("pp", "pp", "pp")
            xlo, xhi = map(float, handle.readline().split()[:2])
            ylo, yhi = map(float, handle.readline().split()[:2])
            zlo, zhi = map(float, handle.readline().split()[:2])
            box = BoxSpec(
                lengths=(xhi - xlo, yhi - ylo, zhi - zlo),
                bounds=((xlo, xhi), (ylo, yhi), (zlo, zhi)),
                boundary_tokens=boundary_tokens,
            )

            line = handle.readline()
            if not line.startswith("ITEM: ATOMS"):
                raise RuntimeError(f"{dump_path}: expected 'ITEM: ATOMS' near timestep {timestep}")
            columns = line.strip().split()[2:]
            col_index = {name: idx for idx, name in enumerate(columns)}

            data = np.empty((n_atoms, len(columns)), dtype=float)
            for i in range(n_atoms):
                parts = handle.readline().split()
                if len(parts) < len(columns):
                    raise RuntimeError(
                        f"{dump_path}: malformed atom row at timestep {timestep}, line {i + 1}: "
                        f"got {len(parts)} fields, expected {len(columns)}"
                    )
                data[i, :] = np.array(parts[: len(columns)], dtype=float)

            yield timestep, box, columns, col_index, data


def strip_lammps_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def parse_lammps_variables(lmp_path: Path) -> Dict[str, str]:
    variables: Dict[str, str] = {}
    for raw_line in lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = strip_lammps_comment(raw_line)
        if not line:
            continue
        parts = line.split(None, 3)
        if len(parts) >= 4 and parts[0] == "variable" and parts[2] in {"equal", "index"}:
            variables[parts[1]] = parts[3].strip().strip("\"'")
    return variables


def resolve_lammps_token(token: str, variables: Dict[str, str]) -> str:
    value = token.strip()
    pattern = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")
    previous = None
    while previous != value:
        previous = value
        value = pattern.sub(lambda match: variables.get(match.group(1), match.group(0)), value)
    return value


def evaluate_lammps_numeric_expression(expression: str, variables: Dict[str, str], stack: Optional[set[str]] = None) -> float:
    stack = set() if stack is None else set(stack)

    def resolve_named_variable(name: str) -> str:
        if name in stack:
            raise ValueError(f"recursive LAMMPS variable reference: {name}")
        if name not in variables:
            raise ValueError(f"unknown LAMMPS variable: {name}")
        value = variables[name]
        try:
            return str(evaluate_lammps_numeric_expression(value, variables, stack | {name}))
        except ValueError:
            return resolve_lammps_token(value, variables)

    resolved = resolve_lammps_token(expression, variables)
    resolved = re.sub(r"\bv_([A-Za-z_][A-Za-z0-9_]*)\b", lambda match: resolve_named_variable(match.group(1)), resolved)
    resolved = resolve_lammps_token(resolved, variables)
    python_expr = resolved.replace("^", "**")
    if not re.fullmatch(r"[0-9eE+\-*/().,_ A-Za-z]+", python_expr):
        raise ValueError(f"unsupported LAMMPS numeric expression: {expression}")
    allowed = {
        "sqrt": math.sqrt,
        "exp": math.exp,
        "log": math.log,
        "sin": math.sin,
        "cos": math.cos,
        "tan": math.tan,
        "abs": abs,
        "min": min,
        "max": max,
        "PI": math.pi,
        "pi": math.pi,
    }
    return float(eval(python_expr, {"__builtins__": {}}, allowed))


def lammps_token_float(token: str, variables: Dict[str, str]) -> float:
    return evaluate_lammps_numeric_expression(token, variables)


def parse_gayberne_params_from_lmp(lmp_path: Path) -> GayBerneParams:
    variables = parse_lammps_variables(lmp_path)
    gamma = upsilon = mu = cutoff = None
    epsilon = sigma = None
    eps_i: Optional[Tuple[float, float, float]] = None
    eps_j: Optional[Tuple[float, float, float]] = None

    for raw_line in lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = strip_lammps_comment(raw_line)
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 6 and parts[0] == "pair_style":
            for idx, token in enumerate(parts):
                if token == "gayberne" and idx + 4 < len(parts):
                    gamma = lammps_token_float(parts[idx + 1], variables)
                    upsilon = lammps_token_float(parts[idx + 2], variables)
                    mu = lammps_token_float(parts[idx + 3], variables)
                    cutoff = lammps_token_float(parts[idx + 4], variables)
                    break
        if len(parts) >= 11 and parts[0] == "pair_coeff" and parts[1] == "1" and parts[2] == "1":
            coeff_offset = 4 if len(parts) >= 12 and parts[3] == "gayberne" else 3
            if len(parts) < coeff_offset + 8:
                continue
            epsilon = lammps_token_float(parts[coeff_offset], variables)
            sigma = lammps_token_float(parts[coeff_offset + 1], variables)
            eps_i = (
                lammps_token_float(parts[coeff_offset + 2], variables),
                lammps_token_float(parts[coeff_offset + 3], variables),
                lammps_token_float(parts[coeff_offset + 4], variables),
            )
            eps_j = (
                lammps_token_float(parts[coeff_offset + 5], variables),
                lammps_token_float(parts[coeff_offset + 6], variables),
                lammps_token_float(parts[coeff_offset + 7], variables),
            )
            if len(parts) > coeff_offset + 8:
                cutoff = lammps_token_float(parts[coeff_offset + 8], variables)

    missing = [
        name
        for name, value in (
            ("gamma", gamma),
            ("upsilon", upsilon),
            ("mu", mu),
            ("cutoff", cutoff),
            ("epsilon", epsilon),
            ("sigma", sigma),
            ("eps_i", eps_i),
            ("eps_j", eps_j),
        )
        if value is None
    ]
    if missing:
        raise RuntimeError(f"{lmp_path}: missing Gay-Berne parameters: {missing}")
    return GayBerneParams(
        gamma=float(gamma),
        upsilon=0.5 * float(upsilon),
        mu=float(mu),
        cutoff=float(cutoff),
        epsilon=float(epsilon),
        sigma=float(sigma),
        eps_i=eps_i,
        eps_j=eps_j,
    )


def sanitize_output_name(path: Path | str) -> str:
    path_obj = Path(path)
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", path_obj.name).strip("_") or "dump"
    digest = hashlib.sha1(str(path_obj).encode("utf-8")).hexdigest()[:8]
    return f"{base}_{digest}"


def sanitize_output_stem(path: Path | str) -> str:
    path_obj = Path(path)
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", path_obj.stem).strip("_")
    return stem or "dump"


def build_output_stems(paths: Sequence[Path]) -> Dict[Path, str]:
    """Build output stems without hashes unless duplicate input stems collide."""
    base_counts: Dict[str, int] = {}
    for path in paths:
        base = sanitize_output_stem(path)
        base_counts[base] = base_counts.get(base, 0) + 1

    stems: Dict[Path, str] = {}
    for path in paths:
        base = sanitize_output_stem(path)
        if base_counts[base] == 1:
            stems[path] = base
            continue
        digest = hashlib.sha1(str(path).encode("utf-8")).hexdigest()[:8]
        stems[path] = f"{base}_{digest}"
    return stems


def dump_sort_key(path: Path | str) -> Tuple[int, object, str]:
    """Sort one-frame dump sequences by the last integer in the filename stem."""
    path_obj = Path(path)
    matches = re.findall(r"\d+", path_obj.stem)
    if matches:
        return 0, int(matches[-1]), str(path_obj)
    return 1, path_obj.stem, str(path_obj)


def resolve_axis(data: np.ndarray, col_index: Dict[str, int], axis_choice: str) -> str:
    if axis_choice != "auto":
        return axis_choice
    if all(key in col_index for key in ("shapex", "shapey", "shapez")):
        sx = float(np.median(data[:, col_index["shapex"]]))
        sy = float(np.median(data[:, col_index["shapey"]]))
        sz = float(np.median(data[:, col_index["shapez"]]))
        if sx >= sy and sx >= sz:
            return "x"
        if sy >= sx and sy >= sz:
            return "y"
        return "z"
    return "z"


def extract_positions(data: np.ndarray, col_index: Dict[str, int]) -> np.ndarray:
    if all(key in col_index for key in ("xu", "yu", "zu")):
        return np.column_stack(
            (
                data[:, col_index["xu"]],
                data[:, col_index["yu"]],
                data[:, col_index["zu"]],
            )
        )

    required = ("x", "y", "z")
    missing = [key for key in required if key not in col_index]
    if missing:
        raise RuntimeError(f"missing required ATOMS columns for positions: {missing}; prefer xu/yu/zu when available")
    return np.column_stack(
        (
            data[:, col_index["x"]],
            data[:, col_index["y"]],
            data[:, col_index["z"]],
        )
    )


def extract_particle_ids(data: np.ndarray, col_index: Dict[str, int]) -> np.ndarray:
    if "id" in col_index:
        return data[:, col_index["id"]].astype(int)
    return np.arange(data.shape[0], dtype=int)


def extract_particle_types(data: np.ndarray, col_index: Dict[str, int]) -> np.ndarray:
    if "type" in col_index:
        return data[:, col_index["type"]].astype(int)
    return np.ones(data.shape[0], dtype=int)


def extract_quaternions(data: np.ndarray, col_index: Dict[str, int]) -> np.ndarray:
    required = ("quatw", "quati", "quatj", "quatk")
    missing = [key for key in required if key not in col_index]
    if missing:
        raise RuntimeError(f"missing required ATOMS columns for quaternions: {missing}")
    return np.column_stack(
        (
            data[:, col_index["quatw"]],
            data[:, col_index["quati"]],
            data[:, col_index["quatj"]],
            data[:, col_index["quatk"]],
        )
    )


def extract_shape_axes(data: np.ndarray, col_index: Dict[str, int]) -> Optional[np.ndarray]:
    if all(key in col_index for key in ("shapex", "shapey", "shapez")):
        return np.column_stack(
            (
                data[:, col_index["shapex"]],
                data[:, col_index["shapey"]],
                data[:, col_index["shapez"]],
            )
        )
    return None


def select_mesogen_indices(types: np.ndarray, mesogen_type: int) -> np.ndarray:
    return np.nonzero(types == int(mesogen_type))[0]


def infer_chain_indices(
    data: np.ndarray,
    col_index: Dict[str, int],
    mesogen_indices: np.ndarray,
    particle_ids: np.ndarray,
    positions: Optional[np.ndarray] = None,
    types: Optional[np.ndarray] = None,
    anchor_types: Tuple[int, ...] = (2, 3),
) -> np.ndarray:
    """Infer sequence coordinate s_i for mesogens.

    Explicit dump columns are preferred. Without topology, the conservative
    fallback is id order within molecule, then global id order.
    """
    for key in ("s", "chain_s", "mesogen_s", "mesogen_index", "seq"):
        if key in col_index:
            return data[mesogen_indices, col_index[key]].astype(int)

    mesogen_ids = particle_ids[mesogen_indices].astype(int)
    s_values = np.zeros(mesogen_indices.size, dtype=int)
    if "mol" in col_index:
        if positions is not None and types is not None:
            all_mol_values = data[:, col_index["mol"]].astype(int)
            type_array = np.asarray(types, dtype=int)
            anchor_type_array = np.array(anchor_types, dtype=int)
            filled = np.zeros(mesogen_indices.size, dtype=bool)
            for mol in sorted(set(int(item) for item in data[mesogen_indices, col_index["mol"]].astype(int))):
                local = np.nonzero(data[mesogen_indices, col_index["mol"]].astype(int) == mol)[0]
                molecule_anchor_mask = (all_mol_values == mol) & np.isin(type_array, anchor_type_array)
                anchor_pos = positions[molecule_anchor_mask, :]
                if anchor_pos.shape[0] < 2:
                    continue
                axis = infer_stretch_axis(anchor_pos, np.full(anchor_pos.shape[0], anchor_types[0], dtype=int), (anchor_types[0],))
                if axis is None:
                    continue
                mesogen_pos = positions[mesogen_indices[local], :]
                origin = anchor_pos[np.argmin(anchor_pos @ axis), :]
                projection = (mesogen_pos - origin) @ axis
                ordered = local[np.argsort(projection)]
                for rank, local_idx in enumerate(ordered):
                    s_values[int(local_idx)] = rank
                    filled[int(local_idx)] = True
            if bool(np.all(filled)):
                return s_values

        mol_values = data[mesogen_indices, col_index["mol"]].astype(int)
        for mol in sorted(set(int(item) for item in mol_values)):
            local = np.nonzero(mol_values == mol)[0]
            ordered = local[np.argsort(mesogen_ids[local])]
            for rank, local_idx in enumerate(ordered):
                s_values[int(local_idx)] = rank
        return s_values

    ordered = np.argsort(mesogen_ids)
    for rank, local_idx in enumerate(ordered):
        s_values[int(local_idx)] = rank
    return s_values


def extract_orientations(data: np.ndarray, col_index: Dict[str, int], axis_choice: str) -> Tuple[np.ndarray, str]:
    required = ("quatw", "quati", "quatj", "quatk")
    missing = [key for key in required if key not in col_index]
    if missing:
        raise RuntimeError(f"missing required ATOMS columns for quaternions: {missing}")

    axis_used = resolve_axis(data, col_index, axis_choice)
    u = np.zeros((data.shape[0], 3), dtype=float)
    w = data[:, col_index["quatw"]]
    x = data[:, col_index["quati"]]
    y = data[:, col_index["quatj"]]
    z = data[:, col_index["quatk"]]

    for i in range(data.shape[0]):
        u[i, :] = axis_from_quat(float(w[i]), float(x[i]), float(y[i]), float(z[i]), axis_used)
    return u, axis_used


def compute_gr_like_curve(
    pos: np.ndarray,
    box: Tuple[float, float, float],
    nbins: int,
    rmax: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute a simple g(r)-like curve for first-shell cutoff estimation."""
    n_particles = pos.shape[0]
    lx, ly, lz = box
    volume = float(lx) * float(ly) * float(lz)
    rho = n_particles / volume
    if rmax is None:
        rmax = 0.5 * min(float(value) for value in box)
    if nbins <= 0:
        raise ValueError("nbins must be positive")

    dr = rmax / nbins
    counts = np.zeros(nbins, dtype=float)
    for i in range(n_particles - 1):
        delta = pos[i + 1 :, :] - pos[i, :]
        delta[:, 0] = minimum_image(delta[:, 0], lx)
        delta[:, 1] = minimum_image(delta[:, 1], ly)
        delta[:, 2] = minimum_image(delta[:, 2], lz)
        distances = np.linalg.norm(delta, axis=1)
        mask = (distances > 1e-12) & (distances < rmax)
        if not np.any(mask):
            continue
        bins = np.floor(distances[mask] / dr).astype(int)
        for bin_idx in bins:
            if 0 <= bin_idx < nbins:
                counts[bin_idx] += 1.0

    r_edges = np.linspace(0.0, rmax, nbins + 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    shell_volume = (4.0 / 3.0) * math.pi * (r_edges[1:] ** 3 - r_edges[:-1] ** 3)
    ideal = 0.5 * n_particles * rho * shell_volume
    g_curve = np.zeros_like(counts)
    valid = ideal > 0.0
    g_curve[valid] = counts[valid] / ideal[valid]
    return r_centers, g_curve


def smooth_curve(values: np.ndarray, window: int = 5) -> np.ndarray:
    if window <= 1 or values.size < window:
        return values.copy()
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values, kernel, mode="same")


def estimate_first_shell_cutoff(r_values: np.ndarray, g_values: np.ndarray) -> Optional[float]:
    """Pick the first minimum after the first visible peak in a smoothed g(r) curve."""
    if r_values.size < 5 or g_values.size != r_values.size:
        return None

    smooth = smooth_curve(g_values, window=5)
    peak_idx: Optional[int] = None
    for idx in range(1, smooth.size - 1):
        if smooth[idx] > smooth[idx - 1] and smooth[idx] >= smooth[idx + 1] and smooth[idx] > 1.0:
            peak_idx = idx
            break
    if peak_idx is None:
        peak_idx = int(np.argmax(smooth))
    if peak_idx >= smooth.size - 2:
        return None

    for idx in range(peak_idx + 1, g_values.size - 1):
        if g_values[idx] <= g_values[idx - 1] and g_values[idx] <= g_values[idx + 1]:
            return float(r_values[idx])

    below_one = np.where(smooth[peak_idx + 1 :] <= 1.0)[0]
    if below_one.size:
        return float(r_values[peak_idx + 1 + int(below_one[0])])
    return None


def auto_r_cut_cap_from_frame(
    data: np.ndarray,
    col_index: Dict[str, int],
    mesogen_indices: np.ndarray,
    config: AggregationConfig,
) -> Optional[float]:
    if config.max_auto_r_cut is not None:
        return float(config.max_auto_r_cut)

    shapes = extract_shape_axes(data, col_index)
    if shapes is None or mesogen_indices.size == 0:
        return None
    mesogen_shapes = shapes[mesogen_indices, :]
    long_diameters = np.max(mesogen_shapes, axis=1)
    long_diameters = long_diameters[np.isfinite(long_diameters) & (long_diameters > 0.0)]
    if long_diameters.size == 0:
        return None
    return float(np.median(long_diameters) * config.auto_r_cut_shape_factor)


def apply_auto_r_cut_cap(candidate: float, cap: Optional[float]) -> float:
    if cap is None or cap <= 0.0:
        return float(candidate)
    return float(min(float(candidate), float(cap)))


def determine_r_cut_for_dump(dump_path: Path, config: AggregationConfig) -> float:
    if config.r_cut_mode == "manual":
        if config.manual_r_cut is None or config.manual_r_cut <= 0.0:
            raise ValueError("manual r_cut must be a positive number")
        return float(config.manual_r_cut)

    cutoff_candidates: List[float] = []
    sampled = 0
    frame_index = 0
    for _timestep, box, _columns, col_index, data in parse_dump_frames(dump_path):
        frame_index += 1
        if (frame_index - 1) % config.every != 0:
            continue
        pos = extract_positions(data, col_index)
        types = extract_particle_types(data, col_index)
        mesogen_indices = select_mesogen_indices(types, config.mesogen_type)
        if mesogen_indices.size:
            pos = pos[mesogen_indices, :]
        r_values, g_values = compute_gr_like_curve(pos, box, nbins=config.cutoff_bins)
        candidate = estimate_first_shell_cutoff(r_values, g_values)
        if candidate is not None and candidate > 0.0:
            cap = auto_r_cut_cap_from_frame(data, col_index, mesogen_indices, config)
            cutoff_candidates.append(apply_auto_r_cut_cap(float(candidate), cap))
        sampled += 1
        if sampled >= config.cutoff_frames:
            break

    if not cutoff_candidates:
        raise RuntimeError(
            f"{dump_path}: automatic r_cut estimation failed. "
            "Please pass a numeric value with --r-cut."
        )
    return float(np.median(np.array(cutoff_candidates, dtype=float)))


def estimate_cutoff_from_frame(
    pos: np.ndarray,
    box: Tuple[float, float, float],
    config: AggregationConfig,
    cap: Optional[float] = None,
) -> Optional[float]:
    r_values, g_values = compute_gr_like_curve(pos, box, nbins=config.cutoff_bins)
    candidate = estimate_first_shell_cutoff(r_values, g_values)
    if candidate is not None and candidate > 0.0:
        return apply_auto_r_cut_cap(float(candidate), cap)
    return None


def determine_shared_r_cut_for_files(files: Sequence[Path], config: AggregationConfig) -> float:
    if config.r_cut_mode == "manual":
        if config.manual_r_cut is None or config.manual_r_cut <= 0.0:
            raise ValueError("manual r_cut must be a positive number")
        return float(config.manual_r_cut)

    cutoff_candidates: List[float] = []
    sampled = 0
    for dump_path in files:
        frame_index = 0
        for _timestep, box, _columns, col_index, data in parse_dump_frames(dump_path):
            frame_index += 1
            if (frame_index - 1) % config.every != 0:
                continue
            pos = extract_positions(data, col_index)
            types = extract_particle_types(data, col_index)
            mesogen_indices = select_mesogen_indices(types, config.mesogen_type)
            if mesogen_indices.size:
                pos = pos[mesogen_indices, :]
            cap = auto_r_cut_cap_from_frame(data, col_index, mesogen_indices, config)
            candidate = estimate_cutoff_from_frame(pos, box, config, cap=cap)
            if candidate is not None:
                cutoff_candidates.append(candidate)
            sampled += 1
            if sampled >= config.cutoff_frames:
                break
        if sampled >= config.cutoff_frames:
            break

    if not cutoff_candidates:
        raise RuntimeError("automatic shared r_cut estimation failed. Please pass a numeric value with --r-cut.")
    return float(np.median(np.array(cutoff_candidates, dtype=float)))


def normalized_pair(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a <= b else (b, a)


def build_explicit_adjacent_pairs(adjacent_pairs: Optional[Iterable[Tuple[int, int]]]) -> set[Tuple[int, int]]:
    if adjacent_pairs is None:
        return set()
    return {normalized_pair(int(a), int(b)) for a, b in adjacent_pairs}


def read_pair_list(path: Optional[str | Path]) -> set[Tuple[int, int]]:
    if path is None:
        return set()
    pair_path = Path(path)
    if not pair_path.exists():
        raise FileNotFoundError(pair_path)
    pairs: set[Tuple[int, int]] = set()
    with pair_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\s,]+", line)
            if len(parts) < 2 or parts[0].lower() in {"atom_i", "i", "id1"}:
                continue
            try:
                atom_i = int(float(parts[0]))
                atom_j = int(float(parts[1]))
            except ValueError:
                continue
            pairs.add(normalized_pair(atom_i, atom_j))
    return pairs


def q_threshold_from_p2(p2_cutoff: float) -> float:
    return float(math.sqrt(max(0.0, (2.0 * float(p2_cutoff) + 1.0) / 3.0)))


def resolved_q_thresholds(config: AggregationConfig, p2_cutoff: float) -> Tuple[float, float]:
    q_on = float(config.q_on) if config.q_on is not None else q_threshold_from_p2(p2_cutoff)
    q_off = float(config.q_off) if config.q_off is not None else max(0.0, q_on - 0.10)
    return q_on, q_off


def is_topologically_adjacent(
    atom_i: int,
    atom_j: int,
    explicit_adjacent_pairs: set[Tuple[int, int]],
    adjacent_id_gap: int,
) -> bool:
    pair = normalized_pair(int(atom_i), int(atom_j))
    if pair in explicit_adjacent_pairs:
        return True
    return adjacent_id_gap > 0 and abs(int(atom_i) - int(atom_j)) == adjacent_id_gap


def contact_scale_from_shape(
    shape_i: Optional[np.ndarray],
    shape_j: Optional[np.ndarray],
    u_i: np.ndarray,
    u_j: np.ndarray,
    r_hat: np.ndarray,
    fallback: float,
) -> float:
    """Approximate ellipsoid contact distance along the center-center vector.

    If shape columns are unavailable this intentionally collapses to the
    geometry-only version, g_ij = r_ij / r_cut.
    """
    if shape_i is None or shape_j is None:
        return float(fallback)

    def radius(shape: np.ndarray, u_vec: np.ndarray, sign: float) -> float:
        long_diameter = float(np.max(shape))
        transverse = float(np.median(np.sort(shape)[:2]))
        a = max(0.5 * long_diameter, 1e-12)
        b = max(0.5 * transverse, 1e-12)
        cos_theta = float(abs(np.dot(u_vec, sign * r_hat)))
        sin2 = max(0.0, 1.0 - cos_theta * cos_theta)
        return 1.0 / math.sqrt((cos_theta * cos_theta) / (a * a) + sin2 / (b * b))

    return max(radius(shape_i, u_i, 1.0) + radius(shape_j, u_j, -1.0), 1e-12)


def gayberne_pair_metrics(
    r12: np.ndarray,
    quat_i: np.ndarray,
    quat_j: np.ndarray,
    shape_i: np.ndarray,
    shape_j: np.ndarray,
    params: GayBerneParams,
) -> Optional[Tuple[float, float]]:
    """Return LAMMPS pair_gayberne energy and orientation-specific well depth."""
    r_vec = np.asarray(r12, dtype=float)
    rsq = float(np.dot(r_vec, r_vec))
    if rsq <= 1e-24 or rsq >= float(params.cutoff) * float(params.cutoff):
        return None
    r = math.sqrt(rsq)
    rhat = r_vec / r

    a1 = quat_to_mat_trans_lammps(np.asarray(quat_i, dtype=float))
    a2 = quat_to_mat_trans_lammps(np.asarray(quat_j, dtype=float))
    # LAMMPS dump shapex/shapey/shapez are diameters; pair_gayberne uses internal semiaxes.
    shape_i = 0.5 * np.asarray(shape_i, dtype=float)
    shape_j = 0.5 * np.asarray(shape_j, dtype=float)
    shape2_i = shape_i * shape_i
    shape2_j = shape_j * shape_j
    well_i = np.power(np.asarray(params.eps_i, dtype=float), -1.0 / float(params.mu))
    well_j = np.power(np.asarray(params.eps_j, dtype=float), -1.0 / float(params.mu))

    g1 = a1.T @ (np.diag(shape2_i) @ a1)
    g2 = a2.T @ (np.diag(shape2_j) @ a2)
    g12 = g1 + g2
    try:
        kappa = np.linalg.solve(g12, r_vec)
    except np.linalg.LinAlgError:
        return None
    sigma_term = float(np.dot(rhat, kappa / r))
    if sigma_term <= 0.0:
        return None
    sigma12 = math.pow(0.5 * sigma_term, -0.5)
    h12 = r - sigma12
    denominator = h12 + float(params.gamma) * float(params.sigma)
    if abs(denominator) <= 1e-14:
        return None
    varrho = float(params.sigma) / denominator
    varrho6 = math.pow(varrho, 6.0)
    varrho12 = varrho6 * varrho6
    u_r = 4.0 * float(params.epsilon) * (varrho12 - varrho6)

    lshape_i = (float(shape_i[0]) * float(shape_i[1]) + float(shape_i[2]) * float(shape_i[2])) * math.sqrt(
        max(float(shape_i[0]) * float(shape_i[1]), 1e-24)
    )
    lshape_j = (float(shape_j[0]) * float(shape_j[1]) + float(shape_j[2]) * float(shape_j[2])) * math.sqrt(
        max(float(shape_j[0]) * float(shape_j[1]), 1e-24)
    )
    det_g12 = float(np.linalg.det(g12))
    if det_g12 <= 0.0:
        return None
    eta = math.pow((2.0 * lshape_i * lshape_j) / det_g12, float(params.upsilon))

    b1 = a1.T @ (np.diag(well_i) @ a1)
    b2 = a2.T @ (np.diag(well_j) @ a2)
    b12 = b1 + b2
    try:
        iota = np.linalg.solve(b12, r_vec)
    except np.linalg.LinAlgError:
        return None
    chi_base = 2.0 * float(np.dot(rhat, iota / r))
    if chi_base < 0.0 and abs(float(params.mu) - round(float(params.mu))) > 1e-12:
        return None
    chi = math.pow(chi_base, float(params.mu))
    well_depth = float(params.epsilon) * eta * chi
    if well_depth <= 0.0:
        return None
    return float(u_r * eta * chi), float(well_depth)


def gayberne_pair_energy(
    r12: np.ndarray,
    quat_i: np.ndarray,
    quat_j: np.ndarray,
    shape_i: np.ndarray,
    shape_j: np.ndarray,
    params: GayBerneParams,
) -> Optional[float]:
    """LAMMPS pair_gayberne ellipsoid-ellipsoid energy for one type 1-1 pair."""
    metrics = gayberne_pair_metrics(
        r12=r12,
        quat_i=quat_i,
        quat_j=quat_j,
        shape_i=shape_i,
        shape_j=shape_j,
        params=params,
    )
    if metrics is None:
        return None
    return metrics[0]


def build_support_graph(
    pos: np.ndarray,
    u: np.ndarray,
    box: Tuple[float, float, float],
    r_cut: float,
    p2_cutoff: float,
    particle_ids: Optional[np.ndarray] = None,
    chain_indices: Optional[np.ndarray] = None,
    shape_axes: Optional[np.ndarray] = None,
    quaternions: Optional[np.ndarray] = None,
    config: Optional[AggregationConfig] = None,
    adjacent_pairs: Optional[Iterable[Tuple[int, int]]] = None,
    excluded_pairs: Optional[Iterable[Tuple[int, int]]] = None,
    adjacent_id_gap: int = 1,
) -> Tuple[List[List[int]], List[SupportEdge]]:
    """Build strong/gray/local support contacts from center-distance and alignment."""
    n_particles = pos.shape[0]
    ids = np.asarray(particle_ids if particle_ids is not None else np.arange(n_particles), dtype=int)
    chain_indices_provided = chain_indices is not None
    s_values = np.asarray(chain_indices if chain_indices is not None else np.arange(n_particles), dtype=int)
    effective_config = config or make_default_config("z", r_cut, p2_cutoff, 1)
    q_on, q_off = resolved_q_thresholds(effective_config, p2_cutoff)
    explicit_pairs = build_explicit_adjacent_pairs(adjacent_pairs)
    excluded_pair_set = build_explicit_adjacent_pairs(excluded_pairs)
    lx, ly, lz = box
    adjacency: List[List[int]] = [[] for _ in range(n_particles)]
    edges: List[SupportEdge] = []
    if effective_config.contact_mode == "gayberne":
        gb_params = effective_config.gayberne_params or GayBerneParams()
        if quaternions is None:
            raise RuntimeError("contact-mode gayberne requires quaternion columns")
        if shape_axes is None:
            raise RuntimeError("contact-mode gayberne requires shapex/shapey/shapez columns")
        broad_cut = float(effective_config.r_energy_cap or gb_params.cutoff)
    else:
        gb_params = None
        broad_cut = max(float(r_cut), float(effective_config.g_off) * float(r_cut))

    for i in range(n_particles - 1):
        delta = pos[i + 1 :, :] - pos[i, :]
        delta[:, 0] = minimum_image(delta[:, 0], lx)
        delta[:, 1] = minimum_image(delta[:, 1], ly)
        delta[:, 2] = minimum_image(delta[:, 2], lz)
        distances = np.linalg.norm(delta, axis=1)
        close_mask = (distances > 1e-12) & (distances <= broad_cut)
        if not np.any(close_mask):
            continue

        neighbors = np.nonzero(close_mask)[0] + i + 1
        dots = np.clip(u[neighbors, :] @ u[i, :], -1.0, 1.0)
        q_scores = np.abs(dots)
        pair_scores = p2(q_scores)
        for neighbor, q_score, score, distance in zip(neighbors, q_scores, pair_scores, distances[close_mask]):
            pair_energy: Optional[float] = None
            well_depth: Optional[float] = None
            attraction_strength: Optional[float] = None
            if effective_config.contact_mode == "gayberne":
                neighbor_idx = int(neighbor)
                delta_vec = minimum_image_vector(pos[neighbor_idx, :] - pos[i, :], box)
                assert gb_params is not None
                pair_metrics = gayberne_pair_metrics(
                    r12=delta_vec,
                    quat_i=quaternions[i, :],
                    quat_j=quaternions[neighbor_idx, :],
                    shape_i=shape_axes[i, :],
                    shape_j=shape_axes[neighbor_idx, :],
                    params=gb_params,
                )
                if pair_metrics is None:
                    continue
                pair_energy, well_depth = pair_metrics
                attraction_strength = max(0.0, -float(pair_energy) / max(float(well_depth), 1e-24))
                if effective_config.gb_threshold_mode == "relative":
                    g_value = -float(attraction_strength)
                    is_strong = attraction_strength >= effective_config.gb_on_strength
                    is_gray = attraction_strength >= effective_config.gb_off_strength
                else:
                    g_value = float(pair_energy)
                    is_strong = pair_energy <= effective_config.u_on
                    is_gray = pair_energy <= effective_config.u_off
                if effective_config.gb_require_pair_orientation:
                    is_strong = bool(is_strong and float(q_score) > q_on)
                    is_gray = bool(is_gray and float(q_score) > q_off)
            elif effective_config.contact_mode == "ellipsoid":
                delta_vec = minimum_image_vector(pos[int(neighbor), :] - pos[i, :], box)
                r_hat = delta_vec / max(float(distance), 1e-12)
                shape_i = shape_axes[i, :] if shape_axes is not None else None
                shape_j = shape_axes[int(neighbor), :] if shape_axes is not None else None
                contact_scale = contact_scale_from_shape(
                    shape_i,
                    shape_j,
                    u[i, :],
                    u[int(neighbor), :],
                    r_hat,
                    fallback=r_cut,
                )
                g_value = float(distance) / contact_scale
                is_strong = g_value < effective_config.g_on and float(q_score) > q_on
                is_gray = g_value < effective_config.g_off and float(q_score) > q_off
            else:
                contact_scale = float(r_cut)
                g_value = float(distance) / contact_scale
                is_strong = g_value < effective_config.g_on and float(q_score) > q_on
                is_gray = g_value < effective_config.g_off and float(q_score) > q_off
            if not (is_strong or is_gray):
                continue
            atom_i = int(ids[i])
            atom_j = int(ids[int(neighbor)])
            if normalized_pair(atom_i, atom_j) in excluded_pair_set:
                continue
            delta_s = int(abs(int(s_values[i]) - int(s_values[int(neighbor)])))
            is_explicit_adjacent = normalized_pair(atom_i, atom_j) in explicit_pairs
            is_id_adjacent_fallback = (
                not chain_indices_provided
                and is_topologically_adjacent(atom_i, atom_j, explicit_pairs, adjacent_id_gap)
            )
            is_local = delta_s <= effective_config.s_excl or is_explicit_adjacent or is_id_adjacent_fallback
            if is_strong and not is_local:
                edge_type = "strong"
                can_seed = True
            elif is_strong and is_local:
                edge_type = "local-support"
                can_seed = False
            else:
                edge_type = "gray"
                can_seed = False
            edge = SupportEdge(
                i=i,
                j=int(neighbor),
                atom_i=atom_i,
                atom_j=atom_j,
                distance=float(distance),
                g_value=float(g_value),
                q_score=float(q_score),
                p2_score=float(score),
                delta_s=delta_s,
                edge_type=edge_type,
                can_seed=can_seed,
                is_local=is_local,
                pair_energy=pair_energy,
                well_depth=well_depth,
                attraction_strength=attraction_strength,
                contact_mode=effective_config.contact_mode,
            )
            adjacency[i].append(int(neighbor))
            adjacency[int(neighbor)].append(i)
            edges.append(edge)
    return adjacency, edges


def build_neighbor_graph(
    pos: np.ndarray,
    u: np.ndarray,
    box: Tuple[float, float, float],
    r_cut: float,
    p2_cutoff: float,
) -> Tuple[List[List[int]], int]:
    """Backward-compatible wrapper returning the untyped support graph and pair count."""
    adjacency, edges = build_support_graph(pos, u, box=box, r_cut=r_cut, p2_cutoff=p2_cutoff)
    return adjacency, len(edges)


def connected_components(nodes: Iterable[int], adjacency: Sequence[Sequence[int]], allowed: Optional[set[int]] = None) -> List[List[int]]:
    allowed_nodes = set(nodes) if allowed is None else set(nodes) & allowed
    seen: set[int] = set()
    components: List[List[int]] = []
    for node in sorted(allowed_nodes):
        if node in seen:
            continue
        stack = [node]
        seen.add(node)
        component: List[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor in allowed_nodes and neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return components


def assign_clusters(adjacency: Sequence[Sequence[int]], min_core_neighbors: int) -> Tuple[List[List[int]], np.ndarray]:
    """Legacy DBSCAN-like core/border assignment kept for compatibility."""
    n_particles = len(adjacency)
    degrees = np.array([len(neighbors) for neighbors in adjacency], dtype=int)

    if min_core_neighbors <= 1:
        active_nodes = {idx for idx, degree in enumerate(degrees) if degree > 0}
        clusters = [component for component in connected_components(active_nodes, adjacency) if len(component) >= 2]
        return clusters, degrees

    core_nodes = {idx for idx, degree in enumerate(degrees) if degree >= min_core_neighbors}
    if not core_nodes:
        return [], degrees

    clusters = connected_components(core_nodes, adjacency)
    cluster_lookup: Dict[int, int] = {}
    cluster_sets: List[set[int]] = []
    for cluster_id, component in enumerate(clusters):
        comp_set = set(component)
        cluster_sets.append(comp_set)
        for node in component:
            cluster_lookup[node] = cluster_id

    for particle in range(n_particles):
        if particle in cluster_lookup or degrees[particle] == 0:
            continue
        votes: Dict[int, int] = {}
        for neighbor in adjacency[particle]:
            if neighbor in cluster_lookup:
                cluster_id = cluster_lookup[neighbor]
                votes[cluster_id] = votes.get(cluster_id, 0) + 1
        if not votes:
            continue
        chosen_cluster = min(votes.items(), key=lambda item: (-item[1], item[0]))[0]
        cluster_sets[chosen_cluster].add(particle)

    final_clusters = [sorted(cluster) for cluster in cluster_sets if len(cluster) >= 2]
    final_clusters.sort(key=lambda cluster: (-len(cluster), cluster))
    return final_clusters, degrees


def resolve_visual_cluster_cut(
    r_cut: float,
    shape_axes: Optional[np.ndarray],
    config: AggregationConfig,
) -> float:
    if config.cluster_cut is not None:
        return float(config.cluster_cut)
    if shape_axes is not None and shape_axes.size:
        long_diameters = np.max(shape_axes, axis=1)
        long_diameters = long_diameters[np.isfinite(long_diameters) & (long_diameters > 0.0)]
        if long_diameters.size:
            shape_cut = float(np.median(long_diameters) * config.cluster_cut_shape_factor)
            return float(min(float(r_cut), shape_cut))
    return float(r_cut)


def build_visual_clusters(
    pos: np.ndarray,
    particle_ids: np.ndarray,
    box: Tuple[float, float, float],
    cluster_cut: float,
    min_size: int = 2,
    chain_indices: Optional[np.ndarray] = None,
    s_excl: int = 0,
) -> List[Dict[str, object]]:
    """Build visual spatial clusters for OVITO coloring.

    With chain indices, nonlocal contacts seed clusters; adjacent chain contacts
    may attach particles but cannot create or bridge clusters by themselves.
    """
    n_particles = pos.shape[0]
    if n_particles == 0 or cluster_cut <= 0.0:
        return []

    s_values = np.asarray(chain_indices if chain_indices is not None else np.arange(n_particles), dtype=int)
    use_local_rule = chain_indices is not None
    lx, ly, lz = box
    seed_adjacency: List[List[int]] = [[] for _ in range(n_particles)]
    local_edges: List[Tuple[int, int, float]] = []
    for i in range(n_particles - 1):
        delta = pos[i + 1 :, :] - pos[i, :]
        delta[:, 0] = minimum_image(delta[:, 0], lx)
        delta[:, 1] = minimum_image(delta[:, 1], ly)
        delta[:, 2] = minimum_image(delta[:, 2], lz)
        distances = np.linalg.norm(delta, axis=1)
        neighbors = np.nonzero((distances > 1e-12) & (distances <= cluster_cut))[0] + i + 1
        for neighbor in neighbors:
            neighbor_idx = int(neighbor)
            distance = float(pbc_distance(pos[i, :], pos[neighbor_idx, :], box))
            is_local = use_local_rule and abs(int(s_values[i]) - int(s_values[neighbor_idx])) <= int(s_excl)
            if is_local:
                local_edges.append((i, neighbor_idx, distance))
                continue
            seed_adjacency[i].append(neighbor_idx)
            seed_adjacency[neighbor_idx].append(i)

    seed_nodes = {idx for idx, neighbors in enumerate(seed_adjacency) if neighbors}
    components = connected_components(seed_nodes, seed_adjacency)
    components.sort(key=lambda component: (-len(component), [int(particle_ids[idx]) for idx in component]))
    cluster_sets: List[set[int]] = [set(component) for component in components]
    owner: Dict[int, int] = {}
    for cluster_idx, members in enumerate(cluster_sets):
        for member in members:
            owner[member] = cluster_idx

    changed = True
    while changed:
        changed = False
        for particle in range(n_particles):
            if particle in owner:
                continue
            candidate_scores: Dict[int, Tuple[int, float]] = {}
            for left, right, distance in local_edges:
                other: Optional[int] = None
                if left == particle and right in owner:
                    other = right
                elif right == particle and left in owner:
                    other = left
                if other is None:
                    continue
                cluster_idx = owner[other]
                count, best_distance = candidate_scores.get(cluster_idx, (0, math.inf))
                candidate_scores[cluster_idx] = (count + 1, min(best_distance, distance))
            if not candidate_scores:
                continue
            chosen = min(candidate_scores.items(), key=lambda item: (-item[1][0], item[1][1], item[0]))[0]
            cluster_sets[chosen].add(particle)
            owner[particle] = chosen
            changed = True

    components = [sorted(cluster_set) for cluster_set in cluster_sets if len(cluster_set) >= int(min_size)]
    components.sort(key=lambda component: (-len(component), [int(particle_ids[idx]) for idx in component]))

    records: List[Dict[str, object]] = []
    for cluster_id, component in enumerate(components, start=1):
        component_array = np.asarray(component, dtype=int)
        points = pos[component_array, :]
        center = pbc_centroid(points, box)
        offsets = pbc_offsets_from_center(points, center, box)
        distances = np.linalg.norm(offsets, axis=1)
        radius_of_gyration = float(math.sqrt(np.mean(distances * distances))) if distances.size else 0.0
        max_radius = float(np.max(distances)) if distances.size else 0.0
        records.append(
            {
                "cluster_id": int(cluster_id),
                "size": int(len(component)),
                "members": [int(item) for item in component],
                "atom_ids": [int(particle_ids[idx]) for idx in component],
                "cluster_cut": float(cluster_cut),
                "center": [float(value) for value in center],
                "radius_of_gyration": radius_of_gyration,
                "max_radius": max_radius,
            }
        )
    return records


def make_default_config(axis_used: str, r_cut: float, p2_cutoff: float, min_core_neighbors: int) -> AggregationConfig:
    return AggregationConfig(
        axis=axis_used,
        every=1,
        p2_cut=float(p2_cutoff),
        min_core_neighbors=int(min_core_neighbors),
        cutoff_bins=120,
        cutoff_frames=5,
        r_cut_mode="manual",
        manual_r_cut=float(r_cut),
    )


def support_edge_maps(edges: Sequence[SupportEdge]) -> Tuple[Dict[Tuple[int, int], SupportEdge], Dict[int, List[SupportEdge]]]:
    pair_map: Dict[Tuple[int, int], SupportEdge] = {}
    by_node: Dict[int, List[SupportEdge]] = {}
    for edge in edges:
        pair_map[normalized_pair(edge.i, edge.j)] = edge
        by_node.setdefault(edge.i, []).append(edge)
        by_node.setdefault(edge.j, []).append(edge)
    return pair_map, by_node


def is_component_stable_under_perturbation(
    component: Sequence[int],
    pair_map: Dict[Tuple[int, int], SupportEdge],
    r_cut: float,
    p2_cutoff: float,
    config: AggregationConfig,
) -> bool:
    if len(component) < 2:
        return False

    q_on, _q_off = resolved_q_thresholds(config, p2_cutoff)
    strict_g_on = max(0.0, config.g_on * config.perturbation_r_cut_scale)
    strict_u_on = config.u_on * (1.0 + max(0.0, 1.0 - config.perturbation_r_cut_scale))
    strict_gb_strength = min(
        1.0,
        config.gb_on_strength * (1.0 + max(0.0, 1.0 - config.perturbation_r_cut_scale)),
    )
    strict_q_on = min(1.0, q_on + config.perturbation_p2_margin)
    local_index = {node: idx for idx, node in enumerate(component)}
    strict_adjacency: List[List[int]] = [[] for _ in component]

    for left_pos, left in enumerate(component[:-1]):
        for right in component[left_pos + 1 :]:
            edge = pair_map.get(normalized_pair(left, right))
            if edge is None:
                continue
            if config.contact_mode == "gayberne":
                if config.gb_threshold_mode == "relative":
                    stable_edge = (
                        edge.attraction_strength is not None
                        and float(edge.attraction_strength) >= strict_gb_strength
                    )
                else:
                    stable_edge = edge.pair_energy is not None and float(edge.pair_energy) <= strict_u_on
            else:
                stable_edge = edge.g_value < strict_g_on and edge.q_score > strict_q_on
            if stable_edge:
                li = local_index[left]
                ri = local_index[right]
                strict_adjacency[li].append(ri)
                strict_adjacency[ri].append(li)

    active = {idx for idx, neighbors in enumerate(strict_adjacency) if neighbors}
    if not active:
        return False
    components = connected_components(active, strict_adjacency)
    largest = max((len(item) for item in components), default=0)
    required = max(2, int(math.ceil(config.stable_overlap_fraction * len(component))))
    return largest >= required


def build_domain_candidates(
    pos: np.ndarray,
    u: np.ndarray,
    particle_ids: np.ndarray,
    adjacency: Sequence[Sequence[int]],
    edges: Sequence[SupportEdge],
    box: Tuple[float, float, float],
    r_cut: float,
    p2_cutoff: float,
    config: AggregationConfig,
) -> List[DomainCandidate]:
    seed_adjacency: List[List[int]] = [[] for _ in range(pos.shape[0])]
    for edge in edges:
        if edge.can_seed:
            seed_adjacency[edge.i].append(edge.j)
            seed_adjacency[edge.j].append(edge.i)

    seed_nodes = {idx for idx, neighbors in enumerate(seed_adjacency) if neighbors}
    seed_components = [component for component in connected_components(seed_nodes, seed_adjacency) if len(component) >= 2]
    seed_components.sort(key=lambda component: (-len(component), component))

    domain_sets: List[set[int]] = [set(component) for component in seed_components]
    owner: Dict[int, int] = {}
    for domain_idx, members in enumerate(domain_sets):
        for member in members:
            owner[member] = domain_idx

    # Attach gray/local border particles to existing seeds without merging seeds.
    changed = True
    while changed:
        changed = False
        for particle in range(pos.shape[0]):
            if particle in owner:
                continue
            candidate_scores: Dict[int, Tuple[int, int, float]] = {}
            for edge in edges:
                if edge.edge_type not in {"gray", "local-support"}:
                    continue
                other: Optional[int] = None
                if edge.i == particle and edge.j in owner:
                    other = edge.j
                elif edge.j == particle and edge.i in owner:
                    other = edge.i
                if other is None:
                    continue
                domain_idx = owner[other]
                count, strong_like, best_g = candidate_scores.get(domain_idx, (0, 0, math.inf))
                candidate_scores[domain_idx] = (
                    count + 1,
                    strong_like + (1 if edge.edge_type == "local-support" else 0),
                    min(best_g, float(edge.g_value)),
                )
            if not candidate_scores:
                continue
            chosen = min(
                candidate_scores.items(),
                key=lambda item: (-item[1][0], -item[1][1], item[1][2], item[0]),
            )[0]
            domain_sets[chosen].add(particle)
            owner[particle] = chosen
            changed = True

    domains: List[DomainCandidate] = []
    pair_map, _by_node = support_edge_maps(edges)
    for members in domain_sets:
        component = sorted(members)
        component_set = set(component)
        component_edges = [edge for edge in edges if edge.i in component_set and edge.j in component_set]
        local_edges = sum(1 for edge in component_edges if edge.is_local)
        nonlocal_edges = len(component_edges) - local_edges
        domain_s2, director = compute_Q_and_S(u[np.asarray(component, dtype=int), :])
        stable = is_component_stable_under_perturbation(component, pair_map, r_cut, p2_cutoff, config)
        atom_ids = [int(particle_ids[idx]) for idx in component]
        domains.append(
            DomainCandidate(
                domain_id=0,
                members=component,
                atom_ids=atom_ids,
                s2=float(domain_s2),
                director=director,
                edge_count=len(component_edges),
                adjacent_edge_count=local_edges,
                nonlocal_edge_count=nonlocal_edges,
                stable_under_perturbation=stable,
                classification="robust",
            )
        )

    assigned = set(owner)
    weak_nodes = {idx for idx, neighbors in enumerate(adjacency) if neighbors and idx not in assigned}
    weak_components = [component for component in connected_components(weak_nodes, adjacency) if len(component) >= 2]
    for component in weak_components:
        component_set = set(component)
        component_edges = [edge for edge in edges if edge.i in component_set and edge.j in component_set]
        local_edges = sum(1 for edge in component_edges if edge.is_local)
        nonlocal_edges = len(component_edges) - local_edges
        domain_s2, director = compute_Q_and_S(u[np.asarray(component, dtype=int), :])
        atom_ids = [int(particle_ids[idx]) for idx in component]
        domains.append(
            DomainCandidate(
                domain_id=0,
                members=list(component),
                atom_ids=atom_ids,
                s2=float(domain_s2),
                director=director,
                edge_count=len(component_edges),
                adjacent_edge_count=local_edges,
                nonlocal_edge_count=nonlocal_edges,
                stable_under_perturbation=is_component_stable_under_perturbation(
                    component,
                    pair_map,
                    r_cut,
                    p2_cutoff,
                    config,
                ),
                classification="weak",
            )
        )

    domains.sort(key=lambda domain: (-len(domain.members), domain.classification, domain.atom_ids))
    for domain_id, domain in enumerate(domains, start=1):
        domain.domain_id = domain_id
    return domains


def classify_domains(domains: Sequence[DomainCandidate], config: AggregationConfig) -> None:
    for domain in domains:
        evidence = {
            "size": len(domain.members) >= max(config.n_min, config.robust_min_size),
            "orientation": domain.s2 >= config.robust_min_s2,
            "persistence": domain.age >= config.domain_min_lifetime,
            "nonlocal_support": domain.nonlocal_edge_count > 0,
            "parameter_stability": domain.stable_under_perturbation,
        }
        evidence_count = int(sum(1 for value in evidence.values() if value))
        is_robust = (
            evidence["size"]
            and (evidence["orientation"] or not config.robust_require_orientation)
            and evidence_count >= config.robust_min_evidence
            and (not config.robust_require_nonlocal or evidence["nonlocal_support"])
            and (evidence["nonlocal_support"] or evidence["persistence"])
        )
        domain.evidence = evidence
        domain.evidence_count = evidence_count
        domain.classification = "robust" if is_robust else "weak"


def circular_mean_coordinate(values: np.ndarray, box_length: float) -> float:
    if values.size == 0:
        return 0.0
    angles = 2.0 * math.pi * np.mod(values, box_length) / box_length
    mean_sin = float(np.mean(np.sin(angles)))
    mean_cos = float(np.mean(np.cos(angles)))
    angle = math.atan2(mean_sin, mean_cos)
    if angle < 0.0:
        angle += 2.0 * math.pi
    return float(angle * box_length / (2.0 * math.pi))


def pbc_centroid(points: np.ndarray, box: Tuple[float, float, float]) -> np.ndarray:
    values = []
    for axis in range(3):
        if getattr(box[axis], "periodic", True):
            values.append(circular_mean_coordinate(points[:, axis], box[axis]))
        else:
            values.append(float(np.mean(points[:, axis])) if points.size else 0.0)
    return np.array(values, dtype=float)


def pbc_offsets_from_center(points: np.ndarray, center: np.ndarray, box: Tuple[float, float, float]) -> np.ndarray:
    offsets = points - center
    offsets[:, 0] = minimum_image(offsets[:, 0], box[0])
    offsets[:, 1] = minimum_image(offsets[:, 1], box[1])
    offsets[:, 2] = minimum_image(offsets[:, 2], box[2])
    return offsets


def spatial_contact_support(
    left_members: Sequence[int],
    right_members: Sequence[int],
    pos: np.ndarray,
    box: Tuple[float, float, float],
    gap_cut: float,
) -> Tuple[float, int, int, int]:
    best = math.inf
    contact_count = 0
    left_supported: set[int] = set()
    right_supported: set[int] = set()
    right_array = np.asarray(right_members, dtype=int)
    for left in left_members:
        delta = pos[right_array, :] - pos[left, :]
        delta[:, 0] = minimum_image(delta[:, 0], box[0])
        delta[:, 1] = minimum_image(delta[:, 1], box[1])
        delta[:, 2] = minimum_image(delta[:, 2], box[2])
        distances = np.linalg.norm(delta, axis=1)
        if distances.size:
            best = min(best, float(np.min(distances)))
        supported = np.nonzero(distances <= gap_cut)[0]
        if supported.size:
            left_supported.add(int(left))
            for idx in supported:
                right_supported.add(int(right_array[int(idx)]))
            contact_count += int(supported.size)
    return best, contact_count, len(left_supported), len(right_supported)


def merged_aspect_ratio(members: Sequence[int], pos: np.ndarray, box: Tuple[float, float, float]) -> float:
    if len(members) < 3:
        return math.inf
    points = pos[np.asarray(members, dtype=int), :]
    center = pbc_centroid(points, box)
    offsets = pbc_offsets_from_center(points.copy(), center, box)
    cov = offsets.T @ offsets / max(offsets.shape[0], 1)
    evals = np.linalg.eigvalsh(cov)
    positive = evals[evals > 1e-12]
    if positive.size < 2:
        return math.inf
    return float(math.sqrt(float(np.max(positive)) / float(np.min(positive))))


def component_members(component: Sequence[int], domains: Sequence[DomainCandidate]) -> List[int]:
    members: set[int] = set()
    for domain_idx in component:
        members.update(domains[int(domain_idx)].members)
    return sorted(members)


def refine_pearl_components(
    raw_components: Sequence[Sequence[int]],
    domain_adjacency: Sequence[Sequence[int]],
    robust_domains: Sequence[DomainCandidate],
    pos: np.ndarray,
    box: Tuple[float, float, float],
    max_aspect_ratio: float,
) -> List[List[int]]:
    """Prevent transitive chain merges from becoming one elongated pearl."""
    refined: List[List[int]] = []
    for raw_component in raw_components:
        component = sorted(int(item) for item in raw_component)
        if len(component) <= 1:
            refined.append(component)
            continue
        whole_members = component_members(component, robust_domains)
        if merged_aspect_ratio(whole_members, pos, box) <= max_aspect_ratio:
            refined.append(component)
            continue

        unassigned = set(component)
        while unassigned:
            seed = max(unassigned, key=lambda idx: len(robust_domains[idx].members))
            group = {seed}
            unassigned.remove(seed)
            changed = True
            while changed:
                changed = False
                candidates = sorted(
                    (
                        idx for idx in unassigned
                        if any(idx in domain_adjacency[group_idx] for group_idx in group)
                    ),
                    key=lambda idx: (-len(robust_domains[idx].members), idx),
                )
                for candidate in candidates:
                    trial_group = sorted(group | {candidate})
                    trial_members = component_members(trial_group, robust_domains)
                    if merged_aspect_ratio(trial_members, pos, box) <= max_aspect_ratio:
                        group.add(candidate)
                        unassigned.remove(candidate)
                        changed = True
                        break
            refined.append(sorted(group))
    return refined


def build_pearls(
    robust_domains: Sequence[DomainCandidate],
    pos: np.ndarray,
    particle_ids: np.ndarray,
    box: Tuple[float, float, float],
    r_cut: float,
    config: AggregationConfig,
) -> List[PearlCandidate]:
    if not robust_domains:
        return []

    if config.pearl_gap_cut is not None:
        pearl_gap = float(config.pearl_gap_cut)
    elif config.contact_mode == "gayberne" and config.gayberne_params is not None:
        pearl_gap = float(config.gayberne_params.cutoff)
    else:
        pearl_gap = float(r_cut)
    n_domains = len(robust_domains)
    domain_adjacency: List[List[int]] = [[] for _ in range(n_domains)]

    for i in range(n_domains - 1):
        for j in range(i + 1, n_domains):
            min_distance, contact_count, left_supported, right_supported = spatial_contact_support(
                robust_domains[i].members,
                robust_domains[j].members,
                pos,
                box,
                pearl_gap,
            )
            merged_members = sorted(set(robust_domains[i].members) | set(robust_domains[j].members))
            aspect_ratio = merged_aspect_ratio(
                merged_members,
                pos,
                box,
            )
            has_boundary_support = (
                contact_count >= config.pearl_min_cross_contacts
                and left_supported >= min(config.pearl_min_boundary_particles, len(robust_domains[i].members))
                and right_supported >= min(config.pearl_min_boundary_particles, len(robust_domains[j].members))
            )
            is_bead_like = aspect_ratio <= config.pearl_max_aspect_ratio
            if min_distance <= pearl_gap and has_boundary_support and is_bead_like:
                domain_adjacency[i].append(j)
                domain_adjacency[j].append(i)

    raw_components = connected_components(range(n_domains), domain_adjacency)
    components = refine_pearl_components(
        raw_components=raw_components,
        domain_adjacency=domain_adjacency,
        robust_domains=robust_domains,
        pos=pos,
        box=box,
        max_aspect_ratio=config.pearl_max_aspect_ratio,
    )
    pearls: List[PearlCandidate] = []
    for pearl_id, component in enumerate(components, start=1):
        member_set: set[int] = set()
        domain_ids: List[int] = []
        for domain_idx in component:
            domain = robust_domains[domain_idx]
            domain_ids.append(int(domain.domain_id))
            member_set.update(domain.members)
        members = sorted(member_set)
        points = pos[np.asarray(members, dtype=int), :]
        center = pbc_centroid(points, box)
        offsets = pbc_offsets_from_center(points.copy(), center, box)
        distances = np.linalg.norm(offsets, axis=1)
        radius_of_gyration = float(math.sqrt(np.mean(distances * distances))) if distances.size else 0.0
        max_radius = float(np.max(distances)) if distances.size else 0.0
        compactness = float(radius_of_gyration / max(r_cut, 1e-12))
        aspect_ratio = merged_aspect_ratio(members, pos, box)
        if aspect_ratio > config.pearl_max_aspect_ratio:
            continue
        atom_ids = [int(particle_ids[idx]) for idx in members]
        pearls.append(
            PearlCandidate(
                pearl_id=pearl_id,
                domain_ids=sorted(domain_ids),
                members=members,
                atom_ids=atom_ids,
                center=center,
                radius_of_gyration=radius_of_gyration,
                max_radius=max_radius,
                compactness=compactness,
                aspect_ratio=aspect_ratio,
            )
        )

    pearls.sort(key=lambda pearl: (-len(pearl.members), pearl.atom_ids))
    for idx, pearl in enumerate(pearls, start=1):
        pearl.pearl_id = idx
    return pearls


def resolved_pearl_gap(config: AggregationConfig, r_cut: float) -> float:
    if config.pearl_gap_cut is not None:
        return float(config.pearl_gap_cut)
    if config.contact_mode == "gayberne" and config.gayberne_params is not None:
        return float(config.gayberne_params.cutoff)
    return float(r_cut)


def evaluate_pearl_domain_pairs(
    robust_domains: Sequence[DomainCandidate],
    pos: np.ndarray,
    particle_ids: np.ndarray,
    box: Tuple[float, float, float],
    r_cut: float,
    config: AggregationConfig,
) -> List[Dict[str, object]]:
    pearl_gap = resolved_pearl_gap(config, r_cut)
    rows: List[Dict[str, object]] = []
    for i in range(len(robust_domains) - 1):
        for j in range(i + 1, len(robust_domains)):
            left = robust_domains[i]
            right = robust_domains[j]
            min_distance, contact_count, left_supported, right_supported = spatial_contact_support(
                left.members,
                right.members,
                pos,
                box,
                pearl_gap,
            )
            merged_members = sorted(set(left.members) | set(right.members))
            aspect_ratio = merged_aspect_ratio(merged_members, pos, box)
            has_contact = min_distance <= pearl_gap
            has_boundary_support = (
                contact_count >= config.pearl_min_cross_contacts
                and left_supported >= min(config.pearl_min_boundary_particles, len(left.members))
                and right_supported >= min(config.pearl_min_boundary_particles, len(right.members))
            )
            is_bead_like = aspect_ratio <= config.pearl_max_aspect_ratio
            adjacency_accepted = bool(has_contact and has_boundary_support and is_bead_like)
            reject_reasons: List[str] = []
            if not has_contact:
                reject_reasons.append("gap")
            if not has_boundary_support:
                reject_reasons.append("boundary_support")
            if not is_bead_like:
                reject_reasons.append("aspect_ratio")
            rows.append(
                {
                    "domain_i": int(left.domain_id),
                    "domain_j": int(right.domain_id),
                    "domain_i_size": int(len(left.members)),
                    "domain_j_size": int(len(right.members)),
                    "domain_i_atom_ids": json.dumps([int(item) for item in left.atom_ids]),
                    "domain_j_atom_ids": json.dumps([int(item) for item in right.atom_ids]),
                    "min_distance": float(min_distance),
                    "pearl_gap_cut_used": float(pearl_gap),
                    "contact_count": int(contact_count),
                    "left_supported": int(left_supported),
                    "right_supported": int(right_supported),
                    "boundary_support_min": int(config.pearl_min_boundary_particles),
                    "aspect_ratio": float(aspect_ratio),
                    "aspect_ratio_max": float(config.pearl_max_aspect_ratio),
                    "adjacency_accepted": adjacency_accepted,
                    "reject_reason": "adjacency_accepted" if adjacency_accepted else ",".join(reject_reasons),
                }
            )
    return rows


def fit_pearl_parameters_from_candidates(
    rows: Sequence[Dict[str, object]],
    *,
    current_gap: float,
    current_min_cross_contacts: int,
    current_min_boundary_particles: int,
    current_max_aspect_ratio: float,
    min_candidates: int = 5,
) -> Dict[str, object]:
    accepted = [row for row in rows if bool(row.get("adjacency_accepted", row.get("accepted")))]
    if len(rows) < min_candidates or len(accepted) < max(1, min_candidates // 2):
        return {
            "status": "low",
            "apply_allowed": False,
            "reason": "Too few pearl domain-pair candidates for automatic pearl parameter recommendation.",
            "current": {
                "pearl_gap_cut": float(current_gap),
                "pearl_min_cross_contacts": int(current_min_cross_contacts),
                "pearl_min_boundary_particles": int(current_min_boundary_particles),
                "pearl_max_aspect_ratio": float(current_max_aspect_ratio),
            },
            "recommended": {
                "pearl_gap_cut": float(current_gap),
                "pearl_min_cross_contacts": int(current_min_cross_contacts),
                "pearl_min_boundary_particles": int(current_min_boundary_particles),
                "pearl_max_aspect_ratio": float(current_max_aspect_ratio),
            },
            "sample_sizes": {"candidate_domain_pairs": int(len(rows)), "adjacency_accepted_domain_pairs": int(len(accepted))},
        }
    distances = np.asarray([float(row.get("min_distance", current_gap)) for row in accepted], dtype=float)
    contacts = np.asarray([float(row.get("contact_count", current_min_cross_contacts)) for row in accepted], dtype=float)
    left_support = np.asarray([float(row.get("left_supported", current_min_boundary_particles)) for row in accepted], dtype=float)
    right_support = np.asarray([float(row.get("right_supported", current_min_boundary_particles)) for row in accepted], dtype=float)
    aspects = np.asarray([float(row.get("aspect_ratio", current_max_aspect_ratio)) for row in accepted], dtype=float)
    gap = float(np.max(distances[np.isfinite(distances)])) if np.isfinite(distances).any() else float(current_gap)
    recommended_contacts = max(1, int(math.floor(np.quantile(contacts[np.isfinite(contacts)], 0.25)))) if np.isfinite(contacts).any() else int(current_min_cross_contacts)
    support_min = np.minimum(left_support, right_support)
    recommended_boundary = max(1, int(math.floor(np.quantile(support_min[np.isfinite(support_min)], 0.25)))) if np.isfinite(support_min).any() else int(current_min_boundary_particles)
    aspect = float(np.quantile(aspects[np.isfinite(aspects)], 0.90)) if np.isfinite(aspects).any() else float(current_max_aspect_ratio)
    return {
        "status": "medium",
        "apply_allowed": False,
        "reason": "Recommend-only pearl parameter fit from adjacency-accepted robust-domain pair candidates; do not auto-apply without rerun/stability review.",
        "current": {
            "pearl_gap_cut": float(current_gap),
            "pearl_min_cross_contacts": int(current_min_cross_contacts),
            "pearl_min_boundary_particles": int(current_min_boundary_particles),
            "pearl_max_aspect_ratio": float(current_max_aspect_ratio),
        },
        "recommended": {
            "pearl_gap_cut": float(max(gap, 0.0)),
            "pearl_min_cross_contacts": int(recommended_contacts),
            "pearl_min_boundary_particles": int(recommended_boundary),
            "pearl_max_aspect_ratio": float(max(aspect, 1.0)),
        },
        "sample_sizes": {"candidate_domain_pairs": int(len(rows)), "adjacency_accepted_domain_pairs": int(len(accepted))},
        "notes": [
            "Pearl parameters are fitted from current robust-domain candidates, so this is a second-pass recommendation, not proof of a unique pearl definition."
        ],
    }


def choose_axial_axis(pearls: Sequence[PearlCandidate], global_director: np.ndarray, box: Tuple[float, float, float]) -> np.ndarray:
    if len(pearls) >= 2:
        centers = np.vstack([pearl.center for pearl in pearls])
        reference = centers[0, :]
        offsets = centers - reference
        offsets[:, 0] = minimum_image(offsets[:, 0], box[0])
        offsets[:, 1] = minimum_image(offsets[:, 1], box[1])
        offsets[:, 2] = minimum_image(offsets[:, 2], box[2])
        cov = offsets.T @ offsets
        evals, evecs = np.linalg.eigh(cov)
        axis = evecs[:, int(np.argmax(evals))]
    else:
        axis = np.asarray(global_director, dtype=float)

    norm = float(np.linalg.norm(axis))
    if norm <= 1e-12:
        return np.array([0.0, 0.0, 1.0], dtype=float)
    axis = axis / norm
    first_nonzero = np.flatnonzero(np.abs(axis) > 1e-12)
    if first_nonzero.size and axis[int(first_nonzero[0])] < 0.0:
        axis = -axis
    return axis


def axial_pearl_statistics(
    pearls: Sequence[PearlCandidate],
    pos: np.ndarray,
    box: Tuple[float, float, float],
    global_director: np.ndarray,
) -> Tuple[np.ndarray, List[float], List[float], List[float]]:
    if not pearls:
        return choose_axial_axis([], global_director, box), [], [], []

    axis = choose_axial_axis(pearls, global_director, box)
    reference = pearls[0].center
    records: List[Tuple[float, float, float]] = []
    for pearl in pearls:
        points = pos[np.asarray(pearl.members, dtype=int), :]
        offsets = points - reference
        offsets[:, 0] = minimum_image(offsets[:, 0], box[0])
        offsets[:, 1] = minimum_image(offsets[:, 1], box[1])
        offsets[:, 2] = minimum_image(offsets[:, 2], box[2])
        projection = offsets @ axis
        center = float(np.mean(projection))
        width = float(np.max(projection) - np.min(projection)) if projection.size else 0.0
        half_width = 0.5 * width
        records.append((center, center - half_width, center + half_width))

    records.sort(key=lambda item: item[0])
    centers = [float(item[0]) for item in records]
    widths = [float(item[2] - item[1]) for item in records]
    connectors: List[float] = []
    for left, right in zip(records[:-1], records[1:]):
        connectors.append(float(max(0.0, right[1] - left[2])))
    return axis, centers, widths, connectors


def infer_stretch_axis(pos: np.ndarray, types: np.ndarray, anchor_types: Tuple[int, ...]) -> Optional[np.ndarray]:
    anchor_mask = np.isin(types, np.array(anchor_types, dtype=int))
    anchor_pos = pos[anchor_mask, :]
    if anchor_pos.shape[0] < 2:
        return None
    best_delta = anchor_pos[1, :] - anchor_pos[0, :]
    best_distance = float(np.linalg.norm(best_delta))
    for i in range(anchor_pos.shape[0] - 1):
        deltas = anchor_pos[i + 1 :, :] - anchor_pos[i, :]
        distances = np.linalg.norm(deltas, axis=1)
        if distances.size == 0:
            continue
        idx = int(np.argmax(distances))
        if float(distances[idx]) > best_distance:
            best_distance = float(distances[idx])
            best_delta = deltas[idx, :]
    norm = float(np.linalg.norm(best_delta))
    if norm <= 1e-12:
        return None
    axis = best_delta / norm
    first_nonzero = np.flatnonzero(np.abs(axis) > 1e-12)
    if first_nonzero.size and axis[int(first_nonzero[0])] < 0.0:
        axis = -axis
    return axis


def compute_mechanics_observables(
    mesogen_pos: np.ndarray,
    orientations: np.ndarray,
    mechanics_positions: np.ndarray,
    stretch_axis: Optional[np.ndarray],
    fallback_axis: np.ndarray,
) -> Tuple[float, float, float, float, np.ndarray]:
    axis = np.asarray(stretch_axis if stretch_axis is not None else fallback_axis, dtype=float)
    norm = float(np.linalg.norm(axis))
    if norm <= 1e-12:
        axis = np.array([0.0, 0.0, 1.0], dtype=float)
    else:
        axis = axis / norm

    reference_positions = mechanics_positions if mechanics_positions.size else mesogen_pos
    projected = reference_positions @ axis
    l_parallel = float(np.max(projected) - np.min(projected)) if projected.size else 0.0

    center = np.mean(mesogen_pos, axis=0) if mesogen_pos.size else np.zeros(3, dtype=float)
    offsets = mesogen_pos - center
    mesogen_proj = offsets @ axis
    parallel_vectors = mesogen_proj[:, None] * axis[None, :] if mesogen_proj.size else np.zeros_like(offsets)
    perp_vectors = offsets - parallel_vectors
    rg_parallel = float(math.sqrt(np.mean(mesogen_proj * mesogen_proj))) if mesogen_proj.size else 0.0
    rg_perp = float(math.sqrt(np.mean(np.sum(perp_vectors * perp_vectors, axis=1)))) if perp_vectors.size else 0.0
    dots = np.clip(np.abs(orientations @ axis), 0.0, 1.0) if orientations.size else np.array([], dtype=float)
    s2_force = float(np.mean(p2(dots))) if dots.size else 0.0
    return l_parallel, rg_parallel, rg_perp, s2_force, axis


def clone_config_for_grid(config: AggregationConfig, g_on: float, q_on: float, s_excl: int, n_min: int) -> AggregationConfig:
    return AggregationConfig(
        axis=config.axis,
        every=config.every,
        p2_cut=config.p2_cut,
        min_core_neighbors=config.min_core_neighbors,
        cutoff_bins=config.cutoff_bins,
        cutoff_frames=config.cutoff_frames,
        r_cut_mode=config.r_cut_mode,
        manual_r_cut=config.manual_r_cut,
        max_auto_r_cut=config.max_auto_r_cut,
        auto_r_cut_shape_factor=config.auto_r_cut_shape_factor,
        mesogen_type=config.mesogen_type,
        anchor_types=config.anchor_types,
        contact_mode=config.contact_mode,
        gayberne_params=config.gayberne_params,
        gb_threshold_mode=config.gb_threshold_mode,
        gb_require_pair_orientation=config.gb_require_pair_orientation,
        u_on=config.u_on,
        u_off=config.u_off,
        gb_on_strength=config.gb_on_strength,
        gb_off_strength=config.gb_off_strength,
        r_energy_cap=config.r_energy_cap,
        g_on=float(g_on),
        g_off=max(config.g_off, float(g_on) + 1e-6),
        q_on=float(q_on),
        q_off=config.q_off,
        s_excl=int(s_excl),
        n_min=int(n_min),
        robust_min_size=config.robust_min_size,
        robust_min_s2=config.robust_min_s2,
        robust_min_evidence=config.robust_min_evidence,
        robust_require_nonlocal=config.robust_require_nonlocal,
        domain_min_lifetime=config.domain_min_lifetime,
        adjacent_id_gap=config.adjacent_id_gap,
        perturbation_r_cut_scale=config.perturbation_r_cut_scale,
        perturbation_p2_margin=config.perturbation_p2_margin,
        stable_overlap_fraction=config.stable_overlap_fraction,
        pearl_gap_cut=config.pearl_gap_cut,
        pearl_min_cross_contacts=config.pearl_min_cross_contacts,
        pearl_min_boundary_particles=config.pearl_min_boundary_particles,
        pearl_max_aspect_ratio=config.pearl_max_aspect_ratio,
        track_jaccard=config.track_jaccard,
        consensus_threshold=config.consensus_threshold,
        enable_robustness_scan=False,
        write_ovito_labels=config.write_ovito_labels,
        write_cluster_envelopes=config.write_cluster_envelopes,
        write_contact_edges=config.write_contact_edges,
        write_contact_segments=config.write_contact_segments,
        write_diagnostics=config.write_diagnostics,
        shared_r_cut=config.shared_r_cut,
        track_across_files=config.track_across_files,
        workers=config.workers,
        cluster_cut=config.cluster_cut,
        cluster_cut_shape_factor=config.cluster_cut_shape_factor,
        cluster_min_size=config.cluster_min_size,
        cluster_envelope_padding=config.cluster_envelope_padding,
    )


def compute_consensus_summary(
    pos: np.ndarray,
    u: np.ndarray,
    particle_ids: np.ndarray,
    chain_indices: np.ndarray,
    shape_axes: Optional[np.ndarray],
    quaternions: Optional[np.ndarray],
    box: Tuple[float, float, float],
    r_cut: float,
    p2_cutoff: float,
    base_config: AggregationConfig,
) -> Tuple[float, float]:
    if not base_config.enable_robustness_scan or pos.shape[0] < 2:
        return 0.0, 0.0

    q_on, _q_off = resolved_q_thresholds(base_config, p2_cutoff)
    g_values = sorted({max(1e-6, base_config.g_on * scale) for scale in (0.95, 1.0, 1.05)})
    q_values = sorted({min(1.0, max(0.0, q_on + delta)) for delta in (-0.05, 0.0, 0.05)})
    s_values = sorted({max(0, base_config.s_excl - 1), base_config.s_excl, base_config.s_excl + 1})
    n_values = sorted({max(2, base_config.n_min - 1), base_config.n_min, base_config.n_min + 1})

    n_particles = pos.shape[0]
    consensus = np.zeros((n_particles, n_particles), dtype=float)
    runs = 0
    for g_on in g_values:
        for q_trial in q_values:
            for s_excl in s_values:
                for n_min in n_values:
                    trial = clone_config_for_grid(base_config, g_on=g_on, q_on=q_trial, s_excl=s_excl, n_min=n_min)
                    adjacency, edges = build_support_graph(
                        pos=pos,
                        u=u,
                        box=box,
                        r_cut=r_cut,
                        p2_cutoff=p2_cutoff,
                        particle_ids=particle_ids,
                        chain_indices=chain_indices,
                        shape_axes=shape_axes,
                        quaternions=quaternions,
                        config=trial,
                    )
                    domains = build_domain_candidates(
                        pos=pos,
                        u=u,
                        particle_ids=particle_ids,
                        adjacency=adjacency,
                        edges=edges,
                        box=box,
                        r_cut=r_cut,
                        p2_cutoff=p2_cutoff,
                        config=trial,
                    )
                    classify_domains(domains, trial)
                    for domain in domains:
                        if domain.classification != "robust":
                            continue
                        members = np.asarray(domain.members, dtype=int)
                        consensus[np.ix_(members, members)] += 1.0
                    runs += 1

    if runs == 0:
        return 0.0, 0.0
    consensus /= float(runs)
    upper = consensus[np.triu_indices(n_particles, k=1)]
    stable_core_pair_fraction = float(np.mean(upper >= base_config.consensus_threshold)) if upper.size else 0.0
    ambiguous = (consensus > 0.0) & (consensus < base_config.consensus_threshold)
    ambiguous_mesogens = np.any(ambiguous, axis=1)
    ambiguous_mesogen_fraction = float(np.mean(ambiguous_mesogens)) if ambiguous_mesogens.size else 0.0
    return stable_core_pair_fraction, ambiguous_mesogen_fraction


def domain_to_record(domain: DomainCandidate) -> Dict[str, object]:
    return {
        "domain_id": int(domain.domain_id),
        "track_id": int(domain.track_id),
        "age": int(domain.age),
        "classification": domain.classification,
        "size": len(domain.members),
        "atom_ids": [int(item) for item in domain.atom_ids],
        "s2": float(domain.s2),
        "edge_count": int(domain.edge_count),
        "adjacent_edge_count": int(domain.adjacent_edge_count),
        "nonlocal_edge_count": int(domain.nonlocal_edge_count),
        "stable_under_perturbation": bool(domain.stable_under_perturbation),
        "evidence_count": int(domain.evidence_count),
        "evidence": dict(domain.evidence or {}),
    }


def pearl_to_record(pearl: PearlCandidate) -> Dict[str, object]:
    return {
        "pearl_id": int(pearl.pearl_id),
        "domain_ids": [int(item) for item in pearl.domain_ids],
        "size": len(pearl.members),
        "atom_ids": [int(item) for item in pearl.atom_ids],
        "center": [float(value) for value in pearl.center],
        "radius_of_gyration": float(pearl.radius_of_gyration),
        "max_radius": float(pearl.max_radius),
        "compactness": float(pearl.compactness),
        "aspect_ratio": float(pearl.aspect_ratio),
    }


def particle_contact_records_from_edges(
    edges: Sequence[SupportEdge],
    particle_ids: np.ndarray,
) -> List[Dict[str, object]]:
    energy_by_atom: Dict[int, List[float]] = {int(atom_id): [] for atom_id in particle_ids}
    strength_by_atom: Dict[int, List[float]] = {int(atom_id): [] for atom_id in particle_ids}
    degree_by_atom: Dict[int, int] = {int(atom_id): 0 for atom_id in particle_ids}
    for edge in edges:
        for atom_id in (edge.atom_i, edge.atom_j):
            degree_by_atom[int(atom_id)] = degree_by_atom.get(int(atom_id), 0) + 1
            if edge.pair_energy is not None:
                energy_by_atom.setdefault(int(atom_id), []).append(float(edge.pair_energy))
            if edge.attraction_strength is not None:
                strength_by_atom.setdefault(int(atom_id), []).append(float(edge.attraction_strength))

    records: List[Dict[str, object]] = []
    for atom_id in particle_ids:
        atom_key = int(atom_id)
        values = energy_by_atom.get(atom_key, [])
        strengths = strength_by_atom.get(atom_key, [])
        records.append(
            {
                "atom_id": atom_key,
                "contact_degree": int(degree_by_atom.get(atom_key, 0)),
                "min_pair_energy": float(min(values)) if values else 0.0,
                "mean_pair_energy": float(np.mean(np.array(values, dtype=float))) if values else 0.0,
                "max_attraction_strength": float(max(strengths)) if strengths else 0.0,
                "mean_attraction_strength": float(np.mean(np.array(strengths, dtype=float))) if strengths else 0.0,
            }
        )
    return records


def contact_edge_records_from_edges(edges: Sequence[SupportEdge]) -> List[Dict[str, object]]:
    edge_type_codes = {"local-support": 1, "gray": 2, "strong": 3}
    records: List[Dict[str, object]] = []
    for edge_id, edge in enumerate(edges, start=1):
        records.append(
            {
                "edge_id": int(edge_id),
                "atom_i": int(edge.atom_i),
                "atom_j": int(edge.atom_j),
                "distance": float(edge.distance),
                "pair_energy": float(edge.pair_energy) if edge.pair_energy is not None else 0.0,
                "well_depth": float(edge.well_depth) if edge.well_depth is not None else 0.0,
                "attraction_strength": float(edge.attraction_strength) if edge.attraction_strength is not None else 0.0,
                "q_score": float(edge.q_score),
                "p2_score": float(edge.p2_score),
                "delta_s": int(edge.delta_s),
                "edge_type": edge.edge_type,
                "edge_type_code": int(edge_type_codes.get(edge.edge_type, 0)),
                "is_local": int(edge.is_local),
            }
        )
    return records


def compute_frame_aggregation(
    source_file: str,
    timestep: int,
    pos: np.ndarray,
    u: np.ndarray,
    axis_used: str,
    box: Tuple[float, float, float],
    r_cut: float,
    p2_cutoff: float,
    min_core_neighbors: int,
    particle_ids: Optional[np.ndarray] = None,
    chain_indices: Optional[np.ndarray] = None,
    shape_axes: Optional[np.ndarray] = None,
    quaternions: Optional[np.ndarray] = None,
    stretch_axis: Optional[np.ndarray] = None,
    mechanics_positions: Optional[np.ndarray] = None,
    config: Optional[AggregationConfig] = None,
    domain_tracker: Optional[DomainTracker] = None,
    adjacent_pairs: Optional[Iterable[Tuple[int, int]]] = None,
    excluded_pairs: Optional[Iterable[Tuple[int, int]]] = None,
) -> FrameAggregationResult:
    n_particles = pos.shape[0]
    ids = np.asarray(particle_ids if particle_ids is not None else np.arange(n_particles), dtype=int)
    s_values = np.asarray(chain_indices if chain_indices is not None else np.arange(n_particles), dtype=int)
    effective_config = config or make_default_config(axis_used, r_cut, p2_cutoff, min_core_neighbors)

    adjacency, edges = build_support_graph(
        pos=pos,
        u=u,
        box=box,
        r_cut=r_cut,
        p2_cutoff=p2_cutoff,
        particle_ids=ids,
        chain_indices=s_values,
        shape_axes=shape_axes,
        quaternions=quaternions,
        config=effective_config,
        adjacent_pairs=adjacent_pairs,
        excluded_pairs=excluded_pairs,
        adjacent_id_gap=effective_config.adjacent_id_gap,
    )
    degrees = np.array([len(neighbors) for neighbors in adjacency], dtype=int)
    global_s2, global_director = compute_Q_and_S(u)
    visual_cluster_cut = resolve_visual_cluster_cut(r_cut, shape_axes, effective_config)
    visual_cluster_records = build_visual_clusters(
        pos=pos,
        particle_ids=ids,
        box=box,
        cluster_cut=visual_cluster_cut,
        min_size=effective_config.cluster_min_size,
        chain_indices=s_values,
        s_excl=effective_config.s_excl,
    )
    visual_cluster_sizes = [int(record["size"]) for record in visual_cluster_records]

    domains = build_domain_candidates(
        pos=pos,
        u=u,
        particle_ids=ids,
        adjacency=adjacency,
        edges=edges,
        box=box,
        r_cut=r_cut,
        p2_cutoff=p2_cutoff,
        config=effective_config,
    )
    if domain_tracker is not None:
        domain_tracker.update(domains)
    classify_domains(domains, effective_config)

    robust_domains = [domain for domain in domains if domain.classification == "robust"]
    weak_domains = [domain for domain in domains if domain.classification == "weak"]

    robust_domains.sort(key=lambda domain: (-len(domain.members), domain.atom_ids))
    weak_domains.sort(key=lambda domain: (-len(domain.members), domain.atom_ids))

    robust_domain_sizes = [len(domain.members) for domain in robust_domains]
    weak_domain_sizes = [len(domain.members) for domain in weak_domains]
    robust_domain_s2 = [float(domain.s2) for domain in robust_domains]

    aggregation_degree = 0.0
    robust_particles = 0
    for domain in robust_domains:
        size = len(domain.members)
        robust_particles += size
        aggregation_degree += (size / n_particles) ** 2 * max(float(domain.s2), 0.0)

    weak_particles = sum(len(domain.members) for domain in weak_domains)
    largest_cluster_fraction = max(robust_domain_sizes) / n_particles if robust_domain_sizes else 0.0
    clustered_fraction = robust_particles / n_particles if n_particles else 0.0
    robust_domain_fraction = clustered_fraction
    weak_domain_fraction = weak_particles / n_particles if n_particles else 0.0
    core_fraction = float(np.sum(degrees >= min_core_neighbors)) / n_particles if n_particles else 0.0
    mean_cluster_s2 = (
        float(np.average(np.array(robust_domain_s2, dtype=float), weights=np.array(robust_domain_sizes, dtype=float)))
        if robust_domain_sizes
        else 0.0
    )
    max_cluster_s2 = max(robust_domain_s2) if robust_domain_s2 else 0.0
    qualified_pairs = len(edges)
    qualified_pair_fraction = (2.0 * qualified_pairs / (n_particles * (n_particles - 1))) if n_particles > 1 else 0.0
    pair_energies = [float(edge.pair_energy) for edge in edges if edge.pair_energy is not None]
    particle_contact_records = particle_contact_records_from_edges(edges, ids)
    contact_edge_records = contact_edge_records_from_edges(edges)
    pearl_candidate_records = evaluate_pearl_domain_pairs(
        robust_domains=robust_domains,
        pos=pos,
        particle_ids=ids,
        box=box,
        r_cut=r_cut,
        config=effective_config,
    )

    pearls = build_pearls(
        robust_domains=robust_domains,
        pos=pos,
        particle_ids=ids,
        box=box,
        r_cut=r_cut,
        config=effective_config,
    )
    axial_axis, pearl_axis_centers, pearl_axis_widths, connector_lengths = axial_pearl_statistics(
        pearls=pearls,
        pos=pos,
        box=box,
        global_director=global_director,
    )
    local_edge_count = sum(1 for edge in edges if edge.is_local)
    nonlocal_edge_count = len(edges) - local_edge_count
    local_edge_fraction = local_edge_count / len(edges) if edges else 0.0
    nonlocal_edge_fraction = nonlocal_edge_count / len(edges) if edges else 0.0
    stable_core_pair_fraction, ambiguous_mesogen_fraction = compute_consensus_summary(
        pos=pos,
        u=u,
        particle_ids=ids,
        chain_indices=s_values,
        shape_axes=shape_axes,
        quaternions=quaternions,
        box=box,
        r_cut=r_cut,
        p2_cutoff=p2_cutoff,
        base_config=effective_config,
    )
    l_parallel, rg_parallel, rg_perp, s2_force, force_axis = compute_mechanics_observables(
        mesogen_pos=pos,
        orientations=u,
        mechanics_positions=mechanics_positions if mechanics_positions is not None else pos,
        stretch_axis=stretch_axis,
        fallback_axis=axial_axis,
    )

    return FrameAggregationResult(
        source_file=source_file,
        timestep=int(timestep),
        axis_used=axis_used,
        n_particles=int(n_particles),
        r_cut=float(r_cut),
        aggregation_degree=float(aggregation_degree),
        largest_cluster_fraction=float(largest_cluster_fraction),
        clustered_fraction=float(clustered_fraction),
        core_fraction=float(core_fraction),
        mean_cluster_s2=float(mean_cluster_s2),
        max_cluster_s2=float(max_cluster_s2),
        global_s2=float(global_s2),
        n_clusters=len(robust_domain_sizes),
        qualified_pairs=int(qualified_pairs),
        qualified_pair_fraction=float(qualified_pair_fraction),
        energy_edge_count=len(pair_energies),
        min_pair_energy=float(min(pair_energies)) if pair_energies else 0.0,
        mean_pair_energy=float(np.mean(np.array(pair_energies, dtype=float))) if pair_energies else 0.0,
        cluster_sizes=robust_domain_sizes,
        cluster_s2=robust_domain_s2,
        visual_cluster_count=len(visual_cluster_records),
        largest_visual_cluster_fraction=(
            max(visual_cluster_sizes, default=0) / n_particles if n_particles else 0.0
        ),
        visual_cluster_sizes=visual_cluster_sizes,
        weak_domain_count=len(weak_domains),
        robust_domain_count=len(robust_domains),
        weak_domain_fraction=float(weak_domain_fraction),
        robust_domain_fraction=float(robust_domain_fraction),
        weak_domain_sizes=weak_domain_sizes,
        robust_domain_sizes=robust_domain_sizes,
        robust_domain_s2=robust_domain_s2,
        pearl_count=len(pearls),
        largest_pearl_fraction=(max((len(pearl.members) for pearl in pearls), default=0) / n_particles if n_particles else 0.0),
        pearl_sizes=[len(pearl.members) for pearl in pearls],
        pearl_domain_counts=[len(pearl.domain_ids) for pearl in pearls],
        pearl_compactness=[float(pearl.compactness) for pearl in pearls],
        pearl_aspect_ratios=[float(pearl.aspect_ratio) for pearl in pearls],
        pearl_axis_centers=[float(value) for value in pearl_axis_centers],
        pearl_axis_widths=[float(value) for value in pearl_axis_widths],
        connector_lengths=[float(value) for value in connector_lengths],
        axial_axis=[float(value) for value in axial_axis],
        local_edge_fraction=float(local_edge_fraction),
        nonlocal_edge_fraction=float(nonlocal_edge_fraction),
        ambiguous_mesogen_fraction=float(ambiguous_mesogen_fraction),
        stable_core_pair_fraction=float(stable_core_pair_fraction),
        l_parallel=float(l_parallel),
        rg_parallel=float(rg_parallel),
        rg_perp=float(rg_perp),
        s2_force=float(s2_force),
        stretch_axis=[float(value) for value in force_axis],
        particle_contact_records=particle_contact_records,
        contact_edge_records=contact_edge_records,
        visual_cluster_records=visual_cluster_records,
        domain_records=[domain_to_record(domain) for domain in domains],
        pearl_records=[pearl_to_record(pearl) for pearl in pearls],
        pearl_candidate_records=pearl_candidate_records,
    )


def process_dump_file(
    dump_path: Path,
    config: AggregationConfig,
    tracker: Optional[DomainTracker] = None,
) -> Tuple[List[FrameAggregationResult], Dict[str, object]]:
    r_cut = determine_r_cut_for_dump(dump_path, config)
    results: List[FrameAggregationResult] = []
    frame_index = 0
    domain_tracker = tracker if tracker is not None else DomainTracker(min_jaccard=config.track_jaccard)
    local_pairs = read_pair_list(config.local_pair_file)
    excluded_pairs = read_pair_list(config.exclude_pair_file)

    for timestep, box, _columns, col_index, data in parse_dump_frames(dump_path):
        frame_index += 1
        if (frame_index - 1) % config.every != 0:
            continue
        all_pos = extract_positions(data, col_index)
        all_ids = extract_particle_ids(data, col_index)
        all_types = extract_particle_types(data, col_index)
        all_shapes = extract_shape_axes(data, col_index)
        mesogen_indices = select_mesogen_indices(all_types, config.mesogen_type)
        if mesogen_indices.size == 0:
            raise RuntimeError(f"{dump_path}: no type {config.mesogen_type} mesogens found in timestep {timestep}")
        pos = all_pos[mesogen_indices, :]
        particle_ids = all_ids[mesogen_indices]
        u, axis_used = extract_orientations(data[mesogen_indices, :], col_index, config.axis)
        all_quaternions = extract_quaternions(data, col_index)
        quaternions = all_quaternions[mesogen_indices, :]
        shape_axes = all_shapes[mesogen_indices, :] if all_shapes is not None else None
        chain_indices = infer_chain_indices(
            data,
            col_index,
            mesogen_indices,
            all_ids,
            positions=all_pos,
            types=all_types,
            anchor_types=config.anchor_types,
        )
        stretch_axis = infer_stretch_axis(all_pos, all_types, config.anchor_types)
        mechanics_positions = all_pos[np.isin(all_types, np.array(config.anchor_types, dtype=int)), :]
        result = compute_frame_aggregation(
            source_file=str(dump_path),
            timestep=timestep,
            pos=pos,
            u=u,
            axis_used=axis_used,
            box=box,
            r_cut=r_cut,
            p2_cutoff=config.p2_cut,
            min_core_neighbors=config.min_core_neighbors,
            particle_ids=particle_ids,
            chain_indices=chain_indices,
            shape_axes=shape_axes,
            quaternions=quaternions,
            stretch_axis=stretch_axis,
            mechanics_positions=mechanics_positions,
            config=config,
            domain_tracker=domain_tracker,
            adjacent_pairs=local_pairs,
            excluded_pairs=excluded_pairs,
        )
        results.append(result)

    if not results:
        raise RuntimeError(f"{dump_path}: no frames were processed. Check --every and dump format.")

    aggregation_values = np.array([item.aggregation_degree for item in results], dtype=float)
    largest_values = np.array([item.largest_cluster_fraction for item in results], dtype=float)
    clustered_values = np.array([item.clustered_fraction for item in results], dtype=float)
    mean_cluster_values = np.array([item.mean_cluster_s2 for item in results], dtype=float)
    global_s2_values = np.array([item.global_s2 for item in results], dtype=float)
    visual_cluster_counts = np.array([item.visual_cluster_count for item in results], dtype=float)
    largest_visual_cluster_values = np.array([item.largest_visual_cluster_fraction for item in results], dtype=float)
    robust_counts = np.array([item.robust_domain_count for item in results], dtype=float)
    weak_counts = np.array([item.weak_domain_count for item in results], dtype=float)
    pearl_counts = np.array([item.pearl_count for item in results], dtype=float)
    largest_pearl_values = np.array([item.largest_pearl_fraction for item in results], dtype=float)
    l_parallel_values = np.array([item.l_parallel for item in results], dtype=float)
    rg_parallel_values = np.array([item.rg_parallel for item in results], dtype=float)
    rg_perp_values = np.array([item.rg_perp for item in results], dtype=float)
    s2_force_values = np.array([item.s2_force for item in results], dtype=float)
    ambiguous_values = np.array([item.ambiguous_mesogen_fraction for item in results], dtype=float)
    min_pair_energy_values = np.array([item.min_pair_energy for item in results], dtype=float)
    mean_pair_energy_values = np.array([item.mean_pair_energy for item in results], dtype=float)

    summary: Dict[str, object] = {
        "source_file": str(dump_path),
        "frames_processed": len(results),
        "r_cut": float(r_cut),
        "mean_aggregation_degree": float(np.mean(aggregation_values)),
        "max_aggregation_degree": float(np.max(aggregation_values)),
        "mean_largest_cluster_fraction": float(np.mean(largest_values)),
        "mean_clustered_fraction": float(np.mean(clustered_values)),
        "mean_cluster_s2": float(np.mean(mean_cluster_values)),
        "mean_global_s2": float(np.mean(global_s2_values)),
        "mean_visual_cluster_count": float(np.mean(visual_cluster_counts)),
        "mean_largest_visual_cluster_fraction": float(np.mean(largest_visual_cluster_values)),
        "min_pair_energy": float(np.min(min_pair_energy_values)) if min_pair_energy_values.size else 0.0,
        "mean_pair_energy": float(np.mean(mean_pair_energy_values)) if mean_pair_energy_values.size else 0.0,
        "mean_robust_domain_count": float(np.mean(robust_counts)),
        "mean_weak_domain_count": float(np.mean(weak_counts)),
        "mean_pearl_count": float(np.mean(pearl_counts)),
        "mean_largest_pearl_fraction": float(np.mean(largest_pearl_values)),
        "mean_L_parallel": float(np.mean(l_parallel_values)),
        "mean_Rg_parallel": float(np.mean(rg_parallel_values)),
        "mean_Rg_perp": float(np.mean(rg_perp_values)),
        "mean_S2_force": float(np.mean(s2_force_values)),
        "mean_ambiguous_mesogen_fraction": float(np.mean(ambiguous_values)),
        "last_timestep": int(results[-1].timestep),
        "last_frame_compact": {
            "aggregation_degree": float(results[-1].aggregation_degree),
            "robust_domain_count": int(results[-1].robust_domain_count),
            "weak_domain_count": int(results[-1].weak_domain_count),
            "pearl_count": int(results[-1].pearl_count),
            "L_parallel": float(results[-1].l_parallel),
        },
    }
    return results, summary


def iter_input_files(inputs: Sequence[Path], pattern: str, recursive: bool) -> List[Path]:
    files: List[Path] = []
    for entry in inputs:
        path = entry.resolve()
        if path.is_file():
            files.append(path)
            continue
        if not path.is_dir():
            raise FileNotFoundError(f"input path does not exist: {entry}")
        if recursive:
            files.extend(sorted(path.rglob(pattern)))
        else:
            files.extend(sorted(path.glob(pattern)))
    return sorted({file.resolve() for file in files if file.is_file()}, key=dump_sort_key)


def process_dump_file_job(
    dump_file: Path,
    config: AggregationConfig,
) -> Tuple[Path, List[FrameAggregationResult], Dict[str, object]]:
    results, summary = process_dump_file(dump_file, config=config)
    return dump_file, results, summary


def iter_processed_dump_file_sequence(
    files: Sequence[Path],
    config: AggregationConfig,
) -> Iterator[Tuple[Path, List[FrameAggregationResult], Dict[str, object]]]:
    """Process a one-frame-per-file sequence with one domain tracker across files."""
    tracker = DomainTracker(min_jaccard=config.track_jaccard)
    for dump_file in sorted(files, key=dump_sort_key):
        results, summary = process_dump_file(dump_file, config=config, tracker=tracker)
        yield dump_file, results, summary


def iter_processed_dump_files(
    files: Sequence[Path],
    config: AggregationConfig,
    worker_count: int,
) -> Iterator[Tuple[Path, List[FrameAggregationResult], Dict[str, object]]]:
    if worker_count <= 1 or len(files) <= 1:
        for dump_file in files:
            yield process_dump_file_job(dump_file, config)
        return

    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        future_to_file = {
            executor.submit(process_dump_file_job, dump_file, config): dump_file
            for dump_file in files
        }
        for future in as_completed(future_to_file):
            yield future.result()


def json_list(values: Sequence[object], precision: Optional[int] = None) -> str:
    if precision is None:
        payload = list(values)
    else:
        payload = [round(float(value), precision) for value in values]
    return json.dumps(payload, ensure_ascii=False)


def result_to_tsv_row(result: FrameAggregationResult) -> str:
    values = [
        result.source_file,
        str(result.timestep),
        result.axis_used,
        str(result.n_particles),
        f"{result.r_cut:.8f}",
        f"{result.aggregation_degree:.8f}",
        f"{result.largest_cluster_fraction:.8f}",
        f"{result.clustered_fraction:.8f}",
        f"{result.core_fraction:.8f}",
        f"{result.mean_cluster_s2:.8f}",
        f"{result.max_cluster_s2:.8f}",
        f"{result.global_s2:.8f}",
        str(result.n_clusters),
        str(result.qualified_pairs),
        f"{result.qualified_pair_fraction:.8f}",
        str(result.energy_edge_count),
        f"{result.min_pair_energy:.8f}",
        f"{result.mean_pair_energy:.8f}",
        json_list(result.cluster_sizes),
        json_list(result.cluster_s2, precision=8),
        str(result.visual_cluster_count),
        f"{result.largest_visual_cluster_fraction:.8f}",
        json_list(result.visual_cluster_sizes),
        str(result.weak_domain_count),
        str(result.robust_domain_count),
        f"{result.weak_domain_fraction:.8f}",
        f"{result.robust_domain_fraction:.8f}",
        json_list(result.weak_domain_sizes),
        json_list(result.robust_domain_sizes),
        json_list(result.robust_domain_s2, precision=8),
        str(result.pearl_count),
        f"{result.largest_pearl_fraction:.8f}",
        json_list(result.pearl_sizes),
        json_list(result.pearl_domain_counts),
        json_list(result.pearl_compactness, precision=8),
        json_list(result.pearl_aspect_ratios, precision=8),
        json_list(result.pearl_axis_centers, precision=8),
        json_list(result.pearl_axis_widths, precision=8),
        json_list(result.connector_lengths, precision=8),
        json_list(result.axial_axis, precision=8),
        f"{result.local_edge_fraction:.8f}",
        f"{result.nonlocal_edge_fraction:.8f}",
        f"{result.ambiguous_mesogen_fraction:.8f}",
        f"{result.stable_core_pair_fraction:.8f}",
        f"{result.l_parallel:.8f}",
        f"{result.rg_parallel:.8f}",
        f"{result.rg_perp:.8f}",
        f"{result.s2_force:.8f}",
        json_list(result.stretch_axis, precision=8),
    ]
    return "\t".join(values) + "\n"


def build_label_maps(
    result: FrameAggregationResult,
) -> Tuple[
    Dict[int, int],
    Dict[int, int],
    Dict[int, int],
    Dict[int, int],
    Dict[int, int],
    Dict[int, int],
    Dict[int, float],
    Dict[int, float],
]:
    atom_to_cluster: Dict[int, int] = {}
    atom_to_cluster_size: Dict[int, int] = {}
    atom_to_domain: Dict[int, int] = {}
    atom_to_state: Dict[int, int] = {}
    atom_to_pearl: Dict[int, int] = {}
    atom_to_contact_degree: Dict[int, int] = {}
    atom_to_min_energy: Dict[int, float] = {}
    atom_to_mean_energy: Dict[int, float] = {}
    atom_to_max_strength: Dict[int, float] = {}
    atom_to_mean_strength: Dict[int, float] = {}

    for record in result.particle_contact_records:
        atom_id = int(record["atom_id"])
        atom_to_contact_degree[atom_id] = int(record["contact_degree"])
        atom_to_min_energy[atom_id] = float(record["min_pair_energy"])
        atom_to_mean_energy[atom_id] = float(record["mean_pair_energy"])
        atom_to_max_strength[atom_id] = float(record.get("max_attraction_strength", 0.0))
        atom_to_mean_strength[atom_id] = float(record.get("mean_attraction_strength", 0.0))

    for cluster in result.visual_cluster_records:
        cluster_id = int(cluster["cluster_id"])
        cluster_size = int(cluster["size"])
        for atom_id in cluster["atom_ids"]:
            atom_to_cluster[int(atom_id)] = cluster_id
            atom_to_cluster_size[int(atom_id)] = cluster_size

    for domain in result.domain_records:
        domain_id = int(domain["domain_id"])
        state = 2 if domain["classification"] == "robust" else 1
        for atom_id in domain["atom_ids"]:
            atom_to_domain[int(atom_id)] = domain_id
            atom_to_state[int(atom_id)] = state

    for pearl in result.pearl_records:
        pearl_id = int(pearl["pearl_id"])
        for atom_id in pearl["atom_ids"]:
            atom_to_pearl[int(atom_id)] = pearl_id

    return (
        atom_to_cluster,
        atom_to_cluster_size,
        atom_to_domain,
        atom_to_pearl,
        atom_to_state,
        atom_to_contact_degree,
        atom_to_min_energy,
        atom_to_mean_energy,
        atom_to_max_strength,
        atom_to_mean_strength,
    )


def format_dump_atom_value(column: str, value: float) -> str:
    """Format original dump values while keeping integer identity columns readable."""
    integer_columns = {"id", "type", "mol", "ix", "iy", "iz"}
    if column in integer_columns:
        return str(int(round(float(value))))
    return f"{float(value):.12g}"


def write_ovito_label_dump(
    dump_path: Path,
    results: Sequence[FrameAggregationResult],
    output_path: Path,
) -> None:
    """Write an OVITO-readable dump preserving atom columns and appending LC labels."""
    by_timestep = {int(result.timestep): result for result in results}
    with output_path.open("w", encoding="utf-8") as handle:
        for timestep, box, columns, col_index, data in parse_dump_frames(dump_path):
            result = by_timestep.get(int(timestep))
            if result is None:
                continue

            ids = extract_particle_ids(data, col_index)
            (
                atom_to_cluster,
                atom_to_cluster_size,
                atom_to_domain,
                atom_to_pearl,
                atom_to_state,
                atom_to_contact_degree,
                atom_to_min_energy,
                atom_to_mean_energy,
                atom_to_max_strength,
                atom_to_mean_strength,
            ) = build_label_maps(result)
            label_columns = [
                "lc_cluster",
                "lc_cluster_size",
                "lc_contact_degree",
                "lc_min_pair_energy",
                "lc_mean_pair_energy",
                "lc_max_gb_strength",
                "lc_mean_gb_strength",
                "lc_domain",
                "lc_pearl",
                "lc_state",
            ]
            base_columns = [column for column in columns if column not in label_columns]
            base_indices = [col_index[column] for column in base_columns]

            handle.write("ITEM: TIMESTEP\n")
            handle.write(f"{int(timestep)}\n")
            handle.write("ITEM: NUMBER OF ATOMS\n")
            handle.write(f"{data.shape[0]}\n")
            write_box_bounds(handle, box)
            handle.write(f"ITEM: ATOMS {' '.join(base_columns + label_columns)}\n")
            for atom_id, row in zip(ids, data):
                atom_key = int(atom_id)
                values = [
                    format_dump_atom_value(column, float(row[index]))
                    for column, index in zip(base_columns, base_indices)
                ]
                values.extend(
                    [
                        str(atom_to_cluster.get(atom_key, 0)),
                        str(atom_to_cluster_size.get(atom_key, 0)),
                        str(atom_to_contact_degree.get(atom_key, 0)),
                        f"{atom_to_min_energy.get(atom_key, 0.0):.12g}",
                        f"{atom_to_mean_energy.get(atom_key, 0.0):.12g}",
                        f"{atom_to_max_strength.get(atom_key, 0.0):.12g}",
                        f"{atom_to_mean_strength.get(atom_key, 0.0):.12g}",
                        str(atom_to_domain.get(atom_key, 0)),
                        str(atom_to_pearl.get(atom_key, 0)),
                        str(atom_to_state.get(atom_key, 0)),
                    ]
                )
                handle.write(" ".join(values) + "\n")


def write_cluster_envelope_dump(
    dump_path: Path,
    results: Sequence[FrameAggregationResult],
    output_path: Path,
    padding: float = 0.4,
) -> None:
    """Write one halo sphere per clustered mesogen for OVITO overlay."""
    by_timestep = {int(result.timestep): result for result in results}
    with output_path.open("w", encoding="utf-8") as handle:
        for timestep, box, _columns, col_index, data in parse_dump_frames(dump_path):
            result = by_timestep.get(int(timestep))
            if result is None:
                continue

            clusters = result.visual_cluster_records
            ids = extract_particle_ids(data, col_index)
            positions = extract_positions(data, col_index)
            shapes = extract_shape_axes(data, col_index)
            atom_lookup = {int(atom_id): idx for idx, atom_id in enumerate(ids)}
            halo_rows: List[Tuple[int, int, np.ndarray, float, int, int, float, float]] = []
            halo_id = 1
            for cluster in clusters:
                cluster_id = int(cluster["cluster_id"])
                cluster_size = int(cluster["size"])
                max_radius = float(cluster.get("max_radius", 0.0))
                rg = float(cluster.get("radius_of_gyration", 0.0))
                for atom_id in cluster["atom_ids"]:
                    atom_key = int(atom_id)
                    atom_index = atom_lookup.get(atom_key)
                    if atom_index is None:
                        continue
                    if shapes is not None:
                        base_radius = 0.5 * float(np.max(shapes[atom_index, :]))
                    else:
                        base_radius = 0.5
                    halo_rows.append(
                        (
                            halo_id,
                            atom_key,
                            positions[atom_index, :],
                            max(0.0, base_radius + float(padding)),
                            cluster_id,
                            cluster_size,
                            rg,
                            max_radius,
                        )
                    )
                    halo_id += 1
            handle.write("ITEM: TIMESTEP\n")
            handle.write(f"{int(timestep)}\n")
            handle.write("ITEM: NUMBER OF ATOMS\n")
            handle.write(f"{len(halo_rows)}\n")
            write_box_bounds(handle, box)
            handle.write(
                "ITEM: ATOMS id type x y z radius lc_cluster lc_cluster_size source_atom_id "
                "lc_cluster_rg lc_cluster_max_radius\n"
            )
            for halo_id, atom_id, position, radius, cluster_id, cluster_size, rg, max_radius in halo_rows:
                handle.write(
                    f"{halo_id} 99 "
                    f"{float(position[0]):.12g} {float(position[1]):.12g} {float(position[2]):.12g} "
                    f"{radius:.12g} {cluster_id} {cluster_size} {atom_id} "
                    f"{rg:.12g} {max_radius:.12g}\n"
                )


def write_contact_edge_dump(
    dump_path: Path,
    results: Sequence[FrameAggregationResult],
    output_path: Path,
) -> None:
    """Write a LAMMPS local-style dump of attractive E-E contact edges."""
    by_timestep = {int(result.timestep): result for result in results}
    with output_path.open("w", encoding="utf-8") as handle:
        for timestep, box, _columns, _col_index, _data in parse_dump_frames(dump_path):
            result = by_timestep.get(int(timestep))
            if result is None:
                continue
            edges = result.contact_edge_records
            handle.write("ITEM: TIMESTEP\n")
            handle.write(f"{int(timestep)}\n")
            handle.write("ITEM: NUMBER OF ENTRIES\n")
            handle.write(f"{len(edges)}\n")
            write_box_bounds(handle, box)
            handle.write(
                "ITEM: ENTRIES edge_id atom_i atom_j distance pair_energy well_depth "
                "attraction_strength q_score p2_score delta_s edge_type_code is_local\n"
            )
            for edge in edges:
                handle.write(
                    f"{int(edge['edge_id'])} {int(edge['atom_i'])} {int(edge['atom_j'])} "
                    f"{float(edge['distance']):.12g} {float(edge['pair_energy']):.12g} "
                    f"{float(edge['well_depth']):.12g} {float(edge['attraction_strength']):.12g} "
                    f"{float(edge['q_score']):.12g} {float(edge['p2_score']):.12g} "
                    f"{int(edge['delta_s'])} {int(edge['edge_type_code'])} {int(edge['is_local'])}\n"
                )


def write_contact_segment_vtks(
    dump_path: Path,
    results: Sequence[FrameAggregationResult],
    output_dir: Path,
    safe_name: str,
) -> None:
    """Write OVITO-loadable VTK line segments for attractive E-E contacts."""
    by_timestep = {int(result.timestep): result for result in results}
    for timestep, box, _columns, col_index, data in parse_dump_frames(dump_path):
        result = by_timestep.get(int(timestep))
        if result is None:
            continue
        ids = extract_particle_ids(data, col_index)
        positions = extract_positions(data, col_index)
        id_to_pos = {int(atom_id): positions[idx, :] for idx, atom_id in enumerate(ids)}

        segment_points: List[Tuple[np.ndarray, np.ndarray]] = []
        edge_rows: List[Dict[str, object]] = []
        for edge in result.contact_edge_records:
            atom_i = int(edge["atom_i"])
            atom_j = int(edge["atom_j"])
            if atom_i not in id_to_pos or atom_j not in id_to_pos:
                continue
            left = np.asarray(id_to_pos[atom_i], dtype=float)
            raw_delta = np.asarray(id_to_pos[atom_j], dtype=float) - left
            right = left + minimum_image_vector(raw_delta, box)
            segment_points.append((left, right))
            edge_rows.append(edge)

        vtk_path = output_dir / f"{safe_name}_t{int(timestep)}_lc_contact_segments.vtk"
        with vtk_path.open("w", encoding="utf-8") as handle:
            n_segments = len(segment_points)
            handle.write("# vtk DataFile Version 3.0\n")
            handle.write(f"LC contact segments timestep {int(timestep)}\n")
            handle.write("ASCII\n")
            handle.write("DATASET POLYDATA\n")
            handle.write(f"POINTS {2 * n_segments} float\n")
            for left, right in segment_points:
                handle.write(f"{left[0]:.12g} {left[1]:.12g} {left[2]:.12g}\n")
                handle.write(f"{right[0]:.12g} {right[1]:.12g} {right[2]:.12g}\n")
            handle.write(f"LINES {n_segments} {3 * n_segments}\n")
            for idx in range(n_segments):
                handle.write(f"2 {2 * idx} {2 * idx + 1}\n")
            if n_segments:
                handle.write(f"CELL_DATA {n_segments}\n")
                for name in ("attraction_strength", "pair_energy", "edge_type_code", "is_local"):
                    handle.write(f"SCALARS {name} float 1\n")
                    handle.write("LOOKUP_TABLE default\n")
                    for edge in edge_rows:
                        handle.write(f"{float(edge[name]):.12g}\n")


def write_diagnostic_outputs(
    output_root: Path,
    config: AggregationConfig,
    all_results: Dict[Path, Tuple[List[FrameAggregationResult], Dict[str, object]]],
) -> None:
    """Write audit tables and plots showing where threshold decisions come from."""
    diagnostics_dir = output_root / "diagnostics"
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    edge_rows: List[Dict[str, object]] = []
    domain_rows: List[Dict[str, object]] = []
    pearl_candidate_rows: List[Dict[str, object]] = []
    for source_file, (results, _summary) in all_results.items():
        for result in results:
            for edge in result.contact_edge_records:
                row = dict(edge)
                row["source_file"] = str(source_file)
                row["timestep"] = int(result.timestep)
                edge_rows.append(row)
            for domain in result.domain_records:
                row = dict(domain)
                row["source_file"] = str(source_file)
                row["timestep"] = int(result.timestep)
                evidence = dict(domain.get("evidence", {}))
                for key, value in evidence.items():
                    row[f"evidence_{key}"] = bool(value)
                domain_rows.append(row)
            for candidate in result.pearl_candidate_records:
                row = dict(candidate)
                row["source_file"] = str(source_file)
                row["timestep"] = int(result.timestep)
                pearl_candidate_rows.append(row)

    edge_table = diagnostics_dir / "edge_diagnostics.tsv"
    edge_columns = [
        "source_file",
        "timestep",
        "edge_id",
        "atom_i",
        "atom_j",
        "distance",
        "pair_energy",
        "well_depth",
        "attraction_strength",
        "q_score",
        "p2_score",
        "delta_s",
        "edge_type",
        "edge_type_code",
        "is_local",
    ]
    edge_table_mode = str(config.edge_diagnostics_table)
    written_edge_rows = 0
    if edge_table_mode == "full":
        rows_to_write = edge_rows
    elif edge_table_mode == "sample" and edge_rows:
        sample_size = max(0, int(config.edge_diagnostics_sample_size))
        if sample_size >= len(edge_rows):
            rows_to_write = edge_rows
        elif sample_size > 0:
            indices = np.linspace(0, len(edge_rows) - 1, sample_size, dtype=int)
            rows_to_write = [edge_rows[int(idx)] for idx in indices]
        else:
            rows_to_write = []
    else:
        rows_to_write = []
    if edge_table_mode in {"full", "sample"}:
        with edge_table.open("w", encoding="utf-8") as handle:
            handle.write("\t".join(edge_columns) + "\n")
            for row in rows_to_write:
                handle.write("\t".join(str(row.get(column, "")) for column in edge_columns) + "\n")
        written_edge_rows = len(rows_to_write)
    else:
        if edge_table.exists():
            edge_table.unlink()
        (diagnostics_dir / "edge_diagnostics_table_disabled.txt").write_text(
            "edge_diagnostics.tsv was not written. Default LC-Pearl V2 diagnostics keep summary JSON and histogram plots only. "
            "Use --edge-diagnostics-table sample or --edge-diagnostics-table full if a TSV edge table is explicitly needed.\n",
            encoding="utf-8",
        )

    domain_table = diagnostics_dir / "domain_diagnostics.tsv"
    with domain_table.open("w", encoding="utf-8") as handle:
        columns = [
            "source_file",
            "timestep",
            "domain_id",
            "track_id",
            "age",
            "classification",
            "size",
            "s2",
            "edge_count",
            "adjacent_edge_count",
            "nonlocal_edge_count",
            "stable_under_perturbation",
            "evidence_count",
            "evidence_size",
            "evidence_orientation",
            "evidence_persistence",
            "evidence_nonlocal_support",
            "evidence_parameter_stability",
        ]
        handle.write("\t".join(columns) + "\n")
        for row in domain_rows:
            handle.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")

    pearl_candidate_table = diagnostics_dir / "pearl_candidate_diagnostics.tsv"
    with pearl_candidate_table.open("w", encoding="utf-8") as handle:
        columns = [
            "source_file",
            "timestep",
            "domain_i",
            "domain_j",
            "domain_i_size",
            "domain_j_size",
            "domain_i_atom_ids",
            "domain_j_atom_ids",
            "min_distance",
            "pearl_gap_cut_used",
            "contact_count",
            "left_supported",
            "right_supported",
            "boundary_support_min",
            "aspect_ratio",
            "aspect_ratio_max",
            "adjacency_accepted",
            "reject_reason",
        ]
        handle.write("\t".join(columns) + "\n")
        for row in pearl_candidate_rows:
            handle.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")

    gap_values = [
        float(row.get("pearl_gap_cut_used", 0.0))
        for row in pearl_candidate_rows
        if math.isfinite(float(row.get("pearl_gap_cut_used", 0.0))) and float(row.get("pearl_gap_cut_used", 0.0)) > 0.0
    ]
    current_pearl_gap = float(np.median(np.asarray(gap_values, dtype=float))) if gap_values else resolved_pearl_gap(config, config.manual_r_cut if config.manual_r_cut is not None else 0.0)
    pearl_recommendation = fit_pearl_parameters_from_candidates(
        pearl_candidate_rows,
        current_gap=current_pearl_gap,
        current_min_cross_contacts=config.pearl_min_cross_contacts,
        current_min_boundary_particles=config.pearl_min_boundary_particles,
        current_max_aspect_ratio=config.pearl_max_aspect_ratio,
        min_candidates=5,
    )
    (diagnostics_dir / "pearl_parameter_recommendations.json").write_text(
        json.dumps(pearl_recommendation, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    sweep_table = diagnostics_dir / "edge_threshold_sweep.tsv"
    strength_values = sorted({0.12, 0.20, 0.30, 0.40, 0.50, config.gb_off_strength, config.gb_on_strength})
    p2_values = sorted({0.0, 0.3, 0.5, config.p2_cut, 0.85})
    with sweep_table.open("w", encoding="utf-8") as handle:
        handle.write("gb_strength_cut\tp2_cut\tedge_count\tnonlocal_edge_count\tlocal_edge_count\n")
        for strength_cut in strength_values:
            for p2_cut in p2_values:
                selected = [
                    row for row in edge_rows
                    if float(row.get("attraction_strength", 0.0)) >= strength_cut
                    and float(row.get("p2_score", -1.0)) >= p2_cut
                ]
                local_count = sum(1 for row in selected if int(row.get("is_local", 0)) == 1)
                handle.write(
                    f"{strength_cut:.6g}\t{p2_cut:.6g}\t{len(selected)}\t"
                    f"{len(selected) - local_count}\t{local_count}\n"
                )

    robust_domains = [row for row in domain_rows if row.get("classification") == "robust"]
    summary = {
        "thresholds": {
            "gb_threshold_mode": config.gb_threshold_mode,
            "gb_require_pair_orientation": config.gb_require_pair_orientation,
            "gb_on_strength": config.gb_on_strength,
            "gb_off_strength": config.gb_off_strength,
            "u_on": config.u_on,
            "u_off": config.u_off,
            "p2_cut": config.p2_cut,
            "robust_min_s2": config.robust_min_s2,
            "robust_require_orientation": config.robust_require_orientation,
            "n_min": config.n_min,
            "robust_min_evidence": config.robust_min_evidence,
            "pearl_gap_cut_used": current_pearl_gap,
            "pearl_min_cross_contacts": config.pearl_min_cross_contacts,
            "pearl_min_boundary_particles": config.pearl_min_boundary_particles,
            "pearl_max_aspect_ratio": config.pearl_max_aspect_ratio,
        },
        "counts": {
            "edge_count": len(edge_rows),
            "strong_edge_count": sum(1 for row in edge_rows if row.get("edge_type") == "strong"),
            "local_edge_count": sum(1 for row in edge_rows if int(row.get("is_local", 0)) == 1),
            "domain_count": len(domain_rows),
            "robust_domain_count": len(robust_domains),
            "robust_below_s2_threshold_count": sum(
                1 for row in robust_domains if float(row.get("s2", 0.0)) < config.robust_min_s2
            ),
            "pearl_candidate_pair_count": len(pearl_candidate_rows),
            "adjacency_accepted_pearl_candidate_pair_count": sum(1 for row in pearl_candidate_rows if bool(row.get("adjacency_accepted"))),
        },
        "edge_diagnostics_table": {
            "mode": edge_table_mode,
            "total_edge_rows": len(edge_rows),
            "written_edge_rows": written_edge_rows,
            "sample_size": int(config.edge_diagnostics_sample_size),
            "path": str(edge_table) if edge_table_mode in {"full", "sample"} else None,
        },
        "pearl_parameter_recommendation": pearl_recommendation,
        "notes": [
            "GB pair energy and S2 are computed quantities.",
            "Main threshold selection should come from lc_threshold_prior.py streaming full candidate-pair histograms, not from screened accepted-edge diagnostics.",
            "special_bonds/topology exclusions are not exact unless an explicit topology workflow is added.",
        ],
    }
    (diagnostics_dir / "diagnostic_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        (diagnostics_dir / "plot_error.txt").write_text(str(exc), encoding="utf-8")
        return

    if edge_rows:
        strengths = np.array([float(row.get("attraction_strength", 0.0)) for row in edge_rows], dtype=float)
        p2_values_array = np.array([float(row.get("p2_score", 0.0)) for row in edge_rows], dtype=float)

        fig, ax = plt.subplots(figsize=(7, 4.2), dpi=160)
        ax.hist(strengths, bins=40, color="#2c7a5b", alpha=0.78)
        annotate_vline(ax, config.gb_off_strength, threshold_label("gb_off", config.gb_off_strength), color="#b7791f", lw=2)
        annotate_vline(ax, config.gb_on_strength, threshold_label("gb_on", config.gb_on_strength), color="#9b2c2c", ymax=0.82, lw=2)
        ax.set_xlabel("GB attraction strength = max(0, -U_GB / U_well)")
        ax.set_ylabel("edge count")
        ax.legend(frameon=False, loc="upper right")
        fig.tight_layout()
        fig.savefig(diagnostics_dir / "gb_strength_hist.png")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 5.0), dpi=160)
        hb = ax.hexbin(strengths, p2_values_array, gridsize=70, bins="log", mincnt=1, cmap="viridis")
        annotate_vline(ax, config.gb_on_strength, threshold_label("gb_on", config.gb_on_strength), color="#9b2c2c", lw=1.5)
        annotate_hline(ax, config.p2_cut, threshold_label("p2_cut", config.p2_cut), color="#b7791f", lw=1.5)
        ax.set_xlabel("GB attraction strength")
        ax.set_ylabel("pair P2")
        ax.set_title("Accepted LC edges under current gates")
        fig.colorbar(hb, ax=ax, label="log10(edge count)")
        ax.legend(frameon=False, loc="upper right")
        fig.tight_layout()
        fig.savefig(diagnostics_dir / "gb_strength_vs_p2.png")
        plt.close(fig)

    if domain_rows:
        sizes = np.array([float(row.get("size", 0.0)) for row in domain_rows], dtype=float)
        s2_values = np.array([float(row.get("s2", 0.0)) for row in domain_rows], dtype=float)
        is_robust = np.array([row.get("classification") == "robust" for row in domain_rows], dtype=bool)
        colors = np.where(is_robust, "#9b2c2c", "#66716b")
        fig, ax = plt.subplots(figsize=(6.4, 5.0), dpi=160)
        ax.scatter(sizes, s2_values, s=28, alpha=0.78, c=colors)
        annotate_hline(ax, config.robust_min_s2, threshold_label("S2_min", config.robust_min_s2), color="#b7791f", lw=1.5)
        size_floor = max(config.n_min, config.robust_min_size)
        annotate_vline(ax, float(size_floor), f"size_floor={size_floor}", color="#315f9f", lw=1.5)
        ax.set_xlabel("domain size")
        ax.set_ylabel("domain S2")
        ax.legend(frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(diagnostics_dir / "domain_size_vs_s2.png")
        plt.close(fig)


def write_results(
    output_root: Path,
    config: AggregationConfig,
    all_results: Dict[Path, Tuple[List[FrameAggregationResult], Dict[str, object]]],
) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    per_file_dir = output_root / "per_file"
    per_file_dir.mkdir(parents=True, exist_ok=True)

    combined_tsv = output_root / "aggregation_timeseries.tsv"
    header = "\t".join(
        [
            "source_file",
            "timestep",
            "axis_used",
            "N",
            "r_cut",
            "aggregation_degree",
            "largest_cluster_fraction",
            "clustered_fraction",
            "core_fraction",
            "mean_cluster_s2",
            "max_cluster_s2",
            "global_s2",
            "n_clusters",
            "qualified_pairs",
            "qualified_pair_fraction",
            "energy_edge_count",
            "min_pair_energy",
            "mean_pair_energy",
            "cluster_sizes",
            "cluster_s2",
            "visual_cluster_count",
            "largest_visual_cluster_fraction",
            "visual_cluster_sizes",
            "weak_domain_count",
            "robust_domain_count",
            "weak_domain_fraction",
            "robust_domain_fraction",
            "weak_domain_sizes",
            "robust_domain_sizes",
            "robust_domain_s2",
            "pearl_count",
            "largest_pearl_fraction",
            "pearl_sizes",
            "pearl_domain_counts",
            "pearl_compactness",
            "pearl_aspect_ratios",
            "pearl_axis_centers",
            "pearl_axis_widths",
            "connector_lengths",
            "axial_axis",
            "local_edge_fraction",
            "nonlocal_edge_fraction",
            "ambiguous_mesogen_fraction",
            "stable_core_pair_fraction",
            "L_parallel",
            "Rg_parallel",
            "Rg_perp",
            "S2_force",
            "stretch_axis",
        ]
    ) + "\n"

    with combined_tsv.open("w", encoding="utf-8") as combined:
        combined.write(header)
        source_file_names: List[str] = []
        total_frames_processed = 0
        output_stems = build_output_stems(list(all_results.keys()))
        for source_file, (results, summary) in sorted(
            all_results.items(),
            key=lambda item: (
                min((result.timestep for result in item[1][0]), default=0),
                dump_sort_key(item[0]),
            ),
        ):
            safe_name = output_stems[source_file]
            file_tsv = per_file_dir / f"{safe_name}_aggregation.tsv"
            file_json = per_file_dir / f"{safe_name}_summary.json"
            frames_jsonl = per_file_dir / f"{safe_name}_frames.jsonl"
            ovito_dump = per_file_dir / f"{safe_name}_lc_labels.dump"
            envelope_dump = per_file_dir / f"{safe_name}_lc_cluster_envelopes.dump"
            contact_edge_dump = per_file_dir / f"{safe_name}_lc_contact_edges.dump"

            with file_tsv.open("w", encoding="utf-8") as handle:
                handle.write(header)
                for result in results:
                    row = result_to_tsv_row(result)
                    handle.write(row)
                    combined.write(row)

            file_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
            if config.write_frame_jsonl:
                with frames_jsonl.open("w", encoding="utf-8") as handle:
                    for result in results:
                        handle.write(json.dumps(asdict(result), ensure_ascii=False) + "\n")
            elif frames_jsonl.exists():
                frames_jsonl.unlink()
            if config.write_ovito_labels:
                write_ovito_label_dump(Path(source_file), results, ovito_dump)
            if config.write_cluster_envelopes:
                write_cluster_envelope_dump(
                    Path(source_file),
                    results,
                    envelope_dump,
                    padding=config.cluster_envelope_padding,
                )
            if config.write_contact_edges:
                write_contact_edge_dump(Path(source_file), results, contact_edge_dump)
            if config.write_contact_segments:
                write_contact_segment_vtks(Path(source_file), results, per_file_dir, safe_name)
            source_file_names.append(str(source_file))
            total_frames_processed += int(summary.get("frames_processed", len(results)))

    run_summary = {
        "summary_mode": "compact",
        "config": asdict(config),
        "n_files": len(all_results),
        "frames_processed": int(total_frames_processed),
        "outputs": {
            "aggregation_timeseries": str(combined_tsv),
            "per_file_directory": str(per_file_dir),
            "per_file_summary_glob": str(per_file_dir / "*_summary.json"),
            "per_file_frames_jsonl_glob": str(per_file_dir / "*_frames.jsonl") if config.write_frame_jsonl else None,
            "ovito_label_dump_glob": str(per_file_dir / "*_lc_labels.dump") if config.write_ovito_labels else None,
            "cluster_envelope_dump_glob": str(per_file_dir / "*_lc_cluster_envelopes.dump") if config.write_cluster_envelopes else None,
            "contact_edge_dump_glob": str(per_file_dir / "*_lc_contact_edges.dump") if config.write_contact_edges else None,
            "contact_segment_vtk_glob": str(per_file_dir / "*_lc_contact_segments.vtk") if config.write_contact_segments else None,
            "diagnostics_directory": str(output_root / "diagnostics") if config.write_diagnostics else None,
        },
        "source_files": {
            "count": len(source_file_names),
            "first": source_file_names[0] if source_file_names else None,
            "last": source_file_names[-1] if source_file_names else None,
            "sample": source_file_names[:5],
        },
    }
    (output_root / "run_summary.json").write_text(
        json.dumps(run_summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    if config.write_diagnostics:
        write_diagnostic_outputs(output_root, config=config, all_results=all_results)


def parse_r_cut(raw_value: str) -> Tuple[str, Optional[float]]:
    lowered = raw_value.strip().lower()
    if lowered == "auto":
        return "auto", None
    try:
        numeric = float(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("r_cut must be a positive float or 'auto'") from exc
    if numeric <= 0.0:
        raise argparse.ArgumentTypeError("r_cut must be positive")
    return "manual", numeric


def parse_worker_count(raw_value: str) -> int:
    lowered = raw_value.strip().lower()
    if lowered == "auto":
        return 0
    try:
        worker_count = int(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("workers must be a positive integer or 'auto'") from exc
    if worker_count <= 0:
        raise argparse.ArgumentTypeError("workers must be a positive integer or 'auto'")
    return worker_count


def resolve_worker_count(worker_count: int, n_files: int, cpu_count: Optional[int] = None) -> int:
    if n_files <= 1:
        return 1
    max_workers = max(1, int(os.environ.get("LC_PEARL_MAX_WORKERS", DEFAULT_MAX_WORKERS)))
    if worker_count == 0:
        available_cpus = int(cpu_count if cpu_count is not None else (os.cpu_count() or 1))
        return max(1, min(int(n_files), max(1, available_cpus), max_workers))
    return max(1, min(int(n_files), int(worker_count), max_workers))


def parse_optional_positive_float(raw_value: str) -> Optional[float]:
    lowered = raw_value.strip().lower()
    if lowered == "auto":
        return None
    try:
        numeric = float(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be a positive float or 'auto'") from exc
    if numeric <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return numeric


def parse_optional_unit_float(raw_value: str) -> Optional[float]:
    lowered = raw_value.strip().lower()
    if lowered == "auto":
        return None
    try:
        numeric = float(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be in [0, 1] or 'auto'") from exc
    if not 0.0 <= numeric <= 1.0:
        raise argparse.ArgumentTypeError("value must be in [0, 1]")
    return numeric


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute LC weak/robust domains and 3D pearls from LAMMPS dump files.")
    parser.add_argument("inputs", nargs="*", type=Path, help="Dump files or directories. Defaults to the current directory.")
    parser.add_argument("--pattern", default="*.dump", help="Glob used when an input is a directory.")
    parser.add_argument("--recursive", action="store_true", help="Recursively scan directory inputs.")
    parser.add_argument("--axis", default="auto", choices=["auto", "x", "y", "z"], help="Particle long axis in body frame.")
    parser.add_argument("--every", type=int, default=1, help="Process every k-th frame.")
    parser.add_argument("--r-cut", default="auto", help="Spatial cutoff for support contacts. Pass a positive float or 'auto'.")
    parser.add_argument("--p2-cut", type=float, default=0.7, help="Nematic pair threshold in P2 space.")
    parser.add_argument("--min-core-neighbors", type=int, default=2, help="Legacy local degree threshold reported as core_fraction.")
    parser.add_argument("--cutoff-bins", type=int, default=120, help="Radial bins used when estimating r_cut automatically.")
    parser.add_argument("--cutoff-frames", type=int, default=5, help="Maximum sampled frames used for automatic cutoff estimation.")
    parser.add_argument("--max-auto-r-cut", default="auto", help="Upper bound for automatically estimated r_cut; 'auto' uses mesogen shape columns when available.")
    parser.add_argument("--auto-r-cut-shape-factor", type=float, default=1.8, help="When shape columns exist, cap auto r_cut at this factor times the mesogen long diameter.")
    parser.add_argument("--per-file-r-cut", dest="shared_r_cut", action="store_false", default=True, help="Estimate auto r_cut separately for each dump instead of sharing one value across the batch.")
    parser.add_argument("--mesogen-type", type=int, default=1, help="Atom type used as mesogen aggregation members.")
    parser.add_argument("--anchor-types", default="2,3", help="Comma-separated atom types used only for endpoints/stretch direction.")
    parser.add_argument("--contact-mode", default="center", choices=["center", "ellipsoid", "gayberne"], help="Contact definition: center uses r_ij/r_cut; ellipsoid uses a shape-based contact estimate; gayberne uses type 1-1 GB pair energy.")
    parser.add_argument("--gb-param-file", type=Path, default=None, help="LAMMPS input file containing pair_style/pair_coeff parameters for --contact-mode gayberne.")
    parser.add_argument("--gb-threshold-mode", default="relative", choices=["absolute", "relative"], help="Gay-Berne contact threshold: relative uses --gb-on-strength/--gb-off-strength on -U/well_depth; absolute uses --u-on/--u-off.")
    parser.add_argument("--gb-pair-orientation", dest="gb_require_pair_orientation", action="store_true", default=True, help="Require pair orientation thresholds in Gay-Berne contact mode. Enabled by default.")
    parser.add_argument("--no-gb-pair-orientation", dest="gb_require_pair_orientation", action="store_false", help="Disable pair orientation filtering in Gay-Berne contact mode.")
    parser.add_argument("--u-on", type=float, default=-0.30, help="Strong attractive edge threshold for Gay-Berne energy.")
    parser.add_argument("--u-off", type=float, default=-0.12, help="Gray attractive edge threshold for Gay-Berne energy; must be less negative than --u-on.")
    parser.add_argument("--gb-on-strength", type=float, default=0.30, help="Strong edge threshold for relative Gay-Berne attraction strength -U/well_depth.")
    parser.add_argument("--gb-off-strength", type=float, default=0.12, help="Gray edge threshold for relative Gay-Berne attraction strength -U/well_depth.")
    parser.add_argument("--r-energy-cap", default="auto", help="Maximum center distance for Gay-Berne candidate pairs; 'auto' uses the GB cutoff from --gb-param-file.")
    parser.add_argument("--cluster-cut", default="auto", help="Visual cluster center-distance cutoff for OVITO labels; 'auto' uses min(r_cut, shape_factor * long diameter).")
    parser.add_argument("--cluster-cut-shape-factor", type=float, default=1.35, help="When shape columns exist, auto cluster_cut is capped at this factor times the mesogen long diameter.")
    parser.add_argument("--cluster-min-size", type=int, default=2, help="Minimum mesogen count assigned a nonzero visual lc_cluster label.")
    parser.add_argument("--g-on", type=float, default=1.0, help="Strong contact threshold for normalized contact g_ij.")
    parser.add_argument("--g-off", type=float, default=1.25, help="Gray contact threshold for normalized contact g_ij.")
    parser.add_argument("--q-on", default="auto", help="Strong orientational threshold for |u_i dot u_j|; auto derives from --p2-cut.")
    parser.add_argument("--q-off", default="auto", help="Gray orientational threshold for |u_i dot u_j|; auto is q_on - 0.10.")
    parser.add_argument("--s-excl", type=int, default=1, help="Sequence separation at or below this value is local support only.")
    parser.add_argument("--n-min", type=int, default=3, help="Minimum mesogen count for robust domains.")
    parser.add_argument("--robust-min-size", type=int, default=3, help="Compatibility size floor; v1 primary parameter is --n-min.")
    parser.add_argument("--robust-min-s2", type=float, default=0.70, help="Internal S2 threshold counted as orientational evidence.")
    parser.add_argument("--robust-min-evidence", type=int, default=4, help="Minimum supporting evidence count retained for robust-domain audit.")
    parser.add_argument("--domain-min-lifetime", type=int, default=2, help="Processed-frame age counted as persistence evidence.")
    parser.add_argument("--require-robust-orientation", dest="robust_require_orientation", action="store_true", default=True, help="Require domain S2 >= robust_min_s2 for robust classification. Enabled by default.")
    parser.add_argument("--no-require-robust-orientation", dest="robust_require_orientation", action="store_false", help="Allow S2 to be evidence only, not a hard robust-domain gate.")
    parser.add_argument("--require-nonlocal-robust", action="store_true", help="Require at least one nonlocal support contact for robust classification.")
    parser.add_argument("--allow-adjacent-only-robust", action="store_true", help="Compatibility flag; adjacent-only robust domains are allowed by default when other evidence passes.")
    parser.add_argument("--adjacent-id-gap", type=int, default=1, help="Atom-id gap treated as topological adjacency; use 0 to disable.")
    parser.add_argument("--local-pair-file", type=Path, default=None, help="Optional two-column atom_i atom_j file treated as explicit local support pairs.")
    parser.add_argument("--exclude-pair-file", type=Path, default=None, help="Optional two-column atom_i atom_j file excluded from contact/domain analysis, e.g. special_bonds-zeroed pairs.")
    parser.add_argument("--perturbation-r-cut-scale", type=float, default=0.95, help="Stricter r_cut scale for parameter-stability evidence.")
    parser.add_argument("--perturbation-p2-margin", type=float, default=0.05, help="Stricter P2 margin for parameter-stability evidence.")
    parser.add_argument("--stable-overlap-fraction", type=float, default=0.75, help="Largest retained fraction needed under stricter parameters.")
    parser.add_argument("--pearl-gap-cut", default="auto", help="3D gap for merging robust domains into pearls; positive float or 'auto'=r_cut.")
    parser.add_argument("--pearl-min-cross-contacts", type=int, default=2, help="Minimum spatial cross-domain contact pairs for pearl merge.")
    parser.add_argument("--pearl-min-boundary-particles", type=int, default=2, help="Minimum supported particles on each domain boundary for pearl merge.")
    parser.add_argument("--pearl-max-aspect-ratio", type=float, default=3.0, help="Maximum merged 3D aspect ratio still treated as bead-like.")
    parser.add_argument("--track-jaccard", type=float, default=0.50, help="Jaccard overlap needed to continue a domain track.")
    parser.add_argument("--consensus-threshold", type=float, default=0.70, help="Pair co-membership threshold for stable core in robustness scan.")
    parser.add_argument("--robustness-scan", action="store_true", help="Enable small parameter-grid consensus audit per frame.")
    parser.add_argument("--write-ovito-labels", dest="write_ovito_labels", action="store_true", default=True, help="Write per-atom LC label dumps for OVITO coloring/filtering. Enabled by default.")
    parser.add_argument("--no-ovito-labels", dest="write_ovito_labels", action="store_false", help="Disable OVITO label dump output.")
    parser.add_argument("--write-cluster-envelopes", dest="write_cluster_envelopes", action="store_true", default=True, help="Write one OVITO overlay sphere per visual cluster. Enabled by default.")
    parser.add_argument("--no-cluster-envelopes", dest="write_cluster_envelopes", action="store_false", help="Disable visual cluster envelope dump output.")
    parser.add_argument("--write-contact-edges", dest="write_contact_edges", action="store_true", default=True, help="Write local-style dumps of attractive E-E contact edges. Enabled by default.")
    parser.add_argument("--no-contact-edges", dest="write_contact_edges", action="store_false", help="Disable contact edge dump output.")
    parser.add_argument("--write-contact-segments", dest="write_contact_segments", action="store_true", default=True, help="Write OVITO-loadable VTK line segments for attractive contacts. Enabled by default.")
    parser.add_argument("--no-contact-segments", dest="write_contact_segments", action="store_false", help="Disable VTK contact segment output.")
    parser.add_argument("--write-frame-jsonl", dest="write_frame_jsonl", action="store_true", default=False, help="Write full per-frame JSONL debug records. Disabled by default because records can be very large.")
    parser.add_argument("--no-frame-jsonl", dest="write_frame_jsonl", action="store_false", help="Disable full per-frame JSONL debug output.")
    parser.add_argument("--write-diagnostics", dest="write_diagnostics", action="store_true", default=True, help="Write threshold audit TSV/JSON/PNG diagnostics. Enabled by default.")
    parser.add_argument("--no-diagnostics", dest="write_diagnostics", action="store_false", help="Disable threshold diagnostic outputs.")
    parser.add_argument("--edge-diagnostics-table", choices=["off", "sample", "full"], default="off", help="Write accepted-edge TSV diagnostics. Default off writes summary JSON and histogram plots only.")
    parser.add_argument("--edge-diagnostics-sample-size", type=int, default=200_000, help="Rows written when --edge-diagnostics-table sample is used.")
    parser.add_argument("--cluster-envelope-padding", type=float, default=0.4, help="Padding added to each visual cluster envelope radius.")
    parser.add_argument("--track-across-files", dest="track_across_files", action="store_true", default=True, help="Track domain persistence across sorted dump files. Enabled by default for one-frame dump sequences.")
    parser.add_argument("--no-track-across-files", dest="track_across_files", action="store_false", help="Process each dump independently; enables full file-level parallelism but disables persistence across one-frame files.")
    parser.add_argument("--workers", type=parse_worker_count, default=0, help="Parallel dump-file workers. Use 'auto' for CPU count capped by file count and LC_PEARL_MAX_WORKERS (default 10), or 1 for serial.")
    parser.add_argument("--output-root", type=Path, default=Path("lc_domain_pearl_v2_output"), help="Output directory for TSV, JSON, and OVITO label dumps.")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.every <= 0:
        raise SystemExit("--every must be positive")
    if args.cutoff_bins <= 0:
        raise SystemExit("--cutoff-bins must be positive")
    if args.cutoff_frames <= 0:
        raise SystemExit("--cutoff-frames must be positive")
    if args.auto_r_cut_shape_factor <= 0.0:
        raise SystemExit("--auto-r-cut-shape-factor must be positive")
    if not -0.5 <= args.p2_cut <= 1.0:
        raise SystemExit("--p2-cut must be in [-0.5, 1.0]")
    if args.min_core_neighbors < 1:
        raise SystemExit("--min-core-neighbors must be at least 1")
    if args.contact_mode == "gayberne" and args.gb_param_file is None:
        raise SystemExit("--contact-mode gayberne requires --gb-param-file")
    if args.contact_mode == "gayberne" and args.gb_threshold_mode == "absolute" and not (args.u_on < args.u_off < 0.0):
        raise SystemExit("--contact-mode gayberne requires --u-on < --u-off < 0")
    if args.contact_mode == "gayberne" and args.gb_threshold_mode == "relative" and not (1.0 >= args.gb_on_strength > args.gb_off_strength > 0.0):
        raise SystemExit("--contact-mode gayberne requires 1 >= --gb-on-strength > --gb-off-strength > 0")
    if args.cluster_cut_shape_factor <= 0.0:
        raise SystemExit("--cluster-cut-shape-factor must be positive")
    if args.cluster_min_size < 2:
        raise SystemExit("--cluster-min-size must be at least 2")
    if args.cluster_envelope_padding < 0.0:
        raise SystemExit("--cluster-envelope-padding must be non-negative")
    if args.g_on <= 0.0 or args.g_off <= 0.0 or args.g_off <= args.g_on:
        raise SystemExit("--g-off must be greater than --g-on, and both must be positive")
    try:
        parsed_anchor_types = [int(item.strip()) for item in args.anchor_types.split(",") if item.strip()]
    except ValueError as exc:
        raise SystemExit("--anchor-types must be comma-separated integers") from exc
    if not parsed_anchor_types:
        raise SystemExit("--anchor-types must contain at least one atom type")
    q_on = parse_optional_unit_float(args.q_on)
    q_off = parse_optional_unit_float(args.q_off)
    if q_on is not None and q_off is not None and q_off >= q_on:
        raise SystemExit("--q-off must be smaller than --q-on")
    if args.s_excl < 0:
        raise SystemExit("--s-excl must be non-negative")
    if args.n_min < 2:
        raise SystemExit("--n-min must be at least 2")
    if args.robust_min_size < 2:
        raise SystemExit("--robust-min-size must be at least 2")
    if not -0.5 <= args.robust_min_s2 <= 1.0:
        raise SystemExit("--robust-min-s2 must be in [-0.5, 1.0]")
    if not 1 <= args.robust_min_evidence <= 5:
        raise SystemExit("--robust-min-evidence must be in [1, 5]")
    if args.domain_min_lifetime < 1:
        raise SystemExit("--domain-min-lifetime must be at least 1")
    if args.adjacent_id_gap < 0:
        raise SystemExit("--adjacent-id-gap must be non-negative")
    if not 0.0 < args.perturbation_r_cut_scale <= 1.0:
        raise SystemExit("--perturbation-r-cut-scale must be in (0, 1]")
    if args.perturbation_p2_margin < 0.0:
        raise SystemExit("--perturbation-p2-margin must be non-negative")
    if not 0.0 < args.stable_overlap_fraction <= 1.0:
        raise SystemExit("--stable-overlap-fraction must be in (0, 1]")
    if args.pearl_min_cross_contacts < 1:
        raise SystemExit("--pearl-min-cross-contacts must be at least 1")
    if args.pearl_min_boundary_particles < 1:
        raise SystemExit("--pearl-min-boundary-particles must be at least 1")
    if args.pearl_max_aspect_ratio < 1.0:
        raise SystemExit("--pearl-max-aspect-ratio must be at least 1")
    if not 0.0 <= args.track_jaccard <= 1.0:
        raise SystemExit("--track-jaccard must be in [0, 1]")
    if not 0.0 <= args.consensus_threshold <= 1.0:
        raise SystemExit("--consensus-threshold must be in [0, 1]")
    if args.edge_diagnostics_sample_size < 0:
        raise SystemExit("--edge-diagnostics-sample-size must be non-negative")


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()
    validate_args(args)

    r_cut_mode, manual_r_cut = parse_r_cut(args.r_cut)
    max_auto_r_cut = parse_optional_positive_float(args.max_auto_r_cut)
    pearl_gap_cut = parse_optional_positive_float(args.pearl_gap_cut)
    cluster_cut = parse_optional_positive_float(args.cluster_cut)
    r_energy_cap = parse_optional_positive_float(args.r_energy_cap)
    gayberne_params = parse_gayberne_params_from_lmp(args.gb_param_file) if args.contact_mode == "gayberne" else None
    q_on = parse_optional_unit_float(args.q_on)
    q_off = parse_optional_unit_float(args.q_off)
    anchor_types = tuple(
        int(item.strip())
        for item in args.anchor_types.split(",")
        if item.strip()
    )
    config = AggregationConfig(
        axis=args.axis,
        every=int(args.every),
        p2_cut=float(args.p2_cut),
        min_core_neighbors=int(args.min_core_neighbors),
        cutoff_bins=int(args.cutoff_bins),
        cutoff_frames=int(args.cutoff_frames),
        r_cut_mode=r_cut_mode,
        manual_r_cut=manual_r_cut,
        max_auto_r_cut=max_auto_r_cut,
        auto_r_cut_shape_factor=float(args.auto_r_cut_shape_factor),
        mesogen_type=int(args.mesogen_type),
        anchor_types=anchor_types,
        contact_mode=args.contact_mode,
        gayberne_params=gayberne_params,
        gb_threshold_mode=args.gb_threshold_mode,
        gb_require_pair_orientation=bool(args.gb_require_pair_orientation),
        u_on=float(args.u_on),
        u_off=float(args.u_off),
        gb_on_strength=float(args.gb_on_strength),
        gb_off_strength=float(args.gb_off_strength),
        r_energy_cap=r_energy_cap,
        cluster_cut=cluster_cut,
        cluster_cut_shape_factor=float(args.cluster_cut_shape_factor),
        cluster_min_size=int(args.cluster_min_size),
        g_on=float(args.g_on),
        g_off=float(args.g_off),
        q_on=q_on,
        q_off=q_off,
        s_excl=int(args.s_excl),
        n_min=int(args.n_min),
        robust_min_size=int(args.robust_min_size),
        robust_min_s2=float(args.robust_min_s2),
        robust_min_evidence=int(args.robust_min_evidence),
        robust_require_orientation=bool(args.robust_require_orientation),
        robust_require_nonlocal=bool(args.require_nonlocal_robust) and not bool(args.allow_adjacent_only_robust),
        domain_min_lifetime=int(args.domain_min_lifetime),
        adjacent_id_gap=int(args.adjacent_id_gap),
        local_pair_file=str(args.local_pair_file.resolve()) if args.local_pair_file is not None else None,
        exclude_pair_file=str(args.exclude_pair_file.resolve()) if args.exclude_pair_file is not None else None,
        perturbation_r_cut_scale=float(args.perturbation_r_cut_scale),
        perturbation_p2_margin=float(args.perturbation_p2_margin),
        stable_overlap_fraction=float(args.stable_overlap_fraction),
        pearl_gap_cut=pearl_gap_cut,
        pearl_min_cross_contacts=int(args.pearl_min_cross_contacts),
        pearl_min_boundary_particles=int(args.pearl_min_boundary_particles),
        pearl_max_aspect_ratio=float(args.pearl_max_aspect_ratio),
        track_jaccard=float(args.track_jaccard),
        consensus_threshold=float(args.consensus_threshold),
        enable_robustness_scan=bool(args.robustness_scan),
        write_ovito_labels=bool(args.write_ovito_labels),
        write_cluster_envelopes=bool(args.write_cluster_envelopes),
        write_contact_edges=bool(args.write_contact_edges),
        write_contact_segments=bool(args.write_contact_segments),
        write_frame_jsonl=bool(args.write_frame_jsonl),
        write_diagnostics=bool(args.write_diagnostics),
        edge_diagnostics_table=str(args.edge_diagnostics_table),
        edge_diagnostics_sample_size=int(args.edge_diagnostics_sample_size),
        shared_r_cut=bool(args.shared_r_cut),
        track_across_files=bool(args.track_across_files),
        workers=int(args.workers),
        cluster_envelope_padding=float(args.cluster_envelope_padding),
    )

    input_paths = args.inputs if args.inputs else [Path.cwd()]
    files = iter_input_files(input_paths, pattern=args.pattern, recursive=bool(args.recursive))
    if not files:
        raise SystemExit("No dump files matched the given inputs.")

    worker_count = resolve_worker_count(config.workers, n_files=len(files))
    config = replace(config, workers=worker_count)

    if config.r_cut_mode == "auto" and config.shared_r_cut:
        shared_r_cut = determine_shared_r_cut_for_files(files, config)
        config = replace(config, r_cut_mode="manual", manual_r_cut=shared_r_cut)
        print(f"[INFO] shared r_cut={shared_r_cut:.6f} estimated from up to {args.cutoff_frames} processed frame(s)")

    if config.track_across_files and len(files) > 1:
        print(f"[INFO] processing {len(files)} dump file(s) sequentially with cross-file domain tracking")
        processed_iter = iter_processed_dump_file_sequence(files, config=config)
    else:
        print(f"[INFO] processing {len(files)} dump file(s) with {worker_count} worker(s)")
        processed_iter = iter_processed_dump_files(files, config=config, worker_count=worker_count)

    all_results: Dict[Path, Tuple[List[FrameAggregationResult], Dict[str, object]]] = {}
    for dump_file, results, summary in processed_iter:
        all_results[dump_file] = (results, summary)
        print(
            f"[OK] {dump_file.name}: frames={summary['frames_processed']}, "
            f"r_cut={summary['r_cut']:.6f}, "
            f"mean_A_lc={summary['mean_aggregation_degree']:.6f}, "
            f"mean_robust_domains={summary['mean_robust_domain_count']:.2f}, "
            f"mean_pearls={summary['mean_pearl_count']:.2f}"
        )

    write_results(args.output_root.resolve(), config=config, all_results=all_results)
    print(f"Results written to {args.output_root.resolve()}")


if __name__ == "__main__":
    main()
