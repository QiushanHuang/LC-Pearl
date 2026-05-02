#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import liquid_crystal_aggregation as lca
import lc_domain_pearl_pipeline as pipeline


SCHEMA_VERSION = 1
MARKER = "LC_RUN0_PAIR_ENERGY"


def lammps_boundary_token(token: str) -> str:
    if token.startswith("p"):
        return "p"
    if token.startswith("f"):
        return "f"
    if token.startswith("s"):
        return "s"
    if token.startswith("m"):
        return "m"
    return "p"


def wrapped_position_for_lammps(position: np.ndarray, box: lca.BoxSpec) -> np.ndarray:
    wrapped = np.asarray(position, dtype=float).copy()
    for axis, periodic in enumerate(box.periodic):
        if not periodic:
            continue
        lo, hi = box.bounds[axis]
        length = hi - lo
        if length <= 0.0:
            continue
        wrapped[axis] = lo + ((wrapped[axis] - lo) % length)
    return wrapped


def quaternion_to_set_quat_args(quaternion: np.ndarray) -> Tuple[float, float, float, float]:
    q = np.asarray(quaternion, dtype=float)
    norm = float(np.linalg.norm(q))
    if norm <= 1e-12:
        return 1.0, 0.0, 0.0, 0.0
    q = q / norm
    if q[0] < 0.0:
        q = -q
    w = max(-1.0, min(1.0, float(q[0])))
    sin_half = math.sqrt(max(0.0, 1.0 - w * w))
    if sin_half <= 1e-12:
        return 1.0, 0.0, 0.0, 0.0
    axis = q[1:4] / sin_half
    angle_degrees = 2.0 * math.acos(w) * 180.0 / math.pi
    return float(axis[0]), float(axis[1]), float(axis[2]), float(angle_degrees)


def frame_by_index(dump_file: Path, frame_index: int) -> Tuple[int, lca.BoxSpec, List[str], Dict[str, int], np.ndarray]:
    for idx, frame in enumerate(lca.parse_dump_frames(dump_file)):
        if idx == int(frame_index):
            return frame
    raise RuntimeError(f"{dump_file}: frame index {frame_index} not found")


def extract_type1_frame(
    dump_file: Path,
    *,
    frame_index: int,
    mesogen_type: int,
) -> Dict[str, object]:
    timestep, box, columns, col_index, data = frame_by_index(dump_file, frame_index)
    positions = lca.extract_positions(data, col_index)
    ids = lca.extract_particle_ids(data, col_index).astype(int)
    types = lca.extract_particle_types(data, col_index)
    shapes = lca.extract_shape_axes(data, col_index)
    quats = lca.extract_quaternions(data, col_index)
    if shapes is None:
        raise RuntimeError(f"{dump_file}: run0 validation requires shapex/shapey/shapez columns")
    mesogen_indices = lca.select_mesogen_indices(types, mesogen_type)
    if mesogen_indices.size < 2:
        raise RuntimeError(f"{dump_file}: need at least two mesogen atoms for run0 validation")
    return {
        "timestep": int(timestep),
        "box": box,
        "columns": columns,
        "atom_ids": ids[mesogen_indices],
        "positions": positions[mesogen_indices, :],
        "shapes": shapes[mesogen_indices, :],
        "quaternions": quats[mesogen_indices, :],
    }


def python_type1_pair_table(
    frame: Dict[str, object],
    params: lca.GayBerneParams,
    output_path: Path,
    records: Optional[Sequence[Dict[str, object]]] = None,
) -> Dict[str, object]:
    pair_records = list(python_type1_pair_records(frame, params) if records is None else records)
    timestep = int(frame["timestep"])
    total = sum(float(record["pair_energy"]) for record in pair_records)
    attractive_count = sum(1 for record in pair_records if float(record["pair_energy"]) < 0.0)
    min_energy = min((float(record["pair_energy"]) for record in pair_records), default=math.inf)
    max_abs_energy = max((abs(float(record["pair_energy"])) for record in pair_records), default=0.0)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("timestep\tatom_i\tatom_j\tdistance\tpair_energy\twell_depth\n")
        for record in pair_records:
            handle.write(
                f"{timestep}\t{int(record['atom_i'])}\t{int(record['atom_j'])}\t{float(record['distance']):.12g}\t"
                f"{float(record['pair_energy']):.12g}\t{float(record['well_depth']):.12g}\n"
            )
    return {
        "timestep": timestep,
        "mesogen_count": int(np.asarray(frame["positions"], dtype=float).shape[0]),
        "pair_count_within_cutoff": int(len(pair_records)),
        "attractive_pair_count": int(attractive_count),
        "python_pair_energy_sum": float(total),
        "min_pair_energy": float(min_energy) if math.isfinite(min_energy) else None,
        "max_abs_pair_energy": float(max_abs_energy),
        "python_pair_table": str(output_path),
    }


def python_type1_pair_records(frame: Dict[str, object], params: lca.GayBerneParams) -> List[Dict[str, object]]:
    timestep = int(frame["timestep"])
    box = frame["box"]
    atom_ids = np.asarray(frame["atom_ids"], dtype=int)
    positions = np.asarray(frame["positions"], dtype=float)
    shapes = np.asarray(frame["shapes"], dtype=float)
    quats = np.asarray(frame["quaternions"], dtype=float)
    records: List[Dict[str, object]] = []
    for i in range(positions.shape[0] - 1):
        for j in range(i + 1, positions.shape[0]):
            delta = lca.minimum_image_vector(positions[j, :] - positions[i, :], box)
            distance = float(np.linalg.norm(delta))
            metrics = lca.gayberne_pair_metrics(
                r12=delta,
                quat_i=quats[i, :],
                quat_j=quats[j, :],
                shape_i=shapes[i, :],
                shape_j=shapes[j, :],
                params=params,
            )
            if metrics is None:
                continue
            pair_energy, well_depth = metrics
            records.append(
                {
                    "timestep": timestep,
                    "index_i": int(i),
                    "index_j": int(j),
                    "atom_i": int(atom_ids[i]),
                    "atom_j": int(atom_ids[j]),
                    "r12": [float(value) for value in delta.tolist()],
                    "distance": float(distance),
                    "pair_energy": float(pair_energy),
                    "well_depth": float(well_depth),
                }
            )
    return records


def select_microstate_records(
    records: Sequence[Dict[str, object]],
    *,
    sample_count: int,
    sample_percent: float,
    seed: int,
) -> List[Dict[str, object]]:
    total = len(records)
    if total == 0:
        return []
    if int(sample_count) < 0:
        raise ValueError("--microstate-sample-count must be non-negative")
    if float(sample_percent) < 0.0 or float(sample_percent) > 100.0:
        raise ValueError("--microstate-sample-percent must be in [0, 100]")
    percent_count = int(math.ceil(total * float(sample_percent) / 100.0)) if float(sample_percent) > 0.0 else 0
    desired = max(int(sample_count), percent_count)
    if desired <= 0:
        return []
    method = "all" if desired >= total else "random_without_replacement"
    if desired >= total:
        selected_indices = list(range(total))
    else:
        rng = np.random.default_rng(int(seed))
        selected_indices = sorted(int(index) for index in rng.choice(total, size=desired, replace=False).tolist())
    selected_total = len(selected_indices)
    sampling_rate = float(selected_total / max(total, 1))
    selected_records: List[Dict[str, object]] = []
    for index in selected_indices:
        row = dict(records[index])
        row["microstate_candidate_index"] = int(index)
        row["microstate_candidate_pair_count"] = int(total)
        row["microstate_selected_pair_count"] = int(selected_total)
        row["microstate_sampling_rate"] = float(sampling_rate)
        row["microstate_sampling_probability"] = float(sampling_rate)
        row["microstate_sampling_method"] = method
        selected_records.append(row)
    return selected_records


def pair_modify_lines(lmp_path: Path) -> List[str]:
    lines: List[str] = []
    variables = lca.parse_lammps_variables(lmp_path)
    for raw_line in lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = lca.strip_lammps_comment(raw_line)
        if not line:
            continue
        parts = line.split()
        if parts and parts[0] == "pair_modify":
            lines.append(" ".join(lca.resolve_lammps_token(part, variables) for part in parts))
    return lines


def write_lammps_run0_input(
    *,
    frame: Dict[str, object],
    params: lca.GayBerneParams,
    gb_param_file: Path,
    output_path: Path,
) -> Dict[str, object]:
    box: lca.BoxSpec = frame["box"]  # type: ignore[assignment]
    positions = np.asarray(frame["positions"], dtype=float)
    shapes = np.asarray(frame["shapes"], dtype=float)
    quats = np.asarray(frame["quaternions"], dtype=float)
    atom_ids = np.asarray(frame["atom_ids"], dtype=int)
    boundaries = [lammps_boundary_token(token) for token in box.boundary_tokens]
    original_upsilon = float(params.upsilon) * 2.0
    pair_modify = pair_modify_lines(gb_param_file)

    lines = [
        "# Auto-generated type-1-only Gay-Berne run0 validation input.",
        "# This validates the Python type1-type1 GB energy reconstruction against LAMMPS total pair energy.",
        "clear",
        "units lj",
        "dimension 3",
        f"boundary {' '.join(boundaries)}",
        "atom_style ellipsoid",
        f"region simbox block {box.bounds[0][0]:.16g} {box.bounds[0][1]:.16g} {box.bounds[1][0]:.16g} {box.bounds[1][1]:.16g} {box.bounds[2][0]:.16g} {box.bounds[2][1]:.16g} units box",
        "create_box 1 simbox",
    ]
    for idx, position in enumerate(positions, start=1):
        wrapped = wrapped_position_for_lammps(position, box)
        lines.append(
            f"create_atoms 1 single {wrapped[0]:.16g} {wrapped[1]:.16g} {wrapped[2]:.16g} units box"
        )
        shape = shapes[idx - 1, :]
        quat = quats[idx - 1, :]
        qx, qy, qz, qangle = quaternion_to_set_quat_args(quat)
        lines.append(f"set atom {idx} mass 1.0")
        lines.append(f"set atom {idx} shape {shape[0]:.16g} {shape[1]:.16g} {shape[2]:.16g}")
        lines.append(f"set atom {idx} quat {qx:.16g} {qy:.16g} {qz:.16g} {qangle:.16g}")
    lines.extend(
        [
            "neighbor 0.3 bin",
            "neigh_modify delay 0 every 1 check yes",
            f"pair_style gayberne {params.gamma:.16g} {original_upsilon:.16g} {params.mu:.16g} {params.cutoff:.16g}",
        ]
    )
    lines.extend(pair_modify)
    lines.append(
        "pair_coeff 1 1 "
        f"{params.epsilon:.16g} {params.sigma:.16g} "
        f"{params.eps_i[0]:.16g} {params.eps_i[1]:.16g} {params.eps_i[2]:.16g} "
        f"{params.eps_j[0]:.16g} {params.eps_j[1]:.16g} {params.eps_j[2]:.16g} {params.cutoff:.16g}"
    )
    lines.extend(
        [
            "compute peatom all pe/atom pair",
            "compute pesum all reduce sum c_peatom",
            "thermo 1",
            "thermo_style custom step pe c_pesum",
            "thermo_modify norm no",
            "run 0",
            "variable pair_energy equal pe",
            f'print "{MARKER} ${{pair_energy}}" file lammps_pair_energy.txt screen yes',
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")
    mapping_path = output_path.with_name("validation_atom_id_map.tsv")
    with mapping_path.open("w", encoding="utf-8") as handle:
        handle.write("validation_atom_id\toriginal_atom_id\n")
        for validation_id, original_id in enumerate(atom_ids.tolist(), start=1):
            handle.write(f"{validation_id}\t{int(original_id)}\n")
    return {
        "input": str(output_path),
        "atom_id_map": str(mapping_path),
        "mesogen_count": int(positions.shape[0]),
        "pair_modify_lines": pair_modify,
    }


def run_lammps(
    *,
    lammps_executable: str,
    mpi_prefix: str,
    input_file: Path,
    work_dir: Path,
) -> Dict[str, object]:
    command = []
    if mpi_prefix.strip():
        command.extend(shlex.split(mpi_prefix))
    command.extend([lammps_executable, "-in", str(input_file)])
    completed = subprocess.run(command, cwd=work_dir, text=True, capture_output=True, check=False)
    (work_dir / "lammps_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (work_dir / "lammps_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return {
        "command": command,
        "returncode": int(completed.returncode),
        "stdout": str(work_dir / "lammps_stdout.txt"),
        "stderr": str(work_dir / "lammps_stderr.txt"),
    }


def parse_lammps_energy(work_dir: Path) -> float:
    candidates = [work_dir / "lammps_pair_energy.txt", work_dir / "lammps_stdout.txt", work_dir / "log.lammps"]
    pattern = re.compile(rf"{re.escape(MARKER)}\s+([-+0-9.eE]+)")
    for path in candidates:
        if not path.exists():
            continue
        match = pattern.search(path.read_text(encoding="utf-8", errors="ignore"))
        if match:
            return float(match.group(1))
    raise RuntimeError(f"LAMMPS output did not contain {MARKER}")


def nontrivial_validation_failure(args: argparse.Namespace, python_summary: Dict[str, object]) -> Optional[str]:
    min_pairs = int(getattr(args, "min_pairs", 1))
    min_attractive_pairs = int(getattr(args, "min_attractive_pairs", 1))
    min_abs_python_total = float(getattr(args, "min_abs_python_total", 1e-12))
    pair_count = int(python_summary.get("pair_count_within_cutoff", 0))
    attractive_count = int(python_summary.get("attractive_pair_count", 0))
    python_total = abs(float(python_summary.get("python_pair_energy_sum", 0.0)))
    if pair_count < min_pairs:
        return f"Python pair table contains only {pair_count} pairs within cutoff; minimum required is {min_pairs}."
    if attractive_count < min_attractive_pairs:
        return f"Python pair table contains only {attractive_count} attractive pairs; minimum required is {min_attractive_pairs}."
    if not math.isfinite(python_total) or python_total < min_abs_python_total:
        return (
            f"Python total pair energy magnitude is {python_total:.6g}; "
            f"minimum required is {min_abs_python_total:.6g}."
        )
    return None


def write_failed_summary(
    *,
    output_root: Path,
    artifact_path: Path,
    fingerprint: Dict[str, object],
    gb_param_file: Path,
    message: str,
    validation_payload: Dict[str, object],
    comparison: Dict[str, object],
) -> None:
    validation_payload["comparison_summary"] = comparison
    (output_root / "run0_validation_summary.json").write_text(
        json.dumps(comparison, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    pipeline.write_validation_artifact(
        artifact_path,
        status="failed",
        fingerprint=fingerprint,
        gb_param_file=gb_param_file,
        message=message,
        validation=validation_payload,
    )


def two_particle_frame(frame: Dict[str, object], record: Dict[str, object]) -> Dict[str, object]:
    i = int(record["index_i"])
    j = int(record["index_j"])
    atom_ids = np.asarray(frame["atom_ids"], dtype=int)
    shapes = np.asarray(frame["shapes"], dtype=float)
    quats = np.asarray(frame["quaternions"], dtype=float)
    r12 = np.asarray(record.get("r12"), dtype=float)
    if r12.shape != (3,):
        positions = np.asarray(frame["positions"], dtype=float)
        r12 = positions[j, :] - positions[i, :]
    margin = max(10.0, float(np.linalg.norm(r12)) + 5.0)
    origin = np.array([margin, margin, margin], dtype=float)
    positions = np.array([origin, origin + r12], dtype=float)
    lo = np.minimum(positions[0, :], positions[1, :]) - margin
    hi = np.maximum(positions[0, :], positions[1, :]) + margin
    box = lca.BoxSpec(
        lengths=(float(hi[0] - lo[0]), float(hi[1] - lo[1]), float(hi[2] - lo[2])),
        bounds=((float(lo[0]), float(hi[0])), (float(lo[1]), float(hi[1])), (float(lo[2]), float(hi[2]))),
        boundary_tokens=("ff", "ff", "ff"),
    )
    return {
        "timestep": int(frame["timestep"]),
        "box": box,
        "atom_ids": np.array([atom_ids[i], atom_ids[j]], dtype=int),
        "positions": positions,
        "shapes": np.array([shapes[i, :], shapes[j, :]], dtype=float),
        "quaternions": np.array([quats[i, :], quats[j, :]], dtype=float),
    }


def run_microstate_checks(
    *,
    frame: Dict[str, object],
    records: Sequence[Dict[str, object]],
    params: lca.GayBerneParams,
    gb_param_file: Path,
    output_root: Path,
    args: argparse.Namespace,
) -> Dict[str, object]:
    sample_count = int(getattr(args, "microstate_sample_count", 12))
    sample_percent = float(getattr(args, "microstate_sample_percent", 0.0))
    seed = int(getattr(args, "microstate_seed", 20260429))
    tolerance_arg = getattr(args, "microstate_max_abs_delta", None)
    tolerance = float(args.max_abs_delta if tolerance_arg is None else tolerance_arg)
    selected = select_microstate_records(
        records,
        sample_count=sample_count,
        sample_percent=sample_percent,
        seed=seed,
    )
    actual_sampling_rate = float(len(selected) / max(len(records), 1)) if records else 0.0
    micro_root = output_root / "microstate_checks"
    micro_root.mkdir(parents=True, exist_ok=True)
    table_path = micro_root / "microstate_pair_checks.tsv"
    rows: List[Dict[str, object]] = []
    abs_deltas: List[float] = []
    failed_count = 0
    with table_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "sample_index\tmicrostate_candidate_index\tatom_i\tatom_j\tdistance\tpython_pair_energy\tlammps_pair_energy\t"
            "delta_python_minus_lammps\tabs_delta\tpassed\tcandidate_pair_count\tselected_pair_count\t"
            "microstate_sampling_rate\tmicrostate_sampling_probability\tmicrostate_sampling_method\twork_dir\n"
        )
        for sample_index, record in enumerate(selected, start=1):
            atom_i = int(record["atom_i"])
            atom_j = int(record["atom_j"])
            work_dir = micro_root / f"microstate_{sample_index:04d}_atoms_{atom_i}_{atom_j}"
            work_dir.mkdir(parents=True, exist_ok=True)
            input_path = work_dir / "validate_pair_run0.in"
            write_lammps_run0_input(
                frame=two_particle_frame(frame, record),
                params=params,
                gb_param_file=gb_param_file,
                output_path=input_path,
            )
            python_energy = float(record["pair_energy"])
            lammps_energy: Optional[float] = None
            failure_reason = ""
            try:
                run_summary = run_lammps(
                    lammps_executable=str(args.lammps_executable),
                    mpi_prefix=str(args.mpi_prefix or ""),
                    input_file=input_path,
                    work_dir=work_dir,
                )
                if run_summary["returncode"] != 0:
                    failure_reason = f"LAMMPS return code {run_summary['returncode']}"
                else:
                    lammps_energy = parse_lammps_energy(work_dir)
            except Exception as exc:
                failure_reason = str(exc)
            if lammps_energy is None:
                delta = float("nan")
                abs_delta = float("nan")
                passed = False
            else:
                delta = python_energy - float(lammps_energy)
                abs_delta = abs(float(delta))
                passed = math.isfinite(abs_delta) and abs_delta <= tolerance
            if not passed:
                failed_count += 1
            if math.isfinite(abs_delta):
                abs_deltas.append(abs_delta)
            row = {
                "sample_index": int(sample_index),
                "microstate_candidate_index": int(record.get("microstate_candidate_index", sample_index - 1)),
                "atom_i": atom_i,
                "atom_j": atom_j,
                "distance": float(record["distance"]),
                "python_pair_energy": python_energy,
                "lammps_pair_energy": lammps_energy,
                "delta_python_minus_lammps": delta,
                "abs_delta": abs_delta,
                "passed": passed,
                "work_dir": str(work_dir),
                "failure_reason": failure_reason,
                "candidate_pair_count": int(len(records)),
                "selected_pair_count": int(len(selected)),
                "microstate_sampling_rate": float(actual_sampling_rate),
                "microstate_sampling_probability": float(record.get("microstate_sampling_probability", actual_sampling_rate)),
                "microstate_sampling_method": str(record.get("microstate_sampling_method", "unknown")),
            }
            rows.append(row)
            handle.write(
                f"{sample_index}\t{int(row['microstate_candidate_index'])}\t{atom_i}\t{atom_j}\t{float(record['distance']):.12g}\t"
                f"{python_energy:.12g}\t{'' if lammps_energy is None else f'{float(lammps_energy):.12g}'}\t"
                f"{delta:.12g}\t{abs_delta:.12g}\t{int(passed)}\t{int(row['candidate_pair_count'])}\t"
                f"{int(row['selected_pair_count'])}\t{float(row['microstate_sampling_rate']):.12g}\t"
                f"{float(row['microstate_sampling_probability']):.12g}\t{row['microstate_sampling_method']}\t{work_dir}\n"
            )
    passed_count = len(selected) - failed_count
    max_abs_delta = max(abs_deltas, default=float("nan"))
    mean_abs_delta = float(sum(abs_deltas) / len(abs_deltas)) if abs_deltas else float("nan")
    rmse = math.sqrt(sum(value * value for value in abs_deltas) / len(abs_deltas)) if abs_deltas else float("nan")
    finite_rows = [row for row in rows if math.isfinite(float(row["abs_delta"]))]
    worst_pair = max(finite_rows, key=lambda row: float(row["abs_delta"])) if finite_rows else None
    summary = {
        "method": "sampled_two_particle_lammps_run0_microstates",
        "candidate_pair_count": int(len(records)),
        "selected_pair_count": int(len(selected)),
        "actual_sampling_rate": float(actual_sampling_rate),
        "passed_pair_count": int(passed_count),
        "sample_count": int(sample_count),
        "sample_percent": float(sample_percent),
        "seed": int(seed),
        "max_abs_tolerance": float(tolerance),
        "max_abs_delta": float(max_abs_delta),
        "mean_abs_delta": float(mean_abs_delta),
        "rmse": float(rmse),
        "failed_pair_count": int(failed_count),
        "passed": int(len(selected)) > 0 and failed_count == 0,
        "table": str(table_path),
        "results_table": str(table_path),
        "output_root": str(micro_root),
        "selection": {
            "candidate_pair_count": int(len(records)),
            "selected_pair_count": int(len(selected)),
            "actual_sampling_rate": float(actual_sampling_rate),
            "sample_count": int(sample_count),
            "sample_percent": float(sample_percent),
            "seed": int(seed),
            "selected_pairs": [
                {
                    "microstate_candidate_index": int(record.get("microstate_candidate_index", idx)),
                    "atom_i": int(record.get("atom_i", -1)),
                    "atom_j": int(record.get("atom_j", -1)),
                    "microstate_sampling_rate": float(record.get("microstate_sampling_rate", actual_sampling_rate)),
                    "microstate_sampling_method": str(record.get("microstate_sampling_method", "unknown")),
                }
                for idx, record in enumerate(selected)
            ],
        },
        "worst_pair": worst_pair,
        "rows": rows,
    }
    (micro_root / "microstate_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary


def validation_forwarded_args(
    *,
    dump_file: Path,
    local_pair_file: Optional[Path],
    exclude_pair_file: Optional[Path],
) -> List[str]:
    forwarded = [str(dump_file), "--contact-mode", "gayberne"]
    if local_pair_file is not None:
        forwarded.extend(["--local-pair-file", str(local_pair_file)])
    if exclude_pair_file is not None:
        forwarded.extend(["--exclude-pair-file", str(exclude_pair_file)])
    return forwarded


def validate_run0(args: argparse.Namespace) -> Dict[str, object]:
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    params = lca.parse_gayberne_params_from_lmp(args.gb_param_file)
    frame = extract_type1_frame(args.dump_file, frame_index=args.frame_index, mesogen_type=args.mesogen_type)
    pair_records = python_type1_pair_records(frame, params)
    python_table = output_root / "python_type1_pair_energies.tsv"
    python_summary = python_type1_pair_table(frame, params, python_table, records=pair_records)
    lammps_input = output_root / "validate_type1_gb_run0.in"
    input_summary = write_lammps_run0_input(
        frame=frame,
        params=params,
        gb_param_file=args.gb_param_file,
        output_path=lammps_input,
    )
    forwarded = validation_forwarded_args(
        dump_file=args.dump_file,
        local_pair_file=args.local_pair_file,
        exclude_pair_file=args.exclude_pair_file,
    )
    fingerprint = pipeline.build_potential_fingerprint(
        gb_param_file=args.gb_param_file,
        forwarded=forwarded,
        lammps_executable=args.lammps_executable,
        mesogen_type=int(args.mesogen_type),
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "type1_only_total_pair_energy",
        "dump_file": str(args.dump_file),
        "frame_index": int(args.frame_index),
        "timestep": int(frame["timestep"]),
        "gb_param_file": str(args.gb_param_file),
        "lammps_input": input_summary,
        "python": python_summary,
        "fingerprint_sha256": fingerprint["sha256"],
        "limitations": [
            "This validates type1-type1 Gay-Berne total pair energy, not sphere/anchor interactions.",
            "It does not use pair/local eng because common pair_gayberne builds do not expose pair-local energies.",
            "Topology exclusion files are included in the fingerprint, but the type1-only validation contains no bonds.",
        ],
    }
    (output_root / "run0_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    artifact_path = args.verified_potential_file or (output_root / "verified_potential.json")
    validation_payload: Dict[str, object] = {
        "method": "standalone_type1_only_lammps_run0_total_pair_energy",
        "manifest": str(output_root / "run0_validation_manifest.json"),
        "python_summary": python_summary,
        "lammps_input": input_summary,
    }
    if args.dry_run:
        pipeline.write_validation_artifact(
            artifact_path,
            status="validation_required",
            fingerprint=fingerprint,
            gb_param_file=args.gb_param_file,
            message="Dry run only: LAMMPS input and Python energy table were generated, but LAMMPS was not executed.",
            validation=validation_payload,
        )
        return {"status": "validation_required", "artifact": str(artifact_path), "manifest": manifest}

    preflight_failure = nontrivial_validation_failure(args, python_summary)
    if preflight_failure is not None:
        comparison = {
            "python_pair_energy_sum": float(python_summary["python_pair_energy_sum"]),
            "pair_count_within_cutoff": int(python_summary["pair_count_within_cutoff"]),
            "attractive_pair_count": int(python_summary["attractive_pair_count"]),
            "passed": False,
            "failure_reason": preflight_failure,
        }
        write_failed_summary(
            output_root=output_root,
            artifact_path=artifact_path,
            fingerprint=fingerprint,
            gb_param_file=args.gb_param_file,
            message="LAMMPS run0 validation was not attempted because the selected frame is not an informative GB validation case.",
            validation_payload=validation_payload,
            comparison=comparison,
        )
        raise SystemExit(f"Potential validation failed preflight. {preflight_failure}")

    try:
        run_summary = run_lammps(
            lammps_executable=str(args.lammps_executable),
            mpi_prefix=str(args.mpi_prefix or ""),
            input_file=lammps_input,
            work_dir=output_root,
        )
    except Exception as exc:
        comparison = {"passed": False, "failure_reason": f"LAMMPS launch failed: {exc}"}
        write_failed_summary(
            output_root=output_root,
            artifact_path=artifact_path,
            fingerprint=fingerprint,
            gb_param_file=args.gb_param_file,
            message="LAMMPS run0 validation failed before a comparison could be made.",
            validation_payload=validation_payload,
            comparison=comparison,
        )
        raise SystemExit(f"LAMMPS run failed before launch completed. See {artifact_path}") from exc
    validation_payload["lammps_run"] = run_summary
    if run_summary["returncode"] != 0:
        write_failed_summary(
            output_root=output_root,
            artifact_path=artifact_path,
            fingerprint=fingerprint,
            gb_param_file=args.gb_param_file,
            message=f"LAMMPS run0 validation failed with return code {run_summary['returncode']}.",
            validation_payload=validation_payload,
            comparison={"passed": False, "failure_reason": f"LAMMPS return code {run_summary['returncode']}"},
        )
        raise SystemExit(f"LAMMPS run failed. See {output_root / 'lammps_stderr.txt'}")

    try:
        lammps_energy = parse_lammps_energy(output_root)
    except Exception as exc:
        write_failed_summary(
            output_root=output_root,
            artifact_path=artifact_path,
            fingerprint=fingerprint,
            gb_param_file=args.gb_param_file,
            message="LAMMPS run0 validation output did not contain a usable pair-energy marker.",
            validation_payload=validation_payload,
            comparison={"passed": False, "failure_reason": str(exc)},
        )
        raise SystemExit(f"LAMMPS run did not produce {MARKER}. See {artifact_path}") from exc
    python_energy = float(python_summary["python_pair_energy_sum"])
    delta = python_energy - lammps_energy
    comparison = {
        "python_pair_energy_sum": python_energy,
        "lammps_pair_energy_sum": float(lammps_energy),
        "delta_python_minus_lammps": float(delta),
        "abs_delta": abs(float(delta)),
        "max_abs_tolerance": float(args.max_abs_delta),
        "passed": abs(float(delta)) <= float(args.max_abs_delta),
    }
    validation_payload["comparison_summary"] = comparison
    (output_root / "run0_validation_summary.json").write_text(
        json.dumps(comparison, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    microstate_summary: Optional[Dict[str, object]] = None
    if comparison["passed"]:
        microstate_summary = run_microstate_checks(
            frame=frame,
            records=pair_records,
            params=params,
            gb_param_file=args.gb_param_file,
            output_root=output_root,
            args=args,
        )
        validation_payload["microstate_summary"] = microstate_summary
    status = "validated" if comparison["passed"] and microstate_summary is not None and bool(microstate_summary.get("passed")) else "failed"
    pipeline.write_validation_artifact(
        artifact_path,
        status=status,
        fingerprint=fingerprint,
        gb_param_file=args.gb_param_file,
        message=(
            "LAMMPS run0 total pair energy and sampled two-particle microstates matched Python type1-type1 GB reconstruction."
            if status == "validated"
            else "LAMMPS run0 validation failed for total energy or sampled two-particle microstates."
        ),
        validation=validation_payload,
    )
    if status != "validated":
        raise SystemExit(f"Potential validation failed. See {artifact_path}")
    return {"status": status, "artifact": str(artifact_path), "comparison": comparison, "microstate_summary": microstate_summary}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Standalone LAMMPS run0 validator for Python type1-type1 Gay-Berne energy reconstruction.")
    parser.add_argument("--dump-file", type=Path, required=True)
    parser.add_argument("--gb-param-file", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("lc_run0_validation"))
    parser.add_argument("--verified-potential-file", type=Path, default=None)
    parser.add_argument("--frame-index", type=int, default=0)
    parser.add_argument("--mesogen-type", type=int, default=1)
    parser.add_argument("--local-pair-file", type=Path, default=None)
    parser.add_argument("--exclude-pair-file", type=Path, default=None)
    parser.add_argument("--lammps-executable", required=True)
    parser.add_argument("--mpi-prefix", default="", help="Optional prefix such as 'mpiexec -np 4'.")
    parser.add_argument("--max-abs-delta", type=float, default=1e-6)
    parser.add_argument("--min-pairs", type=int, default=1, help="Minimum Python type1-type1 pairs within cutoff required before validation can pass.")
    parser.add_argument("--min-attractive-pairs", type=int, default=1, help="Minimum attractive Python type1-type1 pairs required before validation can pass.")
    parser.add_argument("--min-abs-python-total", type=float, default=1e-12, help="Minimum absolute Python total pair energy required before validation can pass.")
    parser.add_argument("--microstate-sample-count", type=int, default=12, help="Number of two-particle pair microstates to sample. Default: 12.")
    parser.add_argument("--microstate-sample-percent", type=float, default=0.0, help="Percent of candidate pairs to sample. Use 100 for full pair sampling.")
    parser.add_argument("--microstate-seed", type=int, default=20260429, help="Deterministic RNG seed for microstate pair sampling.")
    parser.add_argument("--microstate-max-abs-delta", type=float, default=None, help="Optional stricter absolute tolerance for sampled microstates. Defaults to --max-abs-delta.")
    parser.add_argument("--dry-run", action="store_true", help="Generate validation input and Python table without running LAMMPS.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.frame_index < 0:
        raise SystemExit("--frame-index must be non-negative")
    if args.microstate_sample_count < 0:
        raise SystemExit("--microstate-sample-count must be non-negative")
    if args.microstate_sample_percent < 0.0 or args.microstate_sample_percent > 100.0:
        raise SystemExit("--microstate-sample-percent must be in [0, 100]")
    result = validate_run0(args)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
