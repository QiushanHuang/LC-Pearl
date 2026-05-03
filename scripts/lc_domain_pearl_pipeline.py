#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import liquid_crystal_aggregation as lca
import lammps_pairlocal_compare as pair_compare

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_SCRIPT = SCRIPT_DIR / "liquid_crystal_aggregation.py"
RUN0_VALIDATOR_SCRIPT = SCRIPT_DIR / "lc_lammps_run0_validate.py"
THRESHOLD_PRIOR_SCRIPT = SCRIPT_DIR / "lc_threshold_prior.py"
PIPELINE_SCHEMA_VERSION = 1
THRESHOLD_PRIOR_SCHEMA_VERSION = 7
THRESHOLD_PRIOR_METHOD_NAME = "LC-Pearl 2.1.0 core-tier streaming threshold prior"
VALIDATED_STATUSES = {"validated"}
STANDALONE_RUN0_METHODS = {
    "standalone_type1_only_lammps_run0_total_pair_energy",
}
EXTERNAL_PAIR_TABLE_METHODS = {
    "external_pair_table_lammps_pair_energy_comparison",
}
DEFAULT_LAMMPS_EXECUTABLE_CANDIDATES = [
    Path("/Users/joshua/Desktop/Qiushan_Code/lammmps_git/lammps/build-asphere-mpi/lmp"),
]


def find_option_value(args: List[str], option: str) -> Optional[str]:
    for idx, item in enumerate(args):
        if item == option and idx + 1 < len(args):
            return args[idx + 1]
        if item.startswith(option + "="):
            return item.split("=", 1)[1]
    return None


def forwarded_positionals(args: List[str]) -> List[str]:
    values: List[str] = []
    skip_next = False
    options_with_values = {
        "--pattern",
        "--output-root",
        "--axis",
        "--r-cut",
        "--p2-cut",
        "--min-core-neighbors",
        "--cutoff-bins",
        "--cutoff-frames",
        "--max-auto-r-cut",
        "--auto-r-cut-shape-factor",
        "--mesogen-type",
        "--anchor-types",
        "--contact-mode",
        "--s-excl",
        "--exclude-pair-file",
        "--local-pair-file",
        "--gb-param-file",
        "--gb-threshold-mode",
        "--u-on",
        "--u-off",
        "--gb-off-strength",
        "--gb-on-strength",
        "--gb-core-strength",
        "--gb-strict-core-strength",
        "--p2-core-cut",
        "--p2-strict-core-cut",
        "--cluster-cut",
        "--cluster-cut-shape-factor",
        "--cluster-min-size",
        "--g-on",
        "--g-off",
        "--q-on",
        "--q-off",
        "--robust-min-s2",
        "--robust-min-size",
        "--robust-min-evidence",
        "--domain-min-lifetime",
        "--n-min",
        "--adjacent-id-gap",
        "--perturbation-r-cut-scale",
        "--perturbation-p2-margin",
        "--stable-overlap-fraction",
        "--pearl-gap-cut",
        "--pearl-min-cross-contacts",
        "--pearl-min-boundary-particles",
        "--pearl-max-aspect-ratio",
        "--track-jaccard",
        "--consensus-threshold",
        "--cluster-envelope-padding",
        "--every",
        "--workers",
        "--r-energy-cap",
        "--edge-diagnostics-table",
        "--edge-diagnostics-sample-size",
    }
    for item in args:
        if skip_next:
            skip_next = False
            continue
        if item.startswith("--"):
            option = item.split("=", 1)[0]
            if "=" not in item and option in options_with_values:
                skip_next = True
            continue
        values.append(item)
    return values


def discover_lammps_executable(explicit: Optional[object] = None) -> Optional[str]:
    candidates: List[object] = []
    if explicit:
        candidates.append(explicit)
    for env_name in ("LC_LAMMPS_EXECUTABLE", "LAMMPS_EXECUTABLE"):
        value = os.environ.get(env_name)
        if value:
            candidates.append(value)
    candidates.extend(DEFAULT_LAMMPS_EXECUTABLE_CANDIDATES)
    for executable_name in ("lmp", "lmp_mpi", "lammps"):
        found = shutil.which(executable_name)
        if found:
            candidates.append(found)
    for candidate in candidates:
        path = Path(str(candidate)).expanduser()
        if path.exists():
            return str(path.resolve())
    if explicit:
        return str(explicit)
    return None


def replace_or_append_option(args: List[str], option: str, value: object) -> List[str]:
    updated: List[str] = []
    skip_next = False
    replaced = False
    for item in args:
        if skip_next:
            skip_next = False
            continue
        if item == option:
            updated.extend([option, str(value)])
            skip_next = True
            replaced = True
            continue
        if item.startswith(option + "="):
            updated.append(f"{option}={value}")
            replaced = True
            continue
        updated.append(item)
    if not replaced:
        updated.extend([option, str(value)])
    return updated


def has_option(args: Sequence[str], option: str) -> bool:
    return option in args or any(item.startswith(option + "=") for item in args)


def option_value(args: Sequence[str], option: str, default: object) -> object:
    value = find_option_value(list(args), option)
    return default if value is None else value


def bool_flag_enabled(args: Sequence[str], option: str) -> bool:
    return option in args


def run_command(command: List[str]) -> None:
    print("[RUN] " + " ".join(command))
    subprocess.run(command, check=True)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Optional[Path]) -> Optional[str]:
    if path is None:
        return None
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stripped_lammps_lines(path: Path) -> List[str]:
    lines: List[str] = []
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = lca.strip_lammps_comment(raw_line)
        if line:
            lines.append(line)
    return lines


def canonical_potential_lines(lmp_path: Path) -> List[str]:
    variables = lca.parse_lammps_variables(lmp_path)
    interesting = {"units", "atom_style", "pair_style", "pair_modify", "pair_coeff", "special_bonds"}
    canonical: List[str] = []
    for line in stripped_lammps_lines(lmp_path):
        parts = line.split()
        if not parts or parts[0] not in interesting:
            continue
        resolved = " ".join(lca.resolve_lammps_token(part, variables) for part in parts)
        canonical.append(resolved)
    return canonical


def find_input_files_for_schema(forwarded: Sequence[str]) -> List[Path]:
    positionals = [Path(item) for item in forwarded_positionals(list(forwarded))]
    if not positionals:
        positionals = [Path.cwd()]
    pattern = find_option_value(list(forwarded), "--pattern") or "*.dump"
    recursive = has_option(forwarded, "--recursive")
    try:
        return lca.iter_input_files(positionals, pattern=pattern, recursive=recursive)
    except Exception:
        return []


def dump_schema_for_fingerprint(forwarded: Sequence[str]) -> Dict[str, object]:
    files = find_input_files_for_schema(forwarded)
    if not files:
        return {"available": False, "reason": "no_input_dump_found"}
    first = files[0]
    try:
        _timestep, box, columns, _col_index, _data = next(lca.parse_dump_frames(first))
    except Exception as exc:
        return {"available": False, "file": str(first), "reason": str(exc)}
    required_columns = [
        "id",
        "type",
        "xu",
        "yu",
        "zu",
        "x",
        "y",
        "z",
        "quatw",
        "quati",
        "quatj",
        "quatk",
        "shapex",
        "shapey",
        "shapez",
    ]
    present = [column for column in required_columns if column in columns]
    return {
        "available": True,
        "columns": columns,
        "required_columns_present": present,
        "boundary_tokens": list(box.boundary_tokens),
    }


def lammps_executable_version(executable: Optional[str]) -> Dict[str, object]:
    if not executable:
        return {"available": False, "reason": "not_requested"}
    try:
        completed = subprocess.run(
            [executable, "-h"],
            check=False,
            text=True,
            capture_output=True,
            timeout=10,
        )
    except Exception as exc:
        return {"available": False, "executable": executable, "reason": str(exc)}
    first_lines = "\n".join((completed.stdout or completed.stderr).splitlines()[:6])
    return {
        "available": completed.returncode == 0,
        "executable": executable,
        "returncode": completed.returncode,
        "summary": first_lines,
        "sha256": sha256_text(first_lines),
    }


def build_potential_fingerprint(
    *,
    gb_param_file: Path,
    forwarded: Sequence[str],
    lammps_executable: Optional[str] = None,
    mesogen_type: int = 1,
) -> Dict[str, object]:
    params = lca.parse_gayberne_params_from_lmp(gb_param_file)
    local_pair_file = find_option_value(list(forwarded), "--local-pair-file")
    exclude_pair_file = find_option_value(list(forwarded), "--exclude-pair-file")
    acceptance_components: Dict[str, object] = {
        "schema_version": PIPELINE_SCHEMA_VERSION,
        "validation_contract": {
            "scope": "type1_type1_gayberne_energy",
            "accepted_status": "validated",
            "requires_comparison_proof": True,
            "validator": "lc_lammps_run0_validate.py",
        },
        "mesogen_type": int(mesogen_type),
        "run0_validator_code_sha256": sha256_file(RUN0_VALIDATOR_SCRIPT),
        "canonical_potential_lines": canonical_potential_lines(gb_param_file),
        "gayberne_type_1_1": {
            "gamma": params.gamma,
            "upsilon_effective": params.upsilon,
            "mu": params.mu,
            "cutoff": params.cutoff,
            "epsilon": params.epsilon,
            "sigma": params.sigma,
            "eps_i": list(params.eps_i),
            "eps_j": list(params.eps_j),
        },
        "analysis_code_sha256": sha256_file(ANALYSIS_SCRIPT),
        "dump_schema": dump_schema_for_fingerprint(forwarded),
        "topology_pair_files": {
            "local_pair_file_sha256": sha256_file(Path(local_pair_file)) if local_pair_file else None,
            "exclude_pair_file_sha256": sha256_file(Path(exclude_pair_file)) if exclude_pair_file else None,
        },
    }
    components = dict(acceptance_components)
    components["lammps_executable_provenance"] = lammps_executable_version(lammps_executable)
    encoded = json.dumps(acceptance_components, sort_keys=True, ensure_ascii=False)
    provenance_encoded = json.dumps(components, sort_keys=True, ensure_ascii=False)
    return {
        "sha256": sha256_text(encoded),
        "components": components,
        "provenance_sha256": sha256_text(provenance_encoded),
    }


def load_validation_artifact(path: Path) -> Optional[Dict[str, object]]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def validation_artifact_matches(artifact: Dict[str, object], fingerprint: Dict[str, object]) -> bool:
    validation = artifact.get("validation", {})
    if not isinstance(validation, dict):
        return False
    method = str(validation.get("method", ""))
    comparison = validation.get("comparison_summary", {})
    if not isinstance(comparison, dict) or comparison.get("passed") is not True:
        return False

    try:
        if method in EXTERNAL_PAIR_TABLE_METHODS:
            if "matched_pairs" not in comparison:
                return False
            matched = int(comparison.get("matched_pairs", 0))
            missing_in_lammps = int(comparison.get("missing_in_lammps", 0))
            missing_in_python = int(comparison.get("missing_in_python", 0))
            max_abs = float(comparison.get("max_abs_delta", float("nan")))
            tolerance = float(comparison.get("max_abs_tolerance", float("nan")))
            proof_ok = (
                matched > 0
                and missing_in_lammps == 0
                and missing_in_python == 0
                and math.isfinite(max_abs)
                and math.isfinite(tolerance)
                and max_abs <= tolerance
            )
        elif method in STANDALONE_RUN0_METHODS:
            python_summary = validation.get("python_summary", {})
            if not isinstance(python_summary, dict):
                return False
            pair_count = int(python_summary.get("pair_count_within_cutoff", 0))
            attractive_count = int(python_summary.get("attractive_pair_count", 0))
            python_total = abs(float(python_summary.get("python_pair_energy_sum", 0.0)))
            abs_delta = float(comparison.get("abs_delta", float("nan")))
            tolerance = float(comparison.get("max_abs_tolerance", float("nan")))
            proof_ok = (
                pair_count > 0
                and attractive_count > 0
                and math.isfinite(python_total)
                and python_total > 0.0
                and math.isfinite(abs_delta)
                and math.isfinite(tolerance)
                and abs_delta <= tolerance
            )
            microstate = validation.get("microstate_summary", {})
            if not isinstance(microstate, dict):
                return False
            micro_selected = int(microstate.get("selected_pair_count", 0))
            micro_failed = int(microstate.get("failed_pair_count", 0))
            micro_max_abs = float(microstate.get("max_abs_delta", float("nan")))
            micro_tolerance = float(microstate.get("max_abs_tolerance", float("nan")))
            proof_ok = (
                proof_ok
                and microstate.get("passed") is True
                and micro_selected > 0
                and micro_failed == 0
                and math.isfinite(micro_max_abs)
                and math.isfinite(micro_tolerance)
                and micro_max_abs <= micro_tolerance
            )
        else:
            return False
    except (TypeError, ValueError):
        return False

    return (
        int(artifact.get("schema_version", 0)) == PIPELINE_SCHEMA_VERSION
        and artifact.get("status") in VALIDATED_STATUSES
        and artifact.get("fingerprint", {}).get("sha256") == fingerprint.get("sha256")
        and proof_ok
    )


def write_validation_artifact(
    path: Path,
    *,
    status: str,
    fingerprint: Dict[str, object],
    gb_param_file: Path,
    message: str,
    validation: Optional[Dict[str, object]] = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    params = lca.parse_gayberne_params_from_lmp(gb_param_file)
    payload = {
        "schema_version": PIPELINE_SCHEMA_VERSION,
        "status": status,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "fingerprint": fingerprint,
        "gayberne_type_1_1": {
            "gamma": params.gamma,
            "upsilon_effective": params.upsilon,
            "mu": params.mu,
            "cutoff": params.cutoff,
            "epsilon": params.epsilon,
            "sigma": params.sigma,
            "eps_i": list(params.eps_i),
            "eps_j": list(params.eps_j),
        },
        "message": message,
        "validation": validation or {},
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def compare_validation_tables(
    *,
    python_table: Path,
    lammps_table: Path,
    output_root: Path,
    max_abs_tolerance: float,
) -> Tuple[bool, Dict[str, object]]:
    summary = pair_compare.compare_pair_tables(python_table, lammps_table, output_root)
    matched = int(summary.get("matched_pairs", 0))
    max_abs = float(summary.get("max_abs_delta", float("nan")))
    missing_in_lammps = int(summary.get("missing_in_lammps", 0))
    missing_in_python = int(summary.get("missing_in_python", 0))
    passed = (
        matched > 0
        and math.isfinite(max_abs)
        and max_abs <= float(max_abs_tolerance)
        and missing_in_lammps == 0
        and missing_in_python == 0
    )
    summary["passed"] = passed
    summary["max_abs_tolerance"] = float(max_abs_tolerance)
    (output_root / "validation_decision.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return passed, summary


def attach_potential_validation_to_run_summary(output_root: Path, artifact_path: Optional[Path]) -> None:
    if artifact_path is None or not artifact_path.exists():
        return
    summary_path = output_root.resolve() / "run_summary.json"
    if not summary_path.exists():
        return
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return
    summary["potential_validation"] = {
        "artifact_path": str(artifact_path),
        "status": artifact.get("status", "unknown"),
        "fingerprint_sha256": artifact.get("fingerprint", {}).get("sha256"),
        "message": artifact.get("message", ""),
        "validated": artifact.get("status") in VALIDATED_STATUSES,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")


def copy_if_exists(source: Path, destination: Path) -> Optional[Path]:
    if not source.exists() or not source.is_file():
        return None
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def sync_threshold_prior_diagnostics(
    output_root: Path,
    prior_path: Path,
    artifact: Optional[Dict[str, object]],
) -> List[str]:
    """Copy the threshold-prior evidence used by this run into diagnostics."""
    if artifact is None:
        try:
            artifact = load_threshold_prior(prior_path)
        except (OSError, json.JSONDecodeError):
            artifact = {}
    diagnostics_dir = output_root.resolve() / "diagnostics"
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    destination_dir = diagnostics_dir / "threshold_prior"
    copied: List[str] = []

    def record(source: Path, relative_destination: str, *, top_level_alias: Optional[str] = None) -> None:
        destination = copy_if_exists(source, destination_dir / relative_destination)
        if destination is not None:
            copied.append(str(destination))
        if top_level_alias:
            alias = copy_if_exists(source, diagnostics_dir / top_level_alias)
            if alias is not None:
                copied.append(str(alias))

    prior_dir = prior_path.resolve().parent
    record(prior_path.resolve(), "global_thresholds.json")
    record(prior_dir / "threshold_recommendations.json", "threshold_recommendations.json")
    record(prior_dir / "streaming_manifest.json", "streaming_manifest.json")

    outputs = artifact.get("outputs", {}) if isinstance(artifact, dict) else {}
    plot_sources = {
        "gb_strength_hist.png": Path(str(outputs.get("gb_strength_hist_plot", prior_dir / "gb_strength_hist.png"))),
        "gb_strength_vs_p2_stream_hist.png": Path(str(outputs.get("stream_histogram_plot", prior_dir / "gb_strength_vs_p2_stream_hist.png"))),
        "gb_strength_vs_p2_stream_dotgrid.png": Path(str(outputs.get("stream_dotgrid_plot", prior_dir / "gb_strength_vs_p2_stream_dotgrid.png"))),
        "gb_strength_vs_p2_stream_hexbin.png": Path(str(outputs.get("stream_hexbin_plot", prior_dir / "gb_strength_vs_p2_stream_hexbin.png"))),
    }
    for name, source in plot_sources.items():
        record(source, name, top_level_alias=name)

    preview_dir = Path(str(outputs.get("lobe_split_preview_dir", prior_dir / "lobe_split_preview")))
    for name in (
        "gb_strength_vs_p2_stream_lobe_split_dotgrid.png",
        "gb_strength_vs_p2_stream_lobe_split_comparison.png",
        "gb_core_slice_hist.png",
    ):
        record(preview_dir / name, f"lobe_split_preview/{name}", top_level_alias=f"lobe_split_preview/{name}")
        if name == "gb_core_slice_hist.png":
            root_copy = copy_if_exists(preview_dir / name, diagnostics_dir / name)
            if root_copy is not None:
                copied.append(str(root_copy))

    return copied


def ensure_potential_validation(args: argparse.Namespace, forwarded: Sequence[str]) -> Optional[Path]:
    if args.potential_validation == "off":
        return None
    contact_mode = find_option_value(list(forwarded), "--contact-mode")
    if contact_mode is not None and contact_mode != "gayberne":
        return None
    if args.gb_param_file is None:
        return None

    cache_dir = (args.potential_cache_dir or (args.output_root.resolve() / "potential_validation")).resolve()
    explicit_artifact_requested = args.verified_potential_file is not None
    artifact_path = (args.verified_potential_file or (cache_dir / "verified_potential.json")).resolve()
    comparison_dir = cache_dir / "lammps_comparison"
    fingerprint = build_potential_fingerprint(
        gb_param_file=args.gb_param_file,
        forwarded=forwarded,
        lammps_executable=args.lammps_executable,
        mesogen_type=int(float(option_value(forwarded, "--mesogen-type", 1))),
    )

    existing = None if args.potential_validation == "refresh" else load_validation_artifact(artifact_path)
    if existing is not None and validation_artifact_matches(existing, fingerprint):
        print(f"[OK] verified potential cache hit: {artifact_path}")
        return artifact_path
    if existing is not None:
        print(f"[WARN] verified potential cache miss or invalid status: {artifact_path}")

    validation_payload: Dict[str, object] = {
        "method": "pending_standalone_lammps_run0_total_pair_energy",
        "notes": [
            "Python GB reconstruction is not treated as validated until lc_lammps_run0_validate.py writes a matching validated artifact.",
            "compute pair/local eng is not assumed to work for pair_gayberne.",
        ],
    }
    if args.python_pair_energy_table and args.lammps_pair_energy_table:
        passed, summary = compare_validation_tables(
            python_table=args.python_pair_energy_table,
            lammps_table=args.lammps_pair_energy_table,
            output_root=comparison_dir,
            max_abs_tolerance=args.validation_max_abs_delta,
        )
        validation_payload["comparison_summary"] = summary
        if passed:
            validation_payload["method"] = "external_pair_table_lammps_pair_energy_comparison"
            write_validation_artifact(
                artifact_path,
                status="validated",
                fingerprint=fingerprint,
                gb_param_file=args.gb_param_file,
                message="External LAMMPS/Python pair-energy table comparison passed for this potential fingerprint.",
                validation=validation_payload,
            )
            print(f"[OK] verified potential artifact written: {artifact_path}")
            return artifact_path
        write_validation_artifact(
            artifact_path,
            status="failed",
            fingerprint=fingerprint,
            gb_param_file=args.gb_param_file,
            message="LAMMPS energy comparison failed for this potential fingerprint.",
            validation=validation_payload,
        )
        raise SystemExit(f"Potential validation failed. See {artifact_path}")

    input_files = find_input_files_for_schema(forwarded)
    representative_dump = str(input_files[0]) if input_files else "representative.dump"
    suggested_lammps = discover_lammps_executable(args.lammps_executable)
    validation_command = [
        sys.executable,
        str(SCRIPT_DIR / "lc_lammps_run0_validate.py"),
        "--dump-file",
        representative_dump,
        "--gb-param-file",
        str(args.gb_param_file),
        "--output-root",
        str(cache_dir / "run0_validation"),
        "--verified-potential-file",
        str(artifact_path),
        "--lammps-executable",
        suggested_lammps or "REPLACE_WITH_LAMMPS_EXECUTABLE",
    ]
    if suggested_lammps is None:
        validation_payload["standalone_validator_command_requires_lammps_executable"] = True
    local_pair_file = find_option_value(list(forwarded), "--local-pair-file")
    exclude_pair_file = find_option_value(list(forwarded), "--exclude-pair-file")
    if local_pair_file:
        validation_command.extend(["--local-pair-file", local_pair_file])
    if exclude_pair_file:
        validation_command.extend(["--exclude-pair-file", exclude_pair_file])
    validation_payload["standalone_validator_command"] = validation_command
    pending_artifact_path = artifact_path
    if explicit_artifact_requested and existing is not None:
        pending_artifact_path = (cache_dir / "validation_required.json").resolve()
    write_validation_artifact(
        pending_artifact_path,
        status="validation_required",
        fingerprint=fingerprint,
        gb_param_file=args.gb_param_file,
        message=(
            "No matching validated potential artifact was found. Run lc_lammps_run0_validate.py on a representative "
            "dump frame with the same topology pair files to create a validated artifact."
        ),
        validation=validation_payload,
    )
    message = (
        f"Potential validation required. Suggested validator command: {' '.join(validation_command)}. "
        f"Artifact: {pending_artifact_path}"
    )
    if args.potential_validation == "require":
        raise SystemExit(message)
    print(f"[WARN] {message}")
    return pending_artifact_path


def potential_validation_status(artifact_path: Optional[Path]) -> str:
    if artifact_path is None:
        return "off"
    artifact = load_validation_artifact(artifact_path)
    if artifact is None:
        return "unknown"
    return str(artifact.get("status", "unknown"))


def current_thresholds_from_forwarded(forwarded: Sequence[str]) -> Dict[str, object]:
    return {
        "gb_off_strength": float(option_value(forwarded, "--gb-off-strength", 0.12)),
        "gb_on_strength": float(option_value(forwarded, "--gb-on-strength", 0.30)),
        "p2_cut": float(option_value(forwarded, "--p2-cut", 0.70)),
        "gb_core_strength": float(option_value(forwarded, "--gb-core-strength", 0.70)),
        "p2_core_cut": float(option_value(forwarded, "--p2-core-cut", 0.71)),
        "gb_strict_core_strength": float(option_value(forwarded, "--gb-strict-core-strength", 0.90)),
        "p2_strict_core_cut": float(option_value(forwarded, "--p2-strict-core-cut", 0.80)),
        "robust_min_s2": float(option_value(forwarded, "--robust-min-s2", 0.70)),
        "n_min": int(float(option_value(forwarded, "--n-min", 3))),
    }


def append_forwarded_pair_scan_options(command: List[str], forwarded: Sequence[str]) -> None:
    for option in (
        "--pattern",
        "--axis",
        "--mesogen-type",
        "--anchor-types",
        "--s-excl",
        "--local-pair-file",
        "--exclude-pair-file",
        "--r-energy-cap",
    ):
        value = find_option_value(list(forwarded), option)
        if value is not None:
            command.extend([option, value])
    if bool_flag_enabled(forwarded, "--recursive"):
        command.append("--recursive")


def load_threshold_prior(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def threshold_prior_matches(
    artifact: Dict[str, object],
    *,
    args: argparse.Namespace,
    forwarded: Sequence[str],
) -> Tuple[bool, str]:
    if int(artifact.get("schema_version", -1)) != THRESHOLD_PRIOR_SCHEMA_VERSION:
        return False, "threshold prior schema_version is not LC Domain-Pearl V2"
    if str(artifact.get("method_name", "")) != THRESHOLD_PRIOR_METHOD_NAME:
        return False, "threshold prior method_name is not LC Domain-Pearl V2 streaming prior"
    algorithm = artifact.get("algorithm_fingerprint", {})
    if not isinstance(algorithm, dict):
        return False, "threshold prior has no algorithm_fingerprint"
    if int(algorithm.get("schema_version", -1)) != THRESHOLD_PRIOR_SCHEMA_VERSION:
        return False, "threshold prior algorithm schema_version does not match"
    if str(algorithm.get("method_name", "")) != THRESHOLD_PRIOR_METHOD_NAME:
        return False, "threshold prior algorithm method_name does not match"
    for key, expected in (
        ("min_pairs", int(getattr(args, "threshold_prior_min_pairs", 100))),
        ("min_oriented_pairs", int(getattr(args, "threshold_prior_min_oriented_pairs", 60))),
    ):
        saved = algorithm.get(key)
        if saved is not None and int(saved) != int(expected):
            return False, f"threshold prior algorithm {key} does not match current calibration setting"
    build = artifact.get("build_parameters", {})
    if not isinstance(build, dict):
        return False, "threshold prior has no build_parameters"
    for key, expected in (
        ("every", int(getattr(args, "threshold_prior_every", 1))),
        ("global_frame_budget", int(getattr(args, "threshold_prior_global_frame_budget", 0))),
        ("global_frame_stride", int(getattr(args, "threshold_prior_global_frame_stride", 1))),
        ("block_size_frames", int(getattr(args, "threshold_prior_block_size_frames", 100))),
        ("file_chunk_size", int(getattr(args, "threshold_prior_file_chunk_size", 500))),
        ("max_block_histograms", int(getattr(args, "threshold_prior_max_block_histograms", 256))),
    ):
        saved = build.get(key)
        if saved is None:
            return False, f"threshold prior has no build parameter {key}"
        if int(saved) != int(expected):
            return False, f"threshold prior build parameter {key} does not match current setting"
    provenance = artifact.get("input_provenance", {})
    if not isinstance(provenance, dict):
        return False, "threshold prior has no input_provenance"
    if args.gb_param_file is not None:
        expected_sha = sha256_file(args.gb_param_file)
        saved_sha = provenance.get("gb_param_file_sha256")
        if saved_sha is None:
            return False, "threshold prior has no GB parameter file sha256"
        if expected_sha != saved_sha:
            return False, "threshold prior GB parameter file sha256 does not match current input"
    expected_mesogen = int(float(option_value(forwarded, "--mesogen-type", 1)))
    saved_mesogen = provenance.get("mesogen_type")
    if saved_mesogen is not None and int(saved_mesogen) != expected_mesogen:
        return False, "threshold prior mesogen_type does not match current analysis"
    expected_axis = option_value(forwarded, "--axis", "auto")
    saved_axis = provenance.get("axis")
    if saved_axis is not None and str(saved_axis) != str(expected_axis):
        return False, "threshold prior axis setting does not match current analysis"
    expected_anchor_types = option_value(forwarded, "--anchor-types", "2,3")
    saved_anchor_types = provenance.get("anchor_types")
    if saved_anchor_types is not None and str(saved_anchor_types) != str(expected_anchor_types):
        return False, "threshold prior anchor_types setting does not match current analysis"
    expected_s_excl = int(float(option_value(forwarded, "--s-excl", 1)))
    saved_s_excl = provenance.get("s_excl")
    if saved_s_excl is not None and int(saved_s_excl) != expected_s_excl:
        return False, "threshold prior s_excl setting does not match current analysis"
    expected_r_energy_cap = option_value(forwarded, "--r-energy-cap", "auto")
    saved_r_energy_cap = provenance.get("r_energy_cap")
    if saved_r_energy_cap is not None and str(saved_r_energy_cap) != str(expected_r_energy_cap):
        return False, "threshold prior r_energy_cap setting does not match current analysis"
    for option, hash_key, label in (
        ("--local-pair-file", "local_pair_file_sha256", "local pair file"),
        ("--exclude-pair-file", "exclude_pair_file_sha256", "exclude pair file"),
    ):
        saved_sha = provenance.get(hash_key)
        current_path = option_value(forwarded, option, None)
        expected_sha = sha256_file(Path(str(current_path))) if current_path else None
        if saved_sha != expected_sha:
            return False, f"threshold prior {label} sha256 does not match current analysis"
    return True, "threshold prior fingerprint matches current GB parameter file, mesogen type, axis, local/exclude pair files, and pair-scan settings"


def apply_threshold_prior_to_forwarded(
    forwarded: Sequence[str],
    *,
    prior_path: Path,
    gate: str,
) -> Tuple[List[str], Dict[str, object]]:
    _ = gate
    artifact = load_threshold_prior(prior_path)
    status = str(artifact.get("calibration_status", "unknown"))
    apply_allowed = bool(artifact.get("apply_allowed", False))
    updated = list(forwarded)
    recommended = artifact.get("recommended", {})
    if not isinstance(recommended, dict):
        raise SystemExit(f"Threshold prior has no recommended thresholds: {prior_path}")
    if not apply_allowed:
        print(f"[WARN] threshold prior status={status}; applying recommended values anyway under LC-Pearl V2 lobe-prior policy. See {prior_path}")
    for key, option in (
        ("gb_off_strength", "--gb-off-strength"),
        ("gb_on_strength", "--gb-on-strength"),
        ("p2_cut", "--p2-cut"),
        ("gb_core_strength", "--gb-core-strength"),
        ("p2_core_cut", "--p2-core-cut"),
        ("gb_strict_core_strength", "--gb-strict-core-strength"),
        ("p2_strict_core_cut", "--p2-strict-core-cut"),
    ):
        if key in recommended:
            updated = replace_or_append_option(updated, option, recommended[key])
    artifact["applied"] = True
    artifact["applied_thresholds"] = current_thresholds_from_forwarded(updated)
    artifact["applied_policy"] = "LC-Pearl V2 applies available streaming lobe-prior thresholds before main analysis; low status is kept as an audit note, not a gate."
    return updated, artifact


def build_threshold_prior_command(
    *,
    args: argparse.Namespace,
    forwarded: Sequence[str],
    input_args: Sequence[str],
    output_dir: Path,
    prior_file: Path,
) -> List[str]:
    thresholds = current_thresholds_from_forwarded(forwarded)
    command = [
        sys.executable,
        str(THRESHOLD_PRIOR_SCRIPT),
        *(input_args or ["."]),
        "--gb-param-file",
        str(args.gb_param_file),
        "--output-dir",
        str(output_dir),
        "--global-threshold-file",
        str(prior_file),
        "--every",
        str(getattr(args, "threshold_prior_every", 1)),
        "--global-frame-stride",
        str(getattr(args, "threshold_prior_global_frame_stride", 1)),
        "--global-frame-budget",
        str(getattr(args, "threshold_prior_global_frame_budget", 0)),
        "--block-size-frames",
        str(getattr(args, "threshold_prior_block_size_frames", 100)),
        "--file-chunk-size",
        str(getattr(args, "threshold_prior_file_chunk_size", 500)),
        "--max-block-histograms",
        str(getattr(args, "threshold_prior_max_block_histograms", 256)),
        "--audit-example-pairs",
        str(getattr(args, "threshold_prior_audit_example_pairs", 5000)),
        "--workers",
        str(getattr(args, "threshold_prior_workers", "auto")),
        "--sample-seed",
        str(getattr(args, "threshold_prior_seed", 20260429)),
        "--current-gb-off",
        str(thresholds["gb_off_strength"]),
        "--current-gb-on",
        str(thresholds["gb_on_strength"]),
        "--current-p2-cut",
        str(thresholds["p2_cut"]),
        "--current-gb-core",
        str(thresholds["gb_core_strength"]),
        "--current-p2-core-cut",
        str(thresholds["p2_core_cut"]),
        "--current-gb-strict-core",
        str(thresholds["gb_strict_core_strength"]),
        "--current-p2-strict-core-cut",
        str(thresholds["p2_strict_core_cut"]),
        "--current-s2-cut",
        str(thresholds["robust_min_s2"]),
        "--n-min",
        str(thresholds["n_min"]),
        "--min-pairs",
        str(getattr(args, "threshold_prior_min_pairs", 100)),
        "--min-oriented-pairs",
        str(getattr(args, "threshold_prior_min_oriented_pairs", 60)),
    ]
    append_forwarded_pair_scan_options(command, forwarded)
    return command


def ensure_and_apply_threshold_prior(
    *,
    args: argparse.Namespace,
    forwarded: Sequence[str],
    input_args: Sequence[str],
    potential_status: str,
) -> Tuple[List[str], Optional[Path], Optional[Dict[str, object]]]:
    mode = str(getattr(args, "threshold_prior_mode", "auto"))
    if mode == "off":
        return list(forwarded), None, None
    if args.gb_param_file is None:
        if mode == "require":
            raise SystemExit("--threshold-prior-mode require needs --gb-param-file")
        print("[WARN] --gb-param-file not provided; threshold prior was not used.")
        return list(forwarded), None, None
    if potential_status != "validated":
        message = (
            "LC Domain-Pearl V2 threshold prior requires validated Gay-Berne potential before automatic use; "
            f"current potential_status={potential_status}. Run validation first or set --threshold-prior-mode off for exploratory current-threshold analysis."
        )
        if mode in {"require", "refresh"}:
            raise SystemExit(message)
        print(f"[WARN] {message}")
        return list(forwarded), None, None
    prior_file = (getattr(args, "threshold_prior_file", None) or (args.output_root.resolve() / "threshold_prior" / "global_thresholds.json")).resolve()
    output_dir = (getattr(args, "threshold_prior_output_dir", None) or prior_file.parent).resolve()
    use_existing = prior_file.exists() and mode in {"auto", "require"}
    if use_existing:
        artifact = load_threshold_prior(prior_file)
        matches, reason = threshold_prior_matches(artifact, args=args, forwarded=forwarded)
        if matches:
            print(f"[OK] threshold prior cache hit: {prior_file} ({reason})")
            updated, applied_artifact = apply_threshold_prior_to_forwarded(
                forwarded,
                prior_path=prior_file,
                gate="off",
            )
            return updated, prior_file, applied_artifact
        if mode == "require":
            raise SystemExit(f"threshold prior does not match current inputs: {reason}. Use refresh mode or rebuild {prior_file}")
        print(f"[WARN] threshold prior cache miss: {reason}; rebuilding {prior_file}")
    elif mode == "require":
        raise SystemExit(f"required threshold prior not found: {prior_file}")
    command = build_threshold_prior_command(
        args=args,
        forwarded=forwarded,
        input_args=input_args,
        output_dir=output_dir,
        prior_file=prior_file,
    )
    run_command(command)
    updated, artifact = apply_threshold_prior_to_forwarded(
        forwarded,
        prior_path=prior_file,
        gate="off",
    )
    return updated, prior_file, artifact

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Project-level LC domain/pearl pipeline wrapper. Use '--' before arguments "
            "that should be forwarded to liquid_crystal_aggregation.py."
        )
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run LC-Pearl analysis with optional preflight validation and streaming threshold prior.")
    run.add_argument("--output-root", type=Path, required=True, help="Same output root passed to the analysis script.")
    run.add_argument("--gb-param-file", type=Path, default=None, help="LAMMPS input containing Gay-Berne pair parameters.")
    run.add_argument("--threshold-prior-mode", choices=["off", "auto", "require", "refresh"], default="auto", help="LC Domain-Pearl V2: use or build a streaming global threshold prior before main analysis.")
    run.add_argument("--threshold-prior-file", type=Path, default=None, help="Global threshold prior JSON. Defaults to output_root/threshold_prior/global_thresholds.json.")
    run.add_argument("--threshold-prior-output-dir", type=Path, default=None, help="Output directory for streaming threshold prior diagnostics.")
    run.add_argument("--threshold-prior-every", type=int, default=1)
    run.add_argument("--threshold-prior-global-frame-stride", type=int, default=1, help="Explicit global frame stride for threshold prior. 1 = all candidate frames; 10 = every 10th; 100 = every 100th. Takes precedence over --threshold-prior-global-frame-budget when > 1.")
    run.add_argument("--threshold-prior-global-frame-budget", type=int, default=0, help="0 means stream all frames; positive values use deterministic global stride.")
    run.add_argument("--threshold-prior-block-size-frames", type=int, default=100)
    run.add_argument("--threshold-prior-file-chunk-size", type=int, default=500, help="Number of dump files grouped into one threshold-prior worker task.")
    run.add_argument("--threshold-prior-max-block-histograms", type=int, default=256, help="Maximum block histograms saved for diagnostics; full block_count is still recorded.")
    run.add_argument("--threshold-prior-audit-example-pairs", type=int, default=5000)
    run.add_argument("--threshold-prior-workers", default="auto")
    run.add_argument("--threshold-prior-seed", type=int, default=20260429)
    run.add_argument("--threshold-prior-min-pairs", type=int, default=100)
    run.add_argument("--threshold-prior-min-oriented-pairs", type=int, default=60)
    run.add_argument("--potential-validation", choices=["off", "cache", "require", "refresh"], default="cache", help="Validate or reuse a cached Gay-Berne potential artifact before analysis.")
    run.add_argument("--potential-cache-dir", type=Path, default=None, help="Directory for verified potential artifacts and LAMMPS validation templates. Defaults to output_root/potential_validation.")
    run.add_argument("--verified-potential-file", type=Path, default=None, help="Explicit verified potential artifact path. Defaults to potential-cache-dir/verified_potential.json.")
    run.add_argument("--lammps-executable", default=None, help="LAMMPS executable used in suggested standalone validation commands and stored as fingerprint provenance.")
    run.add_argument("--python-pair-energy-table", type=Path, default=None, help="Python pair energy table to compare against a trusted LAMMPS run0/pe table.")
    run.add_argument("--lammps-pair-energy-table", type=Path, default=None, help="Trusted LAMMPS pair energy table generated from run0/pe or equivalent.")
    run.add_argument("--validation-max-abs-delta", type=float, default=1e-6, help="Maximum allowed absolute pair-energy delta for marking a potential artifact validated.")
    run.add_argument("analysis_args", nargs=argparse.REMAINDER, help="Arguments forwarded to liquid_crystal_aggregation.py after '--'.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command != "run":
        raise SystemExit("unknown command")
    if args.threshold_prior_global_frame_budget < 0:
        raise SystemExit("--threshold-prior-global-frame-budget must be non-negative")
    if args.threshold_prior_global_frame_stride <= 0:
        raise SystemExit("--threshold-prior-global-frame-stride must be positive")
    if args.threshold_prior_file_chunk_size <= 0:
        raise SystemExit("--threshold-prior-file-chunk-size must be positive")
    if args.threshold_prior_max_block_histograms < 0:
        raise SystemExit("--threshold-prior-max-block-histograms must be non-negative")
    if args.threshold_prior_min_pairs < 0:
        raise SystemExit("--threshold-prior-min-pairs must be non-negative")
    if args.threshold_prior_min_oriented_pairs < 0:
        raise SystemExit("--threshold-prior-min-oriented-pairs must be non-negative")
    forwarded = list(args.analysis_args)
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    if "--output-root" not in forwarded and not any(item.startswith("--output-root=") for item in forwarded):
        forwarded.extend(["--output-root", str(args.output_root)])
    if args.gb_param_file is not None and "--gb-param-file" not in forwarded and not any(item.startswith("--gb-param-file=") for item in forwarded):
        forwarded.extend(["--gb-param-file", str(args.gb_param_file)])

    input_args = forwarded_positionals(forwarded)
    potential_artifact_path = ensure_potential_validation(args, forwarded)
    potential_status = potential_validation_status(potential_artifact_path)

    threshold_prior_path: Optional[Path] = None
    threshold_prior_artifact: Optional[Dict[str, object]] = None
    forwarded, threshold_prior_path, threshold_prior_artifact = ensure_and_apply_threshold_prior(
        args=args,
        forwarded=forwarded,
        input_args=input_args,
        potential_status=potential_status,
    )

    run_command([sys.executable, str(ANALYSIS_SCRIPT), *forwarded])
    attach_potential_validation_to_run_summary(args.output_root, potential_artifact_path)
    if threshold_prior_path is not None:
        threshold_prior_diagnostic_copies = sync_threshold_prior_diagnostics(
            args.output_root,
            threshold_prior_path,
            threshold_prior_artifact,
        )
        summary_path = args.output_root.resolve() / "run_summary.json"
        if summary_path.exists():
            try:
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
                summary["threshold_prior"] = {
                    "path": str(threshold_prior_path),
                    "calibration_status": threshold_prior_artifact.get("calibration_status") if threshold_prior_artifact else None,
                    "apply_allowed": threshold_prior_artifact.get("apply_allowed") if threshold_prior_artifact else None,
                    "applied": threshold_prior_artifact.get("applied") if threshold_prior_artifact else None,
                    "applied_policy": threshold_prior_artifact.get("applied_policy") if threshold_prior_artifact else None,
                    "applied_thresholds": threshold_prior_artifact.get("applied_thresholds") if threshold_prior_artifact else None,
                    "recommended": threshold_prior_artifact.get("recommended") if threshold_prior_artifact else None,
                    "warnings": threshold_prior_artifact.get("warnings") if threshold_prior_artifact else None,
                    "sample_sizes": threshold_prior_artifact.get("sample_sizes") if threshold_prior_artifact else None,
                    "method_name": threshold_prior_artifact.get("method_name") if threshold_prior_artifact else None,
                    "diagnostic_copies": threshold_prior_diagnostic_copies,
                }
                summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
            except json.JSONDecodeError:
                pass


if __name__ == "__main__":
    main()
