#!/usr/bin/env python3
from __future__ import annotations

import argparse
import site
import shlex
import subprocess
import sys
import sysconfig
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib  # type: ignore[no-redef]


PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"


def script_path(name: str) -> Path:
    candidates = [
        SCRIPTS_DIR / name,
        PROJECT_ROOT / "share" / "lc-pearl" / "scripts" / name,
        Path(sysconfig.get_path("data")) / "share" / "lc-pearl" / "scripts" / name,
        Path(sys.prefix) / "share" / "lc-pearl" / "scripts" / name,
        Path(site.getuserbase()) / "share" / "lc-pearl" / "scripts" / name,
    ]
    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists():
            return candidate
    searched = "\n  - ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"LC-Pearl helper script not found: {name}\nSearched:\n  - {searched}")


def load_config(path: Path) -> Dict[str, Any]:
    with path.open("rb") as handle:
        return dict(tomllib.load(handle))


def table(config: Mapping[str, Any], name: str) -> Dict[str, Any]:
    value = config.get(name, {})
    if isinstance(value, dict):
        return dict(value)
    return {}


def config_path(config: Mapping[str, Any], section: str, key: str, default: Optional[str] = None) -> Optional[str]:
    value = table(config, section).get(key, default)
    if value is None:
        return None
    return str(value)


def add_option(command: List[str], option: str, value: object) -> None:
    if value is None:
        return
    command.extend([option, str(value)])


def add_bool_flag(command: List[str], option: str, value: object) -> None:
    if bool(value):
        command.append(option)


def python_executable(config: Mapping[str, Any]) -> str:
    value = table(config, "runtime").get("python", "auto")
    if value is None or str(value).strip().lower() in {"", "auto", "sys.executable"}:
        return sys.executable
    return str(value)


def resolve_run_path(config: Mapping[str, Any], raw: Optional[str]) -> Optional[str]:
    if raw is None:
        return None
    path = Path(raw).expanduser()
    if path.is_absolute():
        return str(path)
    workdir = config_path(config, "paths", "workdir")
    if workdir:
        return str((Path(workdir).expanduser() / path).resolve())
    return str(path)


def pipeline_output_root(config: Mapping[str, Any]) -> str:
    output_root = config_path(config, "paths", "output_root", "lc_domain_pearl_v2_output")
    resolved = resolve_run_path(config, output_root)
    assert resolved is not None
    return resolved


def pipeline_potential_cache_dir(config: Mapping[str, Any]) -> str:
    explicit = config_path(config, "paths", "potential_cache_dir")
    if explicit:
        resolved = resolve_run_path(config, explicit)
        assert resolved is not None
        return resolved
    return str(Path(pipeline_output_root(config)) / "potential_validation")


def pipeline_verified_potential_file(config: Mapping[str, Any]) -> str:
    explicit = config_path(config, "paths", "verified_potential_file")
    if explicit:
        resolved = resolve_run_path(config, explicit)
        assert resolved is not None
        return resolved
    return str(Path(pipeline_potential_cache_dir(config)) / "verified_potential.json")


def resolve_config_path(config: Mapping[str, Any], raw: Optional[str]) -> Optional[str]:
    if raw is None:
        return None
    return resolve_run_path(config, raw)


def build_pipeline_command(config: Mapping[str, Any], config_path_for_context: Optional[Path] = None) -> List[str]:
    paths = table(config, "paths")
    pipeline = table(config, "pipeline")
    analysis = table(config, "analysis")
    threshold_prior = table(config, "threshold_prior")

    command = [
        python_executable(config),
        str(script_path("lc_domain_pearl_pipeline.py")),
        "run",
        "--output-root",
        pipeline_output_root(config),
    ]
    add_option(command, "--gb-param-file", paths.get("gb_param_file"))
    add_option(command, "--threshold-prior-mode", threshold_prior.get("mode", pipeline.get("threshold_prior_mode", "auto")))
    add_option(command, "--threshold-prior-file", resolve_config_path(config, threshold_prior.get("file", paths.get("threshold_prior_file"))))
    add_option(command, "--threshold-prior-output-dir", resolve_config_path(config, threshold_prior.get("output_dir", paths.get("threshold_prior_output_dir"))))
    add_option(command, "--threshold-prior-every", threshold_prior.get("every", 1))
    add_option(command, "--threshold-prior-global-frame-stride", threshold_prior.get("global_frame_stride", 1))
    add_option(command, "--threshold-prior-global-frame-budget", threshold_prior.get("global_frame_budget", 0))
    add_option(command, "--threshold-prior-block-size-frames", threshold_prior.get("block_size_frames", 100))
    add_option(command, "--threshold-prior-file-chunk-size", threshold_prior.get("file_chunk_size", 500))
    add_option(command, "--threshold-prior-max-block-histograms", threshold_prior.get("max_block_histograms", 256))
    add_option(command, "--threshold-prior-audit-example-pairs", threshold_prior.get("audit_example_pairs", 5000))
    add_option(command, "--threshold-prior-workers", threshold_prior.get("workers", "auto"))
    add_option(command, "--threshold-prior-seed", threshold_prior.get("seed", pipeline.get("threshold_prior_seed", 20260429)))
    add_option(command, "--threshold-prior-min-pairs", threshold_prior.get("min_pairs", 100))
    add_option(command, "--threshold-prior-min-oriented-pairs", threshold_prior.get("min_oriented_pairs", 60))
    add_option(command, "--potential-validation", pipeline.get("potential_validation", "cache"))
    add_option(command, "--potential-cache-dir", pipeline_potential_cache_dir(config))
    add_option(command, "--verified-potential-file", pipeline_verified_potential_file(config))
    add_option(command, "--lammps-executable", paths.get("lammps_executable"))
    command.append("--")
    input_path = resolve_run_path(config, str(paths.get("input", ".")))
    if input_path:
        command.append(input_path)
    add_option(command, "--pattern", paths.get("pattern", "*.dump"))
    if bool(paths.get("recursive", False)):
        command.append("--recursive")
    add_option(command, "--contact-mode", analysis.get("contact_mode", "gayberne"))
    add_option(command, "--mesogen-type", analysis.get("mesogen_type", 1))
    add_option(command, "--anchor-types", analysis.get("anchor_types", "3"))
    add_option(command, "--axis", analysis.get("axis", "auto"))
    add_option(command, "--workers", analysis.get("workers", "auto"))
    for key, option in (
        ("every", "--every"),
        ("r_cut", "--r-cut"),
        ("max_auto_r_cut", "--max-auto-r-cut"),
        ("cluster_cut", "--cluster-cut"),
        ("cluster_cut_shape_factor", "--cluster-cut-shape-factor"),
        ("cluster_min_size", "--cluster-min-size"),
        ("gb_off_strength", "--gb-off-strength"),
        ("gb_on_strength", "--gb-on-strength"),
        ("r_energy_cap", "--r-energy-cap"),
        ("p2_cut", "--p2-cut"),
        ("q_on", "--q-on"),
        ("q_off", "--q-off"),
        ("robust_min_s2", "--robust-min-s2"),
        ("robust_min_evidence", "--robust-min-evidence"),
        ("domain_min_lifetime", "--domain-min-lifetime"),
        ("stable_overlap_fraction", "--stable-overlap-fraction"),
        ("n_min", "--n-min"),
        ("s_excl", "--s-excl"),
        ("pearl_gap_cut", "--pearl-gap-cut"),
        ("pearl_min_cross_contacts", "--pearl-min-cross-contacts"),
        ("pearl_min_boundary_particles", "--pearl-min-boundary-particles"),
        ("pearl_max_aspect_ratio", "--pearl-max-aspect-ratio"),
        ("cluster_envelope_padding", "--cluster-envelope-padding"),
        ("edge_diagnostics_table", "--edge-diagnostics-table"),
        ("edge_diagnostics_sample_size", "--edge-diagnostics-sample-size"),
        ("local_pair_file", "--local-pair-file"),
        ("exclude_pair_file", "--exclude-pair-file"),
    ):
        add_option(command, option, analysis.get(key))
    for key, option in (
        ("require_nonlocal_robust", "--require-nonlocal-robust"),
        ("robustness_scan", "--robustness-scan"),
        ("no_ovito_labels", "--no-ovito-labels"),
        ("no_cluster_envelopes", "--no-cluster-envelopes"),
        ("no_contact_edges", "--no-contact-edges"),
        ("no_contact_segments", "--no-contact-segments"),
        ("write_frame_jsonl", "--write-frame-jsonl"),
        ("no_diagnostics", "--no-diagnostics"),
        ("no_track_across_files", "--no-track-across-files"),
    ):
        add_bool_flag(command, option, analysis.get(key, False))
    return command


def build_validation_command(config: Mapping[str, Any], config_path_for_context: Optional[Path] = None) -> List[str]:
    paths = table(config, "paths")
    validation = table(config, "validation")
    representative = validation.get("representative_dump") or paths.get("representative_dump")
    command = [
        python_executable(config),
        str(script_path("lc_lammps_run0_validate.py")),
        "--dump-file",
        str(resolve_run_path(config, str(representative or "."))),
        "--gb-param-file",
        str(paths.get("gb_param_file")),
        "--output-root",
        str(Path(pipeline_potential_cache_dir(config)) / "run0_validation"),
        "--verified-potential-file",
        pipeline_verified_potential_file(config),
        "--lammps-executable",
        str(paths.get("lammps_executable")),
    ]
    add_option(command, "--mpi-prefix", validation.get("mpi_prefix"))
    add_option(command, "--microstate-sample-count", validation.get("microstate_sample_count", 12))
    add_option(command, "--microstate-sample-percent", validation.get("microstate_sample_percent"))
    add_option(command, "--frame-index", validation.get("frame_index"))
    add_option(command, "--mesogen-type", validation.get("mesogen_type", table(config, "analysis").get("mesogen_type", 1)))
    add_option(command, "--local-pair-file", table(config, "analysis").get("local_pair_file"))
    add_option(command, "--exclude-pair-file", table(config, "analysis").get("exclude_pair_file"))
    return command


def run_command(command: List[str], *, cwd: Optional[str] = None, dry_run: bool = False) -> None:
    print("[LC-Pearl] " + shlex.join(command))
    if dry_run:
        return
    subprocess.run(command, check=True, cwd=cwd)


def main() -> None:
    parser = argparse.ArgumentParser(description="LC-Pearl TOML launcher for LC domain/pearl analysis.")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "validate", "print-run", "print-validate"):
        item = sub.add_parser(name)
        item.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    dry_run = str(args.command).startswith("print-")
    if args.command.endswith("validate"):
        command = build_validation_command(config, args.config)
    else:
        command = build_pipeline_command(config, args.config)
    run_command(command, cwd=config_path(config, "paths", "workdir"), dry_run=dry_run)


if __name__ == "__main__":
    main()
