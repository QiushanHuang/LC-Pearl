#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import shutil
import site
import shlex
import subprocess
import sys
import sysconfig
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib  # type: ignore[no-redef]

import lc_pearl_cli


PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = lc_pearl_cli.script_path("lc_domain_pearl_pipeline.py").parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import lc_domain_pearl_pipeline as pipeline  # noqa: E402
import lc_topology_prepare  # noqa: E402


DEFAULT_PREFLIGHT_DIR = "lc_pearl_preflight"
DEFAULT_CONFIG_NAME = "lc_pearl_config.toml"
CANONICAL_GB_PARAM = "gb_param_source.in"


def load_toml(path: Path) -> Dict[str, Any]:
    with path.open("rb") as handle:
        return dict(tomllib.load(handle))


def find_template_preflight_dir() -> Path:
    candidates = [
        PROJECT_ROOT / "templates" / DEFAULT_PREFLIGHT_DIR,
        PROJECT_ROOT / "share" / "lc-pearl" / "templates" / DEFAULT_PREFLIGHT_DIR,
        Path(sysconfig.get_path("data")) / "share" / "lc-pearl" / "templates" / DEFAULT_PREFLIGHT_DIR,
        Path(sys.prefix) / "share" / "lc-pearl" / "templates" / DEFAULT_PREFLIGHT_DIR,
        Path(site.getuserbase()) / "share" / "lc-pearl" / "templates" / DEFAULT_PREFLIGHT_DIR,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    searched = "\n  - ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"LC-Pearl preflight template not found. Searched:\n  - {searched}")


def create_preflight_from_template(destination: Path) -> None:
    template = find_template_preflight_dir()
    shutil.copytree(template, destination)


def table(config: Mapping[str, Any], name: str) -> Dict[str, Any]:
    value = config.get(name, {})
    return dict(value) if isinstance(value, dict) else {}


def is_blank(value: object) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text.upper().startswith("REPLACE_") or text.upper().startswith("PUT_")


def resolve_path(raw: object, *, base: Path, allow_blank: bool = True) -> Optional[str]:
    if is_blank(raw):
        if allow_blank:
            return None
        raise ValueError("blank path is not allowed")
    path = Path(str(raw)).expanduser()
    if path.is_absolute():
        return str(path)
    return str((base / path).resolve())


def read_lammps_executable(preflight_dir: Path, paths: Mapping[str, Any]) -> Optional[str]:
    explicit = paths.get("lammps_executable")
    if not is_blank(explicit):
        return resolve_path(explicit, base=preflight_dir, allow_blank=False)
    text_file = preflight_dir / "lammps" / "lammps_executable.txt"
    if not text_file.exists():
        return None
    for line in text_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return resolve_path(stripped, base=preflight_dir, allow_blank=False)
    return None


def lammps_input_has_gayberne(path: Path) -> bool:
    text = path.read_text(encoding="utf-8", errors="ignore").lower()
    return "pair_style" in text and "gayberne" in text and "pair_coeff" in text


def discover_gb_param_file(preflight_dir: Path, configured: object) -> Tuple[Optional[Path], List[str], List[str]]:
    notes: List[str] = []
    problems: List[str] = []
    lammps_dir = preflight_dir / "lammps"
    canonical = lammps_dir / CANONICAL_GB_PARAM
    configured_path = resolve_path(configured, base=preflight_dir) if not is_blank(configured) else None
    if configured_path is not None and Path(configured_path).exists():
        path = Path(configured_path)
        if path.resolve() != canonical.resolve():
            canonical.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, canonical)
            notes.append(f"normalized GB parameter input: {path} -> {canonical}")
            return canonical, notes, problems
        return path, notes, problems

    candidates: List[Path] = []
    if lammps_dir.exists():
        for suffix in ("*.in", "*.lmp", "*.include", "*.inc"):
            candidates.extend(path for path in lammps_dir.glob(suffix) if path.is_file() and path.name != CANONICAL_GB_PARAM)
    gb_candidates = [path for path in sorted(set(candidates)) if lammps_input_has_gayberne(path)]
    if canonical.exists():
        return canonical, notes, problems
    if len(gb_candidates) == 1:
        canonical.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(gb_candidates[0], canonical)
        notes.append(f"auto-detected and normalized GB parameter input: {gb_candidates[0].name} -> lammps/{CANONICAL_GB_PARAM}")
        return canonical, notes, problems
    if len(gb_candidates) > 1:
        choices = ", ".join(path.name for path in gb_candidates)
        problems.append(f"multiple Gay-Berne input candidates in lc_pearl_preflight/lammps: {choices}; set [paths].gb_param_file explicitly")
        return None, notes, problems
    return None, notes, problems


def parse_special_lj_from_lammps_input(path: Path) -> Optional[Tuple[float, float, float]]:
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if not parts or parts[0] != "special_bonds":
            continue
        try:
            index = parts.index("lj")
        except ValueError:
            continue
        values = parts[index + 1 : index + 4]
        if len(values) != 3:
            continue
        try:
            return (float(values[0]), float(values[1]), float(values[2]))
        except ValueError:
            continue
    return None


def parse_special_lj(value: object, gb_param_file: Optional[Path]) -> Tuple[float, float, float]:
    text = "auto" if is_blank(value) else str(value).strip()
    if text.lower() == "auto":
        if gb_param_file is not None and gb_param_file.exists():
            parsed = parse_special_lj_from_lammps_input(gb_param_file)
            if parsed is not None:
                return parsed
        return (0.0, 1.0, 1.0)
    parts = [item.strip() for item in text.split(",") if item.strip()]
    if len(parts) != 3:
        raise ValueError("[topology].special_lj must be 'auto' or three comma-separated numbers")
    return (float(parts[0]), float(parts[1]), float(parts[2]))


def discover_topology_data_file(preflight_dir: Path, config: Mapping[str, Any]) -> Optional[Path]:
    topology = table(config, "topology")
    raw = topology.get("data_file")
    if not is_blank(raw) and str(raw).strip().lower() != "auto":
        resolved = resolve_path(raw, base=preflight_dir)
        return Path(resolved) if resolved else None
    topology_dir = preflight_dir / "topology"
    if not topology_dir.exists():
        return None
    candidates: List[Path] = []
    for suffix in ("*.data", "*.dat"):
        candidates.extend(path for path in topology_dir.glob(suffix) if path.is_file())
    candidates = sorted(set(candidates))
    if len(candidates) == 1:
        return candidates[0]
    return None


def maybe_generate_topology_pair_files(
    preflight_dir: Path,
    config: Mapping[str, Any],
    *,
    gb_param_file: Optional[Path],
) -> Tuple[List[str], List[str]]:
    notes: List[str] = []
    problems: List[str] = []
    topology_dir = preflight_dir / "topology"
    local_path = topology_dir / "local_pairs.tsv"
    exclude_path = topology_dir / "exclude_pairs.tsv"
    if local_path.exists() and exclude_path.exists():
        return notes, problems
    data_file = discover_topology_data_file(preflight_dir, config)
    if data_file is None:
        notes.append(
            "no raw topology data file found; local_pairs.tsv/exclude_pairs.tsv were not generated. "
            "Put a LAMMPS data file in lc_pearl_preflight/topology/ if you want automatic topology conversion."
        )
        return notes, problems
    if not data_file.exists():
        problems.append(f"configured [topology].data_file does not exist: {data_file}")
        return notes, problems
    topology = table(config, "topology")
    try:
        special_lj = parse_special_lj(topology.get("special_lj", "auto"), gb_param_file)
    except ValueError as exc:
        problems.append(str(exc))
        return notes, problems
    topology_dir.mkdir(parents=True, exist_ok=True)
    manifest = lc_topology_prepare.prepare_topology(
        data_file=data_file,
        output_root=topology_dir,
        special_lj=special_lj,
        mesogen_type=int(table(config, "analysis").get("mesogen_type", 1)),
        anchor_type=int(str(table(config, "analysis").get("anchor_type", table(config, "analysis").get("anchor_types", "3"))).split(",")[0]),
    )
    generated_local = topology_dir / "local_special_pairs.tsv"
    if generated_local.exists() and not local_path.exists():
        shutil.copy2(generated_local, local_path)
    notes.append(
        "generated topology pair files from raw LAMMPS data: "
        f"excluded={manifest.get('counts', {}).get('excluded_pairs', 0)}, "
        f"local={manifest.get('counts', {}).get('local_nonzero_special_pairs', 0)}"
    )
    return notes, problems


def first_dump_file(workdir: Path, pattern: str, recursive: bool) -> Optional[Path]:
    iterator = workdir.rglob(pattern) if recursive else workdir.glob(pattern)
    candidates = [
        path
        for path in iterator
        if path.is_file()
        and DEFAULT_PREFLIGHT_DIR not in path.parts
        and "lc_domain_pearl_v1_output" not in path.parts
        and "lc_domain_pearl_v2_output" not in path.parts
    ]
    return sorted(candidates)[0] if candidates else None


def resolve_representative_dump(
    *,
    workdir: Path,
    preflight_dir: Path,
    config: Mapping[str, Any],
    override: Optional[str],
) -> Optional[str]:
    validation = table(config, "validation")
    paths = table(config, "paths")
    raw = override or validation.get("representative_dump") or paths.get("representative_dump") or "auto"
    if is_blank(raw) or str(raw).strip().lower() == "auto":
        pattern = str(paths.get("pattern", "*.dump"))
        recursive = bool(paths.get("recursive", False))
        first = first_dump_file(workdir, pattern, recursive)
        return str(first.resolve()) if first else None
    path = Path(str(raw)).expanduser()
    if path.is_absolute():
        return str(path)
    candidate = workdir / path
    if candidate.exists():
        return str(candidate.resolve())
    return str((preflight_dir / path).resolve())


def default_optional_pair_file(preflight_dir: Path, filename: str) -> Optional[str]:
    path = preflight_dir / "topology" / filename
    if path.exists():
        return str(path.resolve())
    return None


def resolve_preflight_relative(preflight_dir: Path, raw: object, default: str) -> str:
    value = default if is_blank(raw) else str(raw)
    path = Path(value).expanduser()
    if path.is_absolute():
        return str(path.resolve())
    return str((preflight_dir / path).resolve())


def build_current_config(
    raw_config: Mapping[str, Any],
    *,
    workdir: Path,
    preflight_dir: Path,
    args: argparse.Namespace,
) -> Tuple[Dict[str, Any], List[str], List[str]]:
    config = copy.deepcopy(dict(raw_config))
    config.setdefault("runtime", {})
    config.setdefault("paths", {})
    config.setdefault("pipeline", {})
    config.setdefault("analysis", {})
    config.setdefault("validation", {})
    config.setdefault("threshold_prior", {})
    paths = config["paths"]
    analysis = config["analysis"]
    validation = config["validation"]
    threshold_prior = config["threshold_prior"]
    missing: List[str] = []
    notes: List[str] = []

    paths["workdir"] = str(workdir.resolve())
    paths["input"] = args.input or str(paths.get("input", ".") or ".")
    paths["output_root"] = args.output_root or str(paths.get("output_root", "lc_domain_pearl_v2_output") or "lc_domain_pearl_v2_output")
    if args.pattern:
        paths["pattern"] = args.pattern
    paths.setdefault("pattern", "*.dump")

    gb_param = args.gb_param_file or paths.get("gb_param_file")
    gb_path, gb_notes, gb_problems = discover_gb_param_file(preflight_dir, gb_param)
    notes.extend(gb_notes)
    missing.extend(gb_problems)
    if gb_path is None or not gb_path.exists():
        missing.append(
            "LAMMPS GB parameter input: put an original .in/.lmp containing pair_style gayberne and type 1-1 pair_coeff in lc_pearl_preflight/lammps/"
        )
    else:
        paths["gb_param_file"] = str(gb_path.resolve())

    lammps_executable = args.lammps_executable or read_lammps_executable(preflight_dir, paths)
    if is_blank(lammps_executable):
        missing.append(
            "lammps_executable: set [paths].lammps_executable or put the executable path in lc_pearl_preflight/lammps/lammps_executable.txt"
        )
    else:
        paths["lammps_executable"] = str(lammps_executable)

    paths["potential_cache_dir"] = str((preflight_dir / "validation").resolve())
    paths["verified_potential_file"] = str((preflight_dir / "validation" / "verified_potential.json").resolve())
    threshold_prior.setdefault("mode", "auto")
    threshold_prior["file"] = resolve_preflight_relative(preflight_dir, threshold_prior.get("file"), "thresholds/global_thresholds.json")
    threshold_prior["output_dir"] = resolve_preflight_relative(preflight_dir, threshold_prior.get("output_dir"), "thresholds")

    representative = resolve_representative_dump(
        workdir=workdir,
        preflight_dir=preflight_dir,
        config=config,
        override=args.representative_dump,
    )
    if representative is None or not Path(representative).exists():
        missing.append(
            "representative dump: set [validation].representative_dump or keep it as 'auto' with at least one dump matching [paths].pattern"
        )
    else:
        validation["representative_dump"] = representative

    topology_notes, topology_problems = maybe_generate_topology_pair_files(
        preflight_dir,
        config,
        gb_param_file=gb_path,
    )
    notes.extend(topology_notes)
    missing.extend(topology_problems)
    topology_data_missing = any("no raw topology data file found" in note for note in topology_notes)

    for key, filename in (("local_pair_file", "local_pairs.tsv"), ("exclude_pair_file", "exclude_pairs.tsv")):
        raw_value = analysis.get(key)
        if is_blank(raw_value):
            default_value = default_optional_pair_file(preflight_dir, filename)
            if default_value:
                analysis[key] = default_value
            else:
                analysis.pop(key, None)
                if not topology_data_missing:
                    notes.append(f"derived topology/{filename} not available; put a raw LAMMPS data file in topology/ to generate it automatically")
        else:
            resolved = resolve_path(raw_value, base=preflight_dir)
            if resolved is None or not Path(resolved).exists():
                missing.append(f"configured analysis.{key} does not exist: {raw_value}")
            else:
                analysis[key] = resolved

    config["runtime"].setdefault("python", "auto")
    return config, missing, notes


def pipeline_forwarded_args(config: Mapping[str, Any]) -> List[str]:
    command = lc_pearl_cli.build_pipeline_command(config)
    if "--" not in command:
        return []
    forwarded = list(command[command.index("--") + 1 :])
    output_root = lc_pearl_cli.pipeline_output_root(config)
    gb_param_file = table(config, "paths").get("gb_param_file")
    if not pipeline.has_option(forwarded, "--output-root"):
        forwarded.extend(["--output-root", output_root])
    if gb_param_file and not pipeline.has_option(forwarded, "--gb-param-file"):
        forwarded.extend(["--gb-param-file", str(gb_param_file)])
    return forwarded


def validation_cache_status(config: Mapping[str, Any]) -> Tuple[bool, str]:
    paths = table(config, "paths")
    artifact_path = Path(lc_pearl_cli.pipeline_verified_potential_file(config))
    artifact = pipeline.load_validation_artifact(artifact_path)
    if artifact is None:
        return False, f"missing or unreadable verified potential: {artifact_path}"
    try:
        fingerprint = pipeline.build_potential_fingerprint(
            gb_param_file=Path(str(paths["gb_param_file"])),
            forwarded=pipeline_forwarded_args(config),
            lammps_executable=str(paths.get("lammps_executable", "")),
            mesogen_type=int(float(table(config, "analysis").get("mesogen_type", 1))),
        )
    except Exception as exc:
        return False, f"could not build validation fingerprint: {exc}"
    if pipeline.validation_artifact_matches(artifact, fingerprint):
        return True, f"validated cache hit: {artifact_path}"
    return False, f"verified potential exists but does not match current fingerprint: {artifact_path}"


def print_preflight_missing(preflight_dir: Path, missing: Sequence[str], notes: Sequence[str]) -> None:
    print("[LC-Pearl preflight]")
    print(f"Directory: {preflight_dir}")
    if missing:
        print("\nMissing required inputs:")
        for item in missing:
            print(f"- {item}")
    if notes:
        print("\nOptional or derived inputs not available:")
        for item in notes:
            print(f"- {item}")
    print("\nRequired file meanings:")
    print("- lammps/*.in or *.lmp: raw LAMMPS input; preflight detects pair_style gayberne/pair_coeff and standardizes it to lammps/gb_param_source.in.")
    print("- lammps/lammps_executable.txt: path used to run LAMMPS run 0 / microstate validation.")
    print("- representative dump: used only for validation; auto selects the first dump matching [paths].pattern.")
    print("- topology/*.data or *.dat: optional raw LAMMPS data; preflight converts Bonds/special_bonds into local_pairs.tsv and exclude_pairs.tsv.")
    print("\nEdit lc_pearl_preflight/lc_pearl_config.toml, then run the same command again.")


def run_command(command: List[str], *, cwd: Path, dry_run: bool) -> None:
    print("[LC-Pearl] " + shlex.join(command))
    if not dry_run:
        subprocess.run(command, cwd=str(cwd), check=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run LC-Pearl from the current dump directory using a portable lc_pearl_preflight folder. "
            "Default mode is auto: validate if needed, then run the full pipeline."
        )
    )
    parser.add_argument(
        "mode",
        nargs="?",
        choices=["auto", "validate", "run", "init", "print-auto", "print-validate", "print-run"],
        default="auto",
    )
    parser.add_argument("--preflight-dir", default=DEFAULT_PREFLIGHT_DIR)
    parser.add_argument("--config", default=None, help="Defaults to <preflight-dir>/lc_pearl_config.toml.")
    parser.add_argument("--input", default=None, help="Input path passed to analysis; default from config, normally '.'.")
    parser.add_argument("--pattern", default=None, help="Dump filename pattern; overrides [paths].pattern.")
    parser.add_argument("--output-root", default=None, help="Output root relative to current dump directory.")
    parser.add_argument("--gb-param-file", default=None, help="Override GB/LAMMPS parameter input file.")
    parser.add_argument("--lammps-executable", default=None, help="Override LAMMPS executable path.")
    parser.add_argument("--representative-dump", default=None, help="Override validation dump file; default is auto.")
    parser.add_argument("--refresh-validation", action="store_true", help="Force run0 validation even if the cache matches.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    workdir = Path.cwd()
    preflight_dir = (workdir / args.preflight_dir).resolve()
    if args.mode == "init":
        if preflight_dir.exists():
            print(f"[LC-Pearl] preflight already exists: {preflight_dir}")
        else:
            create_preflight_from_template(preflight_dir)
        config_path = preflight_dir / DEFAULT_CONFIG_NAME
        config, missing, notes = build_current_config(
            load_toml(config_path),
            workdir=workdir,
            preflight_dir=preflight_dir,
            args=args,
        )
        _ = config
        print_preflight_missing(preflight_dir, missing, notes)
        return
    if not preflight_dir.exists():
        create_preflight_from_template(preflight_dir)
        config_path = preflight_dir / DEFAULT_CONFIG_NAME
        config, missing, notes = build_current_config(
            load_toml(config_path),
            workdir=workdir,
            preflight_dir=preflight_dir,
            args=args,
        )
        _ = config
        print_preflight_missing(preflight_dir, missing, notes)
        raise SystemExit(2)

    config_path = Path(args.config).expanduser() if args.config else preflight_dir / DEFAULT_CONFIG_NAME
    if not config_path.is_absolute():
        config_path = (workdir / config_path).resolve()
    if not config_path.exists():
        print_preflight_missing(preflight_dir, [f"configuration file not found: {config_path}"], [])
        raise SystemExit(2)

    config, missing, notes = build_current_config(
        load_toml(config_path),
        workdir=workdir,
        preflight_dir=preflight_dir,
        args=args,
    )
    if missing:
        print_preflight_missing(preflight_dir, missing, notes)
        raise SystemExit(2)
    for note in notes:
        print(f"[WARN] {note}")

    mode = args.mode[6:] if args.mode.startswith("print-") else args.mode
    dry_run = bool(args.dry_run or args.mode.startswith("print-"))
    valid, reason = validation_cache_status(config)
    if valid:
        print(f"[OK] {reason}")
    else:
        print(f"[WARN] {reason}")

    if mode in {"auto", "validate"} and (args.refresh_validation or not valid):
        run_command(lc_pearl_cli.build_validation_command(config, config_path), cwd=workdir, dry_run=dry_run)
        if not dry_run:
            valid_after, reason_after = validation_cache_status(config)
            if not valid_after:
                raise SystemExit(f"validation did not produce a matching cache: {reason_after}")
            print(f"[OK] {reason_after}")
    elif mode == "validate":
        print("[LC-Pearl] validation cache is already usable; add --refresh-validation to rerun LAMMPS validation.")

    if mode in {"auto", "run"}:
        run_command(lc_pearl_cli.build_pipeline_command(config, config_path), cwd=workdir, dry_run=dry_run)


if __name__ == "__main__":
    main()
