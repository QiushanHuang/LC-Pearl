from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import lc_threshold_prior as prior  # noqa: E402
import lc_domain_pearl_pipeline as pipeline  # noqa: E402


def synthetic_lobe_histogram() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gb_edges = np.linspace(0.0, 1.5, 121)
    p2_edges = np.linspace(-0.5, 1.0, 121)
    gb_centers = 0.5 * (gb_edges[:-1] + gb_edges[1:])
    p2_centers = 0.5 * (p2_edges[:-1] + p2_edges[1:])
    x, y = np.meshgrid(gb_centers, p2_centers, indexing="ij")

    background = 3.0 * np.exp(-((x - 0.12) ** 2) / 0.010) * (0.7 + 0.3 * (y + 0.5))
    left_oriented_lobe = 500.0 * np.exp(-((x - 0.08) ** 2) / 0.002) * np.exp(-((y - 0.92) ** 2) / 0.018)
    right_oriented_lobe = 320.0 * np.exp(-((x - 0.72) ** 2) / 0.035) * np.exp(-((y - 0.91) ** 2) / 0.020)
    hist2d = background + left_oriented_lobe + right_oriented_lobe
    return gb_edges, p2_edges, hist2d


def test_lobe_split_uses_high_p2_2d_valley() -> None:
    gb_edges, p2_edges, hist2d = synthetic_lobe_histogram()
    result = prior.estimate_lobe_thresholds_from_histograms(
        gb_edges=gb_edges,
        p2_edges=p2_edges,
        hist2d=hist2d,
        current={"gb_off_strength": 0.12, "gb_on_strength": 0.30, "p2_cut": 0.70},
        min_pairs=100,
        min_oriented_pairs=60,
    )

    recommended = result["recommended"]
    assert result["calibration_status"] in {"medium", "high"}
    assert result["apply_allowed"] is True
    assert 0.12 <= recommended["gb_off_strength"] < recommended["gb_on_strength"]
    assert 0.22 <= recommended["gb_on_strength"] <= 0.50
    assert 0.45 <= recommended["p2_cut"] <= 0.65
    assert result["parameters"]["gb_on_strength"]["decision"] == "2d_high_p2_lobe_valley"
    assert result["parameters"]["gb_on_strength"]["orientation_gate_for_gb_histogram"] == recommended["p2_cut"]
    assert recommended["p2_core_cut"] == 0.71
    assert recommended["p2_strict_core_cut"] == 0.8
    assert recommended["gb_on_strength"] < recommended["gb_core_strength"] <= recommended["gb_strict_core_strength"]
    assert result["parameters"]["gb_core_strength"]["tier_name"] == "core_shoulder"
    assert result["parameters"]["gb_strict_core_strength"]["tier_name"] == "strict_core"
    assert result["parameters"]["gb_core_strength"]["decision"] in {
        "high_p2_core_shoulder_left_shoulder",
        "fallback_high_p2_core_shoulder_weighted_q25",
        "fallback_current_core_shoulder_threshold",
    }
    assert result["parameters"]["gb_strict_core_strength"]["decision"].startswith("high_p2_strict_core_weighted_q25")


def test_lobe_split_is_independent_of_current_thresholds() -> None:
    gb_edges, p2_edges, hist2d = synthetic_lobe_histogram()
    low_current = prior.estimate_lobe_thresholds_from_histograms(
        gb_edges=gb_edges,
        p2_edges=p2_edges,
        hist2d=hist2d,
        current={
            "gb_off_strength": 0.05,
            "gb_on_strength": 0.20,
            "p2_cut": 0.55,
            "gb_core_strength": 0.60,
            "p2_core_cut": 0.71,
            "gb_strict_core_strength": 0.90,
            "p2_strict_core_cut": 0.80,
        },
        min_pairs=100,
        min_oriented_pairs=60,
    )["recommended"]
    high_current = prior.estimate_lobe_thresholds_from_histograms(
        gb_edges=gb_edges,
        p2_edges=p2_edges,
        hist2d=hist2d,
        current={
            "gb_off_strength": 0.25,
            "gb_on_strength": 0.75,
            "p2_cut": 0.90,
            "gb_core_strength": 0.95,
            "p2_core_cut": 0.85,
            "gb_strict_core_strength": 1.10,
            "p2_strict_core_cut": 0.90,
        },
        min_pairs=100,
        min_oriented_pairs=60,
    )["recommended"]
    assert low_current == high_current


def test_pipeline_positionals_skip_core_threshold_values() -> None:
    args = [
        "input.dump",
        "--gb-core-strength",
        "0.82",
        "--p2-core-cut",
        "0.71",
        "--gb-strict-core-strength",
        "0.94",
        "--p2-strict-core-cut",
        "0.80",
        "--auto-r-cut-shape-factor",
        "1.8",
        "--p2-cut",
        "0.50",
        "--accepted-edge-audit",
        "--pattern",
        "*.dump",
    ]
    assert pipeline.forwarded_positionals(args) == ["input.dump"]


if __name__ == "__main__":
    test_lobe_split_uses_high_p2_2d_valley()
    test_lobe_split_is_independent_of_current_thresholds()
    test_pipeline_positionals_skip_core_threshold_values()
