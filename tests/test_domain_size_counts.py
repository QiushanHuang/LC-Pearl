from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import lc_domain_size_counts as counts  # noqa: E402


def test_domain_size_counts_group_by_frame_size_and_classification(tmp_path: Path) -> None:
    table = tmp_path / "domain_diagnostics.tsv"
    table.write_text(
        "\t".join(["source_file", "timestep", "classification", "size", "s2"]) + "\n"
        + "a.dump\t0\tweak\t2\t0.80\n"
        + "a.dump\t0\tweak\t2\t0.82\n"
        + "a.dump\t0\trobust\t3\t0.91\n"
        + "b.dump\t1\tweak\t2\t0.72\n"
        + "b.dump\t1\trobust\t3\t0.88\n"
        + "b.dump\t1\trobust\t5\t0.94\n"
        + "c.dump\t2\ttransient\t5\t0.51\n",
        encoding="utf-8",
    )

    rows = counts.collect_domain_size_frame_counts(table)

    assert rows == [
        counts.DomainSizeFrameCount(source_file="a.dump", timestep=0, size=2, total_count=2, weak_count=2, robust_count=0, other_count=0),
        counts.DomainSizeFrameCount(source_file="a.dump", timestep=0, size=3, total_count=1, weak_count=0, robust_count=1, other_count=0),
        counts.DomainSizeFrameCount(source_file="b.dump", timestep=1, size=2, total_count=1, weak_count=1, robust_count=0, other_count=0),
        counts.DomainSizeFrameCount(source_file="b.dump", timestep=1, size=3, total_count=1, weak_count=0, robust_count=1, other_count=0),
        counts.DomainSizeFrameCount(source_file="b.dump", timestep=1, size=5, total_count=1, weak_count=0, robust_count=1, other_count=0),
        counts.DomainSizeFrameCount(source_file="c.dump", timestep=2, size=5, total_count=1, weak_count=0, robust_count=0, other_count=1),
    ]
    assert counts.collect_domain_size_frame_counts_from_records(
        [
            {"source_file": "a.dump", "timestep": 0, "classification": "weak", "size": 2},
            {"source_file": "a.dump", "timestep": 0, "classification": "weak", "size": 2},
            {"source_file": "a.dump", "timestep": 0, "classification": "robust", "size": 3},
            {"source_file": "b.dump", "timestep": 1, "classification": "weak", "size": 2},
            {"source_file": "b.dump", "timestep": 1, "classification": "robust", "size": 3},
            {"source_file": "b.dump", "timestep": 1, "classification": "robust", "size": 5},
            {"source_file": "c.dump", "timestep": 2, "classification": "transient", "size": 5},
        ]
    ) == rows


if __name__ == "__main__":
    tmp = Path("/tmp/lc_domain_size_counts_test")
    tmp.mkdir(parents=True, exist_ok=True)
    test_domain_size_counts_group_by_frame_size_and_classification(tmp)
