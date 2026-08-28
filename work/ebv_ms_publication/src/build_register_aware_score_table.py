"""Build auditable register-filtered score tables from frozen benchmark inputs."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any

from register_aware_scoring import score_same_register_alignment


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
BENCHMARK = PROC / "register_aware_benchmark"
OUT = PROC / "register_aware_scoring"
GEOMETRY_CONTEXT_STATUS = "whole_original_local_alignment_only"
GEOMETRY_CONTEXT_NOTE = (
    "Precomputed RMSD spans the original local alignment; no per-residue model "
    "coordinates are available to recompute a same-register-only RMSD."
)
EXPECTED_ASSESSABLE_COUNTS = {"1": 11, "2": 10, "3": 9, "4": 3}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def pair_id(row: dict[str, Any]) -> str:
    return f"{row['ebv_candidate_id']}::{row['human_candidate_id']}"


def build_score_rows(
    universe_rows: list[dict[str, Any]], geometry_rows: list[dict[str, Any]]
) -> list[dict[str, object]]:
    """Join each frozen universe pair to its labelled geometry context."""
    geometry_by_pair: dict[str, dict[str, Any]] = {}
    for geometry_row in geometry_rows:
        identifier = pair_id(geometry_row)
        if identifier in geometry_by_pair:
            raise ValueError(f"duplicate geometry row for {identifier}")
        geometry_by_pair[identifier] = geometry_row

    universe_ids = {str(row["pair_id"]) for row in universe_rows}
    missing = universe_ids - set(geometry_by_pair)
    extras = set(geometry_by_pair) - universe_ids
    if missing or extras:
        raise ValueError(
            "universe and geometry pair IDs differ: "
            f"missing_geometry={sorted(missing)!r}, extra_geometry={sorted(extras)!r}"
        )

    scored: list[dict[str, object]] = []
    for universe_row in universe_rows:
        identifier = str(universe_row["pair_id"])
        geometry_row = geometry_by_pair[identifier]
        score = score_same_register_alignment(universe_row)
        scored.append({
            **universe_row,
            **score,
            "whole_local_alignment_geometry_rmsd_context": geometry_row.get(
                "local_peptide_ca_rmsd_after_hla_fit", ""
            ),
            "geometry_context_status": GEOMETRY_CONTEXT_STATUS,
            "geometry_primary_score_eligible": False,
            "geometry_context_note": GEOMETRY_CONTEXT_NOTE,
        })
    return scored


def assert_expected_assessable_coverage(rows: list[dict[str, object]]) -> None:
    """Reject silent drift in the frozen register-eligibility universe."""
    assessable = [row for row in rows if row["register_assessment"] == "assessable_register_hypothesis"]
    observed = Counter(str(row["same_register_alignment_count"]) for row in assessable)
    if len(assessable) != 33 or dict(observed) != EXPECTED_ASSESSABLE_COUNTS:
        raise ValueError(
            "unexpected frozen assessable coverage: "
            f"n={len(assessable)}, counts={dict(sorted(observed.items()))}"
        )


def coverage_summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter(
        (
            str(row["score_coverage_status"]),
            str(row["register_assessment"]),
            str(row.get("same_register_alignment_count", "")),
        )
        for row in rows
    )
    return [
        {
            "score_coverage_status": status,
            "register_assessment": assessment,
            "same_register_alignment_count": same_count,
            "pair_count": count,
        }
        for (status, assessment, same_count), count in sorted(counts.items())
    ]


def rendered_readme() -> str:
    return "\n".join([
        "# Register-aware scoring outputs",
        "",
        "`register_aware_pair_scores.csv` contains every frozen PASS-universe pair.",
        "Sequence and physicochemical descriptors are recalculated only for existing",
        "local-alignment coordinates that occupy the same P1-P9 register index in both",
        "resolved peptide cores. P1/P4/P6/P9 are labelled HLA-anchor positions and",
        "P2/P3/P5/P7/P8 candidate-exposed positions; these labels do not assert TCR",
        "contacts.",
        "",
        "Rows with one or two retained coordinates are `limited_coverage_report_only`.",
        "Only rows with at least three retained coordinates are eligible for the robust",
        "primary sequence/chemistry ranking. The score is a descriptor, not a molecular",
        "mimicry probability or evidence of presentation, shared-TCR binding,",
        "activation, cross-reactivity, or MS mechanism.",
        "",
        "`whole_local_alignment_geometry_rmsd_context` is included strictly as context:",
        "it was precomputed across the original local alignment. It is not a",
        "same-register structural score and cannot enter the primary score until the",
        "original per-residue pMHC coordinates are recovered and re-analysed.",
        "",
        "A positive computational prioritization result requires an independently",
        "interpretable target, complete frozen strict decoys, robust score coverage,",
        "and persistence under retained register-sensitivity analysis. Otherwise the",
        "result is negative/mixed: the current evidence does not support robust",
        "register-aware prioritization. This is not evidence against all EBV-myelin biology.",
        "",
        "Regenerate:",
        "",
        "```bash",
        "PYTHONPATH=src python3 src/build_register_aware_score_table.py",
        "PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v",
        "```",
        "",
    ])


def main() -> None:
    universe = read_csv(BENCHMARK / "benchmark_pair_universe.csv")
    geometry = read_csv(PROC / "colabfold_tier1_ebv_myelin_geometry_matrix.csv")
    assert_expected_assessable_coverage(universe)
    score_rows = build_score_rows(universe, geometry)
    write_csv(OUT / "register_aware_pair_scores.csv", score_rows)
    write_csv(OUT / "score_coverage_summary.csv", coverage_summary_rows(score_rows))
    write_csv(OUT / "register_aware_geometry_limitations.csv", [{
        "geometry_context_status": GEOMETRY_CONTEXT_STATUS,
        "geometry_primary_score_eligible": False,
        "current_limitation": GEOMETRY_CONTEXT_NOTE,
        "required_coordinate_inputs": "Original pMHC PDB files for both members of every scored pair.",
        "recovery_condition": (
            "Recover original pMHC PDBs, HLA-fit with the prior procedure, and "
            "recompute RMSD using only retained same-register residue pairs."
        ),
    }])
    (OUT / "README.md").write_text(rendered_readme(), encoding="utf-8")
    print(f"wrote {len(score_rows)} score rows to {OUT}")


if __name__ == "__main__":
    main()
