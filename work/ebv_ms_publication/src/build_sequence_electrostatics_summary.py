"""Combine the immutable pilot and expansion electrostatics target summaries."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
PILOT_DIR = ROOT / "processed/pmhc_surface_electrostatics_pilot_2026-08-29"
EXPANSION_DIR = ROOT / "processed/pmhc_surface_electrostatics_sequence_expansion_2026-08-30"
DEFAULT_OUT = ROOT / "processed/pmhc_surface_electrostatics_all_sequence_candidates_2026-08-30"
EXPECTED_TARGET_IDS = (
    "HY03_SEQ_01",
    "HY03_SEQ_02",
    "HY08_SEQ_01",
    "HY08_SEQ_02",
    "HY13_SEQ_01",
    "HY13_SEQ_02",
    "HY15_SEQ_01",
    "HY15_SEQ_02",
)
CLAIM_BOUNDARY = (
    "Descriptive, model-derived HLA-specific pMHC local-field resemblance only; "
    "not evidence of presentation, TCR recognition, activation, specificity, "
    "cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate."
)


def _truth(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def classify_combined_support(
    primary_rank: int,
    *,
    register_robust: bool,
    dielectric_robust: bool,
) -> str:
    if primary_rank > 3:
        return "sequence_supported_electrostatics_not_supportive"
    if not dielectric_robust:
        return "sequence_supported_electrostatics_dielectric_unstable"
    if not register_robust:
        return "sequence_plus_electrostatics_rank_supported_register_unresolved"
    return "sequence_plus_electrostatics_supported"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_summary(
    output_dir: Path = DEFAULT_OUT,
    source_dirs: Sequence[Path] = (PILOT_DIR, EXPANSION_DIR),
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite existing package: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    source_manifest = []
    for source in source_dirs:
        summary = source / "target_electrostatic_summary.csv"
        checksums = source / "SHA256SUMS.csv"
        if not summary.exists() or not checksums.exists():
            raise FileNotFoundError(f"incomplete electrostatics source package: {source}")
        source_manifest.append(
            {
                "source_package": source.name,
                "target_summary_sha256": _sha256(summary),
                "checksums_sha256": _sha256(checksums),
            }
        )
        for row in _read_csv(summary):
            primary_rank = int(row["primary_full_pmhc_rank"])
            register_robust = _truth(row["register_robust"])
            dielectric_robust = _truth(row["dielectric_rank_class_robust"])
            rows.append(
                {
                    **row,
                    "source_package": source.name,
                    "sequence_support_status": "sequence_supported_upstream",
                    "combined_support_status": classify_combined_support(
                        primary_rank,
                        register_robust=register_robust,
                        dielectric_robust=dielectric_robust,
                    ),
                    "claim_boundary": CLAIM_BOUNDARY,
                }
            )
    rows.sort(key=lambda row: (row["allele"], int(row["primary_full_pmhc_rank"]), row["target_id"]))
    observed = tuple(sorted(row["target_id"] for row in rows))
    if observed != tuple(sorted(EXPECTED_TARGET_IDS)):
        raise ValueError(f"expected exactly the eight sequence targets, observed {observed}")
    _write_csv(output_dir / "all_sequence_candidate_electrostatics.csv", rows)
    _write_csv(output_dir / "source_manifest.csv", source_manifest)

    formal = [row for row in rows if row["combined_support_status"] == "sequence_plus_electrostatics_supported"]
    rank_only = [
        row for row in rows
        if row["combined_support_status"] == "sequence_plus_electrostatics_rank_supported_register_unresolved"
    ]
    result = {
        "candidate_count": len(rows),
        "formal_sequence_plus_electrostatics_supported_count": len(formal),
        "rank_context_supported_register_unresolved_count": len(rank_only),
        "electrostatics_not_supportive_or_unstable_count": len(rows) - len(formal) - len(rank_only),
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (output_dir / "support_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    table = "\n".join(
        f"| {row['target_id']} | {row['allele']} | `{row['ebv_core']}` / `{row['self_core']}` | "
        f"{row['primary_full_pmhc_rank']}/26 | {row['full_rank_eps4']}/26 / {row['full_rank_eps8']}/26 | "
        f"{row['hla_normalized_rank_eps2']}/26 | `{row['combined_support_status']}` |"
        for row in rows
    )
    readme = f"""# All Sequence-Supported Candidate Electrostatics

This package joins the immutable two-candidate pilot with the six-candidate expansion. It does not rerank the discovery universe or create a composite score.

| Target | HLA | EBV / self core | Full-pMHC rank eps-in 2 | eps-in 4 / 8 | HLA-subtracted rank | Combined status |
|---|---|---|---:|---:|---:|---|
{table}

Formal sequence-plus-electrostatics support requires a full-pMHC rank of 1-3, stable rank classification across solute dielectrics 2/4/8, and a robust frozen register. Rank-context support with unresolved register is reported separately and is not promoted to formal support.

N3 comparisons have unknown recognition status and are not specificity negatives. {CLAIM_BOUNDARY}
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    checksum_rows = []
    for path in sorted(value for value in output_dir.rglob("*") if value.is_file() and value.name != "SHA256SUMS.csv"):
        checksum_rows.append(
            {
                "relative_path": path.relative_to(output_dir).as_posix(),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
        )
    _write_csv(output_dir / "SHA256SUMS.csv", checksum_rows)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(build_summary(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
