"""Build the additive register-resolution package for eight high-yield pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from high_yield_control_validation import FROZEN_TARGETS
from high_yield_register_resolution import (
    REGISTER_CLAIM_BOUNDARY,
    build_experimental_peptide_panel,
    build_register_resolution_gate,
    evaluate_target_windows,
    prioritize_register_confirmation,
    summarize_target_windows,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V3 = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
DEFAULT_HIGH_YIELD = ROOT / "processed/high_yield_control_validation_2026-08-28"
DEFAULT_OUT = ROOT / "processed/high_yield_register_resolution_2026-08-28"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] = (),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or sorted({key for row in rows for key in row}))
    if not fieldnames:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sequence_targets(
    v3_rows: Sequence[Mapping[str, Any]],
    summary_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    v3_by_pair = {str(row["pair_id"]): row for row in v3_rows}
    summary_by_target = {str(row["target_id"]): row for row in summary_rows}
    targets: list[dict[str, Any]] = []
    for frozen in FROZEN_TARGETS:
        if str(frozen["lane"]) != "sequence":
            continue
        target_id = str(frozen["target_id"])
        pair_id = str(frozen["pair_id"])
        if pair_id not in v3_by_pair or target_id not in summary_by_target:
            raise ValueError(f"missing frozen sequence target {target_id}")
        v3 = v3_by_pair[pair_id]
        summary = summary_by_target[target_id]
        if int(summary["target_primary_rank"]) > 3:
            raise ValueError(f"sequence target {target_id} is not top-three in its frozen N3 panel")
        target = {
            **dict(frozen),
            "ebv_sequence": v3["ebv_sequence"],
            "self_sequence": v3["self_sequence"],
            "ebv_core_p1_p9": v3["ebv_core_p1_p9"],
            "self_core_p1_p9": v3["self_core_p1_p9"],
            "ebv_declared_core_start_1_based": v3["ebv_declared_core_start_1_based"],
            "self_declared_core_start_1_based": v3["self_declared_core_start_1_based"],
            "declared_register_status": v3["declared_register_status"],
            "frozen_n3_panel_rank": summary["target_primary_rank"],
            "frozen_n3_panel_status": "rank_context_supportive",
            "register_robust_before_this_analysis": v3["register_robust"],
            "source_v3_primary_rank": v3["primary_rank"],
        }
        targets.append(target)
    targets.sort(key=lambda row: str(row["target_id"]))
    counts = Counter(str(row["allele"]) for row in targets)
    if len(targets) != 8 or set(counts.values()) != {2}:
        raise ValueError("register-resolution scope must contain two targets per HLA and eight total")
    return targets


def _panel_groups(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row["target_id"])].append(dict(row))
    for target_id, panel in groups.items():
        if len(panel) != 26:
            raise ValueError(f"{target_id} does not retain 26 frozen panel rows")
        if sum(str(row["row_role"]) == "target" for row in panel) != 1:
            raise ValueError(f"{target_id} does not retain exactly one target row")
        if sum(str(row["row_role"]) == "n3" for row in panel) != 25:
            raise ValueError(f"{target_id} does not retain 25 N3 rows")
    return dict(groups)


def _validate_locked_local_windows(
    rows: Sequence[Mapping[str, Any]],
    locked_rows: Sequence[Mapping[str, Any]],
) -> None:
    locked = {
        (
            str(row["pair_id"]),
            int(row["ebv_window_start_1_based"]),
            int(row["self_window_start_1_based"]),
        ): row
        for row in locked_rows
    }
    for row in rows:
        if not bool(row["is_local_shift_window_pair"]):
            continue
        key = (
            str(row["pair_id"]),
            int(row["ebv_window_start_1_based"]),
            int(row["self_window_start_1_based"]),
        )
        if key not in locked:
            raise ValueError(f"local sensitivity row is absent from frozen V3: {key}")
        source = locked[key]
        if (
            str(source["ebv_window_core"]) != str(row["ebv_window_core"])
            or str(source["self_window_core"]) != str(row["self_window_core"])
            or abs(
                float(source["tcr_facing_blosum62_similarity"])
                - float(row["tcr_facing_blosum62_similarity"])
            )
            > 1e-10
        ):
            raise ValueError(f"local sensitivity values drifted from frozen V3: {key}")


def _results_markdown(summaries: Sequence[Mapping[str, Any]], gate: Mapping[str, Any]) -> str:
    status_counts = Counter(str(row["register_resolution_status"]) for row in summaries)
    lines = [
        "# Register-resolution results",
        "",
        f"Overall sensitivity status: `{gate['status']}`.",
        "",
        "This package preserves the eight top-three N3 results, then asks how strongly each result depends on the assumed nine-residue HLA-II register. Alternate windows use sequence features only; their structures were not remodeled.",
        "",
        "## Counts",
        "",
        f"- All-window robust: {status_counts['all_window_robust']}/8",
        f"- Local +/-1 robust only: {status_counts['local_shift_robust_only']}/8",
        f"- Declared-window only: {status_counts['declared_window_only']}/8",
        f"- Not supportive at the declared window: {status_counts['rank_context_not_supportive']}/8",
        "",
        "## Candidate results",
        "",
        "| Target | HLA | Declared rank | Worst local rank | Local capture | All-window capture | Status |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in summaries:
        lines.append(
            "| {target_id} | {allele} | {declared_window_rank} | {worst_local_shift_rank} | "
            "{local_shift_capture_at_3_count}/{local_shift_window_pair_count} | "
            "{all_window_capture_at_3_count}/{all_window_pair_count} | {register_resolution_status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A declared-window result can remain promising even when shifts fail, but it is then register-dependent and should be tested with nested peptides before biological interpretation. The 25 N3 pairs retain unknown recognition status and do not become specificity negatives.",
            "",
            REGISTER_CLAIM_BOUNDARY,
            "",
        ]
    )
    return "\n".join(lines)


def _readme() -> str:
    return """# High-yield register resolution (2026-08-28)

This additive package follows the high-yield N3 analysis. It does not change V1-V3 rankings, replace the prior 12-target package, freeze weights, or unlock discovery.

## Question

For each of the eight sequence-supported candidates, does its top-three rank survive plausible or exhaustive changes to the assumed P1-P9 register when compared with the same frozen 25 N3 pairs?

## Method

- Recompute TCR-facing BLOSUM62, physicochemical mismatch, five-position identity, full-core identity, and full-core BLOSUM62 for every fully contained 9-mer combination.
- Rank each alternative against the candidate's unchanged 25-pair N3 panel with the frozen V3 lexicographic sequence order.
- Abstain from structural tie-breaking for every alternate-register row because alternate registers were not modeled.
- Report both the local +/-1-by-+/-1 neighborhood and the exhaustive window set.
- Provide a proposed nested-peptide design, marked proposed and not ordered.

HLA-II ligands commonly occur as nested peptides of variable length, so a parent peptide alone does not prove its binding register. Experimental work has used nested register peptides and binding measurements to separate alternative MHC-II registers. See [Chicz et al., 1992](https://pubmed.ncbi.nlm.nih.gov/1380674/) and [Mohan et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3256971/).

## Files

- `protocol_lock.json`: frozen scope, inputs, ranking logic, and claim gates.
- `frozen_sequence_target_registry.csv`: exact eight candidates and source sequences.
- `all_window_panel_ranks.csv`: exhaustive register-pair sensitivity matrix.
- `local_shift_panel_ranks.csv`: declared +/-1 window matrix.
- `target_register_summary.csv`: candidate-level robustness and worst ranks.
- `experimental_register_priority.csv`: assay order based on local-shift robustness; not a discovery rerank.
- `experimental_peptide_panel.csv`: proposed parent and nested 9-mer sequences; nothing has been ordered.
- `register_resolution_gate.json`: machine-readable status and permanent lock flags.
- `RESULTS_SUMMARY.md`: concise result table.
- `SHA256SUMS.csv`: deterministic artifact checksums.

## Claim boundary

""" + REGISTER_CLAIM_BOUNDARY + "\n"


def _experimental_next_steps() -> str:
    return """# Experimental next step

1. Start with the parent peptide and the declared, -1, and +1 nested 9-mer cores in `experimental_peptide_panel.csv` for each arm.
2. Measure binding to the exact HLA allele first. Compare affinity and complex stability across the nested cores; do not infer TCR recognition from binding.
3. Prioritize candidates whose declared core binds and whose local alternatives separate clearly. A flat result means the register remains unresolved.
4. Confirm the displayed register structurally or with an allele-appropriate anchor-perturbation design before calling the geometry register-resolved.
5. T-cell cross-recognition requires a shared paired TCR or clone and matched functional assays for both ligands. This project does not currently supply that evidence for the discovery candidates.

The proposed 9-mers are register-discrimination reagents, not guaranteed optimal HLA-II assay peptides. Parent peptides are included as reference ligands because natural HLA-II ligands commonly have flanking residues and variable lengths.

Nothing in this file is an order, submission, specificity claim, or biological validation result.
"""


def build_package(
    *,
    v3_dir: Path = DEFAULT_V3,
    high_yield_dir: Path = DEFAULT_HIGH_YIELD,
    output_dir: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    v3_dir = Path(v3_dir)
    high_yield_dir = Path(high_yield_dir)
    output_dir = Path(output_dir)
    v3_ranked_path = v3_dir / "v3_all_hla_ranked_pairs.csv"
    v3_local_path = v3_dir / "register_sensitivity_allowed_windows_v2.csv"
    high_yield_summary_path = high_yield_dir / "panel_rank_summary.csv"
    high_yield_matrix_path = high_yield_dir / "panel_feature_matrix.csv"
    input_paths = [v3_ranked_path, v3_local_path, high_yield_summary_path, high_yield_matrix_path]
    for path in input_paths:
        if not path.exists():
            raise FileNotFoundError(path)

    targets = _sequence_targets(read_csv(v3_ranked_path), read_csv(high_yield_summary_path))
    panels = _panel_groups(read_csv(high_yield_matrix_path))
    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for target in targets:
        target_id = str(target["target_id"])
        if target_id not in panels:
            raise ValueError(f"missing frozen N3 panel for {target_id}")
        rows = evaluate_target_windows(target, panels[target_id], local_shift=1)
        all_rows.extend(rows)
        summaries.append(summarize_target_windows(rows))
    _validate_locked_local_windows(all_rows, read_csv(v3_local_path))
    local_rows = [row for row in all_rows if bool(row["is_local_shift_window_pair"])]
    if len(local_rows) != 72:
        raise ValueError("eight targets must produce exactly 72 local +/-1 window pairs")
    peptide_rows = build_experimental_peptide_panel(targets, local_shift=1)
    gate = build_register_resolution_gate(summaries)
    priority_rows = prioritize_register_confirmation(summaries)

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    protocol = {
        "package_id": "high_yield_register_resolution_2026-08-28",
        "protocol_date": "2026-08-28",
        "scope": "eight_sequence_supported_high_yield_pairs_only",
        "target_ids": [row["target_id"] for row in targets],
        "input_checksums": {path.name: sha256_file(path) for path in input_paths},
        "frozen_n3_pairs_per_target": 25,
        "window_width": 9,
        "local_shift_definition": "declared start plus or minus one independently on both arms",
        "exhaustive_definition": "all fully contained P1-P9 windows on both exact source peptides",
        "primary_rank_order": [
            "tcr_facing_blosum62_similarity_desc",
            "tcr_face_physicochemical_mismatch_asc",
            "tcr_facing_sequence_identity_desc",
            "structural_tie_break_abstained",
            "pair_id_asc",
        ],
        "geometry_read_for_alternate_registers": False,
        "alternate_registers_structurally_modeled": False,
        "n3_specificity_role": "unknown_recognition_not_a_specificity_negative",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "experimental_peptides_ordered": False,
        "claim_boundary": REGISTER_CLAIM_BOUNDARY,
    }
    write_json(output_dir / "protocol_lock.json", protocol)
    write_csv(output_dir / "frozen_sequence_target_registry.csv", targets)
    write_csv(output_dir / "all_window_panel_ranks.csv", all_rows)
    write_csv(output_dir / "local_shift_panel_ranks.csv", local_rows)
    write_csv(output_dir / "target_register_summary.csv", summaries)
    write_csv(output_dir / "experimental_register_priority.csv", priority_rows)
    write_csv(output_dir / "experimental_peptide_panel.csv", peptide_rows)
    write_json(output_dir / "register_resolution_gate.json", gate)
    (output_dir / "README.md").write_text(_readme(), encoding="utf-8")
    (output_dir / "RESULTS_SUMMARY.md").write_text(
        _results_markdown(summaries, gate), encoding="utf-8"
    )
    (output_dir / "EXPERIMENTAL_NEXT_STEP.md").write_text(
        _experimental_next_steps(), encoding="utf-8"
    )

    artifacts = sorted(path for path in output_dir.iterdir() if path.is_file())
    artifact_checksums = {path.name: sha256_file(path) for path in artifacts}
    manifest = {
        "package_id": protocol["package_id"],
        "target_count": len(targets),
        "n3_pairs_per_target": 25,
        "all_window_row_count": len(all_rows),
        "local_window_row_count": len(local_rows),
        "experimental_peptide_row_count": len(peptide_rows),
        "experimental_priority_row_count": len(priority_rows),
        "gate_status": gate["status"],
        "artifact_checksums": artifact_checksums,
        "deterministic_rebuild_required": True,
    }
    write_json(output_dir / "analysis_manifest.json", manifest)
    checksum_paths = sorted(path for path in output_dir.iterdir() if path.is_file())
    checksum_rows = [
        {"file": path.name, "sha256": sha256_file(path)} for path in checksum_paths
    ]
    write_csv(output_dir / "SHA256SUMS.csv", checksum_rows, ("file", "sha256"))
    return {
        **manifest,
        "file_checksums": {row["file"]: row["sha256"] for row in checksum_rows},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v3-dir", type=Path, default=DEFAULT_V3)
    parser.add_argument("--high-yield-dir", type=Path, default=DEFAULT_HIGH_YIELD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    manifest = build_package(
        v3_dir=args.v3_dir,
        high_yield_dir=args.high_yield_dir,
        output_dir=args.output_dir,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
