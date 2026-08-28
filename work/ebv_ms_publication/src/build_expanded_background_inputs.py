"""Prepare a source-traceable human-background pMHC expansion batch.

The batch is created before prediction scores, structural models, pair geometry,
or register-aware scores are inspected. It therefore expands the eligible
pre-score comparator universe without changing the strict decoy rule.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from premeeting_rigor import iedb_mhcii_eligible


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "expanded_background"
SOURCE = PROC / "human_drb1501_mhc_ii_iedb_enriched.csv"
EXISTING_BATCH = PROC / "pmhc_colabfold_batch.csv"
BENCHMARK = PROC / "register_aware_benchmark"
COMPARATOR_ARM = "Human background comparator"
COMPARATOR_TIER = "Tier 4: source-validated HLA-DRB1*15:01 comparator record"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def build_background_candidates(
    source_rows: list[dict[str, str]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Select direct, source-validated MHC-II comparator inputs conservatively."""
    registry: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    selected_peptides: set[str] = set()
    for source in sorted(source_rows, key=lambda row: (str(row.get("iedb_epitope_id", "")), str(row.get("peptide", "")))):
        peptide = str(source.get("peptide", "")).upper()
        identifier = str(source.get("iedb_epitope_id", ""))
        candidate_id = f"HUMAN_BACKGROUND_{identifier}" if identifier else ""
        base = {
            "candidate_id": candidate_id,
            "iedb_epitope_id": identifier,
            "peptide": peptide,
            "peptide_length": len(peptide),
            "source_antigen_name": source.get("source_antigen_name", ""),
            "mhc_allele": source.get("mhc_allele", ""),
            "candidate_class": source.get("candidate_class", ""),
            "provenance_status": source.get("provenance_status", ""),
        }
        if source.get("candidate_class") != "human_background":
            registry.append({**base, "selection_status": "excluded_not_human_background", "selection_reason": "Not labelled human_background in the frozen enriched source table."})
        elif source.get("provenance_status") != "coordinate_validated":
            registry.append({**base, "selection_status": "excluded_unvalidated_coordinates", "selection_reason": "Human source coordinates were not validated in the frozen source table."})
        elif not iedb_mhcii_eligible(peptide):
            registry.append({**base, "selection_status": "retained_not_modeled_missing_verified_flank", "selection_reason": "Peptide is outside the 11-30 aa direct IEDB MHC-II range and no verified natural flank is supplied."})
        elif peptide in selected_peptides:
            registry.append({**base, "selection_status": "retained_not_modeled_duplicate_peptide", "selection_reason": "Duplicate comparator peptide; only the first stable source record enters the pre-score batch."})
        else:
            selected_peptides.add(peptide)
            registry.append({**base, "selection_status": "selected_for_direct_iedb_and_pmhc_batch", "selection_reason": "Source-validated human-background peptide in the direct IEDB MHC-II length range."})
            candidates.append({
                "candidate_id": candidate_id,
                "arm": COMPARATOR_ARM,
                "evidence_tier": COMPARATOR_TIER,
                "peptide": peptide,
                "peptide_length": len(peptide),
                "source_antigen": source.get("source_antigen_name", ""),
                "source_accession": "see human_epitope_accession_map.csv",
                "iedb_assay_id": "",
                "iedb_epitope_id": identifier,
                "pubmed_id": "",
                "hla": source.get("mhc_allele", ""),
                "mhc_class": source.get("mhc_class", "II"),
                "modeling_status": "prepared_pre_score_background_batch",
            })
    return registry, candidates


def build_colabfold_batch_rows(
    candidates: list[dict[str, object]], dra_sequence: str, drb_sequence: str
) -> list[dict[str, str]]:
    """Return exact three-chain pMHC queries with untrimmed comparator peptides."""
    if not dra_sequence or not drb_sequence:
        raise ValueError("existing ColabFold batch does not provide both HLA chains")
    return [
        {
            "id": str(candidate["candidate_id"]),
            "sequence": f"{dra_sequence}:{drb_sequence}:{candidate['peptide']}",
        }
        for candidate in candidates
    ]


def potential_target_coverage(
    feasibility_rows: list[dict[str, str]],
    universe_rows: list[dict[str, str]],
    candidates: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Report length-only potential before scores, registers, or models are known."""
    universe_by_pair = {row["pair_id"]: row for row in universe_rows}
    output: list[dict[str, object]] = []
    for target in sorted(feasibility_rows, key=lambda row: row["target_pair_id"]):
        if target.get("target_register_assessment") != "assessable_register_hypothesis":
            continue
        pair = universe_by_pair.get(target["target_pair_id"])
        if pair is None:
            raise ValueError(f"missing target pair in benchmark universe: {target['target_pair_id']}")
        target_length = len(str(pair["human_peptide"]))
        matching = [
            candidate for candidate in candidates
            if abs(len(str(candidate["peptide"])) - target_length) <= 1
        ]
        output.append({
            "target_pair_id": target["target_pair_id"],
            "target_readiness_before_expansion": target.get("readiness_status", ""),
            "target_ebv_peptide_length": len(str(pair["ebv_peptide"])),
            "target_human_peptide_length": target_length,
            "direct_background_length_match_count": len(matching),
            "direct_background_candidate_ids": ";".join(
                str(candidate["candidate_id"]) for candidate in matching
            ),
            "coverage_interpretation": (
                "length_feasible_pending_iedb_and_structure" if matching
                else "no_direct_background_length_match_in_current_source_set"
            ),
            "claim_boundary": (
                "Length feasibility only; binding-bin matching, register eligibility, "
                "and pMHC geometry must be established before a comparator can become a decoy."
            ),
        })
    return output


def hla_chains_from_existing_batch(batch_rows: list[dict[str, str]]) -> tuple[str, str]:
    """Reuse the exact DRA/DRB chains already used for the frozen screen."""
    if not batch_rows:
        raise ValueError("existing ColabFold batch is empty")
    fragments = batch_rows[0]["sequence"].split(":")
    if len(fragments) != 3:
        raise ValueError("existing ColabFold sequence is not a three-chain pMHC query")
    return fragments[0], fragments[1]


def render_readme(registry: list[dict[str, object]], candidates: list[dict[str, object]]) -> str:
    retained = sum(row["selection_status"] == "retained_not_modeled_missing_verified_flank" for row in registry)
    background_records = sum(row.get("candidate_class") == "human_background" for row in registry)
    return "\n".join([
        "# Expanded human-background comparator inputs",
        "",
        f"- Source-table records scanned: **{len(registry)}**",
        f"- Human-background records reviewed: **{background_records}**",
        f"- Direct IEDB MHC-II and pMHC batch inputs: **{len(candidates)}**",
        f"- Retained but not modeled for missing verified natural flanks: **{retained}**",
        "",
        "The selected records are a pre-score human-background comparator arm. They",
        "were chosen from source/provenance and direct IEDB MHC-II eligibility only,",
        "before binding predictions, structure results, register scores, or candidate",
        "priority values were inspected. Full IEDB peptide sequences are preserved.",
        "",
        "This batch can expand strict-decoy feasibility only after all three gates are",
        "met: matching IEDB binding-rank bins, same-register eligibility, and passed",
        "pMHC structure QA. It does not create a biological negative control or",
        "evidence of molecular mimicry.",
        "",
        "Reproduce:",
        "",
        "```bash",
        "PYTHONPATH=src python3 src/build_expanded_background_inputs.py",
        "PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v",
        "```",
        "",
    ])


def main() -> None:
    source_rows = read_csv(SOURCE)
    registry, candidates = build_background_candidates(source_rows)
    existing_batch = read_csv(EXISTING_BATCH)
    dra, drb = hla_chains_from_existing_batch(existing_batch)
    batch_rows = build_colabfold_batch_rows(candidates, dra, drb)
    coverage = potential_target_coverage(
        read_csv(BENCHMARK / "target_feasibility.csv"),
        read_csv(BENCHMARK / "benchmark_pair_universe.csv"),
        candidates,
    )
    registry_fields = [
        "candidate_id", "iedb_epitope_id", "peptide", "peptide_length",
        "source_antigen_name", "mhc_allele", "candidate_class", "provenance_status",
        "selection_status", "selection_reason",
    ]
    candidate_fields = [
        "candidate_id", "arm", "evidence_tier", "peptide", "peptide_length",
        "source_antigen", "source_accession", "iedb_assay_id", "iedb_epitope_id",
        "pubmed_id", "hla", "mhc_class", "modeling_status",
    ]
    write_csv(OUT / "background_candidate_registry.csv", registry, registry_fields)
    write_csv(OUT / "background_pmhc_candidate_manifest.csv", candidates, candidate_fields)
    write_csv(OUT / "background_pmhc_colabfold_batch.csv", batch_rows, ["id", "sequence"])
    write_csv(OUT / "background_iedb_submission_manifest.csv", [
        {
            "submission_id": f"{candidate['candidate_id']}__segment_001",
            "candidate_id": candidate["candidate_id"],
            "peptide": candidate["peptide"],
            "source_start_1_based": 1,
            "submission_strategy": "direct_full_peptide",
            "claim_boundary": "Computational binding/register hypothesis only; not experimental presentation evidence.",
        }
        for candidate in candidates
    ], [
        "submission_id", "candidate_id", "peptide", "source_start_1_based",
        "submission_strategy", "claim_boundary",
    ])
    write_csv(OUT / "target_length_coverage_audit.csv", coverage, [
        "target_pair_id", "target_readiness_before_expansion", "target_ebv_peptide_length",
        "target_human_peptide_length", "direct_background_length_match_count",
        "direct_background_candidate_ids", "coverage_interpretation", "claim_boundary",
    ])
    fasta = []
    for candidate, batch in zip(candidates, batch_rows):
        fasta.extend([
            f">{candidate['candidate_id']}|arm=Human_background_comparator|tier=Tier_4|iedb_epitope_id={candidate['iedb_epitope_id']}|peptide={candidate['peptide']}",
            batch["sequence"],
        ])
    (OUT / "background_pmhc_inputs.fasta").write_text("\n".join(fasta) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(render_readme(registry, candidates), encoding="utf-8")
    print(f"prepared {len(candidates)} direct IEDB/pMHC background inputs in {OUT}")


if __name__ == "__main__":
    main()
