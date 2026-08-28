"""Build the additive held-out human HLA-II positive-control benchmark package."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import urllib.request
from typing import Any, Dict, List, Mapping, Sequence, Tuple

from build_tcell_library_v2 import (
    DRA_SEQUENCE,
    DRB1_1501_SEQUENCE,
    DRB5_0101_SEQUENCE,
    fetch_binding_predictions,
)
from hla2_positive_control_benchmark import (
    CLAIM_BOUNDARY,
    build_af3_job_batches,
    build_pdb_oracle_pairings,
    build_trust_gate,
    generate_weight_grid,
    geometry_from_mmcif,
    pair_features,
    rank_feature_percentiles,
    select_score_blind_comparators,
    validate_comparator_registry,
    validate_control_registry,
    validate_af3_job_package,
)


ROOT = Path(__file__).resolve().parents[1]
V2_PACKAGE = ROOT / "processed/tcell_library_v2_2026-08-22"
MODEL_ANALYSIS = ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25"
DEFAULT_OUT = ROOT / "processed/hla2_positive_control_benchmark_2026-08-25"
PANEL_SEEDS = (104729, 104759)
STRUCTURAL_SEARCH_SNAPSHOT_DATE = "2026-08-26"
FEATURES = (
    "exposed_ca_rmsd_A",
    "exposed_sidechain_vector_rmsd_A",
    "tcr_face_physicochemical_mismatch",
    "anchor_ca_rmsd_A",
)
DQA1_0102_SEQUENCE = (
    "EDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEQFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRY"
    "NSTAATNEVPEVTVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSA"
    "DEIYDCKVEHWGLDQPLLKHWEPEIPAPMSELTE"
)
DQB1_0502_SEQUENCE = (
    "EGRDSPEDFVYQFKGLCYFTNGTERVRGVTRHIYNREEYVRFDSDVGVYRAVTPQGRPVAEYWNSQKEVLEGARASVDRV"
    "CRHNYEVAYRGILQRRVEPTVTISPSRTEALNHHNLLICSVTDFYPSQIKVRWFRNDQEETAGVVSTPLIRNGDWTFQIL"
    "VMLEMTPQRGDVYTCHVEHPSLQSPITVEWRAQSESAQSKVD"
)


def curated_registry() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    systems = [
        {
            "system_id": "SYS_BALF5_MBP_HY2E11", "tcr_id": "Hy.2E11",
            "eligibility": "strict", "independent_system_weight": 1,
            "same_paired_tcr_or_clone": True, "functional_evidence_for_both_arms": True,
            "evidence_summary": "Same human clone recognizes DRB5-BALF5 and DRB1-MBP; both pMHC structures resolved.",
            "primary_source": "PMID:12244309", "doi": "10.1038/ni835",
        },
        {
            "system_id": "SYS_ENGA_MBP_OB1A12", "tcr_id": "Ob.1A12",
            "eligibility": "strict", "independent_system_weight": 1,
            "same_paired_tcr_or_clone": True, "functional_evidence_for_both_arms": True,
            "evidence_summary": "Same human TCR recognizes DR15-MBP and processed bacterial EngA; both ternary structures resolved.",
            "primary_source": "PMID:19303388", "doi": "10.1016/j.immuni.2009.01.009",
        },
        {
            "system_id": "SYS_MICROBIAL_MBP_HY1B11", "tcr_id": "Hy.1B11",
            "eligibility": "strict", "independent_system_weight": 1,
            "same_paired_tcr_or_clone": True, "functional_evidence_for_both_arms": True,
            "evidence_summary": "Same human TCR recognizes DQ1-MBP, HSV UL15, and Pseudomonas PMM; all three ternary structures resolved.",
            "primary_source": "Nature Communications 4:2623", "doi": "10.1038/ncomms3623",
        },
        {
            "system_id": "PROSPECTIVE_EBNA1_ANO2_2026", "tcr_id": "study clone set",
            "eligibility": "prospective", "independent_system_weight": 0,
            "same_paired_tcr_or_clone": True, "functional_evidence_for_both_arms": True,
            "evidence_summary": "Human cross-reactive clone evidence, but exact paired peptides, registers, and structures are unresolved.",
            "primary_source": "PMID:41534529", "doi": "10.1016/j.cell.2025.12.032",
        },
    ]
    ligands = [
        {
            "ligand_id": "HY2E11_BALF5", "system_id": "SYS_BALF5_MBP_HY2E11", "ligand_role": "microbial_positive",
            "source_protein": "EBV BALF5", "sequence": "TGGVYHFVKKHVHES", "core": "YHFVKKHVH",
            "core_start_1_based": 5, "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB5*01:01",
            "pdb_id": "1H15", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "C",
        },
        {
            "ligand_id": "HY2E11_MBP", "system_id": "SYS_BALF5_MBP_HY2E11", "ligand_role": "self_positive",
            "source_protein": "human MBP", "sequence": "ENPVVHFFKNIVTPR", "core": "VHFFKNIVT",
            "core_start_1_based": 5, "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01",
            "pdb_id": "1BX2", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "C",
        },
        {
            "ligand_id": "OB1A12_ENGA", "system_id": "SYS_ENGA_MBP_OB1A12", "ligand_role": "microbial_positive",
            "source_protein": "bacterial EngA", "sequence": "DFARVHFISALHGSG", "core": "VHFISALHG",
            "core_start_1_based": 5, "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01",
            "pdb_id": "2WBJ", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "D",
        },
        {
            "ligand_id": "OB1A12_MBP", "system_id": "SYS_ENGA_MBP_OB1A12", "ligand_role": "self_positive",
            "source_protein": "human MBP", "sequence": "ENPVVHFFKNIVTPR", "core": "VHFFKNIVT",
            "core_start_1_based": 5, "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01",
            "pdb_id": "1YMM", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "C",
        },
        {
            "ligand_id": "HY1B11_UL15", "system_id": "SYS_MICROBIAL_MBP_HY1B11", "ligand_role": "microbial_positive",
            "source_protein": "HSV UL15", "sequence": "QLVHFVRDFAQL", "core": "VHFVRDFAQ",
            "core_start_1_based": 3, "mhc_alpha_allele": "HLA-DQA1*01:02", "mhc_beta_allele": "HLA-DQB1*05:02",
            "pdb_id": "4MAY", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "D",
        },
        {
            "ligand_id": "HY1B11_PMM", "system_id": "SYS_MICROBIAL_MBP_HY1B11", "ligand_role": "microbial_positive",
            "source_protein": "Pseudomonas PMM", "sequence": "RLLMLFAKDVVSRN", "core": "MLFAKDVVS",
            "core_start_1_based": 4, "mhc_alpha_allele": "HLA-DQA1*01:02", "mhc_beta_allele": "HLA-DQB1*05:02",
            "pdb_id": "4GRL", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "D",
        },
        {
            "ligand_id": "HY1B11_MBP", "system_id": "SYS_MICROBIAL_MBP_HY1B11", "ligand_role": "self_positive",
            "source_protein": "human MBP", "sequence": "ENPVVHFFKNIVTPR", "core": "VHFFKNIVT",
            "core_start_1_based": 5, "mhc_alpha_allele": "HLA-DQA1*01:02", "mhc_beta_allele": "HLA-DQB1*05:02",
            "pdb_id": "3PL6", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "D",
        },
    ]
    pairs = [
        {"pair_id": "PAIR_HY2E11_BALF5_MBP", "system_id": "SYS_BALF5_MBP_HY2E11", "left_ligand_id": "HY2E11_BALF5", "right_ligand_id": "HY2E11_MBP", "required_for_system_pass": True},
        {"pair_id": "PAIR_OB1A12_ENGA_MBP", "system_id": "SYS_ENGA_MBP_OB1A12", "left_ligand_id": "OB1A12_ENGA", "right_ligand_id": "OB1A12_MBP", "required_for_system_pass": True},
        {"pair_id": "PAIR_HY1B11_UL15_MBP", "system_id": "SYS_MICROBIAL_MBP_HY1B11", "left_ligand_id": "HY1B11_UL15", "right_ligand_id": "HY1B11_MBP", "required_for_system_pass": True},
        {"pair_id": "PAIR_HY1B11_PMM_MBP", "system_id": "SYS_MICROBIAL_MBP_HY1B11", "left_ligand_id": "HY1B11_PMM", "right_ligand_id": "HY1B11_MBP", "required_for_system_pass": True},
    ]
    sources = [
        {"source_id": "HY2_PRIMARY", "identifier": "PMID:12244309", "url": "https://doi.org/10.1038/ni835", "role": "Hy.2E11 function and structural mimicry"},
        {"source_id": "OB1_PRIMARY", "identifier": "PMID:19303388", "url": "https://doi.org/10.1016/j.immuni.2009.01.009", "role": "Ob.1A12 EngA cross-recognition and structure"},
        {"source_id": "HY1_PRIMARY", "identifier": "DOI:10.1038/ncomms3623", "url": "https://www.nature.com/articles/ncomms3623", "role": "Hy.1B11 MBP, UL15, and PMM function and structures"},
        {"source_id": "ANO2_PRIMARY", "identifier": "PMID:41534529", "url": "https://pubmed.ncbi.nlm.nih.gov/41534529/", "role": "prospective EBNA1-ANO2 evidence"},
        {"source_id": "TCR3D", "identifier": "TCR3d 2.0", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11701517/", "role": "systematic human TCR-pMHC structure discovery"},
        {"source_id": "NETTCR_LOEO", "identifier": "NetTCR-2.2", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10942633/", "role": "leave-one-epitope-out leakage control"},
    ]
    for pdb_id in ("1H15", "1BX2", "2WBJ", "1YMM", "3PL6", "4MAY", "4GRL"):
        sources.append({
            "source_id": f"PDB_{pdb_id}", "identifier": f"PDB:{pdb_id}",
            "url": f"https://www.rcsb.org/structure/{pdb_id}", "role": "experimental positive-control structure",
        })
    return systems, ligands, pairs, sources


def build_evaluation_skeleton(
    positive_pairs: Sequence[Mapping[str, Any]], *, panel_seeds: Sequence[int]
) -> List[Dict[str, Any]]:
    rows = []
    for pair in positive_pairs:
        base = {
            "system_id": pair["system_id"],
            "pair_id": pair["pair_id"],
            "required_for_system_pass": pair["required_for_system_pass"],
            "positive_rank": "",
            "comparison_count": 0,
            "capture_at_3": "",
            "evaluation_status": "pending",
            "claim_boundary": CLAIM_BOUNDARY,
        }
        rows.append({**base, "layer": "pdb_oracle", "panel_seed": "pdb"})
        for seed in panel_seeds:
            rows.append({**base, "layer": "af3", "panel_seed": int(seed)})
    return rows


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or sorted({key for row in rows for key in row}))
    if not fieldnames:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
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


def _tcr3d_discovery_audit(out: Path) -> List[Dict[str, Any]]:
    source_url = "https://tcr3d.ibbr.umd.edu/static/download/tcr_complexes_data.tsv"
    source_path = out / "sources/tcr3d/tcr_complexes_data.tsv"
    source_path.parent.mkdir(parents=True, exist_ok=True)
    if not source_path.exists() or source_path.stat().st_size < 1000:
        urllib.request.urlretrieve(source_url, source_path)
    rows = read_csv_tsv(source_path)
    grouped = defaultdict(list)
    for row in rows:
        if row["TCR_complex"] != "CLASSII" or row["TCR_organism"] != "Human":
            continue
        key = (row["CDR3_alpha"], row["CDR3_beta"])
        if all(key):
            grouped[key].append(row)
    exclusion_reasons = {
        "E8": "designed single-substitution peptide variants; not an independent biological cross-source system",
        "MS2-3C8": "same nine-residue ligand core with different peptide flanks",
        "F11": "single C-terminal peptide variant and exact HLA alpha/beta allotypes not resolved in the snapshot",
        "JR5.1": "future candidate; exact registers and both-arm functional source evidence not yet locked in v1",
        "XPA5": "future candidate; exact HLA alpha/beta allotypes and both-arm functional source evidence not yet locked in v1",
        "A2.13": "future candidate; exact HLA alpha/beta allotypes and both-arm functional source evidence not yet locked in v1",
        "ET650-4": "future candidate; exact HLA alpha/beta allotypes and both-arm functional source evidence not yet locked in v1",
        "G9": "future candidate; exact HLA alpha/beta allotypes and both-arm functional source evidence not yet locked in v1",
    }
    output = []
    for (cdr3_alpha, cdr3_beta), values in sorted(grouped.items()):
        epitopes = sorted({row["Epitope"] for row in values})
        if len(epitopes) < 2:
            continue
        names = sorted({row["TCR_name"] for row in values})
        display_name = "Hy.1B11" if "Hy.1B11" in names else (
            "Ob.1A12" if any(name in {"Ob", "Ob.1A12"} for name in names) else names[0]
        )
        admitted = display_name in {"Hy.1B11", "Ob.1A12"}
        output.append({
            "screen_id": hashlib.sha256(f"{cdr3_alpha}|{cdr3_beta}".encode("utf-8")).hexdigest()[:12],
            "tcr_names": ";".join(names),
            "cdr3_alpha": cdr3_alpha,
            "cdr3_beta": cdr3_beta,
            "distinct_ligand_count": len(epitopes),
            "epitopes": ";".join(epitopes),
            "mhc_labels": ";".join(sorted({row["MHC_allele"] for row in values})),
            "pdb_ids": ";".join(sorted({row["PDB_ID"].upper() for row in values})),
            "pubmed_ids": ";".join(sorted({row["Pubmed"] for row in values if row["Pubmed"] not in {"", "null"}})),
            "screen_status": "admitted_strict_v1" if admitted else "not_admitted_v1",
            "screen_reason": (
                "same paired human TCR, exact structures, registers, allotypes, and both-arm evidence locked"
                if admitted else exclusion_reasons.get(
                    display_name,
                    "future candidate; exact allotypes, registers, and both-arm evidence require manual curation",
                )
            ),
            "snapshot_date": STRUCTURAL_SEARCH_SNAPSHOT_DATE,
            "source_file": "sources/tcr3d/tcr_complexes_data.tsv",
        })
    output.append({
        "screen_id": "literature_hy2e11",
        "tcr_names": "Hy.2E11",
        "cdr3_alpha": "not_structure_indexed",
        "cdr3_beta": "not_structure_indexed",
        "distinct_ligand_count": 2,
        "epitopes": "BALF5;MBP",
        "mhc_labels": "HLA-DR2a;HLA-DR2b",
        "pdb_ids": "1BX2;1H15",
        "pubmed_ids": "12244309",
        "screen_status": "admitted_strict_v1_literature_clone",
        "screen_reason": "same explicitly identified human clone; both exact pMHC structures and functional arms locked",
        "snapshot_date": STRUCTURAL_SEARCH_SNAPSHOT_DATE,
        "source_file": "registry/literature_and_structure_sources.csv",
    })
    write_json(out / "sources/tcr3d/tcr3d_snapshot_metadata.json", {
        "source_url": source_url,
        "snapshot_date": STRUCTURAL_SEARCH_SNAPSHOT_DATE,
        "source_sha256": sha256_file(source_path),
        "all_complex_row_count": len(rows),
        "human_class_ii_multiligand_tcr_group_count": len(output) - 1,
        "grouping_rule": "exact paired CDR3 alpha and beta sequence; at least two distinct epitope strings",
    })
    return sorted(output, key=lambda row: (row["screen_status"], row["tcr_names"]))


def read_csv_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _prediction_cache(
    out: Path, allele: str, rows: Sequence[Mapping[str, Any]], slug: str
) -> List[Dict[str, Any]]:
    normalized_path = out / f"raw_responses/{slug}_normalized.csv"
    raw_path = out / f"raw_responses/{slug}_iedb.tsv"
    expected = [(str(row["candidate_id"]), str(row["sequence"])) for row in rows]
    if normalized_path.exists():
        cached = read_csv(normalized_path)
        observed = [(row["candidate_id"], row["sequence"]) for row in cached]
        if observed == expected and all(row["allele"] == allele for row in cached):
            return [dict(row) for row in cached]
    predictions, raw = fetch_binding_predictions(allele, [dict(row) for row in rows])
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text(raw, encoding="utf-8")
    for row in predictions:
        row["raw_response_file"] = str(raw_path.relative_to(out))
    write_csv(normalized_path, predictions)
    return predictions


def _panel_metadata() -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    panel = read_csv(V2_PACKAGE / "frozen_v2_80_peptide_panel.csv")
    return (
        [row for row in panel if row["kingdom"] == "EBV"],
        [row for row in panel if row["kingdom"] == "human_self"],
    )


def _native_human_controls() -> List[Dict[str, Any]]:
    rows = read_csv(V2_PACKAGE / "frozen_native_hla_controls.csv")
    return [
        {
            "candidate_id": row["control_candidate_id"],
            "sequence": row["control_sequence"],
            "source_protein": row["control_source"],
            "source_accession": row["control_accession"],
            "binding_percentile": "",
        }
        for row in rows if row["arm"] == "self"
    ]


def _select_comparator_rows(
    out: Path, ligands: Sequence[Mapping[str, Any]]
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    ligand_by_id = {str(row["ligand_id"]): row for row in ligands}
    ebv_panel, _ = _panel_metadata()
    human_controls = _native_human_controls()
    dr15_predictions = {
        row["candidate_id"]: row
        for row in read_csv(V2_PACKAGE / "allele_register_predictions_320.csv")
        if row["allele"] == "HLA-DRB1*15:01"
    }
    dr_human_predictions = {
        row["candidate_id"]: row
        for row in read_csv(V2_PACKAGE / "calibration_control_binding_predictions.csv")
        if row["allele"] == "HLA-DRB1*15:01"
    }
    inherited_raw_dir = out / "sources/inherited_v2_raw_responses"
    inherited_raw_dir.mkdir(parents=True, exist_ok=True)
    inherited_raw_files = {
        "foreign": "iedb_panel80_hla_drb1_15_01.tsv",
        "self": "iedb_calibration_human_non_cns_drb1_1501.tsv",
    }
    for filename in inherited_raw_files.values():
        shutil.copy2(V2_PACKAGE / "raw_responses" / filename, inherited_raw_dir / filename)

    def provenance(prediction: Mapping[str, Any], *, raw_response_file: str = "") -> Dict[str, Any]:
        return {
            "predicted_core": prediction["predicted_core"],
            "core_start_1_based": int(prediction["core_start"]),
            "register_resolution": prediction["register_resolution"],
            "seq_num": int(prediction["seq_num"]),
            "raw_response_file": raw_response_file or prediction["raw_response_file"],
            "prediction_method": prediction["prediction_method"],
            "prediction_status": prediction["prediction_status"],
        }

    for row in human_controls:
        prediction = dr_human_predictions[row["candidate_id"]]
        row["binding_percentile"] = prediction["percentile_rank"]
        row.update(provenance(
            prediction,
            raw_response_file=f"sources/inherited_v2_raw_responses/{inherited_raw_files['self']}",
        ))

    enga = ligand_by_id["OB1A12_ENGA"]
    enga_prediction = _prediction_cache(
        out, "HLA-DRB1*15:01",
        [{"candidate_id": "OB1A12_ENGA", "sequence": enga["sequence"]}],
        "ob1a12_enga_dr15",
    )[0]
    dr_foreign_pool = [
        {
            "candidate_id": row["candidate_id"], "sequence": row["sequence"],
            "binding_percentile": dr15_predictions[row["candidate_id"]]["percentile_rank"],
            "source_protein": row["protein_symbol"], "source_accession": row["accession"],
        }
        for row in ebv_panel
        if row["candidate_id"] in dr15_predictions
        and dr15_predictions[row["candidate_id"]]["register_resolution"] == "resolved_unique_fully_contained"
    ]
    for row in dr_foreign_pool:
        row.update(provenance(
            dr15_predictions[row["candidate_id"]],
            raw_response_file=f"sources/inherited_v2_raw_responses/{inherited_raw_files['foreign']}",
        ))
    ob_foreign = select_score_blind_comparators(
        {"sequence": enga["sequence"], "binding_percentile": enga_prediction["percentile_rank"]},
        dr_foreign_pool, count=5, seed=104759,
    )
    dr_source = {row["candidate_id"]: row for row in dr_foreign_pool}
    for row in ob_foreign:
        row.update({
            "system_id": "SYS_ENGA_MBP_OB1A12", "positive_pair_id": "PAIR_OB1A12_ENGA_MBP",
            "comparator_arm": "microbial", "mhc_alpha_allele": "HLA-DRA*01:01",
            "mhc_beta_allele": "HLA-DRB1*15:01", "negative_tier": "N3",
            "recognition_status": "unknown_not_specificity_negative",
            "source_protein": dr_source[row["candidate_id"]]["source_protein"],
            "source_accession": dr_source[row["candidate_id"]]["source_accession"],
            **provenance(
                dr15_predictions[row["candidate_id"]],
                raw_response_file=f"sources/inherited_v2_raw_responses/{inherited_raw_files['foreign']}",
            ),
        })
    ob_self = []
    for selected in select_score_blind_comparators(
        {"sequence": ligand_by_id["OB1A12_MBP"]["sequence"], "binding_percentile": 0.08},
        human_controls, count=5, seed=104759,
    ):
        source = next(row for row in human_controls if row["candidate_id"] == selected["candidate_id"])
        selected.update({
            "system_id": "SYS_ENGA_MBP_OB1A12", "positive_pair_id": "PAIR_OB1A12_ENGA_MBP",
            "comparator_arm": "self", "mhc_alpha_allele": "HLA-DRA*01:01",
            "mhc_beta_allele": "HLA-DRB1*15:01", "negative_tier": "N3",
            "recognition_status": "unknown_not_specificity_negative",
            "source_protein": source["source_protein"], "source_accession": source["source_accession"],
            **{key: source[key] for key in (
                "predicted_core", "core_start_1_based", "register_resolution", "seq_num",
                "raw_response_file", "prediction_method", "prediction_status",
            )},
        })
        ob_self.append(selected)

    dq_targets = [ligand_by_id[identifier] for identifier in ("HY1B11_UL15", "HY1B11_PMM", "HY1B11_MBP")]
    dq_prediction_inputs = [
        {"candidate_id": row["candidate_id"], "sequence": row["sequence"]}
        for row in ebv_panel
    ] + [
        {"candidate_id": row["candidate_id"], "sequence": row["sequence"]}
        for row in human_controls
    ] + [
        {"candidate_id": row["ligand_id"], "sequence": row["sequence"]}
        for row in dq_targets
    ]
    dq_predictions = _prediction_cache(
        out, "HLA-DQA1*01:02/DQB1*05:02", dq_prediction_inputs, "hy1b11_dq1",
    )
    dq_by_id = {row["candidate_id"]: row for row in dq_predictions}
    dq_foreign_pool = [
        {
            "candidate_id": row["candidate_id"], "sequence": row["sequence"],
            "binding_percentile": dq_by_id[row["candidate_id"]]["percentile_rank"],
            "source_protein": row["protein_symbol"], "source_accession": row["accession"],
        }
        for row in ebv_panel
        if dq_by_id[row["candidate_id"]]["register_resolution"] == "resolved_unique_fully_contained"
    ]
    for row in dq_foreign_pool:
        row.update(provenance(dq_by_id[row["candidate_id"]]))
    dq_human_pool = [
        {
            **row,
            "binding_percentile": dq_by_id[row["candidate_id"]]["percentile_rank"],
            **provenance(dq_by_id[row["candidate_id"]]),
        }
        for row in human_controls
        if dq_by_id[row["candidate_id"]]["register_resolution"] == "resolved_unique_fully_contained"
    ]
    hy_rows = []
    for target_id, pair_id in (
        ("HY1B11_UL15", "PAIR_HY1B11_UL15_MBP"),
        ("HY1B11_PMM", "PAIR_HY1B11_PMM_MBP"),
    ):
        target = ligand_by_id[target_id]
        for selected in select_score_blind_comparators(
            {"sequence": target["sequence"], "binding_percentile": dq_by_id[target_id]["percentile_rank"]},
            dq_foreign_pool, count=5, seed=104759,
        ):
            source = next(row for row in dq_foreign_pool if row["candidate_id"] == selected["candidate_id"])
            selected.update({
                "system_id": "SYS_MICROBIAL_MBP_HY1B11", "positive_pair_id": pair_id,
                "comparator_arm": "microbial", "mhc_alpha_allele": "HLA-DQA1*01:02",
                "mhc_beta_allele": "HLA-DQB1*05:02", "negative_tier": "N3",
                "recognition_status": "unknown_not_specificity_negative",
                "source_protein": source["source_protein"], "source_accession": source["source_accession"],
                **{key: source[key] for key in (
                    "predicted_core", "core_start_1_based", "register_resolution", "seq_num",
                    "raw_response_file", "prediction_method", "prediction_status",
                )},
            })
            hy_rows.append(selected)
    for pair_id in ("PAIR_HY1B11_UL15_MBP", "PAIR_HY1B11_PMM_MBP"):
        for selected in select_score_blind_comparators(
            {"sequence": ligand_by_id["HY1B11_MBP"]["sequence"], "binding_percentile": dq_by_id["HY1B11_MBP"]["percentile_rank"]},
            dq_human_pool, count=5, seed=104759,
        ):
            source = next(row for row in dq_human_pool if row["candidate_id"] == selected["candidate_id"])
            selected.update({
                "system_id": "SYS_MICROBIAL_MBP_HY1B11", "positive_pair_id": pair_id,
                "comparator_arm": "self", "mhc_alpha_allele": "HLA-DQA1*01:02",
                "mhc_beta_allele": "HLA-DQB1*05:02", "negative_tier": "N3",
                "recognition_status": "unknown_not_specificity_negative",
                "source_protein": source["source_protein"], "source_accession": source["source_accession"],
                **{key: source[key] for key in (
                    "predicted_core", "core_start_1_based", "register_resolution", "seq_num",
                    "raw_response_file", "prediction_method", "prediction_status",
                )},
            })
            hy_rows.append(selected)
    return ob_foreign + ob_self + hy_rows, [enga_prediction], dq_predictions


def _hla_sequences() -> Dict[Tuple[str, str], Dict[str, str]]:
    return {
        ("HLA-DRA*01:01", "HLA-DRB1*15:01"): {
            "mhc_alpha_sequence": DRA_SEQUENCE, "mhc_beta_sequence": DRB1_1501_SEQUENCE,
        },
        ("HLA-DRA*01:01", "HLA-DRB5*01:01"): {
            "mhc_alpha_sequence": DRA_SEQUENCE, "mhc_beta_sequence": DRB5_0101_SEQUENCE,
        },
        ("HLA-DQA1*01:02", "HLA-DQB1*05:02"): {
            "mhc_alpha_sequence": DQA1_0102_SEQUENCE, "mhc_beta_sequence": DQB1_0502_SEQUENCE,
        },
    }


def _new_af3_ligands(
    ligands: Sequence[Mapping[str, Any]], comparators: Sequence[Mapping[str, Any]]
) -> List[Dict[str, Any]]:
    ligand_by_id = {str(row["ligand_id"]): row for row in ligands}
    rows = [dict(ligand_by_id["OB1A12_ENGA"])]
    rows.extend(dict(ligand_by_id[identifier]) for identifier in ("HY1B11_UL15", "HY1B11_PMM", "HY1B11_MBP"))
    for row in rows:
        row.update({
            "core_sequence": row["core"],
            "register_resolution": "experimentally_resolved",
            "register_source": f"experimental_structure_PDB_{row['pdb_id']}",
        })
    seen = {(row["mhc_alpha_allele"], row["mhc_beta_allele"], row["sequence"]) for row in rows}
    for comparator in comparators:
        # DR15 self arms already have complete/pending seed-specific jobs in the original frozen calibration.
        if comparator["system_id"] == "SYS_ENGA_MBP_OB1A12" and comparator["comparator_arm"] == "self":
            continue
        key = (comparator["mhc_alpha_allele"], comparator["mhc_beta_allele"], comparator["sequence"])
        if key in seen:
            continue
        prefix = "OB1_DR15" if comparator["system_id"] == "SYS_ENGA_MBP_OB1A12" else "HY1_DQ1"
        rows.append({
            "ligand_id": f"{prefix}_{comparator['candidate_id']}",
            "system_id": comparator["system_id"], "ligand_role": "N3_comparator",
            "sequence": comparator["sequence"], "mhc_alpha_allele": comparator["mhc_alpha_allele"],
            "mhc_beta_allele": comparator["mhc_beta_allele"],
            "core_sequence": comparator["predicted_core"],
            "core_start_1_based": comparator["core_start_1_based"],
            "register_resolution": comparator["register_resolution"],
            "register_source": "IEDB_binding_prediction",
            "seq_num": comparator["seq_num"],
            "raw_response_file": comparator["raw_response_file"],
        })
        seen.add(key)
    return rows


def _comparison_universe(
    pairs: Sequence[Mapping[str, Any]],
    ligands: Sequence[Mapping[str, Any]],
    comparators: Sequence[Mapping[str, Any]],
) -> List[Dict[str, Any]]:
    ligand_by_id = {str(row["ligand_id"]): row for row in ligands}
    output = []
    for pair in pairs:
        if pair["system_id"] == "SYS_BALF5_MBP_HY2E11":
            continue
        left = ligand_by_id[str(pair["left_ligand_id"])]
        right = ligand_by_id[str(pair["right_ligand_id"])]
        left_controls = [row for row in comparators if row["positive_pair_id"] == pair["pair_id"] and row["comparator_arm"] == "microbial"]
        right_controls = [row for row in comparators if row["positive_pair_id"] == pair["pair_id"] and row["comparator_arm"] == "self"]
        if len(left_controls) != 5 or len(right_controls) != 5:
            raise ValueError(f"comparison universe requires five controls per arm for {pair['pair_id']}")
        for seed in PANEL_SEEDS:
            output.append({
                "system_id": pair["system_id"], "positive_pair_id": pair["pair_id"], "panel_seed": seed,
                "comparison_role": "positive", "left_id": left["ligand_id"], "right_id": right["ligand_id"],
                "pair_id": f"{pair['pair_id']}|s{seed}|positive", "negative_tier": "positive",
                "geometry_status": "pending_models", "claim_boundary": CLAIM_BOUNDARY,
            })
            for left_control in left_controls:
                for right_control in right_controls:
                    output.append({
                        "system_id": pair["system_id"], "positive_pair_id": pair["pair_id"], "panel_seed": seed,
                        "comparison_role": "N3", "left_id": left_control["candidate_id"],
                        "right_id": right_control["candidate_id"],
                        "pair_id": f"{pair['pair_id']}|s{seed}|{left_control['candidate_id']}|{right_control['candidate_id']}",
                        "negative_tier": "N3_unknown_recognition_ranking_only",
                        "geometry_status": "pending_models", "claim_boundary": CLAIM_BOUNDARY,
                    })
    return output


def _download_pdb_sources(out: Path, ligands: Sequence[Mapping[str, Any]]) -> Dict[str, Path]:
    sources = {}
    for ligand in ligands:
        pdb_id = str(ligand["pdb_id"])
        path = out / f"sources/pdb/{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists() or path.stat().st_size < 1000:
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.cif", path)
        sources[pdb_id] = path
    return sources


def _mhc2_structure_snapshot(out: Path) -> List[Dict[str, Any]]:
    path = out / "sources/tcr3d/mhc2_structures_snapshot.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    with urllib.request.urlopen("https://tcr3d.ibbr.umd.edu/mhc2_chains", timeout=120) as response:
        page = response.read().decode("utf-8")
    match = re.search(r"var data = (\[.*?\]);\s*\n", page, re.DOTALL)
    if not match:
        raise ValueError("TCR3d MHC-II structural table was not found")
    rows = json.loads(match.group(1))
    write_json(path, rows)
    return rows


def _structural_ligand_registry(
    out: Path, ligands: Sequence[Mapping[str, Any]]
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    snapshot = _mhc2_structure_snapshot(out)
    snapshot_by_pdb = {str(row["pdbid"]).upper(): row for row in snapshot}
    structural = []
    for ligand in ligands:
        source = snapshot_by_pdb[str(ligand["pdb_id"])]
        structural.append({
            **dict(ligand),
            "core_sequence": ligand["core"],
            "peptide_sequence": ligand["sequence"],
            "structural_role": "declared_positive",
            "resolution_A": source["resolution"],
            "register_source": f"experimental_structure_PDB_{ligand['pdb_id']}",
            "selected_for_oracle_pool": True,
        })
    structural.extend([
        {
            "ligand_id": "PDB_DRB5_MBP_ALT_1FV1", "system_id": "PDB_DECOY_LIBRARY",
            "ligand_role": "noncognate_structural_ligand", "source_protein": "human MBP",
            "sequence": "NPVVHFFKNIVTPRTPPPSQ", "peptide_sequence": "NPVVHFFKNIVTPRTPPPSQ",
            "core": "FKNIVTPRT", "core_sequence": "FKNIVTPRT", "core_start_1_based": 7,
            "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB5*01:01",
            "pdb_id": "1FV1", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "C",
            "structural_role": "exact_hla_decoy_ligand", "resolution_A": 1.9,
            "register_source": "experimental_structure_PDB_1FV1", "selected_for_oracle_pool": True,
        },
        {
            "ligand_id": "PDB_DR15_ALPHA3_5V4M", "system_id": "PDB_DECOY_LIBRARY",
            "ligand_role": "noncognate_structural_ligand", "source_protein": "human alpha-3 chain",
            "sequence": "GWISLWKGFSF", "peptide_sequence": "GWISLWKGFSF",
            "core": "ISLWKGFSF", "core_sequence": "ISLWKGFSF", "core_start_1_based": 3,
            "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01",
            "pdb_id": "5V4M", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "B",
            "structural_role": "exact_hla_decoy_ligand", "resolution_A": 2.1,
            "register_source": "experimental_structure_PDB_5V4M", "selected_for_oracle_pool": True,
        },
        {
            "ligand_id": "PDB_DR15_GDP_6CPO", "system_id": "PDB_DECOY_LIBRARY",
            "ligand_role": "noncognate_structural_ligand", "source_protein": "human GDP-L-fucose synthase",
            "sequence": "RFYKTLRAEQASQ", "peptide_sequence": "RFYKTLRAEQASQ",
            "core": "YKTLRAEQA", "core_sequence": "YKTLRAEQA", "core_start_1_based": 3,
            "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01",
            "pdb_id": "6CPO", "mhc_alpha_chain": "A", "mhc_beta_chain": "B", "peptide_chain": "C",
            "structural_role": "exact_hla_decoy_ligand", "resolution_A": 2.4,
            "register_source": "experimental_structure_PDB_6CPO", "selected_for_oracle_pool": True,
        },
    ])
    relevant_labels = {"HLA-DR2a", "HLA-DR2b", "HLA-DR15", "HLA-DQ1"}
    selected_pdbs = {str(row["pdb_id"]) for row in structural}
    exclusions = {
        "1HQR": "incomplete exact nine-residue peptide core in coordinates",
        "1ZGL": "technical pMHC duplicate; 1FV1 selected by better resolution before geometry",
        "6CQQ": "technical pMHC duplicate; 6CPO selected by better resolution before geometry",
        "8VSP": "HLA-DQ1 shorthand but alpha/beta sequences do not match DQA1*01:02/DQB1*05:02",
    }
    screen = []
    for row in snapshot:
        pdb_id = str(row["pdbid"]).upper()
        if row["mhc_allele_name"] not in relevant_labels:
            continue
        screen.append({
            "pdb_id": pdb_id,
            "tcr3d_mhc_label": row["mhc_allele_name"],
            "peptide": row["peptide"],
            "core_peptide": row["core_peptide"],
            "resolution_A": row["resolution"],
            "receptor_bound": bool(row["is_bound"]),
            "oracle_pool_status": "selected" if pdb_id in selected_pdbs else "excluded",
            "selection_reason": (
                "declared positive or unique exact-HLA ligand selected before geometry"
                if pdb_id in selected_pdbs else exclusions.get(
                    pdb_id,
                    "technical duplicate or exact alpha/beta allotype not established for the locked pool",
                )
            ),
            "snapshot_date": STRUCTURAL_SEARCH_SNAPSHOT_DATE,
        })
    return structural, screen


def _experimental_pdb_oracle(
    out: Path, ligands: Sequence[Mapping[str, Any]], pairs: Sequence[Mapping[str, Any]]
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    structural, screen = _structural_ligand_registry(out, ligands)
    paths = _download_pdb_sources(out, structural)
    hla = _hla_sequences()
    geometries = {}
    for ligand in structural:
        alpha = str(ligand["mhc_alpha_allele"])
        beta = str(ligand["mhc_beta_allele"])
        starts = (2, 4) if alpha.startswith("HLA-DQA") else (0, 1)
        sequences = hla[(alpha, beta)]
        geometries[ligand["ligand_id"]] = geometry_from_mmcif(
            paths[str(ligand["pdb_id"])], ligand_id=str(ligand["ligand_id"]),
            peptide_sequence=str(ligand["sequence"]), core_sequence=str(ligand["core"]),
            mhc_alpha_chain=str(ligand["mhc_alpha_chain"]), mhc_beta_chain=str(ligand["mhc_beta_chain"]),
            peptide_chain=str(ligand["peptide_chain"]),
            mhc_alpha_reference_sequence=sequences["mhc_alpha_sequence"],
            mhc_beta_reference_sequence=sequences["mhc_beta_sequence"],
            mhc_alpha_reference_start=starts[0], mhc_beta_reference_start=starts[1],
        )
    matrix = []
    summaries = []
    ligand_by_id = {str(row["ligand_id"]): row for row in structural}
    for pair in pairs:
        left_id, right_id = str(pair["left_ligand_id"]), str(pair["right_ligand_id"])
        left, right = ligand_by_id[left_id], ligand_by_id[right_id]
        needed_hlas = {
            (left["mhc_alpha_allele"], left["mhc_beta_allele"]),
            (right["mhc_alpha_allele"], right["mhc_beta_allele"]),
        }
        eligible = [
            row for row in structural
            if (row["mhc_alpha_allele"], row["mhc_beta_allele"]) in needed_hlas
        ]
        deduplicated = []
        by_equivalence = defaultdict(list)
        for row in eligible:
            by_equivalence[(
                row["mhc_alpha_allele"], row["mhc_beta_allele"], row["core_sequence"]
            )].append(row)
        for values in by_equivalence.values():
            values = sorted(values, key=lambda row: (
                0 if row["ligand_id"] in {left_id, right_id} else 1,
                float(row["resolution_A"]), str(row["pdb_id"]), str(row["ligand_id"]),
            ))
            deduplicated.append(values[0])
        pairing_rows, summary = build_pdb_oracle_pairings({
            **dict(pair),
            "left_mhc_alpha_allele": left["mhc_alpha_allele"],
            "left_mhc_beta_allele": left["mhc_beta_allele"],
            "right_mhc_alpha_allele": right["mhc_alpha_allele"],
            "right_mhc_beta_allele": right["mhc_beta_allele"],
        }, deduplicated)
        panel_rows = []
        for pairing in pairing_rows:
            pair_left = ligand_by_id[str(pairing["left_ligand_id"])]
            pair_right = ligand_by_id[str(pairing["right_ligand_id"])]
            panel_rows.append({
                **pairing,
                "layer": "pdb_oracle", "panel_seed": "pdb",
                "left_pdb_id": pair_left["pdb_id"], "right_pdb_id": pair_right["pdb_id"],
                **{
                    key: round(value, 9)
                    for key, value in pair_features(
                        geometries[str(pairing["left_ligand_id"])],
                        geometries[str(pairing["right_ligand_id"])],
                    ).items()
                },
                "decoy_count": summary["decoy_count"],
                "evaluation_status": summary["evaluation_status"],
                "ranking_endpoint": "frozen_exposed_ca_rmsd_A",
                "claim_boundary": CLAIM_BOUNDARY,
            })
        panel_rows = rank_feature_percentiles(panel_rows, FEATURES)
        ordered = sorted(panel_rows, key=lambda row: (
            float(row["exposed_ca_rmsd_A"]), str(row["pair_id"])
        ))
        for rank, row in enumerate(ordered, start=1):
            row["exposed_ca_rank"] = rank
        positive = next(row for row in panel_rows if row["pair_role"] == "positive")
        positive_rank = int(positive["exposed_ca_rank"])
        positive.update({
            "available_positive_rank": positive_rank,
            "positive_rank": positive_rank if summary["evaluation_status"] == "complete" else "",
            "capture_at_3": positive_rank <= 3 if summary["evaluation_status"] == "complete" else "",
        })
        summaries.append({
            **summary,
            "system_id": pair["system_id"],
            "positive_pair_id": pair["pair_id"],
            "available_positive_rank": positive_rank,
            "positive_rank": positive_rank if summary["evaluation_status"] == "complete" else "",
            "capture_at_3": positive_rank <= 3 if summary["evaluation_status"] == "complete" else "",
            **{feature: positive[feature] for feature in FEATURES},
        })
        matrix.extend(panel_rows)
    positives = [row for row in matrix if row["pair_role"] == "positive"]
    return structural, screen, matrix, summaries


def _overlay_evaluation_status(
    skeleton: List[Dict[str, Any]], pdb_summaries: Sequence[Mapping[str, Any]]
) -> List[Dict[str, Any]]:
    pdb_by_pair = {str(row["positive_pair_id"]): row for row in pdb_summaries}
    seed_recovery = {
        int(row["seed"]): row
        for row in read_csv(MODEL_ANALYSIS / "validation/gold_standard_seed_recovery.csv")
    }
    output = []
    for row in skeleton:
        value = dict(row)
        if row["layer"] == "pdb_oracle":
            source = pdb_by_pair[str(row["pair_id"])]
            value.update({
                "evaluation_status": source["evaluation_status"],
                "comparison_count": source["comparison_count"],
                "positive_rank": source["positive_rank"],
                "available_positive_rank": source["available_positive_rank"],
                "capture_at_3": source["capture_at_3"],
                "decoy_count": source["decoy_count"],
                "ranking_endpoint": "frozen_exposed_ca_rmsd_A",
                **{feature: source[feature] for feature in FEATURES},
            })
        elif row["system_id"] == "SYS_BALF5_MBP_HY2E11":
            seed = int(row["panel_seed"])
            source = seed_recovery[seed]
            formal = str(source["formal_seed_evaluable"]).lower() == "true"
            value.update({
                "evaluation_status": "complete" if formal else "missing_required_comparisons",
                "comparison_count": int(source["available_primary_count"]),
                "positive_rank": int(source["available_rank"]) if formal else "",
                "available_positive_rank": int(source["available_rank"]),
                "capture_at_3": int(source["available_rank"]) <= 3 if formal else "",
            })
        else:
            value.update({
                "evaluation_status": "pending_af3_models", "comparison_count": 0,
                "available_positive_rank": "",
            })
        output.append(value)
    return output


def _copy_frozen_hy2_context(out: Path) -> None:
    target = out / "frozen_hy2e11_context"
    target.mkdir(parents=True, exist_ok=True)
    for source in (
        V2_PACKAGE / "frozen_native_hla_controls.csv",
        V2_PACKAGE / "native_hla_calibration_manifest_24.csv",
        V2_PACKAGE / "calibration_comparison_universe_72.csv",
        MODEL_ANALYSIS / "validation/gold_standard_seed_recovery.csv",
        MODEL_ANALYSIS / "validation/gold_standard_capture_summary.json",
    ):
        shutil.copy2(source, target / source.name)


def _write_documentation(
    out: Path, validation: Mapping[str, Any], trust_gate: Mapping[str, Any], job_count: int,
    batch_sizes: Sequence[int], positive_features: Sequence[Mapping[str, Any]],
    pdb_summaries: Sequence[Mapping[str, Any]],
) -> None:
    feature_lines = "\n".join(
        f"- {row['positive_pair_id']}: exposed CA {float(row['exposed_ca_rmsd_A']):.3f} A; "
        f"side-chain vector {float(row['exposed_sidechain_vector_rmsd_A']):.3f} A"
        for row in positive_features
    )
    pdb_lines = "\n".join(
        f"- {row['positive_pair_id']}: {row['decoy_count']} exact-HLA decoys; "
        f"available exposed-CA rank {row['available_positive_rank']}; {row['evaluation_status']}"
        for row in pdb_summaries
    )
    readme = f"""# Held-out human HLA-II positive-control benchmark

This additive package expands the strict denominator from one to **{validation['strict_independent_system_count']} independent human TCR systems** and **{validation['strict_positive_pair_count']} required positive comparisons**. The original discovery rankings are unchanged.

## Current status

- Overall trust status: **{trust_gate['overall_trust_status']}**.
- Discovery reranking allowed: **{str(trust_gate['discovery_reranking_allowed']).lower()}**.
- New AlphaFold jobs: **{job_count}**, split as **{' + '.join(map(str, batch_sizes))}**.
- AlphaFold state: **prepared, not submitted**.
- The two original missing Hy.2E11 jobs were not retried, replaced, or silently excluded.

## Experimental positive geometry

{feature_lines}

## PDB oracle ranks

{pdb_lines}

Hy.2E11 and Ob.1A12 have evaluable PDB oracle panels under the locked five-decoy minimum. Hy.1B11 remains PDB-not-evaluable because only two exact-HLA structural pair decoys are available for each positive. The new Ob.1A12 and Hy.1B11 AlphaFold panels are pending submission and download. These missing layers intentionally block a formal pass.

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")
    methods = f"""# Methods

## Eligibility

One biological system equals one paired human alpha-beta TCR or explicitly identified clone. Every strict system has functional recognition evidence for all declared arms, exact HLA alpha/beta allotypes, exact peptide sequences, resolved nine-residue registers, and an experimental structure for each pMHC arm. Hy.1B11 has two required self-microbial comparisons but one system vote. EBNA1-ANO2 remains prospective.

## Controls and leakage prevention

N3 comparators are selected before geometry is read by peptide-length difference, binding-percentile difference, a fixed seeded SHA-256 tie-break, and candidate ID. N3 rows are unknown recognition and are never used as specificity negatives. All future weight selection must hold out the complete biological system across PDB and AlphaFold layers.

## Features

The locked feature family contains exposed-position C-alpha RMSD, C-alpha-to-side-chain-centroid vector RMSD, a five-property physicochemical mismatch, and anchor-position C-alpha RMSD. Each feature becomes a within-panel average-tie rank percentile. Candidate weights are nonnegative quarter increments summing to one. The exposed-CA baseline remains frozen until held-out testing is complete.

## Trust rule

Every required pair must rank at most 3 in an evaluable PDB oracle panel and in both fixed AlphaFold seeds. A completed rank above 3 fails; a missing layer or seed is not evaluable. Only an overall pass permits separate within-HLA discovery reranking. Cross-allele consensus is prohibited.

{CLAIM_BOUNDARY}
"""
    (out / "METHODS.md").write_text(methods, encoding="utf-8")


def build_package(out: Path = DEFAULT_OUT) -> Dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    systems, ligands, pairs, sources = curated_registry()
    validation = validate_control_registry(systems, ligands, pairs)
    write_csv(out / "registry/control_system_registry.csv", systems)
    write_csv(out / "registry/control_ligand_registry.csv", ligands)
    write_csv(out / "registry/positive_pair_registry.csv", pairs)
    write_csv(out / "registry/literature_and_structure_sources.csv", sources)
    tcr3d_audit = _tcr3d_discovery_audit(out)
    write_csv(out / "registry/tcr3d_human_hla2_multiligand_screen.csv", tcr3d_audit)

    comparators, ob_predictions, dq_predictions = _select_comparator_rows(out, ligands)
    comparator_validation = validate_comparator_registry(
        comparators,
        expected_pair_ids=tuple(
            str(pair["pair_id"]) for pair in pairs
            if str(pair["system_id"]) != "SYS_BALF5_MBP_HY2E11"
        ),
    )
    write_csv(out / "controls/control_decoy_registry.csv", comparators)
    write_csv(out / "controls/ob1a12_positive_binding_prediction.csv", ob_predictions)
    write_csv(out / "controls/hy1b11_dq1_binding_predictions.csv", dq_predictions)
    comparisons = _comparison_universe(pairs, ligands, comparators)
    write_csv(out / "controls/new_control_comparison_universe.csv", comparisons)

    new_ligands = _new_af3_ligands(ligands, comparators)
    jobs, job_manifest, batches = build_af3_job_batches(
        new_ligands, _hla_sequences(), panel_seeds=PANEL_SEEDS, batch_size=30,
    )
    job_validation = validate_af3_job_package(batches, job_manifest, _hla_sequences())
    write_csv(out / "alphafold_jobs/job_manifest.csv", job_manifest)
    for index, batch in enumerate(batches, start=1):
        write_json(out / f"alphafold_jobs/hla2_controls_batch_{index:02d}_{len(batch)}_jobs.json", batch)

    structural_ligands, pdb_screen, pdb_matrix, pdb_summaries = _experimental_pdb_oracle(
        out, ligands, pairs
    )
    positive_features = [row for row in pdb_matrix if row["pair_role"] == "positive"]
    write_csv(out / "benchmark/pdb_structural_ligand_registry.csv", structural_ligands)
    write_csv(out / "benchmark/pdb_exact_hla_candidate_screen.csv", pdb_screen)
    write_csv(out / "benchmark/pdb_oracle_feature_matrix.csv", pdb_matrix)
    write_csv(out / "benchmark/pdb_positive_feature_matrix.csv", positive_features)
    skeleton = build_evaluation_skeleton(pairs, panel_seeds=PANEL_SEEDS)
    evaluations = _overlay_evaluation_status(skeleton, pdb_summaries)
    write_csv(out / "benchmark/evaluation_status.csv", evaluations)
    required_system_ids = validation["strict_system_ids"]
    trust_gate = build_trust_gate(evaluations, required_system_ids=required_system_ids)
    trust_gate.update({
        "frozen_exposed_ca_hy2_available_rank_1_both_seeds": True,
        "hy2_seed_104759_formal_pass": True,
        "hy2_seed_104729_status": "not_evaluable_incomplete_original_calibration",
        "original_missing_jobs_preserved": [
            "ebvms_native_human_background_2258538_s104729",
            "ebvms_native_ebv_iedb_f5efcef02cf6_s104729",
        ],
        "cross_allele_consensus_allowed": False,
    })
    write_json(out / "benchmark/trust_gate.json", trust_gate)

    weight_grid = [
        {"weight_id": index, **weights}
        for index, weights in enumerate(generate_weight_grid(FEATURES), start=1)
    ]
    write_csv(out / "benchmark/weight_grid.csv", weight_grid)
    write_csv(out / "benchmark/selected_weights.csv", [{
        "selection_status": "not_frozen_trust_gate_not_evaluable",
        "exposed_ca_rmsd_A": 1.0,
        "exposed_sidechain_vector_rmsd_A": 0.0,
        "tcr_face_physicochemical_mismatch": 0.0,
        "anchor_ca_rmsd_A": 0.0,
        "endpoint_role": "frozen_baseline_only",
        "discovery_application_allowed": False,
    }])
    outer_rows = [{
        "held_out_system_id": system_id,
        "evaluation_status": "not_evaluable_required_layers_incomplete",
        "training_system_ids": ";".join(value for value in required_system_ids if value != system_id),
        "positive_rank": "", "capture_at_3": "", "weights_frozen": False,
    } for system_id in required_system_ids]
    write_csv(out / "benchmark/outer_fold_results.csv", outer_rows)
    baseline_rows = [{
        "system_id": row["system_id"], "positive_pair_id": row["positive_pair_id"],
        "baseline": "experimental_pdb_exposed_ca_rmsd_A", "value": row["exposed_ca_rmsd_A"],
        "positive_rank": row["positive_rank"],
        "available_positive_rank": row["available_positive_rank"],
        "rank_status": row["evaluation_status"],
    } for row in positive_features]
    baseline_rows.append({
        "system_id": "SYS_BALF5_MBP_HY2E11", "positive_pair_id": "PAIR_HY2E11_BALF5_MBP",
        "baseline": "frozen_af3_exposed_ca_available_rank", "value": "1 in both available seed sets",
        "rank_status": "one_formal_pass_one_incomplete_seed",
    })
    write_csv(out / "benchmark/baseline_comparisons.csv", baseline_rows)
    write_csv(out / "benchmark/permutation_results.csv", [{
        "system_id": system_id, "permutation_status": "not_evaluable_required_rank_panels_incomplete",
        "permutation_count": 0, "p_value": "",
    } for system_id in required_system_ids])

    grouped_comparisons = defaultdict(int)
    for row in comparisons:
        grouped_comparisons[(str(row["positive_pair_id"]), int(row["panel_seed"]))] += 1
    if set(grouped_comparisons.values()) != {26}:
        raise ValueError("every new AlphaFold panel must contain one positive and 25 N3 pair decoys")
    missing_raw_links = sorted({
        str(row["raw_response_file"]) for row in comparators
        if not (out / str(row["raw_response_file"])).is_file()
    })
    if missing_raw_links:
        raise ValueError(f"missing comparator raw-response files: {missing_raw_links}")
    original_missing = {
        "ebvms_native_human_background_2258538_s104729",
        "ebvms_native_ebv_iedb_f5efcef02cf6_s104729",
    }
    if original_missing.intersection(str(job["name"]) for job in jobs):
        raise ValueError("an original missing Hy.2E11 job was silently retried")
    verification = {
        "package_integrity_status": "pass",
        "scientific_trust_gate_status": trust_gate["overall_trust_status"],
        "control_registry_validation": validation,
        "comparator_registry_validation": comparator_validation,
        "alphafold_job_validation": job_validation,
        "new_panel_comparison_counts": {
            f"{pair_id}|s{seed}": count
            for (pair_id, seed), count in sorted(grouped_comparisons.items())
        },
        "pdb_oracle_pair_summaries": pdb_summaries,
        "original_missing_hy2_jobs_retried": False,
        "original_missing_hy2_jobs": sorted(original_missing),
        "discovery_ranking_inputs_used": [],
        "discovery_rankings_modified": False,
        "cross_allele_consensus_generated": False,
        "n3_used_for_specificity_claims": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "validation/package_verification_summary.json", verification)

    _copy_frozen_hy2_context(out)
    _write_documentation(
        out, validation, trust_gate, len(jobs), [len(batch) for batch in batches],
        positive_features, pdb_summaries,
    )
    manifest = {
        "benchmark_version": "EBV_MS_HLA2_HELD_OUT_CONTROLS_V1",
        "strict_independent_system_count": validation["strict_independent_system_count"],
        "strict_positive_pair_count": validation["strict_positive_pair_count"],
        "new_comparator_pair_count": comparator_validation["comparison_pair_count"],
        "new_comparator_row_count": comparator_validation["comparator_row_count"],
        "new_unique_comparator_count": comparator_validation["unique_comparator_count"],
        "tcr3d_human_hla2_multiligand_group_count": len(tcr3d_audit) - 1,
        "pdb_oracle_structural_ligand_count": len(structural_ligands),
        "pdb_oracle_evaluable_positive_pair_count": sum(
            row["evaluation_status"] == "complete" for row in pdb_summaries
        ),
        "new_alphafold_job_count": len(jobs),
        "alphafold_batch_sizes": [len(batch) for batch in batches],
        "alphafold_state": "prepared_not_submitted",
        "overall_trust_status": trust_gate["overall_trust_status"],
        "discovery_rankings_modified": False,
        "cross_allele_consensus_generated": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "analysis_manifest.json", manifest)
    checksum_rows = []
    for path in sorted(value for value in out.rglob("*") if value.is_file() and value.name != "SHA256SUMS.csv"):
        checksum_rows.append({"relative_path": str(path.relative_to(out)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    write_csv(out / "SHA256SUMS.csv", checksum_rows)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    manifest = build_package(args.out)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
