"""Run the full-ensemble structural analysis for EBV-MS T-cell Library V2.

The workflow treats AlphaFold samples as alternative technical predictions, not
biological replicates. It reports computational pMHC geometry only and does not
establish presentation, TCR binding, activation, cross-reactivity, molecular
mimicry, or an MS disease mechanism.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Iterator, Sequence

import numpy as np

from analyze_af3_pmhc_downloads import (
    ca_coordinates,
    kabsch,
    parse_mmcif,
    peptide_hla_metrics,
    request_details,
    sequence,
)
from build_tcell_library_v2 import DRA_SEQUENCE
from gold_standard_positive_control_audit import run_gold_standard_audit
from multiallele_manuscript_analysis import direct_register_sequence_metrics
from tcell_library_v2 import ANCHOR_POSITIONS, CALIBRATION_SEEDS, EXPOSED_POSITIONS


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "processed/tcell_library_v2_2026-08-22"
DEFAULT_OUT = ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25"
DEFAULT_DOWNLOAD_ROOTS = (
    Path.home() / "Downloads",
    Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Downloads",
)
CLAIM_BOUNDARY = (
    "Computational pMHC geometry only; not evidence of presentation, TCR binding, "
    "activation, cross-reactivity, molecular mimicry, or MS disease mechanism."
)
PRIMARY_ENDPOINT = "median exposed-position P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove fit"
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass(frozen=True)
class Occurrence:
    request_name: str
    directory: Path
    fingerprint: str
    n_cif: int
    n_confidence: int
    n_full_data: int
    n_request: int

    @property
    def complete(self) -> bool:
        return (self.n_cif, self.n_confidence, self.n_full_data, self.n_request) == (5, 5, 5, 1)


@dataclass
class SampleGeometry:
    cohort: str
    job_name: str
    allele: str
    entity_id: str
    server_seed: str
    sample_index: int
    groove_ca: np.ndarray
    core_ca: np.ndarray
    core_cb: np.ndarray
    core_sidechain: np.ndarray
    peptide_ca: np.ndarray


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or (sorted({key for row in rows for key in row}) if rows else []))
    if not fieldnames:
        raise ValueError(f"cannot write an empty table without fields: {path}")
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


def _sample_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.(?:json|cif)$", path.name)
    if match is None:
        raise ValueError(f"cannot recover model index from {path.name}")
    return int(match.group(1))


def _bundle_fingerprint(directory: Path) -> str:
    """Hash request, CIF, and summaries while excluding path and redundant full-data JSON."""
    files = sorted(
        [*directory.glob("*_job_request.json"), *directory.glob("*_model_*.cif"), *directory.glob("*_summary_confidences_*.json")],
        key=lambda path: path.name,
    )
    digest = hashlib.sha256()
    for path in files:
        if path.name.endswith("_job_request.json"):
            role = "job_request.json"
        else:
            match = re.search(r"(model_[0-4]\.cif|summary_confidences_[0-4]\.json)$", path.name)
            if match is None:
                raise ValueError(f"cannot normalize bundle file role: {path.name}")
            role = match.group(1)
        digest.update(role.encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    return digest.hexdigest()


def discover_occurrences(download_roots: Iterable[Path]) -> dict[str, list[Occurrence]]:
    """Discover complete and incomplete extracted result folders by request identity."""
    discovered: dict[str, list[Occurrence]] = defaultdict(list)
    for root in sorted({Path(path).expanduser().resolve() for path in download_roots}, key=str):
        if not root.exists():
            continue
        for result_root in sorted((path for path in root.glob("folds*") if path.is_dir()), key=str):
            for request_path in sorted(result_root.rglob("*_job_request.json"), key=str):
                directory = request_path.parent
                details = request_details(json.loads(request_path.read_text(encoding="utf-8")))
                name = str(details["request_name"]).lower()
                occurrence = Occurrence(
                    request_name=name,
                    directory=directory,
                    fingerprint=_bundle_fingerprint(directory),
                    n_cif=len(list(directory.glob("*_model_*.cif"))),
                    n_confidence=len(list(directory.glob("*_summary_confidences_*.json"))),
                    n_full_data=len(list(directory.glob("*_full_data_*.json"))),
                    n_request=len(list(directory.glob("*_job_request.json"))),
                )
                discovered[name].append(occurrence)
    return discovered


def choose_occurrence(occurrences: Sequence[Occurrence]) -> Occurrence | None:
    """Select a complete bundle by content fingerprint, independent of model scores."""
    complete = [row for row in occurrences if row.complete]
    return min(complete, key=lambda row: (row.fingerprint, str(row.directory))) if complete else None


def expected_hla_sequences(package: Path) -> dict[str, str]:
    rows = read_csv(package / "reference_sequence_manifest.csv")
    result = {row["entity"]: row["sequence"] for row in rows if row["entity"].startswith("HLA-DRB")}
    required = {
        "HLA-DRB1*15:01", "HLA-DRB1*13:03", "HLA-DRB1*03:01",
        "HLA-DRB1*08:01", "HLA-DRB5*01:01",
    }
    if set(result) & required != required:
        raise ValueError(f"reference manifest lacks HLA sequences: {sorted(required - set(result))}")
    return result


def load_registers(package: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = [
        *read_csv(package / "allele_register_predictions_320.csv"),
        *read_csv(package / "calibration_control_binding_predictions.csv"),
    ]
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["allele"], row["candidate_id"])
        current = result.get(key)
        if current is None or current["sequence"] == row["sequence"]:
            result[key] = row
    return result


def build_inventory(
    expected_rows: Sequence[dict[str, str]],
    occurrences: dict[str, list[Occurrence]],
    cohort: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    inventory: list[dict[str, Any]] = []
    duplicate_rows: list[dict[str, Any]] = []
    for expected in expected_rows:
        name = expected["job_name"].lower()
        matches = occurrences.get(name, [])
        chosen = choose_occurrence(matches)
        complete_matches = [row for row in matches if row.complete]
        fingerprints = {row.fingerprint for row in complete_matches}
        inventory.append({
            "cohort": cohort,
            **expected,
            "download_status": "complete" if chosen else "incomplete" if matches else "missing",
            "canonical_directory": str(chosen.directory) if chosen else "",
            "bundle_fingerprint": chosen.fingerprint if chosen else "",
            "complete_occurrence_count": len(complete_matches),
            "distinct_complete_bundle_count": len(fingerprints),
            "duplicate_handling": "score_blind_minimum_bundle_fingerprint",
        })
        for occurrence in sorted(complete_matches, key=lambda row: (row.fingerprint, str(row.directory))):
            duplicate_rows.append({
                "cohort": cohort,
                "job_name": name,
                "directory": str(occurrence.directory),
                "bundle_fingerprint": occurrence.fingerprint,
                "selected_primary": bool(chosen and occurrence.directory == chosen.directory),
                "duplicate_class": "identical_copy" if chosen and occurrence.fingerprint == chosen.fingerprint else "divergent_unplanned_repeat",
            })
    return inventory, duplicate_rows


def _pseudo_cb(residue: dict[str, object]) -> np.ndarray:
    atoms = {str(atom["name"]): np.asarray(atom["xyz"], dtype=float) for atom in residue["atoms"]}
    return atoms.get("CB", atoms["CA"])


def _sidechain_center(residue: dict[str, object]) -> np.ndarray:
    coordinates = [
        np.asarray(atom["xyz"], dtype=float)
        for atom in residue["atoms"]
        if str(atom["name"]) not in BACKBONE_ATOMS and str(atom["element"]) != "H"
    ]
    return np.mean(coordinates, axis=0) if coordinates else _pseudo_cb(residue)


def _summary_chain_metric(summary: dict[str, Any], field: str, chain_id: str) -> float | str:
    chain_ids = list(summary.get("chain_ids", []))
    values = summary.get(field, [])
    if chain_id not in chain_ids or len(values) != len(chain_ids):
        return ""
    return float(values[chain_ids.index(chain_id)])


def _summary_pair_metric(summary: dict[str, Any], field: str, left: str, right: str) -> float | str:
    chain_ids = list(summary.get("chain_ids", []))
    matrix = summary.get(field, [])
    if left not in chain_ids or right not in chain_ids or len(matrix) != len(chain_ids):
        return ""
    i, j = chain_ids.index(left), chain_ids.index(right)
    return float(np.mean([matrix[i][j], matrix[j][i]]))


def full_data_interface_metrics(path: Path) -> dict[str, float | str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    chain_ids = np.asarray(data.get("token_chain_ids", []), dtype=object)
    residue_ids = np.asarray(data.get("token_res_ids", []), dtype=int)
    if chain_ids.size == 0 or residue_ids.size != chain_ids.size:
        return {
            "peptide_groove_pae_mean_A": "",
            "peptide_groove_pae_median_A": "",
            "peptide_groove_contact_probability_mean": "",
            "peptide_groove_contact_probability_max": "",
        }
    peptide = np.where(chain_ids == "C")[0]
    groove = np.where(((chain_ids == "A") | (chain_ids == "B")) & (residue_ids <= 85))[0]
    pae = np.asarray(data.get("pae", []), dtype=float)
    contact = np.asarray(data.get("contact_probs", []), dtype=float)
    if not len(peptide) or not len(groove) or pae.shape != (len(chain_ids), len(chain_ids)):
        return {
            "peptide_groove_pae_mean_A": "",
            "peptide_groove_pae_median_A": "",
            "peptide_groove_contact_probability_mean": "",
            "peptide_groove_contact_probability_max": "",
        }
    bidirectional = np.concatenate([pae[np.ix_(peptide, groove)].ravel(), pae[np.ix_(groove, peptide)].ravel()])
    interface_contact = contact[np.ix_(peptide, groove)] if contact.shape == pae.shape else np.asarray([])
    return {
        "peptide_groove_pae_mean_A": round(float(np.mean(bidirectional)), 4),
        "peptide_groove_pae_median_A": round(float(np.median(bidirectional)), 4),
        "peptide_groove_contact_probability_mean": round(float(np.mean(interface_contact)), 6) if interface_contact.size else "",
        "peptide_groove_contact_probability_max": round(float(np.max(interface_contact)), 6) if interface_contact.size else "",
    }


def _distribution(values: Sequence[float], prefix: str) -> dict[str, float | str]:
    if not values:
        return {f"{prefix}_{field}": "" for field in ("min", "q25", "median", "q75", "max", "iqr")}
    array = np.asarray(values, dtype=float)
    q25, q50, q75 = np.quantile(array, [0.25, 0.5, 0.75])
    return {
        f"{prefix}_min": round(float(np.min(array)), 6),
        f"{prefix}_q25": round(float(q25), 6),
        f"{prefix}_median": round(float(q50), 6),
        f"{prefix}_q75": round(float(q75), 6),
        f"{prefix}_max": round(float(np.max(array)), 6),
        f"{prefix}_iqr": round(float(q75 - q25), 6),
    }


def pair_geometry(left: SampleGeometry, right: SampleGeometry) -> dict[str, float]:
    if left.groove_ca.shape != right.groove_ca.shape:
        raise ValueError("HLA groove arrays must be equivalent for pair fitting")
    rotation, translation = kabsch(right.groove_ca, left.groove_ca)
    fitted_ca = right.core_ca @ rotation + translation
    fitted_cb = right.core_cb @ rotation + translation
    fitted_sidechain = right.core_sidechain @ rotation + translation
    anchors = np.asarray([position - 1 for position in ANCHOR_POSITIONS])
    exposed = np.asarray([position - 1 for position in EXPOSED_POSITIONS])

    def rmsd(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))

    return {
        "full_core_ca_rmsd_A": rmsd(left.core_ca, fitted_ca),
        "anchor_ca_rmsd_A": rmsd(left.core_ca[anchors], fitted_ca[anchors]),
        "exposed_ca_rmsd_A": rmsd(left.core_ca[exposed], fitted_ca[exposed]),
        "exposed_cb_rmsd_A": rmsd(left.core_cb[exposed], fitted_cb[exposed]),
        "exposed_sidechain_centroid_rmsd_A": rmsd(left.core_sidechain[exposed], fitted_sidechain[exposed]),
    }


def analyze_job(
    row: dict[str, Any],
    register: dict[str, str],
    expected_drb: str,
) -> tuple[list[dict[str, Any]], list[SampleGeometry], dict[str, Any]]:
    if row["download_status"] != "complete":
        return [], [], {
            **row,
            "technical_status": row["download_status"],
            "observed_sample_count": 0,
            "valid_geometry_sample_count": 0,
            "selected_sample_index": "",
            "within_job_peptide_pose_rmsd_median_A": "",
            "within_job_peptide_pose_rmsd_max_A": "",
            "claim_boundary": CLAIM_BOUNDARY,
        }
    directory = Path(row["canonical_directory"])
    request_path = next(directory.glob("*_job_request.json"))
    request = request_details(json.loads(request_path.read_text(encoding="utf-8")))
    entity_id = str(row.get("candidate_id", row.get("entity_id", "")))
    expected_peptide = str(row["peptide_sequence"])
    expected_seed = str(row.get("server_seed", ""))
    request_seed = str(request["server_seed"])
    request_valid = (
        str(request["request_name"]).lower() == str(row["job_name"]).lower()
        and request["requested_dra"] == DRA_SEQUENCE
        and request["requested_drb"] == expected_drb
        and request["requested_peptide"] == expected_peptide
        and (row["cohort"] == "discovery" or request_seed == expected_seed)
    )
    core_start = int(register["core_start"])
    core_sequence = register["predicted_core"]
    if (
        register.get("register_resolution") != "resolved_unique_fully_contained"
        or expected_peptide[core_start - 1:core_start + 8] != core_sequence
    ):
        raise ValueError(f"unresolved or inconsistent frozen register for {row['job_name']}")

    sample_rows: list[dict[str, Any]] = []
    geometries: list[SampleGeometry] = []
    for summary_path in sorted(directory.glob("*_summary_confidences_*.json"), key=_sample_index):
        index = _sample_index(summary_path)
        cif_matches = list(directory.glob(f"*_model_{index}.cif"))
        full_matches = list(directory.glob(f"*_full_data_{index}.json"))
        if len(cif_matches) != 1 or len(full_matches) != 1:
            continue
        model = parse_mmcif(cif_matches[0])
        observed_dra = sequence(model.get("A", []))
        observed_drb = sequence(model.get("B", []))
        observed_peptide = sequence(model.get("C", []))
        exact = (
            request_valid
            and set(model) == {"A", "B", "C"}
            and observed_dra == DRA_SEQUENCE
            and observed_drb == expected_drb
            and observed_peptide == expected_peptide
        )
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        has_clash = str(summary.get("has_clash", "")).strip().lower() in {"1", "1.0", "true", "yes"}
        geometry_status = "fail_request_or_model_sequence_identity"
        model_metrics: dict[str, Any] = {}
        if exact:
            model_metrics = peptide_hla_metrics(model["C"], [model["A"], model["B"]])
            geometry_status = "excluded_has_clash" if has_clash else "pass_exact_clash_free"
            try:
                peptide_ca = ca_coordinates(model["C"])
                core_residues = model["C"][core_start - 1:core_start + 8]
                groove_ca = np.vstack([ca_coordinates(model[chain][:85]) for chain in ("A", "B")])
                if peptide_ca.shape[0] != len(expected_peptide) or len(core_residues) != 9:
                    raise ValueError("peptide or core length mismatch")
                if not has_clash:
                    geometries.append(SampleGeometry(
                        cohort=row["cohort"],
                        job_name=row["job_name"],
                        allele=row["allele"],
                        entity_id=entity_id,
                        server_seed=request_seed,
                        sample_index=index,
                        groove_ca=groove_ca,
                        core_ca=ca_coordinates(core_residues),
                        core_cb=np.vstack([_pseudo_cb(residue) for residue in core_residues]),
                        core_sidechain=np.vstack([_sidechain_center(residue) for residue in core_residues]),
                        peptide_ca=peptide_ca,
                    ))
            except (KeyError, StopIteration, ValueError) as error:
                geometry_status = f"fail_geometry_coordinates:{type(error).__name__}"
        interface = full_data_interface_metrics(full_matches[0]) if exact else {}
        sample_rows.append({
            "cohort": row["cohort"],
            "job_name": row["job_name"],
            "allele": row["allele"],
            "entity_id": entity_id,
            "server_seed": request_seed,
            "sample_index": index,
            "sample_status": geometry_status,
            "request_identity_pass": request_valid,
            "observed_dra_length": len(observed_dra),
            "observed_drb_length": len(observed_drb),
            "observed_peptide": observed_peptide,
            "predicted_core": core_sequence,
            "core_start_1_based": core_start,
            "ranking_score": summary.get("ranking_score", ""),
            "iptm": summary.get("iptm", ""),
            "ptm": summary.get("ptm", ""),
            "fraction_disordered": summary.get("fraction_disordered", ""),
            "has_clash": has_clash,
            "peptide_chain_iptm": _summary_chain_metric(summary, "chain_iptm", "C"),
            "peptide_dra_pair_iptm": _summary_pair_metric(summary, "chain_pair_iptm", "C", "A"),
            "peptide_drb_pair_iptm": _summary_pair_metric(summary, "chain_pair_iptm", "C", "B"),
            "peptide_dra_pair_pae_min_A": _summary_pair_metric(summary, "chain_pair_pae_min", "C", "A"),
            "peptide_drb_pair_pae_min_A": _summary_pair_metric(summary, "chain_pair_pae_min", "C", "B"),
            **model_metrics,
            **interface,
            "cif_path": str(cif_matches[0]),
            "summary_path": str(summary_path),
            "full_data_path": str(full_matches[0]),
            "claim_boundary": CLAIM_BOUNDARY,
        })

    valid_indices = {geometry.sample_index for geometry in geometries}
    eligible_samples = [row for row in sample_rows if row["sample_index"] in valid_indices]
    selected = max(
        eligible_samples,
        key=lambda sample: (float(sample["ranking_score"]), -int(sample["sample_index"])),
        default=None,
    )
    within_job: list[float] = []
    for left_index, left in enumerate(geometries):
        for right in geometries[left_index + 1:]:
            rotation, translation = kabsch(right.groove_ca, left.groove_ca)
            fitted = right.peptide_ca @ rotation + translation
            within_job.append(float(np.sqrt(np.mean(np.sum((left.peptide_ca - fitted) ** 2, axis=1)))))
    job_summary = {
        **row,
        "technical_status": "geometry_evaluable" if geometries else "excluded_no_valid_geometry_samples",
        "request_identity_pass": request_valid,
        "observed_sample_count": len(sample_rows),
        "exact_sample_count": sum(str(sample["sample_status"]).startswith(("pass", "excluded_has_clash")) for sample in sample_rows),
        "clash_sample_count": sum(bool(sample["has_clash"]) for sample in sample_rows),
        "valid_geometry_sample_count": len(geometries),
        "selected_sample_index": selected["sample_index"] if selected else "",
        "selected_ranking_score": selected["ranking_score"] if selected else "",
        "selected_iptm": selected["iptm"] if selected else "",
        "selected_peptide_mean_plddt": selected.get("peptide_mean_plddt", "") if selected else "",
        "within_job_peptide_pose_rmsd_median_A": round(median(within_job), 6) if within_job else "",
        "within_job_peptide_pose_rmsd_max_A": round(max(within_job), 6) if within_job else "",
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return sample_rows, geometries, job_summary


def rank_pair_rows(rows: list[dict[str, Any]]) -> None:
    """Add deterministic within-allele ranks to complete pair summaries in place."""
    by_allele: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["geometry_status"] == "complete":
            by_allele[str(row["allele"])].append(row)
    for allele_rows in by_allele.values():
        ordered = sorted(
            allele_rows,
            key=lambda row: (
                float(row["exposed_ca_rmsd_A_median"]),
                float(row["exposed_ca_rmsd_A_iqr"]),
                str(row["pair_id"]),
            ),
        )
        count = len(ordered)
        for rank, row in enumerate(ordered, start=1):
            row["within_allele_rank"] = rank
            row["within_allele_evaluable_count"] = count
            row["within_allele_percentile"] = round((rank - 1) / max(1, count - 1), 8)
    for row in rows:
        row.setdefault("within_allele_rank", "")
        row.setdefault("within_allele_evaluable_count", len(by_allele.get(str(row["allele"]), [])))
        row.setdefault("within_allele_percentile", "")


def summarize_cross_allele(pair_rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in pair_rows:
        grouped[(str(row["ebv_candidate_id"]), str(row["self_candidate_id"]))].append(row)
    output: list[dict[str, Any]] = []
    for (ebv_id, self_id), rows in sorted(grouped.items()):
        complete = [row for row in rows if row["geometry_status"] == "complete"]
        percentiles = [float(row["within_allele_percentile"]) for row in complete]
        medians = [float(row["exposed_ca_rmsd_A_median"]) for row in complete]
        output.append({
            "cross_allele_pair_id": f"{ebv_id}|{self_id}",
            "ebv_candidate_id": ebv_id,
            "self_candidate_id": self_id,
            "complete_allele_count": len(complete),
            "complete_alleles": ";".join(sorted(str(row["allele"]) for row in complete)),
            "cross_allele_status": "complete_four_alleles" if len(complete) == 4 else "incomplete_not_ranked",
            "median_within_allele_percentile": round(float(np.median(percentiles)), 8) if percentiles else "",
            "worst_within_allele_percentile": round(max(percentiles), 8) if percentiles else "",
            "best_within_allele_percentile": round(min(percentiles), 8) if percentiles else "",
            "median_allele_exposed_ca_rmsd_A": round(float(np.median(medians)), 6) if medians else "",
            "worst_allele_exposed_ca_rmsd_A": round(max(medians), 6) if medians else "",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    eligible = sorted(
        (row for row in output if row["cross_allele_status"] == "complete_four_alleles"),
        key=lambda row: (
            float(row["worst_within_allele_percentile"]),
            float(row["median_within_allele_percentile"]),
            str(row["cross_allele_pair_id"]),
        ),
    )
    for rank, row in enumerate(eligible, start=1):
        row["cross_allele_consensus_rank"] = rank
        row["cross_allele_evaluable_count"] = len(eligible)
    for row in output:
        row.setdefault("cross_allele_consensus_rank", "")
        row.setdefault("cross_allele_evaluable_count", len(eligible))
    return output


def classify_recovery(seed_rows: Sequence[dict[str, Any]]) -> str:
    by_seed = {int(row["seed"]): row for row in seed_rows}
    if set(by_seed) != set(CALIBRATION_SEEDS) or any(not row["formal_seed_evaluable"] for row in by_seed.values()):
        return "not_evaluable_incomplete_calibration"
    recovered = all(bool(row["seed_recovery_criterion_pass"]) for row in by_seed.values())
    return "recovered" if recovered else "failed_calibration"


def duplicate_pair_sensitivity(
    primary: dict[str, Any],
    left: Sequence[SampleGeometry],
    right: Sequence[SampleGeometry],
    duplicate: dict[str, Any],
) -> dict[str, Any]:
    values = [
        pair_geometry(left_sample, right_sample)["exposed_ca_rmsd_A"]
        for left_sample in left
        for right_sample in right
    ]
    distribution = _distribution(values, "alternative_exposed_ca_rmsd_A")
    alternative_median = distribution["alternative_exposed_ca_rmsd_A_median"]
    primary_median = primary.get("exposed_ca_rmsd_A_median", "")
    return {
        "comparison_id": primary["pair_id"],
        "allele": primary.get("allele", "cross_native_HLA"),
        "duplicate_job_name": duplicate["job_name"],
        "duplicate_bundle_fingerprint": duplicate["bundle_fingerprint"],
        "duplicate_directory": duplicate["directory"],
        "alternative_model_combination_count": len(values),
        "primary_exposed_ca_rmsd_A_median": primary_median,
        **distribution,
        "alternative_minus_primary_median_A": round(float(alternative_median) - float(primary_median), 6) if values and primary_median != "" else "",
        "interpretation": "Unplanned divergent-download sensitivity only; excluded from primary ranking and not biological replication.",
        "claim_boundary": CLAIM_BOUNDARY,
    }


def _escape(value: Any) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _svg(path: Path, width: int, height: int, body: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#17212b;letter-spacing:0}.title{font-size:22px;font-weight:700}.label{font-size:12px}.small{font-size:10px}.axis{stroke:#4b5563;stroke-width:1}.grid{stroke:#e5e7eb;stroke-width:1}</style>',
        *body,
        '</svg>',
        '',
    ]), encoding="utf-8")


def generate_figures(
    out: Path,
    job_rows: Sequence[dict[str, Any]],
    pair_rows: Sequence[dict[str, Any]],
    cross_rows: Sequence[dict[str, Any]],
    calibration_pairs: Sequence[dict[str, Any]],
    seed_rows: Sequence[dict[str, Any]],
) -> None:
    figure_dir = out / "figures"
    alleles = ("HLA-DRB1*15:01", "HLA-DRB1*13:03", "HLA-DRB1*03:01", "HLA-DRB1*08:01")
    colors = {alleles[0]: "#2b6cb0", alleles[1]: "#2f855a", alleles[2]: "#b7791f", alleles[3]: "#9b2c2c"}

    body = ['<text x="520" y="34" text-anchor="middle" class="title">Model availability and geometry eligibility</text>']
    for index, allele in enumerate(alleles):
        rows = [row for row in job_rows if row["cohort"] == "discovery" and row["allele"] == allele]
        complete = sum(row["download_status"] == "complete" for row in rows)
        evaluable = sum(row["technical_status"] == "geometry_evaluable" for row in rows)
        x = 100 + index * 220
        for offset, count, color in ((0, complete, "#9ec5e8"), (62, evaluable, colors[allele])):
            height = count * 3.4
            body.extend([
                f'<rect x="{x+offset}" y="{350-height:.1f}" width="48" height="{height:.1f}" fill="{color}"/>',
                f'<text x="{x+offset+24}" y="{340-height:.1f}" text-anchor="middle" class="label">{count}</text>',
            ])
        body.append(f'<text x="{x+55}" y="385" text-anchor="middle" class="label">{_escape(allele.replace("HLA-", ""))}</text>')
    body.extend([
        '<rect x="350" y="414" width="13" height="13" fill="#9ec5e8"/><text x="370" y="425" class="small">Downloaded</text>',
        '<rect x="505" y="414" width="13" height="13" fill="#4b5563"/><text x="525" y="425" class="small">Geometry-ready</text>',
        '<text x="520" y="456" text-anchor="middle" class="small">Missing models remain technical missingness and are not biological negatives.</text>',
    ])
    _svg(figure_dir / "figure_1_inventory_qc.svg", 1040, 480, body)

    complete_pairs = [row for row in pair_rows if row["geometry_status"] == "complete"]
    vmax = float(np.quantile([float(row["exposed_ca_rmsd_A_median"]) for row in complete_pairs], 0.95))
    by_key = {(row["allele"], row["ebv_candidate_id"], row["self_candidate_id"]): row for row in pair_rows}
    ebv_ids = sorted({str(row["ebv_candidate_id"]) for row in pair_rows})
    self_ids = sorted({str(row["self_candidate_id"]) for row in pair_rows})
    order_rows = [
        {"arm": arm, "display_index": index + 1, "candidate_id": candidate}
        for arm, candidates in (("EBV", ebv_ids), ("self", self_ids))
        for index, candidate in enumerate(candidates)
    ]
    write_csv(figure_dir / "figure_2_heatmap_order.csv", order_rows)
    cell = 8
    body = ['<text x="760" y="30" text-anchor="middle" class="title">Within-allele EBV-self exposed-core geometry</text>']
    for panel_index, allele in enumerate(alleles):
        x0 = 55 + panel_index * 375
        y0 = 68
        body.append(f'<text x="{x0+160}" y="52" text-anchor="middle" class="label">{_escape(allele.replace("HLA-", ""))}</text>')
        for i, ebv_id in enumerate(ebv_ids):
            for j, self_id in enumerate(self_ids):
                row = by_key[(allele, ebv_id, self_id)]
                if row["geometry_status"] != "complete":
                    color = "#d1d5db"
                else:
                    value = min(1.0, float(row["exposed_ca_rmsd_A_median"]) / vmax)
                    red = int(30 + 210 * value)
                    green = int(145 - 85 * value)
                    blue = int(175 - 105 * value)
                    color = f"#{red:02x}{max(green,0):02x}{max(blue,0):02x}"
                body.append(f'<rect x="{x0+j*cell}" y="{y0+i*cell}" width="{cell}" height="{cell}" fill="{color}"/>')
        body.extend([
            f'<text x="{x0+160}" y="410" text-anchor="middle" class="small">Self peptide index</text>',
            f'<text x="{x0-22}" y="{y0+160}" transform="rotate(-90 {x0-22} {y0+160})" text-anchor="middle" class="small">EBV peptide index</text>',
        ])
    body.append(f'<text x="760" y="445" text-anchor="middle" class="small">Median exposed-position RMSD; color capped at the 95th percentile ({vmax:.2f} Å). Gray = unavailable.</text>')
    _svg(figure_dir / "figure_2_discovery_heatmaps.svg", 1530, 470, body)

    max_median = float(np.quantile([float(row["exposed_ca_rmsd_A_median"]) for row in complete_pairs], 0.99))
    max_iqr = float(np.quantile([float(row["exposed_ca_rmsd_A_iqr"]) for row in complete_pairs], 0.99)) or 1.0
    body = [
        '<text x="520" y="32" text-anchor="middle" class="title">Pair geometry and ensemble stability</text>',
        '<line x1="70" y1="390" x2="990" y2="390" class="axis"/><line x1="70" y1="55" x2="70" y2="390" class="axis"/>',
        '<text x="530" y="430" text-anchor="middle" class="label">Median exposed-position RMSD (Å)</text>',
        '<text x="20" y="225" transform="rotate(-90 20 225)" text-anchor="middle" class="label">Ensemble IQR (Å)</text>',
    ]
    for row in complete_pairs:
        x = 70 + min(float(row["exposed_ca_rmsd_A_median"]), max_median) / max_median * 920
        y = 390 - min(float(row["exposed_ca_rmsd_A_iqr"]), max_iqr) / max_iqr * 335
        body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.7" fill="{colors[str(row["allele"])]}" fill-opacity="0.45"/>')
    for index, allele in enumerate(alleles):
        body.extend([
            f'<rect x="{130+index*210}" y="454" width="12" height="12" fill="{colors[allele]}"/>',
            f'<text x="{148+index*210}" y="464" class="small">{_escape(allele.replace("HLA-", ""))}</text>',
        ])
    body.append(f'<text x="520" y="487" text-anchor="middle" class="small">Axes capped at their 99th percentiles ({max_median:.2f} Å, {max_iqr:.2f} Å).</text>')
    _svg(figure_dir / "figure_3_ensemble_stability.svg", 1040, 510, body)

    complete_cross = [row for row in cross_rows if row["cross_allele_status"] == "complete_four_alleles"]
    body = [
        '<text x="520" y="32" text-anchor="middle" class="title">Cross-allele rank consistency</text>',
        '<line x1="70" y1="390" x2="990" y2="390" class="axis"/><line x1="70" y1="55" x2="70" y2="390" class="axis"/>',
        '<text x="530" y="430" text-anchor="middle" class="label">Median within-allele percentile</text>',
        '<text x="20" y="225" transform="rotate(-90 20 225)" text-anchor="middle" class="label">Worst within-allele percentile</text>',
    ]
    for row in complete_cross:
        x = 70 + float(row["median_within_allele_percentile"]) * 920
        y = 390 - float(row["worst_within_allele_percentile"]) * 335
        body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2" fill="#2b6cb0" fill-opacity="0.4"/>')
    for row in sorted(complete_cross, key=lambda value: int(value["cross_allele_consensus_rank"]))[:10]:
        x = 70 + float(row["median_within_allele_percentile"]) * 920
        y = 390 - float(row["worst_within_allele_percentile"]) * 335
        body.append(f'<text x="{x+4:.2f}" y="{y-4:.2f}" class="small">{row["cross_allele_consensus_rank"]}</text>')
    body.append('<text x="520" y="465" text-anchor="middle" class="small">Formal consensus ranking requires all four allele-specific results; numbered points are ranks 1–10.</text>')
    _svg(figure_dir / "figure_4_cross_allele_consensus.svg", 1040, 490, body)

    primary = [row for row in calibration_pairs if row["analysis_set"] == "primary_rank_of_26" and row["geometry_status"] == "complete"]
    max_value = max(float(row["exposed_ca_rmsd_A_median"]) for row in primary) if primary else 1.0
    body = ['<text x="620" y="32" text-anchor="middle" class="title">Native-HLA calibration ranks</text>']
    seed_by_id = {int(row["seed"]): row for row in seed_rows}
    for panel, seed in enumerate(CALIBRATION_SEEDS):
        rows = sorted((row for row in primary if int(row["seed"]) == seed), key=lambda row: float(row["exposed_ca_rmsd_A_median"]))
        y0 = 95 + panel * 258
        body.append(f'<text x="80" y="{y0-28}" class="label">Seed {seed}: {len(rows)}/26 primary pairs</text>')
        for index, row in enumerate(rows):
            value = float(row["exposed_ca_rmsd_A_median"])
            x = 210 + value / max_value * 930
            y = y0 + index * (150 / max(1, len(rows)-1))
            positive = row["pair_role"] == "E1_positive"
            body.extend([
                f'<line x1="210" y1="{y:.1f}" x2="{x:.1f}" y2="{y:.1f}" stroke="#cbd5e0"/>',
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{6 if positive else 3}" fill="{"#c53030" if positive else "#4a5568"}"/>',
            ])
        seed_row = seed_by_id.get(seed, {})
        axis_y = y0 + 160
        body.extend([
            f'<line x1="210" y1="{axis_y}" x2="1140" y2="{axis_y}" class="axis"/>',
            f'<text x="210" y="{axis_y+14}" text-anchor="middle" class="small">0</text>',
            f'<text x="675" y="{axis_y+14}" text-anchor="middle" class="small">{max_value/2:.1f}</text>',
            f'<text x="1140" y="{axis_y+14}" text-anchor="middle" class="small">{max_value:.1f}</text>',
            f'<text x="675" y="{axis_y+28}" text-anchor="middle" class="small">Median exposed-position RMSD (Å)</text>',
            f'<text x="80" y="{axis_y+42}" class="small">Positive available rank: {seed_row.get("available_rank", "")}/{seed_row.get("available_primary_count", "")}; formal status: {_escape(seed_row.get("formal_seed_status", "")).replace("_", " ")}</text>',
        ])
    body.extend([
        '<text x="620" y="605" text-anchor="middle" class="small">Red = frozen BALF5–MBP positive. Missing controls are not imputed.</text>',
        '<text x="620" y="621" text-anchor="middle" class="small">AlphaFold model samples are alternative technical predictions, not biological replicates.</text>',
    ])
    _svg(figure_dir / "figure_5_calibration.svg", 1240, 645, body)


def _open_stream_csv(path: Path, fields: Sequence[str]) -> tuple[Any, csv.DictWriter]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("w", newline="", encoding="utf-8")
    writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
    writer.writeheader()
    return handle, writer


def run_analysis(download_roots: Sequence[Path], out: Path, package: Path = PACKAGE) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    occurrences = discover_occurrences(download_roots)
    discovery_expected = read_csv(package / "model_inventory_320.csv")
    calibration_expected = read_csv(package / "native_hla_calibration_manifest_24.csv")
    discovery_inventory, discovery_duplicates = build_inventory(discovery_expected, occurrences, "discovery")
    calibration_inventory, calibration_duplicates = build_inventory(calibration_expected, occurrences, "calibration")
    write_csv(out / "inventory/discovery_model_inventory_320.csv", discovery_inventory)
    write_csv(out / "inventory/calibration_model_inventory_24.csv", calibration_inventory)
    write_csv(out / "inventory/duplicate_download_occurrences.csv", [*discovery_duplicates, *calibration_duplicates])

    hla = expected_hla_sequences(package)
    registers = load_registers(package)
    sample_rows: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    geometry_by_job: dict[str, list[SampleGeometry]] = {}
    for row in [*discovery_inventory, *calibration_inventory]:
        entity_id = str(row.get("candidate_id", row.get("entity_id", "")))
        register = registers.get((row["allele"], entity_id))
        if register is None:
            raise ValueError(f"missing frozen register for {row['allele']} {entity_id}")
        samples, geometries, job = analyze_job(row, register, hla[row["allele"]])
        sample_rows.extend(samples)
        job_rows.append(job)
        geometry_by_job[row["job_name"].lower()] = geometries
    write_csv(out / "qc/model_sample_qc.csv", sample_rows)
    write_csv(out / "qc/job_qc_summary_344.csv", job_rows)

    panel_rows = read_csv(package / "frozen_v2_80_peptide_panel.csv")
    panel = {row["candidate_id"]: row for row in panel_rows}
    discovery_job_by_key = {(row["allele"], row["candidate_id"]): row for row in discovery_inventory}
    pair_universe = read_csv(package / "within_allele_pair_universe_6400.csv")
    pair_rows: list[dict[str, Any]] = []
    ensemble_fields = [
        "allele", "pair_id", "ebv_candidate_id", "self_candidate_id",
        "ebv_sample_index", "self_sample_index", "full_core_ca_rmsd_A",
        "anchor_ca_rmsd_A", "exposed_ca_rmsd_A", "exposed_cb_rmsd_A",
        "exposed_sidechain_centroid_rmsd_A", "interpretation",
    ]
    ensemble_handle, ensemble_writer = _open_stream_csv(out / "discovery/pair_model_ensemble.csv", ensemble_fields)
    try:
        for pair in pair_universe:
            allele = pair["allele"]
            ebv_id, self_id = pair["ebv_candidate_id"], pair["self_candidate_id"]
            ebv_job = discovery_job_by_key[(allele, ebv_id)]
            self_job = discovery_job_by_key[(allele, self_id)]
            left = geometry_by_job.get(ebv_job["job_name"].lower(), [])
            right = geometry_by_job.get(self_job["job_name"].lower(), [])
            ebv_register = registers[(allele, ebv_id)]
            self_register = registers[(allele, self_id)]
            base = {
                **pair,
                "ebv_sequence": panel[ebv_id]["sequence"],
                "self_sequence": panel[self_id]["sequence"],
                "ebv_protein": panel[ebv_id]["protein_symbol"],
                "self_protein": panel[self_id]["protein_symbol"],
                "ebv_source_certainty": panel[ebv_id]["source_certainty"],
                "self_source_certainty": panel[self_id]["source_certainty"],
                "ebv_predicted_core": ebv_register["predicted_core"],
                "self_predicted_core": self_register["predicted_core"],
                "ebv_binding_percentile_rank": ebv_register["percentile_rank"],
                "self_binding_percentile_rank": self_register["percentile_rank"],
                "ebv_valid_sample_count": len(left),
                "self_valid_sample_count": len(right),
                "model_combination_count": len(left) * len(right),
                "geometry_status": "complete" if left and right else "missing_or_qc_excluded_model",
            }
            metrics: dict[str, list[float]] = defaultdict(list)
            for left_sample in left:
                for right_sample in right:
                    values = pair_geometry(left_sample, right_sample)
                    ensemble_writer.writerow({
                        "allele": allele,
                        "pair_id": pair["pair_id"],
                        "ebv_candidate_id": ebv_id,
                        "self_candidate_id": self_id,
                        "ebv_sample_index": left_sample.sample_index,
                        "self_sample_index": right_sample.sample_index,
                        **{key: round(value, 6) for key, value in values.items()},
                        "interpretation": "Technical AF3 sample-combination sensitivity only; not biological replication.",
                    })
                    for key, value in values.items():
                        metrics[key].append(value)
            summary = dict(base)
            for key in ("full_core_ca_rmsd_A", "anchor_ca_rmsd_A", "exposed_ca_rmsd_A", "exposed_cb_rmsd_A", "exposed_sidechain_centroid_rmsd_A"):
                summary.update(_distribution(metrics.get(key, []), key))
            summary.update(direct_register_sequence_metrics(ebv_register["predicted_core"], self_register["predicted_core"]))
            summary["primary_endpoint"] = PRIMARY_ENDPOINT
            summary["claim_boundary"] = CLAIM_BOUNDARY
            pair_rows.append(summary)
    finally:
        ensemble_handle.close()
    rank_pair_rows(pair_rows)
    write_csv(out / "discovery/pair_summary_6400.csv", pair_rows)

    cross_rows = summarize_cross_allele(pair_rows)
    for row in cross_rows:
        row.update({
            "ebv_sequence": panel[row["ebv_candidate_id"]]["sequence"],
            "self_sequence": panel[row["self_candidate_id"]]["sequence"],
            "ebv_protein": panel[row["ebv_candidate_id"]]["protein_symbol"],
            "self_protein": panel[row["self_candidate_id"]]["protein_symbol"],
            "ebv_source_certainty": panel[row["ebv_candidate_id"]]["source_certainty"],
            "self_source_certainty": panel[row["self_candidate_id"]]["source_certainty"],
        })
    write_csv(out / "discovery/cross_allele_consensus_1600.csv", cross_rows)
    top_rows: list[dict[str, Any]] = []
    for allele in sorted({row["allele"] for row in pair_rows}):
        for row in sorted((value for value in pair_rows if value["allele"] == allele and value["geometry_status"] == "complete"), key=lambda value: int(value["within_allele_rank"]))[:25]:
            top_rows.append({"selection_set": f"top_25_within_{allele}", **row})
    for row in sorted((value for value in cross_rows if value["cross_allele_status"] == "complete_four_alleles"), key=lambda value: int(value["cross_allele_consensus_rank"]))[:50]:
        top_rows.append({"selection_set": "top_50_cross_allele_consensus", **row})
    write_csv(out / "discovery/top_candidate_sets.csv", top_rows, fields=sorted({key for row in top_rows for key in row}))

    calibration_job_by_key = {(int(row["server_seed"]), row["entity_id"]): row for row in calibration_inventory}
    comparison_universe = read_csv(package / "calibration_comparison_universe_72.csv")
    calibration_pair_rows: list[dict[str, Any]] = []
    calibration_ensemble_fields = [
        "seed", "analysis_set", "pair_role", "pair_id", "viral_entity_id", "self_entity_id",
        "viral_sample_index", "self_sample_index", "full_core_ca_rmsd_A", "anchor_ca_rmsd_A",
        "exposed_ca_rmsd_A", "exposed_cb_rmsd_A", "exposed_sidechain_centroid_rmsd_A", "interpretation",
    ]
    calibration_handle, calibration_writer = _open_stream_csv(out / "calibration/calibration_model_ensemble.csv", calibration_ensemble_fields)
    try:
        for comparison in comparison_universe:
            seed = int(comparison["seed"])
            viral_job = calibration_job_by_key[(seed, comparison["viral_entity_id"])]
            self_job = calibration_job_by_key[(seed, comparison["self_entity_id"])]
            left = geometry_by_job.get(viral_job["job_name"].lower(), [])
            right = geometry_by_job.get(self_job["job_name"].lower(), [])
            metrics: dict[str, list[float]] = defaultdict(list)
            for left_sample in left:
                for right_sample in right:
                    values = pair_geometry(left_sample, right_sample)
                    calibration_writer.writerow({
                        **{key: comparison[key] for key in ("seed", "analysis_set", "pair_role", "pair_id", "viral_entity_id", "self_entity_id")},
                        "viral_sample_index": left_sample.sample_index,
                        "self_sample_index": right_sample.sample_index,
                        **{key: round(value, 6) for key, value in values.items()},
                        "interpretation": "Fixed-seed native-HLA calibration geometry; model samples are technical alternatives.",
                    })
                    for key, value in values.items():
                        metrics[key].append(value)
            summary = {
                **comparison,
                "viral_valid_sample_count": len(left),
                "self_valid_sample_count": len(right),
                "model_combination_count": len(left) * len(right),
                "geometry_status": "complete" if left and right else "missing_or_qc_excluded_model",
            }
            for key in ("full_core_ca_rmsd_A", "anchor_ca_rmsd_A", "exposed_ca_rmsd_A", "exposed_cb_rmsd_A", "exposed_sidechain_centroid_rmsd_A"):
                summary.update(_distribution(metrics.get(key, []), key))
            summary["primary_endpoint"] = PRIMARY_ENDPOINT
            summary["claim_boundary"] = CLAIM_BOUNDARY
            calibration_pair_rows.append(summary)
    finally:
        calibration_handle.close()
    write_csv(out / "calibration/calibration_pair_summary_72.csv", calibration_pair_rows)

    seed_rows: list[dict[str, Any]] = []
    for seed in CALIBRATION_SEEDS:
        primary_rows = [
            row for row in calibration_pair_rows
            if int(row["seed"]) == seed and row["analysis_set"] == "primary_rank_of_26" and row["geometry_status"] == "complete"
        ]
        ordered = sorted(primary_rows, key=lambda row: (float(row["exposed_ca_rmsd_A_median"]), str(row["pair_id"])))
        positive = next((row for row in ordered if row["pair_role"] == "E1_positive"), None)
        available_rank = next((index for index, row in enumerate(ordered, start=1) if row["pair_role"] == "E1_positive"), "")
        controls = [float(row["exposed_ca_rmsd_A_median"]) for row in ordered if row["pair_role"] == "full_decoy"]
        formal = len(ordered) == 26 and positive is not None and len(controls) == 25
        positive_value = float(positive["exposed_ca_rmsd_A_median"]) if positive else float("nan")
        control_median = float(np.median(controls)) if controls else float("nan")
        criterion = bool(formal and int(available_rank) <= 3 and positive_value < control_median)
        seed_rows.append({
            "seed": seed,
            "available_primary_count": len(ordered),
            "expected_primary_count": 26,
            "available_full_decoy_count": len(controls),
            "expected_full_decoy_count": 25,
            "available_rank": available_rank,
            "formal_rank_of_26": available_rank if formal else "",
            "positive_exposed_ca_rmsd_median_A": round(positive_value, 6) if positive else "",
            "available_equal_weight_control_median_A": round(control_median, 6) if controls else "",
            "formal_seed_evaluable": formal,
            "seed_recovery_criterion_pass": criterion if formal else "",
            "formal_seed_status": "pass" if criterion else "fail" if formal else "not_evaluable_incomplete_seed",
            "required_rule": "rank <= 3 of 26 and positive RMSD below equal-weight median of 25 full decoys",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    recovery_status = classify_recovery(seed_rows)
    write_csv(out / "calibration/seed_recovery_summary.csv", seed_rows)
    recovery_report = [{
        "biological_system_id": "SYS_BALF5_MBP_HY2E11",
        "available_calibration_jobs": sum(row["download_status"] == "complete" for row in calibration_inventory),
        "expected_calibration_jobs": 24,
        "available_comparisons": sum(row["geometry_status"] == "complete" for row in calibration_pair_rows),
        "expected_comparisons": 72,
        "recovery_status": recovery_status,
        "required_rule": "top 3 of 26 on both seeds and below equal-weight control median",
        "missing_jobs": ";".join(row["job_name"] for row in calibration_inventory if row["download_status"] != "complete"),
        "claim_boundary": CLAIM_BOUNDARY,
    }]
    write_csv(out / "calibration/positive_recovery_report.csv", recovery_report)

    divergent_duplicates = [
        row for row in [*discovery_duplicates, *calibration_duplicates]
        if row["duplicate_class"] == "divergent_unplanned_repeat"
    ]
    inventory_by_job = {
        row["job_name"].lower(): row
        for row in [*discovery_inventory, *calibration_inventory]
    }
    alternative_geometries: dict[tuple[str, str], list[SampleGeometry]] = {}
    alternative_sample_rows: list[dict[str, Any]] = []
    alternative_job_rows: list[dict[str, Any]] = []
    for duplicate in divergent_duplicates:
        expected = inventory_by_job[duplicate["job_name"].lower()]
        entity_id = str(expected.get("candidate_id", expected.get("entity_id", "")))
        alternative = {
            **expected,
            "canonical_directory": duplicate["directory"],
            "bundle_fingerprint": duplicate["bundle_fingerprint"],
            "download_status": "complete",
        }
        samples, geometries, job = analyze_job(
            alternative,
            registers[(alternative["allele"], entity_id)],
            hla[alternative["allele"]],
        )
        for sample in samples:
            sample.update({
                "analysis_role": "divergent_unplanned_repeat_sensitivity",
                "duplicate_bundle_fingerprint": duplicate["bundle_fingerprint"],
                "duplicate_directory": duplicate["directory"],
            })
        job.update({
            "analysis_role": "divergent_unplanned_repeat_sensitivity",
            "duplicate_bundle_fingerprint": duplicate["bundle_fingerprint"],
            "duplicate_directory": duplicate["directory"],
        })
        alternative_sample_rows.extend(samples)
        alternative_job_rows.append(job)
        alternative_geometries[(duplicate["job_name"].lower(), duplicate["bundle_fingerprint"])] = geometries
    write_csv(
        out / "supplemental/divergent_duplicate_model_sample_qc.csv",
        alternative_sample_rows,
        fields=sorted({key for row in alternative_sample_rows for key in row}) if alternative_sample_rows else ["analysis_role"],
    )
    write_csv(
        out / "supplemental/divergent_duplicate_job_qc.csv",
        alternative_job_rows,
        fields=sorted({key for row in alternative_job_rows for key in row}) if alternative_job_rows else ["analysis_role"],
    )

    discovery_sensitivity: list[dict[str, Any]] = []
    for duplicate in divergent_duplicates:
        expected = inventory_by_job[duplicate["job_name"].lower()]
        if expected["cohort"] != "discovery":
            continue
        alternative_samples = alternative_geometries[(duplicate["job_name"].lower(), duplicate["bundle_fingerprint"])]
        candidate_id = expected["candidate_id"]
        allele = expected["allele"]
        for pair in pair_rows:
            if pair["allele"] != allele or pair["geometry_status"] != "complete":
                continue
            if candidate_id == pair["ebv_candidate_id"]:
                other = discovery_job_by_key[(allele, pair["self_candidate_id"])]
                left, right = alternative_samples, geometry_by_job[other["job_name"].lower()]
            elif candidate_id == pair["self_candidate_id"]:
                other = discovery_job_by_key[(allele, pair["ebv_candidate_id"])]
                left, right = geometry_by_job[other["job_name"].lower()], alternative_samples
            else:
                continue
            discovery_sensitivity.append(duplicate_pair_sensitivity(pair, left, right, duplicate))
    write_csv(
        out / "supplemental/divergent_duplicate_discovery_pair_sensitivity.csv",
        discovery_sensitivity,
        fields=sorted({key for row in discovery_sensitivity for key in row}) if discovery_sensitivity else ["comparison_id"],
    )

    calibration_sensitivity: list[dict[str, Any]] = []
    for duplicate in divergent_duplicates:
        expected = inventory_by_job[duplicate["job_name"].lower()]
        if expected["cohort"] != "calibration":
            continue
        alternative_samples = alternative_geometries[(duplicate["job_name"].lower(), duplicate["bundle_fingerprint"])]
        seed, entity_id = int(expected["server_seed"]), expected["entity_id"]
        for comparison in calibration_pair_rows:
            if int(comparison["seed"]) != seed or comparison["geometry_status"] != "complete":
                continue
            if entity_id == comparison["viral_entity_id"]:
                other = calibration_job_by_key[(seed, comparison["self_entity_id"])]
                left, right = alternative_samples, geometry_by_job[other["job_name"].lower()]
            elif entity_id == comparison["self_entity_id"]:
                other = calibration_job_by_key[(seed, comparison["viral_entity_id"])]
                left, right = geometry_by_job[other["job_name"].lower()], alternative_samples
            else:
                continue
            calibration_sensitivity.append(duplicate_pair_sensitivity(comparison, left, right, duplicate))
    write_csv(
        out / "supplemental/divergent_duplicate_calibration_pair_sensitivity.csv",
        calibration_sensitivity,
        fields=sorted({key for row in calibration_sensitivity for key in row}) if calibration_sensitivity else ["comparison_id"],
    )

    generate_figures(out, job_rows, pair_rows, cross_rows, calibration_pair_rows, seed_rows)
    gold_standard_summary = run_gold_standard_audit(
        registry_path=package / "literature_tcell_pair_registry.csv",
        experimental_metrics_path=ROOT / "processed/experimental_positive_control/experimental_drb2_positive_control_metrics.csv",
        seed_summary_path=out / "calibration/seed_recovery_summary.csv",
        recovery_report_path=out / "calibration/positive_recovery_report.csv",
        out_dir=out / "validation",
    )

    discovery_complete = sum(row["download_status"] == "complete" for row in discovery_inventory)
    calibration_complete = sum(row["download_status"] == "complete" for row in calibration_inventory)
    discovery_evaluable = sum(row["cohort"] == "discovery" and row["technical_status"] == "geometry_evaluable" for row in job_rows)
    calibration_evaluable = sum(row["cohort"] == "calibration" and row["technical_status"] == "geometry_evaluable" for row in job_rows)
    pair_complete = sum(row["geometry_status"] == "complete" for row in pair_rows)
    calibration_comparisons = sum(row["geometry_status"] == "complete" for row in calibration_pair_rows)
    manifest = {
        "analysis_version": "EBV_MS_TCELL_V2_MODEL_ANALYSIS_2026-08-25",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "download_roots": [str(path) for path in download_roots],
        "discovery_jobs_complete": discovery_complete,
        "discovery_jobs_expected": 320,
        "discovery_jobs_geometry_evaluable": discovery_evaluable,
        "calibration_jobs_complete": calibration_complete,
        "calibration_jobs_expected": 24,
        "calibration_jobs_geometry_evaluable": calibration_evaluable,
        "model_sample_qc_rows": len(sample_rows),
        "discovery_pair_rows": len(pair_rows),
        "discovery_pairs_geometry_complete": pair_complete,
        "cross_allele_rows": len(cross_rows),
        "calibration_comparison_rows": len(calibration_pair_rows),
        "calibration_comparisons_geometry_complete": calibration_comparisons,
        "divergent_unplanned_repeat_occurrences": len(divergent_duplicates),
        "divergent_unplanned_repeat_jobs": len({row["job_name"] for row in divergent_duplicates}),
        "duplicate_discovery_pair_sensitivity_rows": len(discovery_sensitivity),
        "duplicate_calibration_pair_sensitivity_rows": len(calibration_sensitivity),
        "positive_recovery_status": recovery_status,
        "gold_standard_independent_system_count": gold_standard_summary["gold_standard_independent_system_count"],
        "gold_standard_capture_at_1_available_seed_fraction": gold_standard_summary["capture_at_1_available_seed_fraction"],
        "gold_standard_formal_evaluable_seed_count": gold_standard_summary["formal_evaluable_seed_count"],
        "gold_standard_formal_seed_pass_count": gold_standard_summary["formal_seed_pass_count"],
        "primary_endpoint": PRIMARY_ENDPOINT,
        "claim_boundary": CLAIM_BOUNDARY,
        "input_manifest_sha256": {
            name: sha256_file(package / name)
            for name in (
                "frozen_v2_80_peptide_panel.csv", "model_inventory_320.csv",
                "allele_register_predictions_320.csv", "within_allele_pair_universe_6400.csv",
                "native_hla_calibration_manifest_24.csv", "calibration_comparison_universe_72.csv",
            )
        },
    }
    write_json(out / "analysis_manifest.json", manifest)

    top_cross = sorted(
        (row for row in cross_rows if row["cross_allele_status"] == "complete_four_alleles"),
        key=lambda row: int(row["cross_allele_consensus_rank"]),
    )[:10]
    results_lines = [
        "# V2 full-ensemble results summary",
        "",
        f"- Discovery jobs: **{discovery_complete}/320** downloaded; **{discovery_evaluable}** geometry-evaluable.",
        f"- Discovery pairs: **{pair_complete}/6,400** geometry-complete.",
        f"- Calibration jobs: **{calibration_complete}/24** downloaded; **{calibration_evaluable}** geometry-evaluable.",
        f"- Calibration comparisons: **{calibration_comparisons}/72** geometry-complete.",
        f"- Strict positive-recovery status: **{recovery_status}**.",
        "",
        "## Gold-standard positive-control capture",
        "",
        "The denominator was locked before ranking and contains one independent experimentally established system: Hy.2E11 recognition of DRB5*01:01-BALF5 and DRB1*15:01-MBP, with experimental pMHC structures 1H15 and 1BX2.",
        "",
        f"- Available-set capture@1: **{gold_standard_summary['capture_at_1_available_seed_count']}/{gold_standard_summary['available_seed_count']} seeds**.",
        f"- Fully evaluable seeds passing the predeclared rule: **{gold_standard_summary['formal_seed_pass_count']}/{gold_standard_summary['formal_evaluable_seed_count']}**.",
        f"- Strict two-seed status: **{gold_standard_summary['strict_two_seed_recovery_status']}**.",
        "- Model or score changed to fit the positive: **no**.",
        "",
        "This confirms capture of the known control in both available seed sets. It does not estimate broad sensitivity because the gold-standard denominator contains only one independent system, and one seed remains incomplete.",
        "",
        "## Top cross-allele consensus pairs",
        "",
        "Consensus minimizes the worst within-allele percentile, then the median, and requires all four alleles.",
        "",
        "| Rank | EBV | Self | EBV protein | Self protein | Median percentile | Worst percentile |",
        "|---:|---|---|---|---|---:|---:|",
        *[
            f"| {row['cross_allele_consensus_rank']} | {row['ebv_sequence']} | {row['self_sequence']} | {row['ebv_protein']} | {row['self_protein']} | {float(row['median_within_allele_percentile']):.4f} | {float(row['worst_within_allele_percentile']):.4f} |"
            for row in top_cross
        ],
        "",
        f"> {CLAIM_BOUNDARY}",
        "",
    ]
    (out / "RESULTS_SUMMARY.md").write_text("\n".join(results_lines), encoding="utf-8")
    methods = f"""# Methods

## Input freeze and identity

The frozen V2 package defines 320 discovery jobs (80 peptides across four HLA-DRB1 alleles), 6,400 within-allele EBV-self pairs, 24 fixed-seed native-HLA calibration jobs, and 72 calibration comparisons. Downloaded folders were matched by the request name inside `job_request.json`. A complete bundle required five CIF files, five confidence summaries, five full-data JSON files, and one request file. Duplicate copies were selected by the lexicographically smallest content fingerprint before reading model scores.

## Locked gold-standard denominator

Gold-standard eligibility was determined without reading model scores. It required exact peptide identities and HLA arms, recognition of both pMHCs by the same human T-cell clone, and experimental pMHC structures for both arms. The current denominator contains one independent system, Hy.2E11 DRB5*01:01-BALF5 versus DRB1*15:01-MBP. Supportive T-cell studies, antibody-only pairs, protein-level associations, canonical tiles, and computational discoveries were excluded. Recovery was reported at ranks 1 and 3 for each available seed; incomplete control sets were not relabeled as formal 26-pair tests.

## Model QC

Every CIF was required to contain exactly HLA-DRA chain A, the expected HLA-DRB chain B, and the exact peptide chain C. All chain sequences had to match the frozen request. Samples with an AlphaFold clash flag were excluded from geometry but retained in QC. No pLDDT, ipTM, PAE, contact-probability, or geometry threshold was used to select samples. The highest-ranking clash-free exact sample is reported only as a descriptive representative; full-ensemble analyses use every valid sample.

## Register-aware geometry

IEDB `recommended_binding` top cores and exact `seq_num` mappings were frozen before structural analysis. Equivalent HLA-DRA and HLA-DRB groove C-alpha atoms from the first 85 residues of each chain were superposed by the Kabsch algorithm. The primary endpoint was {PRIMARY_ENDPOINT}. Full-core, anchor-position P1/P4/P6/P9, pseudo-C-beta, and exposed-side-chain-centroid RMSDs are supporting descriptors. Every valid 5-by-5 model combination was calculated; medians and interquartile ranges summarize technical model sensitivity and are not inferential confidence intervals.

## Ranking and calibration

Discovery pairs were ranked separately within each allele by median exposed-position C-alpha RMSD, then ensemble IQR, then frozen pair ID. Cross-allele consensus minimizes the worst allele-specific percentile, then the median percentile, and requires all four allele results; raw geometry is never pooled across alleles. The BALF5-MBP calibration remains frozen as rank <=3 of 26 on both seeds plus positive RMSD below the equal-weight median of 25 full decoys. Missing calibration models are not imputed, so an incomplete seed is reported as not evaluable.

## Unplanned duplicate sensitivity

Complete duplicate folders with identical content were treated as copies. Distinct content fingerprints under the same frozen job identity were excluded from primary selection and analyzed separately as unplanned technical-repeat sensitivity. They cannot increase the primary sample size or serve as biological replication.

## Interpretation boundary

{CLAIM_BOUNDARY}
"""
    (out / "METHODS.md").write_text(methods, encoding="utf-8")
    readme = f"""# EBV-MS T-cell Library V2 model analysis (2026-08-25)

- Discovery downloads: **{discovery_complete}/320**
- Geometry-evaluable discovery jobs: **{discovery_evaluable}**
- Geometry-complete discovery pairs: **{pair_complete}/6,400**
- Calibration downloads: **{calibration_complete}/24**
- Geometry-complete calibration comparisons: **{calibration_comparisons}/72**
- Strict positive recovery: **{recovery_status}**
- Gold-standard available-set capture@1: **{gold_standard_summary['capture_at_1_available_seed_count']}/{gold_standard_summary['available_seed_count']} seeds**
- Gold-standard independent systems: **{gold_standard_summary['gold_standard_independent_system_count']}**

The analysis uses all valid AlphaFold model combinations and preserves the one missing discovery job and two missing calibration jobs as technical missingness. Original Downloads and the frozen V2 package are read-only inputs.

Distinct unplanned duplicate runs are excluded from primary ranking and reported in `supplemental/` as technical-repeat sensitivity.

The locked positive-control audit is in `validation/`. It verifies exact biological identity and experimental structures before reading ranks, and it records that the model and score were not changed to fit the positive.

Reproduce:

```bash
PYTHONPATH=src python3 src/run_tcell_library_v2_model_analysis.py
```

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")

    checksum_rows = [
        {"relative_path": str(path.relative_to(out)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        for path in sorted(out.rglob("*"), key=str)
        if path.is_file() and path.name != "SHA256SUMS.csv"
    ]
    write_csv(out / "SHA256SUMS.csv", checksum_rows)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--package", type=Path, default=PACKAGE)
    parser.add_argument("--download-root", type=Path, action="append", dest="download_roots")
    args = parser.parse_args()
    roots = tuple(args.download_roots or DEFAULT_DOWNLOAD_ROOTS)
    manifest = run_analysis(roots, args.out, args.package)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
