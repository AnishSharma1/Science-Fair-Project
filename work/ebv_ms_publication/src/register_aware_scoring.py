"""Deterministic sequence/chemistry scores at matched HLA-II register positions.

These helpers score an existing local alignment only after retaining residue
pairs that occupy the same P1--P9 index in the recorded peptide cores.  They
are computational pMHC descriptors, not predictions of TCR recognition or
cross-reactivity.
"""

from __future__ import annotations

import re
from statistics import mean
from typing import Any


# Same simple descriptor conventions used by biochemical_similarity.py.
KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5,
    "N": -3.5, "K": -3.9, "R": -4.5,
}
KD_MIN, KD_MAX = min(KD.values()), max(KD.values())
AROMATIC = set("FWY")
CHARGE = {aa: (1 if aa in "KR" else -1 if aa in "DE" else 0) for aa in KD}
SIZE = {aa: (0 if aa in "AGSTCP" else 2 if aa in "FWYH" else 1) for aa in KD}

ANCHOR_POSITIONS = {1, 4, 6, 9}
CANDIDATE_EXPOSED_POSITIONS = {2, 3, 5, 7, 8}
_COORDINATE = re.compile(r"(\d+)([A-Z]?):(\d+)([A-Z]?)$")


def parse_local_alignment_positions(text: str) -> list[tuple[int, str, int, str]]:
    """Parse ``4Y:3H;7F:6F`` alignment coordinates with strict validation."""
    if not text:
        return []
    pairs: list[tuple[int, str, int, str]] = []
    for field in text.split(";"):
        match = _COORDINATE.fullmatch(field.strip())
        if not match:
            raise ValueError(f"invalid local-alignment coordinate: {field!r}")
        ebv_index, ebv_residue, human_index, human_residue = match.groups()
        pairs.append((int(ebv_index), ebv_residue, int(human_index), human_residue))
    return pairs


def register_position(peptide_index_1_based: int, core_start_1_based: int) -> int | None:
    """Return P1--P9 index for a peptide coordinate, or ``None`` outside the core."""
    position = peptide_index_1_based - core_start_1_based + 1
    return position if 1 <= position <= 9 else None


def property_similarity(ebv_residue: str, human_residue: str) -> dict[str, float]:
    """Return four descriptor similarities, each on the inclusive 0--1 scale.

    Hydropathy is ``1 - abs(normalized KD difference)``; charge and coarse size
    are one minus their category distance divided by two; aromaticity is an
    exact shared-category indicator.
    """
    for residue in (ebv_residue, human_residue):
        if residue not in KD:
            raise ValueError(f"unsupported amino-acid residue: {residue!r}")
    ebv_hydro = (KD[ebv_residue] - KD_MIN) / (KD_MAX - KD_MIN)
    human_hydro = (KD[human_residue] - KD_MIN) / (KD_MAX - KD_MIN)
    return {
        "hydrophobicity_similarity": 1 - abs(ebv_hydro - human_hydro),
        "charge_similarity": 1 - abs(CHARGE[ebv_residue] - CHARGE[human_residue]) / 2,
        "aromatic_similarity": float(
            (ebv_residue in AROMATIC) == (human_residue in AROMATIC)
        ),
        "size_similarity": 1 - abs(SIZE[ebv_residue] - SIZE[human_residue]) / 2,
    }


def _class_summary(pairs: list[tuple[int, str, str]]) -> dict[str, object]:
    if not pairs:
        return {
            "positions": "",
            "count": 0,
            "identity_count": 0,
            "identity_fraction": "",
            "hydrophobicity_similarity": "",
            "charge_similarity": "",
            "aromatic_similarity": "",
            "size_similarity": "",
            "property_similarity": "",
        }
    properties = [property_similarity(ebv, human) for _, ebv, human in pairs]
    component_means = {
        name: round(mean(item[name] for item in properties), 6)
        for name in properties[0]
    }
    return {
        "positions": ";".join(f"P{position}" for position, _, _ in pairs),
        "count": len(pairs),
        "identity_count": sum(ebv == human for _, ebv, human in pairs),
        "identity_fraction": round(sum(ebv == human for _, ebv, human in pairs) / len(pairs), 6),
        **component_means,
        "property_similarity": round(mean(component_means.values()), 6),
    }


def _core_start(row: dict[str, Any], key: str) -> int:
    try:
        start = int(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing valid {key}") from exc
    if start < 1:
        raise ValueError(f"{key} must be one-based positive")
    return start


def score_same_register_alignment(row: dict[str, Any]) -> dict[str, object]:
    """Score only aligned pairs occupying the same P1--P9 register position."""
    if row.get("register_assessment") != "assessable_register_hypothesis":
        return {
            "score_coverage_status": "excluded_nonprimary_register_status",
            "same_register_alignment_count": 0,
            "same_register_positions": "",
            "all_same_register_positions": "",
            "anchor_same_register_positions": "",
            "candidate_exposed_same_register_positions": "",
        }

    ebv_peptide = str(row["ebv_peptide"])
    human_peptide = str(row["human_peptide"])
    ebv_core_start = _core_start(row, "ebv_top_core_start_1_based")
    human_core_start = _core_start(row, "human_top_core_start_1_based")
    retained: list[tuple[int, str, str]] = []
    for ebv_index, ebv_letter, human_index, human_letter in parse_local_alignment_positions(
        str(row.get("original_local_alignment_coordinates", ""))
    ):
        if not (1 <= ebv_index <= len(ebv_peptide) and 1 <= human_index <= len(human_peptide)):
            raise ValueError("local-alignment coordinate is outside its peptide")
        actual_ebv, actual_human = ebv_peptide[ebv_index - 1], human_peptide[human_index - 1]
        if ebv_letter and ebv_letter != actual_ebv:
            raise ValueError("EBV coordinate residue does not match peptide")
        if human_letter and human_letter != actual_human:
            raise ValueError("human coordinate residue does not match peptide")
        ebv_position = register_position(ebv_index, ebv_core_start)
        human_position = register_position(human_index, human_core_start)
        if ebv_position is not None and ebv_position == human_position:
            retained.append((ebv_position, actual_ebv, actual_human))

    retained.sort(key=lambda item: item[0])
    all_summary = _class_summary(retained)
    anchor_summary = _class_summary([item for item in retained if item[0] in ANCHOR_POSITIONS])
    exposed_summary = _class_summary(
        [item for item in retained if item[0] in CANDIDATE_EXPOSED_POSITIONS]
    )
    count = len(retained)
    status = (
        "no_same_register_alignment" if count == 0
        else "limited_coverage_report_only" if count < 3
        else "robust_primary_ranking_eligible"
    )
    return {
        "score_coverage_status": status,
        "same_register_alignment_count": count,
        "same_register_positions": all_summary["positions"],
        "all_same_register_positions": all_summary["positions"],
        "same_register_identity_count": all_summary["identity_count"],
        "same_register_identity_fraction": all_summary["identity_fraction"],
        "same_register_hydrophobicity_similarity": all_summary["hydrophobicity_similarity"],
        "same_register_charge_similarity": all_summary["charge_similarity"],
        "same_register_aromatic_similarity": all_summary["aromatic_similarity"],
        "same_register_size_similarity": all_summary["size_similarity"],
        "same_register_property_similarity": all_summary["property_similarity"],
        "anchor_same_register_positions": anchor_summary["positions"],
        "anchor_same_register_alignment_count": anchor_summary["count"],
        "anchor_same_register_identity_fraction": anchor_summary["identity_fraction"],
        "anchor_same_register_property_similarity": anchor_summary["property_similarity"],
        "candidate_exposed_same_register_positions": exposed_summary["positions"],
        "candidate_exposed_same_register_alignment_count": exposed_summary["count"],
        "candidate_exposed_same_register_identity_fraction": exposed_summary["identity_fraction"],
        "candidate_exposed_same_register_property_similarity": exposed_summary["property_similarity"],
    }
