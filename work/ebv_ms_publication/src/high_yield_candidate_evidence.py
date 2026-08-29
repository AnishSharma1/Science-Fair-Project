"""Pure evidence logic for the eight-candidate HLA-II dossier.

The functions in this module never fetch remote data and never convert missing
evidence into support. They classify cached records, compare independent
predictors, compute sequence-rarity diagnostics, and apply the locked assay
gates.
"""

from __future__ import annotations

import re
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np

from build_same_register_hla_rankings_v2 import sequence_metrics
from hla2_positive_control_benchmark import (
    AROMATIC,
    FORMAL_CHARGE,
    H_BOND_ACCEPTOR,
    H_BOND_DONOR,
    HYDROPATHY,
)
from hla2_positive_control_benchmark_v2 import (
    TCR_FACING_INDICES,
    _BLOSUM62,
    _BLOSUM62_ALPHABET,
)


CLAIM_BOUNDARY = (
    "This dossier prioritizes peptide-HLA binding and register experiments. It does "
    "not establish natural presentation, TCR recognition, specificity, "
    "cross-reactivity, molecular mimicry, probability, false-discovery rate, or an "
    "MS mechanism."
)

INDEPENDENT_NET_PREDICTOR = "netmhciipan_4_3_el"
INDEPENDENT_MIX_PREDICTOR = "mixmhc2pred_2_1_context"


def _text(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item) for item in value if item is not None)
    return "" if value is None else str(value)


def normalize_hla(value: Any) -> str:
    """Normalize common HLA names without changing the biological allotype."""
    text = _text(value).strip().upper().replace(" ", "")
    if not text:
        return ""
    text = text.replace("_", "")
    text = text.removeprefix("HLA-")
    match = re.fullmatch(r"(DRB[1-9]|DQA1|DQB1|DPA1|DPB1|[ABC])\*?(\d{2})(?::?)(\d{2})", text)
    if match:
        locus, group, protein = match.groups()
        return f"HLA-{locus}*{group}:{protein}"
    return f"HLA-{text}" if text.startswith(("DR", "DQ", "DP", "A*", "B*", "C*")) else text


def sequence_relation(query_sequence: Any, observed_sequence: Any) -> str:
    query = _text(query_sequence).strip().upper()
    observed = _text(observed_sequence).strip().upper()
    if not query or not observed:
        return "none"
    if query == observed:
        return "exact"
    if query in observed or observed in query:
        return "overlap"
    return "none"


def _is_human(value: Any) -> bool:
    text = _text(value).lower()
    return not text or "homo sapiens" in text or text.strip() in {"human", "9606"}


def _is_class_i(hla: str, mhc_class: str = "") -> bool:
    normalized_class = mhc_class.strip().upper().replace(" ", "")
    return normalized_class in {"I", "CLASSI", "MHC-I"} or bool(
        re.match(r"^HLA-[ABC]\*", hla)
    )


def classify_assay_evidence(
    record: Mapping[str, Any],
    target_sequence: str,
    target_hla: str,
) -> str:
    """Classify an IEDB-like assay record against one exact modeled arm."""
    observed = (
        record.get("epitope_sequence")
        or record.get("linear_sequence")
        or record.get("epitope__name")
        or ""
    )
    relation = sequence_relation(target_sequence, observed)
    hla = normalize_hla(
        record.get("mhc_allele")
        or record.get("mhc_allele_name")
        or record.get("mhc_restriction__name")
        or ""
    )
    mhc_class = _text(record.get("mhc_class") or record.get("mhc_restriction__class"))
    host = record.get("host_organism") or record.get("host_organism_name") or record.get("host__name")
    if not _is_human(host) or hla.startswith("H-2"):
        return "nonhuman"
    if _is_class_i(hla, mhc_class):
        return "class_i"
    if not hla:
        return "untyped"
    if hla == normalize_hla(target_hla):
        if relation == "exact":
            return "exact_sequence_exact_hla"
        if relation == "overlap":
            return "overlap_exact_hla"
        return "none"
    if hla.startswith("HLA-DR") or hla.startswith("HLA-DQ") or hla.startswith("HLA-DP"):
        return "other_human_hla"
    return "untyped"


def _number(value: Any) -> Optional[float]:
    try:
        if value in (None, "", "NA", "N/A", "nan"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def summarize_predictor_evidence(
    records: Sequence[Mapping[str, Any]],
    declared_core: str,
) -> dict[str, Any]:
    """Summarize the locked independent predictor pair for one peptide arm."""
    by_name = {str(row.get("predictor", "")): row for row in records}
    net = by_name.get(INDEPENDENT_NET_PREDICTOR)
    mix = by_name.get(INDEPENDENT_MIX_PREDICTOR)
    if not net or not mix:
        return {
            "predictor_status": "not_evaluable",
            "binding_consensus": False,
            "binding_supported": False,
            "both_predictors_above_20": False,
            "register_consensus": False,
            "register_consensus_core": "",
            "register_consensus_matches_declared": False,
            "netmhciipan_el_percentile": _number(net.get("percentile_rank")) if net else None,
            "mixmhc2pred_context_percentile": _number(mix.get("percentile_rank")) if mix else None,
            "netmhciipan_el_core": _text(net.get("core")).upper() if net else "",
            "mixmhc2pred_context_core": _text(mix.get("core")).upper() if mix else "",
            "mixmhc2pred_reverse_orientation": bool(mix and str(mix.get("orientation", "1")) == "-1"),
        }
    net_rank = _number(net.get("percentile_rank"))
    mix_rank = _number(mix.get("percentile_rank"))
    net_core = _text(net.get("core")).strip().upper()
    mix_core = _text(mix.get("core")).strip().upper()
    if net_rank is None or mix_rank is None or len(net_core) != 9 or len(mix_core) != 9:
        status = "not_evaluable"
    else:
        status = "complete"
    binding_consensus = status == "complete" and net_rank <= 5.0 and mix_rank <= 5.0
    binding_supported = status == "complete" and (
        min(net_rank, mix_rank) <= 5.0 or max(net_rank, mix_rank) <= 20.0
    )
    both_above_20 = status == "complete" and net_rank > 20.0 and mix_rank > 20.0
    register_consensus = status == "complete" and net_core == mix_core
    consensus_core = net_core if register_consensus else ""
    return {
        "predictor_status": status,
        "binding_consensus": binding_consensus,
        "binding_supported": binding_supported,
        "both_predictors_above_20": both_above_20,
        "register_consensus": register_consensus,
        "register_consensus_core": consensus_core,
        "register_consensus_matches_declared": (
            register_consensus and consensus_core == declared_core.strip().upper()
        ),
        "netmhciipan_el_percentile": net_rank,
        "mixmhc2pred_context_percentile": mix_rank,
        "netmhciipan_el_core": net_core,
        "mixmhc2pred_context_core": mix_core,
        "mixmhc2pred_reverse_orientation": str(mix.get("orientation", "1")) == "-1",
    }


def classify_ligand_hit(
    target_sequence: str,
    observed_sequence: str,
    *,
    exact_hla: bool,
    monoallelic: bool,
    target_core: str = "",
) -> str:
    relation = sequence_relation(target_sequence, observed_sequence)
    if relation == "none" and target_core and sequence_relation(target_core, observed_sequence) != "none":
        relation = "core_overlap"
    if relation == "none":
        return "not_a_hit"
    if relation == "core_overlap":
        if exact_hla and monoallelic:
            return "core_overlap_monoallelic_exact_hla"
        if exact_hla:
            return "core_overlap_multiallelic_compatible"
        return "core_overlap_other_or_untyped_hla"
    if exact_hla and monoallelic and relation == "exact":
        return "exact_sequence_monoallelic_exact_hla"
    if exact_hla and monoallelic:
        return "nested_overlap_monoallelic_exact_hla"
    if exact_hla and relation == "exact":
        return "exact_sequence_multiallelic_compatible"
    if exact_hla:
        return "nested_overlap_multiallelic_compatible"
    return f"{relation}_sequence_other_or_untyped_hla"


def iter_nine_mer_windows(records: Iterable[Mapping[str, Any]]) -> Iterable[dict[str, Any]]:
    for record in records:
        sequence = _text(record.get("sequence")).strip().upper()
        for start in range(max(0, len(sequence) - 8)):
            core = sequence[start : start + 9]
            if len(core) == 9 and all(residue in "ACDEFGHIKLMNPQRSTVWY" for residue in core):
                yield {
                    "accession": _text(record.get("accession")),
                    "protein": _text(record.get("protein")),
                    "start_1_based": start + 1,
                    "core": core,
                }


def _feature_key(metrics: Mapping[str, Any]) -> tuple[float, float, float]:
    return (
        -float(metrics["tcr_facing_blosum62_similarity"]),
        float(metrics["tcr_face_physicochemical_mismatch"]),
        -float(metrics["tcr_facing_sequence_identity"]),
    )


def scan_similarity_rarity(
    *,
    query_core: str,
    paired_core: str,
    database_records: Sequence[Mapping[str, Any]],
    exclude_accession: str = "",
    exclude_core: str = "",
    top_n: int = 100,
) -> dict[str, Any]:
    """Compare a frozen pair with all unrelated nine-mers in a sequence database."""
    target_metrics = sequence_metrics(query_core, paired_core)
    target_key = _feature_key(target_metrics)
    scored: list[dict[str, Any]] = []
    excluded = 0
    for window in iter_nine_mer_windows(database_records):
        if (
            exclude_accession
            and window["accession"] == exclude_accession
            and window["core"] == exclude_core
        ):
            excluded += 1
            continue
        metrics = sequence_metrics(query_core, window["core"])
        scored.append({**window, **metrics})
    scored.sort(
        key=lambda row: (
            *_feature_key(row),
            str(row["accession"]),
            int(row["start_1_based"]),
            str(row["core"]),
        )
    )
    at_least_as_good = sum(_feature_key(row) <= target_key for row in scored)
    count = len(scored)
    percentile = round(100.0 * at_least_as_good / count, 12) if count else None
    neighbors = []
    for rank, row in enumerate(scored[:top_n], start=1):
        neighbors.append({**row, "neighbor_rank": rank})
    return {
        "query_core": query_core,
        "paired_core": paired_core,
        "evaluated_window_count": count,
        "excluded_target_window_count": excluded,
        "at_least_as_good_count": at_least_as_good,
        "empirical_percentile": percentile,
        "target_metrics": target_metrics,
        "nearest_neighbors": neighbors,
    }


def scan_similarity_rarity_fast(
    *,
    query_core: str,
    paired_core: str,
    database_records: Sequence[Mapping[str, Any]],
    exclude_accession: str = "",
    exclude_core: str = "",
    top_n: int = 100,
) -> dict[str, Any]:
    """Vectorized equivalent of :func:`scan_similarity_rarity` for proteomes."""
    alphabet = _BLOSUM62_ALPHABET
    aa_to_index = {aa: index for index, aa in enumerate(alphabet)}
    blosum = np.asarray(
        [[_BLOSUM62[(left, right)] for right in alphabet] for left in alphabet],
        dtype=np.float64,
    )
    diagonal = np.diag(blosum)
    descriptors = np.asarray(
        [
            [
                (FORMAL_CHARGE[aa] + 1.0) / 2.0,
                (HYDROPATHY[aa] + 4.5) / 9.0,
                H_BOND_DONOR[aa],
                H_BOND_ACCEPTOR[aa],
                AROMATIC[aa],
            ]
            for aa in alphabet
        ],
        dtype=np.float64,
    )
    exposed = np.asarray(TCR_FACING_INDICES, dtype=np.int64)
    query_indices = np.asarray([aa_to_index[query_core[index]] for index in exposed])
    target_metrics = sequence_metrics(query_core, paired_core)
    target_blosum = float(target_metrics["tcr_facing_blosum62_similarity"])
    target_physiochemical = float(target_metrics["tcr_face_physicochemical_mismatch"])
    target_identity = float(target_metrics["tcr_facing_sequence_identity"])
    evaluated = 0
    excluded = 0
    at_least_as_good = 0
    nearest: list[dict[str, Any]] = []

    for record in database_records:
        sequence = _text(record.get("sequence")).strip().upper()
        if len(sequence) < 9:
            continue
        try:
            encoded = np.fromiter((aa_to_index[aa] for aa in sequence), dtype=np.int16)
        except KeyError:
            valid = "".join(aa if aa in aa_to_index else "X" for aa in sequence)
            starts = [
                start for start in range(len(valid) - 8)
                if "X" not in valid[start : start + 9]
            ]
            if not starts:
                continue
            windows = np.asarray(
                [[aa_to_index[valid[start + index]] for index in exposed] for start in starts],
                dtype=np.int16,
            )
            start_values = np.asarray(starts, dtype=np.int64)
        else:
            start_values = np.arange(len(sequence) - 8, dtype=np.int64)
            windows = np.stack([encoded[start_values + index] for index in exposed], axis=1)

        numerator = blosum[query_indices[None, :], windows].sum(axis=1)
        denominator = np.maximum(
            diagonal[query_indices][None, :], diagonal[windows]
        ).sum(axis=1)
        blosum_values = numerator / denominator
        phys_values = np.abs(
            descriptors[query_indices][None, :, :] - descriptors[windows]
        ).mean(axis=(1, 2))
        identity_values = (windows == query_indices[None, :]).mean(axis=1)

        keep = np.ones(len(start_values), dtype=bool)
        if exclude_accession and _text(record.get("accession")) == exclude_accession:
            for row_index, start in enumerate(start_values):
                if sequence[int(start) : int(start) + 9] == exclude_core:
                    keep[row_index] = False
                    excluded += 1
        blosum_values = blosum_values[keep]
        phys_values = phys_values[keep]
        identity_values = identity_values[keep]
        kept_starts = start_values[keep]
        evaluated += len(kept_starts)

        equal_blosum = np.isclose(blosum_values, target_blosum, rtol=0.0, atol=1e-12)
        equal_phys = np.isclose(phys_values, target_physiochemical, rtol=0.0, atol=1e-12)
        better = (
            (blosum_values > target_blosum + 1e-12)
            | (equal_blosum & (phys_values < target_physiochemical - 1e-12))
            | (
                equal_blosum
                & equal_phys
                & (identity_values >= target_identity - 1e-12)
            )
        )
        at_least_as_good += int(better.sum())

        if len(kept_starts):
            order = np.lexsort((-identity_values, phys_values, -blosum_values))[:top_n]
            for row_index in order:
                start = int(kept_starts[row_index])
                core = sequence[start : start + 9]
                nearest.append(
                    {
                        "accession": _text(record.get("accession")),
                        "protein": _text(record.get("protein")),
                        "start_1_based": start + 1,
                        "core": core,
                        "tcr_facing_blosum62_similarity": float(blosum_values[row_index]),
                        "tcr_face_physicochemical_mismatch": float(phys_values[row_index]),
                        "tcr_facing_sequence_identity": float(identity_values[row_index]),
                    }
                )
            nearest.sort(
                key=lambda row: (
                    -float(row["tcr_facing_blosum62_similarity"]),
                    float(row["tcr_face_physicochemical_mismatch"]),
                    -float(row["tcr_facing_sequence_identity"]),
                    str(row["accession"]),
                    int(row["start_1_based"]),
                    str(row["core"]),
                )
            )
            del nearest[top_n:]

    percentile = round(100.0 * at_least_as_good / evaluated, 12) if evaluated else None
    for rank, row in enumerate(nearest, start=1):
        row["neighbor_rank"] = rank
    return {
        "query_core": query_core,
        "paired_core": paired_core,
        "evaluated_window_count": evaluated,
        "excluded_target_window_count": excluded,
        "at_least_as_good_count": at_least_as_good,
        "empirical_percentile": percentile,
        "target_metrics": target_metrics,
        "nearest_neighbors": nearest,
    }


def presentation_conditioned_rarity(
    *,
    target_pair_id: str,
    target_allele: str,
    candidate_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Locate one frozen pair within its exact-HLA V3 candidate universe."""
    normalized_allele = normalize_hla(target_allele)
    eligible = [
        row for row in candidate_rows
        if normalize_hla(row.get("allele")) == normalized_allele
    ]
    matches = [row for row in eligible if _text(row.get("pair_id")) == target_pair_id]
    if len(matches) != 1:
        return {
            "status": "not_evaluable_target_pair_not_uniquely_resolved",
            "candidate_pair_count": len(eligible),
            "at_least_as_good_count": None,
            "empirical_percentile": None,
        }

    def ranking_key(row: Mapping[str, Any]) -> tuple[float, float, float, str]:
        return (
            -float(row["tcr_facing_blosum62_similarity"]),
            float(row["tcr_face_physicochemical_mismatch"]),
            -float(row["tcr_facing_sequence_identity"]),
            _text(row.get("pair_id")),
        )

    try:
        target_key = ranking_key(matches[0])
        unrelated = [row for row in eligible if _text(row.get("pair_id")) != target_pair_id]
        at_least_as_good = sum(ranking_key(row) < target_key for row in unrelated)
    except (KeyError, TypeError, ValueError):
        return {
            "status": "not_evaluable_missing_sequence_metric",
            "candidate_pair_count": len(eligible),
            "at_least_as_good_count": None,
            "empirical_percentile": None,
        }
    count = len(unrelated)
    return {
        "status": "evaluable" if count else "not_evaluable_empty_comparator_library",
        "candidate_pair_count": len(eligible),
        "unrelated_pair_count": count,
        "at_least_as_good_count": at_least_as_good,
        "empirical_percentile": round(100.0 * at_least_as_good / count, 12) if count else None,
        "target_sequence_metrics": {
            "tcr_facing_blosum62_similarity": float(matches[0]["tcr_facing_blosum62_similarity"]),
            "tcr_face_physicochemical_mismatch": float(matches[0]["tcr_face_physicochemical_mismatch"]),
            "tcr_facing_sequence_identity": float(matches[0]["tcr_facing_sequence_identity"]),
        },
    }


def summarize_conservation(
    peptide: str,
    core: str,
    homologous_sequences: Sequence[str],
) -> dict[str, Any]:
    normalized = [str(item).strip().upper() for item in homologous_sequences if item]
    exact_peptide = sum(peptide in item for item in normalized)
    exact_core = sum(core in item for item in normalized)
    total = len(normalized)
    return {
        "homologous_sequence_count": total,
        "exact_peptide_count": exact_peptide,
        "exact_core_count": exact_core,
        "exact_peptide_fraction": round(exact_peptide / total, 12) if total else None,
        "exact_core_fraction": round(exact_core / total, 12) if total else None,
        "conservation_status": "evaluable" if total else "not_evaluable_missing_sequence_coverage",
    }


def classify_stage1(
    ebv_arm: Mapping[str, Any],
    self_arm: Mapping[str, Any],
    *,
    rarity_percentile: Optional[float],
) -> str:
    arms = (ebv_arm, self_arm)
    if any(str(arm.get("predictor_status")) != "complete" for arm in arms):
        return "stage1_hold"
    if any(bool(arm.get("identity_or_hla_conflict")) for arm in arms):
        return "stage1_hold"
    if any(bool(arm.get("both_predictors_above_20")) for arm in arms):
        return "stage1_hold"
    if (
        all(bool(arm.get("binding_consensus")) for arm in arms)
        and all(bool(arm.get("register_consensus_matches_declared")) for arm in arms)
        and rarity_percentile is not None
        and rarity_percentile <= 1.0
    ):
        return "stage1_high_priority"
    if all(bool(arm.get("binding_supported")) for arm in arms):
        return "stage1_medium_priority"
    return "stage1_hold"


def build_stage2_gate(candidates: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    return {
        "gate_name": "experimental_binding_and_register_before_tcell_assays",
        "status": "not_evaluable_pending_experimental_binding_and_register",
        "candidate_count": len(candidates),
        "required_experimental_results_per_candidate": [
            "ebv_exact_hla_binding",
            "ebv_experimental_p1_p9_register",
            "self_exact_hla_binding",
            "self_experimental_p1_p9_register",
        ],
        "tcell_assay_recommendation_allowed": False,
        "specificity_claim_allowed": False,
        "cross_reactivity_claim_allowed": False,
        "molecular_mimicry_claim_allowed": False,
        "discovery_unlock_allowed": False,
        "weights_frozen": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
