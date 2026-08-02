"""Transparent physicochemical similarity and matched-control analysis.

This is a hypothesis-generating descriptor model, not an immunogenicity or
TCR-cross-reactivity predictor.
"""

from __future__ import annotations

import csv
import json
import random
from pathlib import Path

import pandas as pd

from explore_similarity import smith_waterman


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
SEED = 20260801
N_PERM = 10000

# Kyte-Doolittle hydropathy, normalized to [0, 1]. The remaining descriptors
# are deliberately simple categorical physicochemical features.
KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5,
    "N": -3.5, "K": -3.9, "R": -4.5,
}
KD_MIN, KD_MAX = min(KD.values()), max(KD.values())
HYDROPHOBIC = set("AVILMFCWY")
AROMATIC = set("FWY")
SMALL = set("AGST")
CHARGE = {aa: (1 if aa in "KR" else -1 if aa in "DE" else 0) for aa in KD}
SIZE = {aa: (0 if aa in "AGSTCP" else 2 if aa in "FWYH" else 1) for aa in KD}


def smith_waterman_pairs(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2):
    """Return the best local aligned residue pairs under the documented score."""
    m, n = len(a), len(b)
    score = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[None] * (n + 1) for _ in range(m + 1)]
    best = (0, 0, 0)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            choices = [
                (0, None),
                (score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch), "diag"),
                (score[i - 1][j] + gap, "up"),
                (score[i][j - 1] + gap, "left"),
            ]
            score[i][j], trace[i][j] = max(choices, key=lambda x: x[0])
            if score[i][j] > best[0]:
                best = (score[i][j], i, j)
    s, i, j = best
    pairs = []
    while i > 0 and j > 0 and score[i][j] > 0:
        step = trace[i][j]
        if step == "diag":
            pairs.append((a[i - 1], b[j - 1]))
            i -= 1
            j -= 1
        elif step == "up":
            i -= 1
        elif step == "left":
            j -= 1
        else:
            break
    pairs.reverse()
    return int(s), pairs


def pair_similarity(x: str, y: str):
    score, pairs = smith_waterman_pairs(x, y)
    if len(pairs) < 5:
        return {"local_score": score, "aligned_length": len(pairs), "property_similarity": 0.0}
    hydro = [1 - abs((KD[a] - KD_MIN) / (KD_MAX - KD_MIN) - (KD[b] - KD_MIN) / (KD_MAX - KD_MIN)) for a, b in pairs]
    charge = [1 - abs(CHARGE[a] - CHARGE[b]) / 2 for a, b in pairs]
    aromatic = [float((a in AROMATIC) == (b in AROMATIC)) for a, b in pairs]
    size = [1 - abs(SIZE[a] - SIZE[b]) / 2 for a, b in pairs]
    return {
        "local_score": score,
        "aligned_length": len(pairs),
        "property_similarity": round(sum(hydro + charge + aromatic + size) / (4 * len(pairs)), 4),
        "hydrophobicity_similarity": round(sum(hydro) / len(hydro), 4),
        "charge_similarity": round(sum(charge) / len(charge), 4),
        "aromatic_similarity": round(sum(aromatic) / len(aromatic), 4),
        "size_similarity": round(sum(size) / len(size), 4),
    }


def max_similarity(peptide: str, targets: list[str]) -> float:
    vals = [pair_similarity(peptide, target)["property_similarity"] for target in targets]
    return max(vals, default=0.0)


def nearest_length_controls(positives: pd.DataFrame, negatives: pd.DataFrame) -> pd.DataFrame:
    """Deterministically match each positive to a closest-length negative."""
    rows = []
    for _, p in positives.sort_values(["peptide_length", "peptide"]).iterrows():
        pool = negatives.assign(distance=(negatives.peptide_length - p.peptide_length).abs())
        best_distance = pool.distance.min()
        chosen = pool[pool.distance == best_distance].sort_values("peptide").iloc[0]
        rows.append({
            "positive_peptide": p.peptide,
            "positive_length": int(p.peptide_length),
            "negative_control_peptide": chosen.peptide,
            "negative_length": int(chosen.peptide_length),
            "absolute_length_difference": int(best_distance),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ebv = pd.read_csv(PROC / "tcell_ebv_drb1501.csv").drop_duplicates("peptide")
    human = pd.read_csv(PROC / "human_drb1501_mhc_ii_iedb_enriched.csv")
    human = human[human.provenance_status == "coordinate_validated"].drop_duplicates("peptide")
    positive_labels = {"Positive", "Positive-Low", "Positive-Intermediate", "Positive-High"}
    positives = ebv[ebv.outcome.isin(positive_labels)].copy()
    negatives = ebv[ebv.outcome == "Negative"].copy()
    myelin = human[human.candidate_class == "myelin_candidate"].peptide.tolist()
    background = human[human.candidate_class == "human_background"].peptide.tolist()

    controls = nearest_length_controls(positives, negatives)
    controls.to_csv(PROC / "ebv_length_matched_negative_controls.csv", index=False)

    rows = []
    for _, e in ebv.iterrows():
        for group, targets in [("myelin", myelin), ("human_background", background)]:
            for target in targets:
                metrics = pair_similarity(e.peptide, target)
                rows.append({"ebv_peptide": e.peptide, "ebv_outcome": e.outcome, "target_group": group, "human_peptide": target, **metrics})
    pair_df = pd.DataFrame(rows)
    pair_df.to_csv(PROC / "ebv_physicochemical_pairwise.csv", index=False)

    peptide_rows = []
    for _, e in ebv.iterrows():
        peptide_rows.append({
            "ebv_peptide": e.peptide,
            "ebv_outcome": e.outcome,
            "ebv_length": int(e.peptide_length),
            "max_myelin_property_similarity": max_similarity(e.peptide, myelin),
            "max_background_property_similarity": max_similarity(e.peptide, background),
        })
    pep_df = pd.DataFrame(peptide_rows)
    pep_df.to_csv(PROC / "ebv_physicochemical_summary.csv", index=False)

    # Label-permutation test for the difference in max-myelin similarity.
    vals = pep_df.max_myelin_property_similarity.tolist()
    labels = [x in positive_labels for x in pep_df.ebv_outcome]
    observed = pep_df[pep_df.ebv_outcome.isin(positive_labels)].max_myelin_property_similarity.mean() - pep_df[pep_df.ebv_outcome == "Negative"].max_myelin_property_similarity.mean()
    rng = random.Random(SEED)
    null = []
    n_pos = sum(labels)
    for _ in range(N_PERM):
        shuffled = labels[:]
        rng.shuffle(shuffled)
        a = [v for v, lab in zip(vals, shuffled) if lab]
        b = [v for v, lab in zip(vals, shuffled) if not lab]
        null.append(sum(a) / len(a) - sum(b) / len(b))
    p_value = (1 + sum(v >= observed for v in null)) / (N_PERM + 1)
    summary = {
        "metric": "maximum local-alignment physicochemical similarity across aligned residues",
        "descriptor_components": ["normalized Kyte-Doolittle hydropathy", "charge compatibility", "aromatic compatibility", "coarse size compatibility"],
        "n_positive_unique_ebv": int(sum(labels)),
        "n_negative_unique_ebv": int(len(labels) - sum(labels)),
        "n_coordinate_validated_myelin_targets": len(myelin),
        "n_coordinate_validated_human_background_targets": len(background),
        "observed_positive_minus_negative_mean": float(observed),
        "label_permutation_p_value": float(p_value),
        "permutations": N_PERM,
        "seed": SEED,
        "interpretation": "exploratory physicochemical resemblance only; not evidence of mimicry or pathogenicity",
    }
    (PROC / "physicochemical_similarity_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
