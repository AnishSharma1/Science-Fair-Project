"""Exploratory sequence comparison; not a pathogenicity or binding model."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"


def smith_waterman(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2):
    # Simple ungapped-chemistry-neutral local alignment. We report the scoring
    # rule explicitly so this exploratory screen cannot be mistaken for an
    # established immunogenicity predictor.
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
    matches = aligned = 0
    while i > 0 and j > 0 and score[i][j] > 0:
        step = trace[i][j]
        if step == "diag":
            aligned += 1
            matches += a[i - 1] == b[j - 1]
            i -= 1
            j -= 1
        elif step == "up":
            aligned += 1
            i -= 1
        elif step == "left":
            aligned += 1
            j -= 1
        else:
            break
    return int(s), int(matches), int(aligned), float(matches / aligned) if aligned else 0.0


def main() -> None:
    ebv = pd.read_csv(PROC / "tcell_ebv_drb1501.csv")
    human = pd.read_csv(PROC / "human_drb1501_mhc_ii_iedb.csv")
    positive_labels = {"Positive", "Positive-Low", "Positive-Intermediate", "Positive-High"}
    ebv = ebv[ebv["outcome"].isin(positive_labels)].copy()
    ebv = ebv.drop_duplicates("peptide")
    human = human[human["candidate_class"] == "myelin_candidate"].drop_duplicates("peptide")

    rows = []
    for _, e in ebv.iterrows():
        for _, h in human.iterrows():
            score, matches, aligned, identity = smith_waterman(e.peptide, h.peptide)
            rows.append(
                {
                    "ebv_peptide": e.peptide,
                    "ebv_source_antigen": e.source_antigen_name,
                    "ebv_iedb_assay_id": e.iedb_assay_id,
                    "human_peptide": h.peptide,
                    "human_source_antigen": h.source_antigen_name,
                    "human_iedb_epitope_id": h.iedb_epitope_id,
                    "local_score": score,
                    "local_matches": matches,
                    "local_aligned_length": aligned,
                    "local_identity": round(identity, 4),
                    "alignment_fraction_of_shorter": round(aligned / min(len(e.peptide), len(h.peptide)), 4),
                    "exact_substring": h.peptide in e.peptide or e.peptide in h.peptide,
                    "interpretation": "exploratory sequence resemblance only",
                }
            )
    # Do not let a one- or two-residue perfect match dominate the screen.
    # The shortlist requires at least five aligned residues; all pairs remain
    # in the full table for transparency.
    rows.sort(
        key=lambda r: (
            r["local_aligned_length"] >= 5,
            r["local_matches"],
            r["local_identity"],
            r["local_score"],
        ),
        reverse=True,
    )
    out = PROC / "exploratory_ebv_myelin_similarity.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    top = pd.DataFrame(rows).head(25)
    top.to_csv(PROC / "exploratory_top25_similarity.csv", index=False)
    print(f"positive EBV peptides: {len(ebv)}")
    print(f"myelin peptides: {len(human)}")
    print(f"pairwise comparisons: {len(rows)}")
    print(top[["ebv_peptide", "human_peptide", "local_score", "local_identity", "exact_substring"]].to_string(index=False))


if __name__ == "__main__":
    main()
