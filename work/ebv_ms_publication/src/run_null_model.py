"""Run the pre-specified composition-preserving sequence null model."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
SEED = 20260801
N_PERM = 10000


def peptide_max_5mer_overlap(peptide: str, human_5mers: list[set[str]]) -> float:
    """Fast exact 5-mer overlap sensitivity metric."""
    if len(peptide) < 5:
        return 0.0
    kmers = {peptide[i : i + 5] for i in range(len(peptide) - 4)}
    return max((len(kmers & target) / len(kmers) for target in human_5mers), default=0.0)


def main() -> None:
    ebv = pd.read_csv(PROC / "tcell_ebv_drb1501.csv")
    human = pd.read_csv(PROC / "human_drb1501_mhc_ii_iedb_enriched.csv")
    positive = {"Positive", "Positive-Low", "Positive-Intermediate", "Positive-High"}
    ebv = ebv[ebv.outcome.isin(positive)].drop_duplicates("peptide")
    human = human[
        (human.candidate_class == "myelin_candidate")
        & (human.provenance_status == "coordinate_validated")
    ].drop_duplicates("peptide")
    ebv_peptides = ebv.peptide.tolist()
    human_peptides = human.peptide.tolist()
    human_5mers = [{target[i : i + 5] for i in range(len(target) - 4)} for target in human_peptides]
    observed_by_peptide = {p: peptide_max_5mer_overlap(p, human_5mers) for p in ebv_peptides}
    observed = sum(observed_by_peptide.values()) / len(observed_by_peptide)

    rng = random.Random(SEED)
    null_means = []
    for _ in range(N_PERM):
        shuffled = ["".join(rng.sample(list(p), len(p))) for p in ebv_peptides]
        null_means.append(sum(peptide_max_5mer_overlap(p, human_5mers) for p in shuffled) / len(shuffled))
    p_value = (1 + sum(v >= observed for v in null_means)) / (N_PERM + 1)
    result = {
        "seed": SEED,
        "permutations": N_PERM,
        "n_ebv_positive_peptides": len(ebv_peptides),
        "n_human_coordinate_validated_myelin_peptides": len(human_peptides),
        "metric": "maximum fraction of EBV contiguous 5-mers shared with any human myelin peptide",
        "observed_mean_max_5mer_overlap": observed,
        "empirical_p_value": p_value,
        "null_mean": sum(null_means) / len(null_means),
        "null_95_percent_interval": [
            sorted(null_means)[int(0.025 * N_PERM)],
            sorted(null_means)[int(0.975 * N_PERM) - 1],
        ],
        "per_peptide_observed_max_5mer_overlap": observed_by_peptide,
        "interpretation": "sequence-resemblance screen only; not evidence of pathogenic molecular mimicry",
    }
    (PROC / "null_model_5mer_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    pd.DataFrame(
        [{"ebv_peptide": p, "max_5mer_overlap": s} for p, s in observed_by_peptide.items()]
    ).sort_values("max_5mer_overlap", ascending=False).to_csv(PROC / "ebv_peptide_5mer_summary.csv", index=False)
    print(json.dumps({k: v for k, v in result.items() if k != "per_peptide_observed_max_5mer_overlap"}, indent=2))


if __name__ == "__main__":
    main()
