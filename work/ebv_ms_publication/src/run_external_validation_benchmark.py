"""External validation overlay for the pMHC candidate-prioritization screen.

This benchmark asks whether independent literature-positive EBV or MBP-region
candidates appear unusually high in the predeclared modeled-pMHC shortlist.
It is a rank-recovery test only; it does not test TCR binding or activation.
"""

from __future__ import annotations

import csv
import random
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "external_validation_benchmark"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def empirical_ge(observed: float, null_values: list[float]) -> float:
    return (1 + sum(value >= observed for value in null_values)) / (len(null_values) + 1)


def empirical_le(observed: float, null_values: list[float]) -> float:
    return (1 + sum(value <= observed for value in null_values)) / (len(null_values) + 1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    random.seed(1501)
    panel = read_csv(PROC / "external_validation_panel.csv")
    shortlist = read_csv(PROC / "fullscreen_tier1_ebv_myelin_shortlist.csv")
    manifest = {row["candidate_id"]: row for row in read_csv(PROC / "pmhc_candidate_manifest.csv")}

    candidate_groups: dict[str, set[str]] = {}
    candidate_basis: dict[str, list[str]] = {}
    for row in panel:
        candidate_groups.setdefault(row["candidate_id"], set()).add(row["validation_group"])
        candidate_basis.setdefault(row["candidate_id"], []).append(row["validation_basis"])

    annotated = []
    for rank, row in enumerate(shortlist, start=1):
        ebv_id = row["ebv_candidate_id"]
        human_id = row["human_candidate_id"]
        groups = sorted(candidate_groups.get(ebv_id, set()) | candidate_groups.get(human_id, set()))
        pair_validation = "background"
        if (
            ebv_id == "EBV_TCELL_63843"
            and "classic_BALF5_MBP_structural_positive" in candidate_groups.get(human_id, set())
        ):
            pair_validation = "classic_BALF5_MBP_pair"
        elif "drosu_2024_DRB1501_EBV_glycoprotein" in groups and "wang_2026_MBP90_region" in groups:
            pair_validation = "combined_new_literature_overlay"
        elif "drosu_2024_DRB1501_EBV_glycoprotein" in groups:
            pair_validation = "drosu_2024_EBV_glycoprotein"
        elif "wang_2026_MBP90_region" in groups:
            pair_validation = "wang_2026_MBP90_region"
        elif "classic_BALF5_MBP_structural_positive" in groups:
            pair_validation = "classic_component_only"
        score = float(row["review_priority_heuristic"])
        annotated.append({
            **row,
            "rank": rank,
            "rank_percentile_top_is_high": round(1 - (rank - 1) / (len(shortlist) - 1), 6),
            "pair_validation": pair_validation,
            "is_external_overlay": pair_validation != "background",
            "is_strict_external_new_literature": pair_validation in {
                "drosu_2024_EBV_glycoprotein",
                "wang_2026_MBP90_region",
                "combined_new_literature_overlay",
            },
            "ebv_source_antigen": manifest.get(ebv_id, {}).get("source_antigen", ""),
            "human_source_antigen": manifest.get(human_id, {}).get("source_antigen", ""),
            "ebv_peptide": manifest.get(ebv_id, {}).get("peptide", ""),
            "human_peptide": manifest.get(human_id, {}).get("peptide", ""),
            "validation_basis": " | ".join(candidate_basis.get(ebv_id, []) + candidate_basis.get(human_id, [])),
            "score_numeric": score,
        })

    total = len(annotated)
    overlay = [row for row in annotated if row["is_external_overlay"]]
    strict = [row for row in annotated if row["is_strict_external_new_literature"]]
    top10 = [row for row in annotated if int(row["rank"]) <= min(10, total)]

    scores = [row["score_numeric"] for row in annotated]
    ranks = [int(row["rank"]) for row in annotated]
    overlay_count = len(overlay)
    strict_count = len(strict)

    def permuted_stats(count: int, iterations: int = 10000) -> tuple[list[float], list[float], list[int]]:
        null_mean_scores = []
        null_mean_ranks = []
        null_top10_counts = []
        indices = list(range(total))
        for _ in range(iterations):
            selected = set(random.sample(indices, count))
            null_mean_scores.append(statistics.mean(scores[index] for index in selected))
            null_mean_ranks.append(statistics.mean(ranks[index] for index in selected))
            null_top10_counts.append(sum(1 for index in selected if ranks[index] <= min(10, total)))
        return null_mean_scores, null_mean_ranks, null_top10_counts

    summary_rows = []
    classic_pairs = [row for row in annotated if row["pair_validation"] == "classic_BALF5_MBP_pair"]
    for label, rows, count in (
        ("classic_BALF5_MBP_pair_recovery", classic_pairs, len(classic_pairs)),
        ("any_external_overlay", overlay, overlay_count),
        ("strict_new_literature_overlay", strict, strict_count),
    ):
        if not rows:
            continue
        null_scores, null_ranks, null_top10 = permuted_stats(count)
        observed_mean_score = statistics.mean(row["score_numeric"] for row in rows)
        observed_mean_rank = statistics.mean(int(row["rank"]) for row in rows)
        observed_top10_count = sum(1 for row in rows if int(row["rank"]) <= min(10, total))
        summary_rows.append({
            "test": label,
            "scored_pair_universe_n": total,
            "positive_pair_count": count,
            "observed_mean_priority_score": round(observed_mean_score, 6),
            "empirical_p_mean_score_ge_observed": empirical_ge(observed_mean_score, null_scores),
            "observed_mean_rank": round(observed_mean_rank, 3),
            "empirical_p_mean_rank_le_observed": empirical_le(observed_mean_rank, null_ranks),
            "observed_top10_positive_count": observed_top10_count,
            "top10_size": min(10, total),
            "empirical_p_top10_count_ge_observed": empirical_ge(observed_top10_count, null_top10),
            "interpretation": "Rank-recovery benchmark only; does not test TCR binding, activation, affinity, or patient pathogenicity.",
        })

    with (OUT / "external_validation_pair_annotations.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [key for key in annotated[0] if key != "score_numeric"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in annotated:
            clean = dict(row)
            clean.pop("score_numeric")
            writer.writerow(clean)

    with (OUT / "external_validation_rank_recovery_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    with (OUT / "README.md").open("w", encoding="utf-8") as handle:
        handle.write("# External validation rank-recovery benchmark\n\n")
        handle.write("This analysis overlays independently reported EBV/MS-positive candidates onto the predeclared pMHC shortlist.\n\n")
        handle.write("It is a prioritization benchmark only. It does not test TCR binding, T-cell activation, affinity, or patient pathogenicity.\n\n")
        for row in summary_rows:
            handle.write(f"- {row['test']}: {row['positive_pair_count']} positive pairs in a {row['scored_pair_universe_n']}-pair universe; ")
            handle.write(f"mean rank {row['observed_mean_rank']}; top-10 positives {row['observed_top10_positive_count']}; ")
            handle.write(f"top-10 empirical p={row['empirical_p_top10_count_ge_observed']:.4f}\n")

    print(f"Wrote {OUT}")
    for row in summary_rows:
        print(
            row["test"],
            "n=", row["positive_pair_count"],
            "mean_rank=", row["observed_mean_rank"],
            "top10=", row["observed_top10_positive_count"],
            "p_top10=", round(row["empirical_p_top10_count_ge_observed"], 4),
        )


if __name__ == "__main__":
    main()
