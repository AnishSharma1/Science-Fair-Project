"""PyTorch expression classifier for the GSE190847 antigen-presentation panel.

This is a bounded transcriptomic-context analysis, not a pMHC, EBV infection,
or TCR-recognition assay. It asks whether the pre-registered seven-gene
HLA-II/APC panel carries out-of-sample signal separating untreated PPMS B-cell
samples from healthy controls.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
from pathlib import Path

import torch


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "raw" / "geo" / "GSE190847_series_matrix.txt.gz"
OUT = ROOT / "processed" / "geo_gse190847" / "pytorch_expression_classifier"
GENES = {
    "HLA-DRA": "TC0600007650.hg.1",
    "HLA-DRB1": "TC0600014273.hg.1",
    "CD74": "TC0500012470.hg.1",
    "CIITA": "TC1600006888.hg.1",
    "HLA-DMA": "TC0600014277.hg.1",
    "HLA-DMB": "TC0600014276.hg.1",
    "CTSS": "TC0100015752.hg.1",
}


def quoted_fields(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))


def load_panel() -> tuple[list[str], list[str], torch.Tensor, torch.Tensor]:
    sample_titles: list[str] | None = None
    sample_ids: list[str] | None = None
    values: dict[str, list[float]] = {}
    in_table = False
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                sample_titles = quoted_fields(line)[1:]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif in_table and line.startswith('"ID_REF"'):
                sample_ids = quoted_fields(line)[1:]
            elif in_table and line.startswith('"'):
                fields = quoted_fields(line)
                if fields[0] in GENES.values():
                    values[fields[0]] = [float(value) for value in fields[1:]]
    if sample_titles is None or sample_ids is None or set(values) != set(GENES.values()):
        raise RuntimeError("Missing expected sample metadata or registered probes.")
    rows = []
    labels = []
    kept_ids = []
    kept_titles = []
    for index, title in enumerate(sample_titles):
        lower = title.lower()
        if "healthy control" in lower:
            labels.append(0.0)
        elif "PPMS" in title:
            labels.append(1.0)
        else:
            continue
        kept_ids.append(sample_ids[index])
        kept_titles.append(title)
        rows.append([values[probe][index] for probe in GENES.values()])
    if len(labels) != 41 or labels.count(0.0) != 28 or labels.count(1.0) != 13:
        raise RuntimeError(f"Unexpected retained cohort: healthy={labels.count(0.0)} PPMS={labels.count(1.0)}")
    return kept_ids, kept_titles, torch.tensor(rows, dtype=torch.float32), torch.tensor(labels, dtype=torch.float32)


def auc_score(labels: list[float], scores: list[float]) -> float:
    pairs = sorted(zip(scores, labels), key=lambda pair: pair[0])
    rank_sum = 0.0
    index = 0
    while index < len(pairs):
        end = index + 1
        while end < len(pairs) and pairs[end][0] == pairs[index][0]:
            end += 1
        average_rank = (index + 1 + end) / 2
        for tied_index in range(index, end):
            if pairs[tied_index][1] == 1.0:
                rank_sum += average_rank
        index = end
    positives = sum(1 for label in labels if label == 1.0)
    negatives = len(labels) - positives
    return (rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def train_fold(x_train: torch.Tensor, y_train: torch.Tensor, x_test: torch.Tensor, epochs: int, seed: int) -> float:
    torch.manual_seed(seed)
    model = torch.nn.Linear(x_train.shape[1], 1)
    positives = float(y_train.sum().item())
    negatives = float(len(y_train) - positives)
    loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor([negatives / positives]))
    optimizer = torch.optim.Adam(model.parameters(), lr=0.04, weight_decay=0.02)
    for _ in range(epochs):
        optimizer.zero_grad()
        loss = loss_fn(model(x_train).squeeze(1), y_train)
        loss.backward()
        optimizer.step()
    with torch.no_grad():
        return float(torch.sigmoid(model(x_test).squeeze(1))[0].item())


def leave_one_out_probabilities(features: torch.Tensor, labels: torch.Tensor, epochs: int, seed: int) -> list[float]:
    probabilities = []
    for held_out in range(len(labels)):
        train_mask = torch.ones(len(labels), dtype=torch.bool)
        train_mask[held_out] = False
        x_train = features[train_mask]
        y_train = labels[train_mask]
        means = x_train.mean(dim=0)
        stdevs = x_train.std(dim=0).clamp_min(1e-6)
        x_train = (x_train - means) / stdevs
        x_test = (features[held_out : held_out + 1] - means) / stdevs
        probabilities.append(train_fold(x_train, y_train, x_test, epochs, seed + held_out))
    return probabilities


def summarize(labels: list[float], probabilities: list[float]) -> dict[str, float]:
    predictions = [1.0 if value >= 0.5 else 0.0 for value in probabilities]
    positives = [index for index, label in enumerate(labels) if label == 1.0]
    negatives = [index for index, label in enumerate(labels) if label == 0.0]
    accuracy = sum(1 for label, pred in zip(labels, predictions) if label == pred) / len(labels)
    sensitivity = sum(1 for index in positives if predictions[index] == 1.0) / len(positives)
    specificity = sum(1 for index in negatives if predictions[index] == 0.0) / len(negatives)
    return {
        "auc": auc_score(labels, probabilities),
        "accuracy_at_0.5": accuracy,
        "balanced_accuracy_at_0.5": (sensitivity + specificity) / 2,
        "sensitivity_at_0.5": sensitivity,
        "specificity_at_0.5": specificity,
        "mean_ppms_probability": sum(probabilities[index] for index in positives) / len(positives),
        "mean_control_probability": sum(probabilities[index] for index in negatives) / len(negatives),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=350)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=1501)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    random.seed(args.seed)
    torch.set_num_threads(1)
    sample_ids, titles, features, labels_tensor = load_panel()
    labels = [float(value) for value in labels_tensor.tolist()]
    probabilities = leave_one_out_probabilities(features, labels_tensor, args.epochs, args.seed)
    metrics = summarize(labels, probabilities)

    permuted_aucs = []
    for iteration in range(args.permutations):
        shuffled = labels[:]
        random.shuffle(shuffled)
        shuffled_tensor = torch.tensor(shuffled, dtype=torch.float32)
        perm_probs = leave_one_out_probabilities(features, shuffled_tensor, max(120, args.epochs // 2), args.seed + 1000 + iteration)
        permuted_aucs.append(auc_score(shuffled, perm_probs))
    empirical_p = (1 + sum(auc >= metrics["auc"] for auc in permuted_aucs)) / (args.permutations + 1)

    sample_rows = []
    for sample_id, title, label, probability in zip(sample_ids, titles, labels, probabilities):
        sample_rows.append({
            "sample_id": sample_id,
            "sample_title": title,
            "group": "PPMS" if label == 1.0 else "healthy_control",
            "loocv_ppms_probability": round(probability, 6),
        })
    with (OUT / "gse190847_pytorch_loocv_predictions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_rows[0].keys())
        writer.writeheader()
        writer.writerows(sample_rows)

    metric_rows = [{
        "analysis": "seven_gene_hla_ii_apc_panel_pytorch_logistic_loocv",
        "healthy_n": int(labels.count(0.0)),
        "ppms_n": int(labels.count(1.0)),
        "features": ";".join(GENES.keys()),
        "epochs": args.epochs,
        "permutations": args.permutations,
        "empirical_auc_p_ge_observed": empirical_p,
        "interpretation": "Transcriptomic context only; not EBV, pMHC, TCR, activation, or disease-mechanism evidence.",
        **{key: round(value, 6) for key, value in metrics.items()},
    }]
    with (OUT / "gse190847_pytorch_classifier_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metric_rows[0].keys())
        writer.writeheader()
        writer.writerows(metric_rows)

    with (OUT / "README.md").open("w") as handle:
        handle.write("# PyTorch expression classifier for GSE190847\n\n")
        handle.write("This analysis uses a seven-gene HLA-II/APC expression panel from RMA log2 microarray B-cell data.\n\n")
        handle.write("It is intentionally bounded: it tests out-of-sample group signal in peripheral B cells and does not measure EBV infection, pMHC presentation, TCR binding, T-cell activation, or MS causality.\n\n")
        handle.write(f"- LOOCV AUC: {metrics['auc']:.3f}\n")
        handle.write(f"- Balanced accuracy at 0.5: {metrics['balanced_accuracy_at_0.5']:.3f}\n")
        handle.write(f"- Empirical permutation p(AUC >= observed): {empirical_p:.3f} using {args.permutations} permutations\n")
        handle.write(f"- PyTorch version: {torch.__version__}\n")

    print(f"AUC={metrics['auc']:.3f} balanced_accuracy={metrics['balanced_accuracy_at_0.5']:.3f} empirical_p={empirical_p:.3f}")


if __name__ == "__main__":
    main()
