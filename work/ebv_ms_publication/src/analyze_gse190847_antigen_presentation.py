"""Pre-registered antigen-presentation analysis of GSE190847 B cells.

Compares one peripheral-blood sample per individual: untreated PPMS (n=13)
versus healthy controls (n=28). Values are RMA log2 microarray intensities.
This provides B-cell expression context only, not EBV infection or pMHC/TCR
evidence.
"""

from __future__ import annotations

import csv
import gzip
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "raw" / "geo" / "GSE190847_series_matrix.txt.gz"
OUT = ROOT / "processed" / "geo_gse190847"
GENES = {
    "HLA-DRA": "TC0600007650.hg.1",
    "HLA-DRB1": "TC0600014273.hg.1",
    "CD74": "TC0500012470.hg.1",
    "CIITA": "TC1600006888.hg.1",
    "HLA-DMA": "TC0600014277.hg.1",
    "HLA-DMB": "TC0600014276.hg.1",
    "CTSS": "TC0100015752.hg.1",
}


def normal_cdf(value: float) -> float:
    return (1 + math.erf(value / math.sqrt(2))) / 2


def welch_pvalue(left: list[float], right: list[float]) -> tuple[float, float]:
    """Normal-approximation p value; transparently conservative for this context screen."""
    mean_difference = statistics.mean(right) - statistics.mean(left)
    variance = statistics.variance(left) / len(left) + statistics.variance(right) / len(right)
    z = mean_difference / math.sqrt(variance)
    return mean_difference, 2 * (1 - normal_cdf(abs(z)))


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    count = len(pvalues)
    order = sorted(range(count), key=lambda index: pvalues[index])
    adjusted = [1.0] * count
    running = 1.0
    for rank, index in reversed(list(enumerate(order, start=1))):
        running = min(running, pvalues[index] * count / rank)
        adjusted[index] = running
    return adjusted


def quoted_fields(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
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
        raise RuntimeError("Missing expected sample metadata or one of the registered probes.")
    healthy = [index for index, title in enumerate(sample_titles) if "healthy control" in title.lower()]
    ppms = [index for index, title in enumerate(sample_titles) if "PPMS" in title]
    if len(healthy) != 28 or len(ppms) != 13:
        raise RuntimeError(f"Unexpected group sizes: healthy={len(healthy)}, PPMS={len(ppms)}")
    rows = []
    raw_pvalues = []
    for gene, probe in GENES.items():
        control = [values[probe][index] for index in healthy]
        case = [values[probe][index] for index in ppms]
        difference, pvalue = welch_pvalue(control, case)
        raw_pvalues.append(pvalue)
        pooled_sd = math.sqrt(((len(control) - 1) * statistics.variance(control) + (len(case) - 1) * statistics.variance(case)) / (len(control) + len(case) - 2))
        rows.append({
            "gene": gene,
            "official_clariom_d_probe": probe,
            "healthy_n": len(control),
            "ppms_n": len(case),
            "healthy_mean_log2_rma": round(statistics.mean(control), 4),
            "ppms_mean_log2_rma": round(statistics.mean(case), 4),
            "ppms_minus_healthy_log2_intensity": round(difference, 4),
            "cohens_d": round(difference / pooled_sd, 4),
            "raw_p_normal_approximation": pvalue,
        })
    adjusted = benjamini_hochberg(raw_pvalues)
    for row, value in zip(rows, adjusted):
        row["bh_fdr_7_gene_panel"] = value
    # Module: mean healthy-standardized value across the seven pre-registered genes.
    control_z, case_z = [], []
    for group, store in ((healthy, control_z), (ppms, case_z)):
        for index in group:
            zscores = []
            for probe in GENES.values():
                control_values = [values[probe][i] for i in healthy]
                zscores.append((values[probe][index] - statistics.mean(control_values)) / statistics.stdev(control_values))
            store.append(statistics.mean(zscores))
    module_difference, module_p = welch_pvalue(control_z, case_z)
    module = [{
        "module": "HLA-II_antigen_presentation",
        "genes": ";".join(GENES),
        "healthy_n": len(control_z),
        "ppms_n": len(case_z),
        "healthy_mean_z": round(statistics.mean(control_z), 4),
        "ppms_mean_z": round(statistics.mean(case_z), 4),
        "ppms_minus_healthy_z": round(module_difference, 4),
        "raw_p_normal_approximation": module_p,
        "interpretation": "Peripheral blood B-cell expression context only; not antigen presentation or EBV evidence.",
    }]
    with (OUT / "gse190847_antigen_presentation_genes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "gse190847_antigen_presentation_module.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=module[0].keys())
        writer.writeheader()
        writer.writerows(module)


if __name__ == "__main__":
    main()
