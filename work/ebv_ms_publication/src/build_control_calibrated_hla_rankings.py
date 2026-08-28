"""Build control-referenced rankings separately for each discovery HLA.

The known BALF5-MBP system validates the structural endpoint. It is not used to
pool alleles or tune discovery scores. The formal control reference comes from
the complete fixed-seed calibration set; incomplete seeds remain sensitivity
analyses.
"""

from __future__ import annotations

import csv
import hashlib
import html
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25"
DEFAULT_OUT = ROOT / "processed/tcell_library_v2_control_calibrated_hla_rankings_2026-08-25"
ALLELES = (
    "HLA-DRB1*15:01",
    "HLA-DRB1*13:03",
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
)
ALLELE_SLUGS = {
    "HLA-DRB1*15:01": "hla_drb1_15_01",
    "HLA-DRB1*13:03": "hla_drb1_13_03",
    "HLA-DRB1*03:01": "hla_drb1_03_01",
    "HLA-DRB1*08:01": "hla_drb1_08_01",
}
FORMAL_CONTROL_SEED = 104759
CLAIM_BOUNDARY = (
    "Computational pMHC geometry prioritization only; control-reference metrics are not "
    "probabilities and do not establish presentation, TCR binding, activation, "
    "cross-reactivity, molecular mimicry, or MS disease mechanism."
)


def _truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str] | None = None) -> None:
    fieldnames = list(fields or (list(rows[0]) if rows else []))
    if not fieldnames:
        raise ValueError(f"cannot write empty table without fields: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize_calibration_seed(
    calibration_rows: Sequence[dict[str, Any]],
    seed_rows: Sequence[dict[str, Any]],
    seed: int,
) -> dict[str, Any]:
    matching_seed_rows = [row for row in seed_rows if int(row["seed"]) == seed]
    if len(matching_seed_rows) != 1:
        raise ValueError(f"expected one seed summary for {seed}")
    seed_summary = matching_seed_rows[0]
    primary = [
        row for row in calibration_rows
        if int(row["seed"]) == seed
        and row["analysis_set"] == "primary_rank_of_26"
        and row["geometry_status"] == "complete"
    ]
    positives = [row for row in primary if row["pair_role"] == "E1_positive"]
    decoys = [row for row in primary if row["pair_role"] == "full_decoy"]
    if len(positives) != 1:
        raise ValueError(f"seed {seed} must contain exactly one complete E1 positive")
    decoy_values = sorted(float(row["exposed_ca_rmsd_A_median"]) for row in decoys)
    if not decoy_values:
        raise ValueError(f"seed {seed} has no complete full decoys")
    positive_value = float(positives[0]["exposed_ca_rmsd_A_median"])
    rank = 1 + sum(float(row["exposed_ca_rmsd_A_median"]) < positive_value for row in primary)
    formal = _truth(seed_summary.get("formal_seed_evaluable"))
    return {
        "seed": seed,
        "calibration_role": "formal_control_reference" if formal else "incomplete_seed_sensitivity_only",
        "formal_seed_evaluable": formal,
        "available_primary_count": len(primary),
        "expected_primary_count": 26,
        "available_decoy_count": len(decoy_values),
        "expected_decoy_count": 25,
        "positive_available_rank": rank,
        "positive_exposed_ca_rmsd_median_A": positive_value,
        "decoy_exposed_ca_rmsd_min_A": min(decoy_values),
        "decoy_exposed_ca_rmsd_median_A": median(decoy_values),
        "decoy_exposed_ca_rmsd_max_A": max(decoy_values),
        "decoy_values": decoy_values,
    }


def derive_formal_control_reference(
    calibration_rows: Sequence[dict[str, Any]],
    seed_rows: Sequence[dict[str, Any]],
    seed: int = FORMAL_CONTROL_SEED,
) -> dict[str, Any]:
    reference = summarize_calibration_seed(calibration_rows, seed_rows, seed)
    if not reference["formal_seed_evaluable"]:
        raise ValueError(f"formal control seed {seed} is incomplete")
    if reference["available_primary_count"] != 26 or reference["available_decoy_count"] != 25:
        raise ValueError(f"formal control seed {seed} does not have 1 positive plus 25 decoys")
    if reference["positive_available_rank"] > 3:
        raise ValueError(f"gold-standard positive was not recovered in the top three for seed {seed}")
    if reference["positive_exposed_ca_rmsd_median_A"] >= reference["decoy_exposed_ca_rmsd_median_A"]:
        raise ValueError(f"gold-standard positive did not beat the decoy median for seed {seed}")
    return reference


def control_reference_metrics(value: float, reference: dict[str, Any]) -> dict[str, Any]:
    positive = float(reference["positive_exposed_ca_rmsd_median_A"])
    decoy_median = float(reference["decoy_exposed_ca_rmsd_median_A"])
    decoys = [float(item) for item in reference["decoy_values"]]
    scale = decoy_median - positive
    if scale <= 0:
        raise ValueError("control reference requires positive RMSD below decoy median")
    index = (decoy_median - value) / scale
    if value <= positive:
        band = "at_or_below_gold_positive_median"
    elif value < decoy_median:
        band = "between_gold_positive_and_decoy_medians"
    else:
        band = "at_or_above_decoy_median"
    return {
        "formal_control_seed": int(reference["seed"]),
        "gold_positive_reference_median_A": round(positive, 6),
        "full_decoy_reference_median_A": round(decoy_median, 6),
        "control_separation_index": round(index, 8),
        "positive_reference_delta_A": round(value - positive, 6),
        "decoy_median_delta_A": round(value - decoy_median, 6),
        "conservative_decoy_tail_fraction": round((1 + sum(item <= value for item in decoys)) / (len(decoys) + 1), 8),
        "control_geometry_band": band,
        "control_metric_interpretation": "descriptive_method_reference_not_probability",
    }


def rank_within_hla(
    pair_rows: Sequence[dict[str, Any]],
    reference: dict[str, Any],
) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    by_allele: dict[str, list[dict[str, Any]]] = defaultdict(list)
    missing_rows: list[dict[str, Any]] = []
    for row in pair_rows:
        allele = str(row["allele"])
        if allele not in ALLELES:
            raise ValueError(f"unexpected discovery allele: {allele}")
        if row["geometry_status"] == "complete":
            by_allele[allele].append(dict(row))
        else:
            missing_rows.append({
                "allele": allele,
                "pair_id": row["pair_id"],
                "ebv_candidate_id": row["ebv_candidate_id"],
                "self_candidate_id": row["self_candidate_id"],
                "geometry_status": row["geometry_status"],
                "ranking_status": "not_ranked_missing_geometry",
            })

    ranked: dict[str, list[dict[str, Any]]] = {}
    for allele in ALLELES:
        ordered = sorted(
            by_allele.get(allele, []),
            key=lambda row: (
                float(row["exposed_ca_rmsd_A_median"]),
                float(row["exposed_ca_rmsd_A_iqr"]),
                str(row["pair_id"]),
            ),
        )
        count = len(ordered)
        output_rows: list[dict[str, Any]] = []
        for rank, row in enumerate(ordered, start=1):
            value = float(row["exposed_ca_rmsd_A_median"])
            output_rows.append({
                "allele": allele,
                "hla_rank": rank,
                "hla_evaluable_pair_count": count,
                "hla_percentile": round((rank - 1) / max(1, count - 1), 8),
                "rank_scope": "within_hla_only",
                "pair_id": row["pair_id"],
                "ebv_candidate_id": row["ebv_candidate_id"],
                "ebv_protein": row["ebv_protein"],
                "ebv_sequence": row["ebv_sequence"],
                "ebv_predicted_core": row["ebv_predicted_core"],
                "ebv_binding_percentile_rank": row["ebv_binding_percentile_rank"],
                "ebv_source_certainty": row["ebv_source_certainty"],
                "self_candidate_id": row["self_candidate_id"],
                "self_protein": row["self_protein"],
                "self_sequence": row["self_sequence"],
                "self_predicted_core": row["self_predicted_core"],
                "self_binding_percentile_rank": row["self_binding_percentile_rank"],
                "self_source_certainty": row["self_source_certainty"],
                "exposed_ca_rmsd_A_median": round(value, 6),
                "exposed_ca_rmsd_A_iqr": round(float(row["exposed_ca_rmsd_A_iqr"]), 6),
                "exposed_ca_rmsd_A_q25": row["exposed_ca_rmsd_A_q25"],
                "exposed_ca_rmsd_A_q75": row["exposed_ca_rmsd_A_q75"],
                "model_combination_count": row["model_combination_count"],
                **control_reference_metrics(value, reference),
                "primary_endpoint": row["primary_endpoint"],
                "claim_boundary": CLAIM_BOUNDARY,
            })
        ranked[allele] = output_rows
    return ranked, sorted(missing_rows, key=lambda row: (row["allele"], row["pair_id"]))


def _write_figure(path: Path, ranked: dict[str, list[dict[str, Any]]], reference: dict[str, Any]) -> None:
    width, height = 1200, 760
    panel_width, panel_height = 500, 250
    panel_origins = ((70, 90), (650, 90), (70, 410), (650, 410))
    all_values = [float(row["exposed_ca_rmsd_A_median"]) for rows in ranked.values() for row in rows]
    y_max = max(all_values) if all_values else 1.0
    positive = float(reference["positive_exposed_ca_rmsd_median_A"])
    decoy = float(reference["decoy_exposed_ca_rmsd_median_A"])
    body = [
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#17212b;letter-spacing:0}.title{font-size:22px;font-weight:700}.panel{font-size:14px;font-weight:700}.label{font-size:11px}.small{font-size:9px}.axis{stroke:#4b5563;stroke-width:1}</style>',
        '<text x="600" y="34" text-anchor="middle" class="title">Control-referenced rankings remain separate by HLA</text>',
        '<line x1="765" y1="56" x2="800" y2="56" stroke="#c53030" stroke-width="2"/><text x="808" y="60" class="small">gold-positive median</text>',
        '<line x1="965" y1="56" x2="1000" y2="56" stroke="#805ad5" stroke-width="2"/><text x="1008" y="60" class="small">decoy median</text>',
    ]
    for allele, (x0, y0) in zip(ALLELES, panel_origins):
        rows = ranked[allele]
        body.extend([
            f'<text x="{x0}" y="{y0-12}" class="panel">{html.escape(allele)}</text>',
            f'<line x1="{x0}" y1="{y0+panel_height}" x2="{x0+panel_width}" y2="{y0+panel_height}" class="axis"/>',
            f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+panel_height}" class="axis"/>',
            f'<text x="{x0+panel_width/2}" y="{y0+panel_height+28}" text-anchor="middle" class="label">Within-HLA percentile</text>',
            f'<text x="{x0-42}" y="{y0+panel_height/2}" transform="rotate(-90 {x0-42} {y0+panel_height/2})" text-anchor="middle" class="label">Exposed RMSD (A)</text>',
        ])
        for value, color in ((positive, "#c53030"), (decoy, "#805ad5")):
            y = y0 + min(value, y_max) / y_max * panel_height
            body.append(f'<line x1="{x0}" y1="{y:.2f}" x2="{x0+panel_width}" y2="{y:.2f}" stroke="{color}" stroke-width="1.5"/>')
        for row in rows:
            x = x0 + float(row["hla_percentile"]) * panel_width
            y = y0 + min(float(row["exposed_ca_rmsd_A_median"]), y_max) / y_max * panel_height
            top = int(row["hla_rank"]) <= 10
            body.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{3 if top else 1.5}" fill="{"#1f77b4" if top else "#9ecae1"}" fill-opacity="{1 if top else 0.65}"/>'
            )
        body.extend([
            f'<text x="{x0}" y="{y0+panel_height+13}" class="small">0</text>',
            f'<text x="{x0+panel_width}" y="{y0+panel_height+13}" text-anchor="end" class="small">1</text>',
            f'<text x="{x0-6}" y="{y0+4}" text-anchor="end" class="small">0</text>',
            f'<text x="{x0-6}" y="{y0+panel_height}" text-anchor="end" class="small">{y_max:.1f}</text>',
        ])
    body.append('<text x="600" y="746" text-anchor="middle" class="small">Blue points are ranked only within their own HLA. Control lines are method references, not biological probabilities.</text>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        + "\n".join(body) + "\n</svg>\n",
        encoding="utf-8",
    )


def run(source: Path = SOURCE, out: Path = DEFAULT_OUT) -> dict[str, Any]:
    pair_path = source / "discovery/pair_summary_6400.csv"
    calibration_path = source / "calibration/calibration_pair_summary_72.csv"
    seed_path = source / "calibration/seed_recovery_summary.csv"
    gold_path = source / "validation/gold_standard_capture_summary.json"
    pair_rows = read_csv(pair_path)
    calibration_rows = read_csv(calibration_path)
    seed_rows = read_csv(seed_path)
    gold_summary = json.loads(gold_path.read_text(encoding="utf-8"))
    if gold_summary["gold_standard_independent_system_count"] != 1:
        raise ValueError("control-calibrated ranking requires one locked gold-standard system")
    if gold_summary["capture_at_1_available_seed_fraction"] != 1.0:
        raise ValueError("gold-standard positive was not captured at rank 1 in all available seeds")

    formal_reference = derive_formal_control_reference(calibration_rows, seed_rows)
    calibration_parameters = [
        summarize_calibration_seed(calibration_rows, seed_rows, seed)
        for seed in sorted({int(row["seed"]) for row in seed_rows})
    ]
    for row in calibration_parameters:
        row.pop("decoy_values")
        row["use_in_control_index"] = row["seed"] == FORMAL_CONTROL_SEED
        row["claim_boundary"] = CLAIM_BOUNDARY

    ranked, missing = rank_within_hla(pair_rows, formal_reference)
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "calibration_parameters.csv", calibration_parameters)
    all_ranked = [row for allele in ALLELES for row in ranked[allele]]
    write_csv(out / "all_hla_ranked_pairs.csv", all_ranked)
    for allele in ALLELES:
        write_csv(out / f"rankings/{ALLELE_SLUGS[allele]}_ranked_pairs.csv", ranked[allele])
    top_rows = [row for allele in ALLELES for row in ranked[allele][:25]]
    write_csv(out / "top_25_by_hla.csv", top_rows)
    write_csv(
        out / "missing_unranked_pairs.csv",
        missing,
        fields=("allele", "pair_id", "ebv_candidate_id", "self_candidate_id", "geometry_status", "ranking_status"),
    )
    _write_figure(out / "figures/separate_hla_rankings.svg", ranked, formal_reference)

    results = [
        "# Control-calibrated HLA-specific results",
        "",
        "The four HLA alleles are ranked independently. There is no cross-allele consensus rank in this package.",
        "",
        f"Formal control seed: **{FORMAL_CONTROL_SEED}**; gold-positive median: **{formal_reference['positive_exposed_ca_rmsd_median_A']:.6f} A**; 25-decoy median: **{formal_reference['decoy_exposed_ca_rmsd_median_A']:.6f} A**.",
        "",
    ]
    for allele in ALLELES:
        results.extend([
            f"## {allele}",
            "",
            f"Evaluable pairs: **{len(ranked[allele])}**.",
            "",
            "| Rank | EBV peptide | Self peptide | EBV protein | Self protein | RMSD (A) | Control index |",
            "|---:|---|---|---|---|---:|---:|",
            *[
                f"| {row['hla_rank']} | {row['ebv_sequence']} | {row['self_sequence']} | {row['ebv_protein']} | {row['self_protein']} | {float(row['exposed_ca_rmsd_A_median']):.3f} | {float(row['control_separation_index']):.3f} |"
                for row in ranked[allele][:10]
            ],
            "",
        ])
    results.extend([f"> {CLAIM_BOUNDARY}", ""])
    (out / "RESULTS_SUMMARY.md").write_text("\n".join(results), encoding="utf-8")

    methods = f"""# Methods

## Ranking scope

Each HLA allele is a separate 40-by-40 EBV-self screen. Complete pairs are ranked within that HLA by median exposed-position P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove fit, then ensemble IQR, then frozen pair ID. No geometry, rank, percentile, or score is pooled across alleles.

## Control recalibration

The endpoint was retained because it captured the locked Hy.2E11 BALF5-MBP gold-standard system at rank 1 in both available fixed-seed sets. Seed {FORMAL_CONTROL_SEED} is the formal reference because it contains the complete positive plus all 25 predeclared score-blind full decoys. Seed 104729 remains an incomplete sensitivity analysis and does not set the control index.

The control separation index is `(decoy median - pair RMSD) / (decoy median - positive median)`. A value of 1 equals the modeled gold-positive median and 0 equals the full-decoy median. Values outside 0 to 1 are retained. This index is a descriptive method reference, not a calibrated probability, false-discovery rate, or biological validation score. Within-HLA rank remains the primary result.

## Missingness

Pairs without complete geometry are retained in `missing_unranked_pairs.csv` and excluded from the relevant HLA rank denominator. Missing results are not imputed.

## Interpretation

{CLAIM_BOUNDARY}
"""
    (out / "METHODS.md").write_text(methods, encoding="utf-8")

    manifest = {
        "analysis_version": "EBV_MS_TCELL_V2_CONTROL_CALIBRATED_HLA_RANKINGS_2026-08-25",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_analysis": str(source),
        "cross_allele_consensus_used": False,
        "ranking_scope": "within_hla_only",
        "formal_control_seed": FORMAL_CONTROL_SEED,
        "gold_standard_independent_system_count": gold_summary["gold_standard_independent_system_count"],
        "gold_standard_available_seed_capture_at_1_fraction": gold_summary["capture_at_1_available_seed_fraction"],
        "formal_gold_positive_median_A": formal_reference["positive_exposed_ca_rmsd_median_A"],
        "formal_full_decoy_median_A": formal_reference["decoy_exposed_ca_rmsd_median_A"],
        "evaluable_pairs_by_hla": {allele: len(ranked[allele]) for allele in ALLELES},
        "missing_unranked_pair_count": len(missing),
        "rank_1_pair_by_hla": {allele: ranked[allele][0]["pair_id"] for allele in ALLELES},
        "input_sha256": {
            str(path.relative_to(source)): sha256_file(path)
            for path in (pair_path, calibration_path, seed_path, gold_path)
        },
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (out / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    readme = f"""# Control-calibrated HLA-specific rankings

- Ranking scope: **within each HLA only**
- Cross-allele consensus: **not used**
- Formal control reference: **seed {FORMAL_CONTROL_SEED}**
- Gold-standard available-set capture@1: **2/2 seeds**
- Independent gold-standard systems: **1**

Open `RESULTS_SUMMARY.md` for four separate top-10 tables. Full rankings are under `rankings/`; `top_25_by_hla.csv` is a compact combined view that preserves the HLA boundary.

The earlier full-ensemble package is preserved unchanged. This additive version reuses its frozen pair geometries and control results, changes the reporting and calibration layer, and does not rerank across HLA alleles.

Reproduce:

```bash
PYTHONPATH=src python3 src/build_control_calibrated_hla_rankings.py
```

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")

    checksums = [
        {"relative_path": str(path.relative_to(out)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        for path in sorted(out.rglob("*"), key=str)
        if path.is_file() and path.name != "SHA256SUMS.csv"
    ]
    write_csv(out / "SHA256SUMS.csv", checksums)
    return manifest


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))

