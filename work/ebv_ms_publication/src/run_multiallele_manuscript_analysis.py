"""Generate the dated multi-allele EBV--MS manuscript analysis package."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from statistics import median

import numpy as np

from analyze_af3_pmhc_downloads import parse_mmcif
from multiallele_manuscript_analysis import (
    ALLELES,
    ALLELE_CODES,
    ANCHOR_EBV,
    ANCHOR_HUMAN,
    CLAIM_BOUNDARY,
    analyze_job,
    build_pair_universe,
    build_prediction_submissions,
    build_robustness_jobs,
    direct_register_sequence_metrics,
    discover_download_jobs,
    fetch_iedb_predictions,
    inventory_against_manifest,
    prediction_records_from_tsv,
    read_csv,
    same_register_geometry_from_models,
    select_score_blind_controls,
    sha256_file,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[1]
DOWNLOADS = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Downloads"
INPUT_PACKAGE = DOWNLOADS / "alphafold_multiallele_5x30_2026-08-20"
OUT = ROOT / "processed/multiallele_analysis_2026-08-21"
CONTROL_SOURCE = ROOT / "processed/structural_control_expansion_2026-08-15/frozen_control_universe.csv"
FLANK_SOURCE = ROOT / "raw/iedb_natural_flank_extensions.csv"
RETRY_SOURCE = Path.home() / ".codex/.chatgpt-projects/g-p-6a795d2b4b0081919da21a9ed96c5ed5/outputs/alphafold_131_download_audit_2026-08-21/alphafold_server_retry_20_jobs.json"
RESULT_GLOB = "folds_2026_08_*"


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _bootstrap_ci(values: list[float], seed: int, iterations: int = 2000) -> tuple[float, float]:
    if not values:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    array = np.asarray(values, dtype=float)
    draws = rng.choice(array, size=(iterations, len(array)), replace=True)
    medians = np.median(draws, axis=1)
    return float(np.quantile(medians, 0.025)), float(np.quantile(medians, 0.975))


def _stable_seed(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def _saved_raw_prediction_records(
    allele: str, submissions: list[dict[str, object]], raw_path: Path
) -> list[dict[str, object]]:
    if raw_path.exists():
        raw_text = raw_path.read_text(encoding="utf-8")
    else:
        raw_text = fetch_iedb_predictions(allele, submissions)
        raw_path.write_text(raw_text, encoding="utf-8")
    retrieved = datetime.fromtimestamp(raw_path.stat().st_mtime, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return prediction_records_from_tsv(allele, submissions, raw_text, str(raw_path), retrieved_utc=retrieved)


def _model_cache(sample_rows: list[dict[str, object]]) -> tuple[dict[Path, object], dict[tuple[str, str], list[dict[str, object]]]]:
    valid = [
        row for row in sample_rows
        if str(row["sequence_layout_status"]).startswith("pass") and not bool(row["has_clash"])
    ]
    by_candidate: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in valid:
        by_candidate[(str(row["allele"]), str(row["candidate_id"]))].append(row)
    return {}, by_candidate


def _write_svg(path: Path, width: int, height: int, body: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#1f2933}.title{font-size:20px;font-weight:bold}.label{font-size:12px}.small{font-size:10px}</style>',
        *body,
        "</svg>",
        "",
    ]), encoding="utf-8")


def _generate_svg_figures(
    inventory: list[dict[str, object]],
    registers: list[dict[str, object]],
    pair_summaries: list[dict[str, object]],
    controls: list[dict[str, object]],
) -> None:
    """Dependency-free vector figures for the manuscript package."""
    figure_dir = OUT / "figures"
    labels = ["DRB1*15:01 pilot", "Fixed 50-peptide panel", "Three new alleles", "IEDB registers", "Anchor + controls", "Within-allele analysis"]
    body = ['<text x="600" y="36" text-anchor="middle" class="title">Figure 1. Fixed multi-allele transfer workflow</text>']
    for index, label in enumerate(labels):
        x = 24 + index * 196
        body.extend([
            f'<rect x="{x}" y="92" width="164" height="66" rx="9" fill="#e8f1fb" stroke="#2c5f8a"/>',
            f'<text x="{x+82}" y="130" text-anchor="middle" class="label">{escape(label)}</text>',
        ])
        if index < len(labels) - 1:
            body.append(f'<path d="M {x+164} 125 L {x+192} 125" stroke="#555" stroke-width="2"/><path d="M {x+192} 125 l -8 -5 v 10 z" fill="#555"/>')
    _write_svg(figure_dir / "figure_1_project_flow.svg", 1200, 230, body)

    body = ['<text x="550" y="34" text-anchor="middle" class="title">Figure 2. Technical and register eligibility coverage</text>']
    for index, allele in enumerate(ALLELES):
        complete = sum(row["allele"] == allele and row["download_status"] == "complete" for row in inventory)
        resolved = sum(row["allele"] == allele and row["register_status"] == "resolved_unique_fully_contained" for row in registers)
        x = 90 + index * 150
        body.extend([
            f'<rect x="{x}" y="{390-complete*6}" width="44" height="{complete*6}" fill="#3977a8"/>',
            f'<rect x="{x}" y="90" width="44" height="{(50-complete)*6}" fill="#d7dce2"/>',
            f'<text x="{x+22}" y="418" text-anchor="middle" class="small">{escape(allele.replace("HLA-", ""))}</text>',
            f'<text x="{x+22}" y="{382-complete*6}" text-anchor="middle" class="small">{complete}/50</text>',
        ])
        x2 = 620 + index * 150
        body.extend([
            f'<rect x="{x2}" y="{390-resolved*6}" width="44" height="{resolved*6}" fill="#4d9c74"/>',
            f'<rect x="{x2}" y="90" width="44" height="{(50-resolved)*6}" fill="#e5b66b"/>',
            f'<text x="{x2+22}" y="418" text-anchor="middle" class="small">{escape(allele.replace("HLA-", ""))}</text>',
            f'<text x="{x2+22}" y="{382-resolved*6}" text-anchor="middle" class="small">{resolved}/50</text>',
        ])
    body.extend(['<text x="270" y="70" text-anchor="middle" class="label">AF3 canonical jobs</text>', '<text x="800" y="70" text-anchor="middle" class="label">Unique fully contained registers</text>'])
    _write_svg(figure_dir / "figure_2_qc_register_coverage.svg", 1100, 450, body)

    anchors = {row["allele"]: row for row in pair_summaries if row["analysis_role"] == "primary_anchor"}
    body = ['<text x="475" y="34" text-anchor="middle" class="title">Figure 3. Anchor geometry and frozen-control run status</text>',
            '<line x1="70" y1="390" x2="920" y2="390" stroke="#333"/><line x1="70" y1="70" x2="70" y2="390" stroke="#333"/>',
            '<text x="20" y="230" transform="rotate(-90 20 230)" class="label">Exposed-position Cα RMSD (Å)</text>']
    for index, allele in enumerate(ALLELES):
        x = 210 + index * 270
        row = anchors.get(allele)
        body.append(f'<text x="{x}" y="418" text-anchor="middle" class="label">{escape(allele.replace("HLA-", ""))}</text>')
        if row and row["geometry_status"] == "complete":
            value = float(row["exposed_ca_rmsd_median_A"])
            y = 390 - min(value, 3.0) / 3.0 * 300
            body.extend([f'<circle cx="{x-45}" cy="{y:.1f}" r="8" fill="#bd3f32"/>', f'<text x="{x-45}" y="{y-13:.1f}" text-anchor="middle" class="small">{value:.2f}</text>'])
        for offset in (0, 30, 60):
            body.append(f'<path d="M {x+offset-5} 380 l 10 10 M {x+offset+5} 380 l -10 10" stroke="#777" stroke-width="2"/>')
    body.append('<text x="475" y="455" text-anchor="middle" class="small">Gray crosses: frozen controls awaiting fixed-seed AF3 geometry</text>')
    _write_svg(figure_dir / "figure_3_anchor_controls_status.svg", 950, 480, body)

    panel = read_csv(INPUT_PACKAGE / "peptide_panel_manifest.csv")
    ebv_ids = sorted(row["candidate_id"] for row in panel if row["arm_group"] == "EBV")
    human_ids = sorted(row["candidate_id"] for row in panel if row["arm_group"] == "CNS/self")
    by_key = {(str(row["allele"]), str(row["ebv_candidate_id"]), str(row["human_candidate_id"])): row for row in pair_summaries}
    finite = [float(row["exposed_ca_rmsd_median_A"]) for row in pair_summaries if row["geometry_status"] == "complete"]
    vmax = float(np.quantile(finite, 0.95)) if finite else 1.0
    body = ['<text x="760" y="34" text-anchor="middle" class="title">Figure 4. Exploratory within-allele EBV–CNS geometry</text>']
    cell = 16
    for panel_index, allele in enumerate(ALLELES):
        x0 = 65 + panel_index * 495
        y0 = 85
        body.append(f'<text x="{x0+200}" y="66" text-anchor="middle" class="label">{escape(allele.replace("HLA-", ""))}</text>')
        for i, ebv_id in enumerate(ebv_ids):
            for j, human_id in enumerate(human_ids):
                row = by_key.get((allele, ebv_id, human_id))
                if not row or row["geometry_status"] != "complete":
                    color = "#e7e7e7"
                else:
                    scaled = min(1.0, float(row["exposed_ca_rmsd_median_A"]) / vmax)
                    red = int(30 + 210 * scaled)
                    green = int(145 - 90 * scaled)
                    blue = int(175 - 120 * scaled)
                    color = f"#{red:02x}{max(green,0):02x}{max(blue,0):02x}"
                body.append(f'<rect x="{x0+j*cell}" y="{y0+i*cell}" width="{cell}" height="{cell}" fill="{color}"/>')
        body.extend([f'<text x="{x0+200}" y="{y0+425}" text-anchor="middle" class="small">CNS/self peptide index</text>',
                     f'<text x="{x0-30}" y="{y0+200}" transform="rotate(-90 {x0-30} {y0+200})" text-anchor="middle" class="small">EBV peptide index</text>'])
    body.append(f'<text x="760" y="540" text-anchor="middle" class="small">Color scale capped at the 95th percentile ({vmax:.2f} Å); gray denotes unresolved register or missing model.</text>')
    _write_svg(figure_dir / "figure_4_exploratory_heatmaps.svg", 1530, 565, body)


def _generate_figures(
    inventory: list[dict[str, object]],
    registers: list[dict[str, object]],
    pair_summaries: list[dict[str, object]],
    controls: list[dict[str, object]],
) -> None:
    try:
        import matplotlib
    except ModuleNotFoundError:
        _generate_svg_figures(inventory, registers, pair_summaries, controls)
        return
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = OUT / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.axis("off")
    labels = [
        "DRB1*15:01\npilot",
        "Fixed 50-peptide\ntransfer panel",
        "3 new alleles\n150 AF3 jobs",
        "Allele-specific\nIEDB registers",
        "Anchor + frozen\ncontrols",
        "Within-allele\nmanuscript analysis",
    ]
    xs = np.linspace(0.06, 0.94, len(labels))
    for index, (x, label) in enumerate(zip(xs, labels)):
        ax.text(x, 0.55, label, ha="center", va="center", fontsize=9,
                bbox={"boxstyle": "round,pad=0.5", "facecolor": "#e8f1fb", "edgecolor": "#2c5f8a"})
        if index < len(labels) - 1:
            ax.annotate("", xy=(xs[index + 1] - 0.075, 0.55), xytext=(x + 0.075, 0.55),
                        arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.5})
    ax.set_title("Figure 1. Fixed multi-allele transfer workflow", fontsize=13, weight="bold")
    fig.tight_layout()
    fig.savefig(figure_dir / "figure_1_project_flow.png", dpi=220)
    fig.savefig(figure_dir / "figure_1_project_flow.svg")
    plt.close(fig)

    completed = [sum(row["allele"] == allele and row["download_status"] == "complete" for row in inventory) for allele in ALLELES]
    missing = [50 - value for value in completed]
    resolved = [sum(row["allele"] == allele and row["register_status"] == "resolved_unique_fully_contained" for row in registers) for allele in ALLELES]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    x = np.arange(3)
    labels_short = [a.replace("HLA-", "") for a in ALLELES]
    axes[0].bar(x, completed, label="Complete", color="#3977a8")
    axes[0].bar(x, missing, bottom=completed, label="Missing", color="#d7dce2")
    axes[0].set_xticks(x, labels_short, rotation=15)
    axes[0].set_ylim(0, 52)
    axes[0].set_ylabel("Canonical jobs")
    axes[0].legend(frameon=False)
    axes[0].set_title("AF3 coverage")
    axes[1].bar(x, resolved, color="#4d9c74")
    axes[1].bar(x, [50-r for r in resolved], bottom=resolved, color="#e5b66b")
    axes[1].set_xticks(x, labels_short, rotation=15)
    axes[1].set_ylim(0, 52)
    axes[1].set_ylabel("Panel peptides")
    axes[1].set_title("Unique contained register")
    fig.suptitle("Figure 2. Technical and register eligibility coverage", weight="bold")
    fig.tight_layout()
    fig.savefig(figure_dir / "figure_2_qc_register_coverage.png", dpi=220)
    fig.savefig(figure_dir / "figure_2_qc_register_coverage.svg")
    plt.close(fig)

    anchors = {
        row["allele"]: row for row in pair_summaries
        if row["ebv_candidate_id"] == ANCHOR_EBV and row["human_candidate_id"] == ANCHOR_HUMAN
    }
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for x_index, allele in enumerate(ALLELES):
        row = anchors.get(allele)
        if row and row["geometry_status"] == "complete":
            ax.scatter([x_index - 0.18], [float(row["exposed_ca_rmsd_median_A"])], s=75, color="#bd3f32", label="Anchor" if x_index == 0 else None, zorder=3)
        selected = [r for r in controls if r["allele"] == allele]
        for offset, _ in zip((-0.03, 0.12, 0.27), selected):
            ax.scatter([x_index + offset], [0], marker="x", s=55, color="#777777", label="Frozen control; AF3 pending" if x_index == 0 and offset == -0.03 else None)
    ax.set_xticks(range(3), labels_short)
    ax.set_ylabel("Exposed-position Cα RMSD (Å)")
    ax.set_title("Figure 3. Anchor geometry and frozen-control run status")
    ax.text(0.5, 0.02, "Control geometry remains pending the fixed-seed robustness batch", transform=ax.transAxes, ha="center", fontsize=9)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(figure_dir / "figure_3_anchor_controls_status.png", dpi=220)
    fig.savefig(figure_dir / "figure_3_anchor_controls_status.svg")
    plt.close(fig)

    panel = read_csv(INPUT_PACKAGE / "peptide_panel_manifest.csv")
    ebv_ids = sorted(row["candidate_id"] for row in panel if row["arm_group"] == "EBV")
    human_ids = sorted(row["candidate_id"] for row in panel if row["arm_group"] == "CNS/self")
    by_key = {(str(row["allele"]), str(row["ebv_candidate_id"]), str(row["human_candidate_id"])): row for row in pair_summaries}
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharex=True, sharey=True)
    finite_values = [float(row["exposed_ca_rmsd_median_A"]) for row in pair_summaries if row["geometry_status"] == "complete"]
    vmax = float(np.quantile(finite_values, 0.95)) if finite_values else 1.0
    image = None
    for ax, allele in zip(axes, ALLELES):
        matrix = np.full((25, 25), np.nan)
        for i, ebv_id in enumerate(ebv_ids):
            for j, human_id in enumerate(human_ids):
                row = by_key.get((allele, ebv_id, human_id))
                if row and row["geometry_status"] == "complete":
                    matrix[i, j] = float(row["exposed_ca_rmsd_median_A"])
        image = ax.imshow(matrix, aspect="auto", cmap="viridis_r", vmin=0, vmax=vmax)
        ax.set_title(allele.replace("HLA-", ""))
        ax.set_xlabel("CNS/self peptide index")
    axes[0].set_ylabel("EBV peptide index")
    if image is not None:
        fig.colorbar(image, ax=axes, label="Median exposed-position Cα RMSD (Å)", shrink=0.82)
    fig.suptitle("Figure 4. Exploratory within-allele EBV–CNS geometry", weight="bold")
    fig.subplots_adjust(left=0.06, right=0.91, bottom=0.12, top=0.86, wspace=0.08)
    fig.savefig(figure_dir / "figure_4_exploratory_heatmaps.png", dpi=220)
    fig.savefig(figure_dir / "figure_4_exploratory_heatmaps.svg")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expected = read_csv(INPUT_PACKAGE / "job_manifest_150.csv")
    panel = read_csv(INPUT_PACKAGE / "peptide_panel_manifest.csv")
    hla_rows = read_csv(INPUT_PACKAGE / "hla_sequence_manifest.csv")
    dra = next(row["sequence"] for row in hla_rows if row["chain"] == "HLA-DRA")
    drb_by_allele = {
        f"HLA-{row['allele_or_name'][:10]}": row["sequence"]
        for row in hla_rows if row["chain"] == "HLA-DRB"
    }
    result_roots = sorted(path for path in DOWNLOADS.glob(RESULT_GLOB) if path.is_dir())
    discovered = discover_download_jobs(result_roots)
    inventory, duplicates = inventory_against_manifest(expected, discovered)
    write_csv(OUT / "model_inventory_150.csv", inventory)
    write_csv(OUT / "duplicate_run_sensitivity_manifest.csv", duplicates, ["canonical_request_name", "request_name", "job_directory_path", "handling"])
    missing = [row for row in inventory if row["download_status"] != "complete"]
    write_csv(OUT / "unresolved_missingness.csv", missing, list(inventory[0]))

    sample_rows: list[dict[str, object]] = []
    job_rows: list[dict[str, object]] = []
    for row in inventory:
        if row["download_status"] != "complete":
            continue
        samples, job = analyze_job(row, dra, drb_by_allele[str(row["allele"])])
        sample_rows.extend(samples)
        job_rows.append(job)
    write_csv(OUT / "af3_sample_qc_750.csv", sample_rows)
    write_csv(OUT / "af3_job_qc_150.csv", job_rows)

    flanks = {row["candidate_id"]: row for row in read_csv(FLANK_SOURCE)}
    panel_submissions = build_prediction_submissions(panel, flanks)
    write_csv(OUT / "iedb_panel_submission_manifest.csv", panel_submissions)
    panel_predictions: list[dict[str, object]] = []
    for allele in ALLELES:
        code = ALLELE_CODES[allele]
        raw_path = OUT / f"iedb_panel_{code}_raw.tsv"
        panel_predictions.extend(_saved_raw_prediction_records(allele, panel_submissions, raw_path))
    panel_metadata = {row["candidate_id"]: row for row in panel}
    panel_registers = [{**row, "arm_group": panel_metadata[str(row["candidate_id"])]["arm_group"], "source_antigen": panel_metadata[str(row["candidate_id"])]["source_antigen"]} for row in panel_predictions]
    write_csv(OUT / "allele_specific_register_table_150.csv", panel_registers)

    study_peptides = {row["peptide_sequence"] for row in panel}
    forbidden = ("myelin basic", "proteolipid", "oligodendrocyte glycoprotein", "glialcam")
    control_pool = [
        row for row in read_csv(CONTROL_SOURCE)
        if row["selection_status"] == "eligible_pre_prediction"
        and abs(int(row["peptide_length"]) - len(panel_metadata[ANCHOR_HUMAN]["peptide_sequence"])) <= 1
        and row["peptide"] not in study_peptides
        and not any(term in row["source_antigen_name"].lower() for term in forbidden)
    ]
    control_submissions = build_prediction_submissions(control_pool, {})
    write_csv(OUT / "iedb_control_submission_manifest.csv", control_submissions)
    control_predictions: list[dict[str, object]] = []
    for allele in ALLELES:
        code = ALLELE_CODES[allele]
        raw_path = OUT / f"iedb_controls_{code}_raw.tsv"
        predicted = _saved_raw_prediction_records(allele, control_submissions, raw_path)
        source_by_id = {row["candidate_id"]: row for row in control_pool}
        control_predictions.extend([{**source_by_id[str(row["candidate_id"])], **row} for row in predicted])
    write_csv(OUT / "control_binding_register_predictions.csv", control_predictions)

    register_by_key = {(str(row["allele"]), str(row["candidate_id"])): row for row in panel_predictions}
    controls: list[dict[str, object]] = []
    anchor_status: list[dict[str, object]] = []
    for allele in ALLELES:
        ebv_register = register_by_key[(allele, ANCHOR_EBV)]
        human_register = register_by_key[(allele, ANCHOR_HUMAN)]
        evaluable = ebv_register["register_status"] == human_register["register_status"] == "resolved_unique_fully_contained"
        selected: list[dict[str, object]] = []
        if evaluable:
            rows = [row for row in control_predictions if row["allele"] == allele]
            selected = select_score_blind_controls(
                panel_metadata[ANCHOR_HUMAN]["peptide_sequence"],
                str(human_register["binding_rank_bin"]), rows, limit=3,
            )
            selected = [{"allele": allele, **row} for row in selected]
            controls.extend(selected)
        anchor_status.append({
            "allele": allele,
            "anchor_ebv_register_status": ebv_register["register_status"],
            "anchor_mbp_register_status": human_register["register_status"],
            "anchor_primary_endpoint_status": "evaluable" if evaluable and len(selected) == 3 else "not_evaluable",
            "frozen_control_count": len(selected),
            "non_evaluable_reason": "" if evaluable and len(selected) == 3 else "unresolved anchor register" if not evaluable else "fewer than three exact-bin controls",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    write_csv(OUT / "anchor_evaluability.csv", anchor_status)
    write_csv(OUT / "frozen_control_manifest.csv", controls, list(controls[0]) if controls else ["allele"])

    entities: list[dict[str, object]] = []
    drb_sequences = {allele: drb_by_allele[allele] for allele in ALLELES}
    for status in anchor_status:
        allele = str(status["allele"])
        if status["anchor_primary_endpoint_status"] != "evaluable":
            continue
        entities.extend([
            {"allele": allele, "entity_id": ANCHOR_EBV, "entity_role": "anchor_ebv", "peptide": panel_metadata[ANCHOR_EBV]["peptide_sequence"], "dra_sequence": dra, "drb_sequence": drb_sequences[allele]},
            {"allele": allele, "entity_id": ANCHOR_HUMAN, "entity_role": "anchor_mbp", "peptide": panel_metadata[ANCHOR_HUMAN]["peptide_sequence"], "dra_sequence": dra, "drb_sequence": drb_sequences[allele]},
        ])
        for control in (row for row in controls if row["allele"] == allele):
            entities.append({"allele": allele, "entity_id": control["candidate_id"], "entity_role": "frozen_non_cns_human_control", "peptide": control["peptide"], "dra_sequence": dra, "drb_sequence": drb_sequences[allele]})
    robustness_jobs, robustness_manifest = build_robustness_jobs(entities)
    _write_json(OUT / "alphafold_robustness_fixed_seed_jobs.json", robustness_jobs)
    write_csv(OUT / "alphafold_robustness_job_manifest.csv", robustness_manifest, list(robustness_manifest[0]) if robustness_manifest else ["job_name"])

    pair_rows = [row for allele in ALLELES for row in build_pair_universe(allele, panel, register_by_key)]
    _, samples_by_candidate = _model_cache(sample_rows)
    parsed_models: dict[Path, object] = {}
    geometry_rows: list[dict[str, object]] = []
    pair_summaries: list[dict[str, object]] = []
    for pair in pair_rows:
        allele = str(pair["allele"])
        left_id = str(pair["ebv_candidate_id"])
        right_id = str(pair["human_candidate_id"])
        values: list[dict[str, float]] = []
        left_samples = samples_by_candidate.get((allele, left_id), [])
        right_samples = samples_by_candidate.get((allele, right_id), [])
        left_register = register_by_key[(allele, left_id)]
        right_register = register_by_key[(allele, right_id)]
        if pair["register_eligible"] and left_samples and right_samples:
            for left_sample in left_samples:
                left_path = Path(str(left_sample["cif_path"]))
                if left_path not in parsed_models:
                    parsed_models[left_path] = parse_mmcif(left_path)
                for right_sample in right_samples:
                    right_path = Path(str(right_sample["cif_path"]))
                    if right_path not in parsed_models:
                        parsed_models[right_path] = parse_mmcif(right_path)
                    metrics = same_register_geometry_from_models(
                        parsed_models[left_path], parsed_models[right_path],
                        int(left_register["core_start_1_based"]), int(right_register["core_start_1_based"]),
                    )
                    values.append(metrics)
                    geometry_rows.append({
                        "allele": allele,
                        "pair_id": pair["pair_id"],
                        "ebv_candidate_id": left_id,
                        "human_candidate_id": right_id,
                        "ebv_server_seed": left_sample["server_seed"],
                        "ebv_sample_index": left_sample["sample_index"],
                        "human_server_seed": right_sample["server_seed"],
                        "human_sample_index": right_sample["sample_index"],
                        **{key: round(value, 6) for key, value in metrics.items()},
                        "claim_boundary": CLAIM_BOUNDARY,
                    })
        summary = {**pair}
        if values:
            exposed = [row["exposed_ca_rmsd_A"] for row in values]
            low, high = _bootstrap_ci(exposed, _stable_seed(str(pair["pair_id"])))
            summary.update({
                "geometry_status": "complete",
                "valid_sample_combination_count": len(values),
                "exposed_ca_rmsd_median_A": round(median(exposed), 6),
                "exposed_ca_rmsd_bootstrap_95_low_A": round(low, 6),
                "exposed_ca_rmsd_bootstrap_95_high_A": round(high, 6),
                "full_core_ca_rmsd_median_A": round(median(row["full_core_ca_rmsd_A"] for row in values), 6),
                "anchor_ca_rmsd_median_A": round(median(row["anchor_ca_rmsd_A"] for row in values), 6),
                **direct_register_sequence_metrics(str(left_register["predicted_core_peptide"]), str(right_register["predicted_core_peptide"])),
            })
        else:
            summary.update({
                "geometry_status": "missing_model" if pair["register_eligible"] else "excluded_unresolved_register",
                "valid_sample_combination_count": 0,
                "exposed_ca_rmsd_median_A": "",
                "exposed_ca_rmsd_bootstrap_95_low_A": "",
                "exposed_ca_rmsd_bootstrap_95_high_A": "",
                "full_core_ca_rmsd_median_A": "",
                "anchor_ca_rmsd_median_A": "",
                "full_core_property_similarity_mean": "",
                "anchor_property_similarity_mean": "",
                "exposed_property_similarity_mean": "",
                "exposed_exact_identity_count": "",
            })
        pair_summaries.append(summary)
    for allele in ALLELES:
        completed = sorted(
            (row for row in pair_summaries if row["allele"] == allele and row["geometry_status"] == "complete"),
            key=lambda row: (float(row["exposed_ca_rmsd_median_A"]), str(row["pair_id"])),
        )
        for rank, row in enumerate(completed, start=1):
            row["within_allele_exploratory_rank"] = rank
    for row in pair_summaries:
        row.setdefault("within_allele_exploratory_rank", "")
    write_csv(OUT / "pair_universe_1875.csv", pair_summaries)
    write_csv(OUT / "geometry_sample_combinations.csv", geometry_rows, list(geometry_rows[0]) if geometry_rows else ["allele"])

    anchor_rows = [row for row in pair_summaries if row["analysis_role"] == "primary_anchor"]
    write_csv(OUT / "primary_anchor_geometry_summary.csv", anchor_rows)
    robustness_status = [{
        "allele": row["allele"],
        "analysis_status": "awaiting_fixed_seed_alphafold_downloads" if row["anchor_primary_endpoint_status"] == "evaluable" else "not_evaluable",
        "prepared_robustness_job_count": sum(job["allele"] == row["allele"] for job in robustness_manifest),
        "equal_control_weighting": True,
        "leave_one_control_out": True,
        "technical_bootstrap": True,
        "empirical_tail_fraction": "pending",
        "p_value": "",
        "claim_boundary": CLAIM_BOUNDARY,
    } for row in anchor_status]
    write_csv(OUT / "robustness_summary.csv", robustness_status)

    _generate_figures(inventory, panel_registers, pair_summaries, controls)

    complete_count = sum(row["download_status"] == "complete" for row in inventory)
    register_counts = {allele: sum(row["allele"] == allele and row["register_status"] == "resolved_unique_fully_contained" for row in panel_registers) for allele in ALLELES}
    anchor_lines = []
    for row in anchor_rows:
        metric = f"{row['exposed_ca_rmsd_median_A']} Å" if row["geometry_status"] == "complete" else row["geometry_status"]
        anchor_lines.append(f"- {row['allele']}: {metric}; primary robustness status is pending fixed-seed controls.")
    (OUT / "MANUSCRIPT_METHODS_RESULTS.md").write_text("\n".join([
        "# Multi-allele EBV–MS pMHC analysis: manuscript-ready working text",
        "",
        "## Methods",
        "",
        "A fixed panel of 25 EBV and 25 CNS/self peptides was modeled independently with HLA-DRB1*13:03, HLA-DRB1*03:01, and HLA-DRB1*08:01. The prespecified primary endpoint was transfer of the EBV_TCELL_950–HUMAN_MYELIN_112214 EBNA1–MBP pair; all other 625 within-allele combinations were exploratory. AlphaFold Server outputs were accepted only when all five CIF, confidence-summary, and full-data files were present and exact DRA, DRB, and peptide sequences agreed with the request. The highest-ranked clash-free sample represented each job for descriptive QC; geometry used all valid clash-free sample combinations without post-hoc confidence thresholds.",
        "",
        "IEDB recommended-binding predictions were generated separately for each allele. Exact seq_num values, raw responses, percentile ranks, predicted cores, and binding-rank bins were retained. The two 10-residue EBV peptides used source-verified natural flanks, but a structural register was eligible only when the predicted P1–P9 core was unique and fully contained in the modeled peptide.",
        "",
        "Structures were superposed within allele using equivalent HLA-DRA and HLA-DRB groove Cα atoms. P1–P9 geometry was summarized for the full core, anchor positions P1/P4/P6/P9, and candidate TCR-facing positions P2/P3/P5/P7/P8. The primary metric was the median exposed-position Cα RMSD across valid sample combinations. No raw metric was pooled across alleles.",
        "",
        "Three non-CNS human controls per evaluable allele were chosen before inspecting any AlphaFold geometry. Controls matched MBP length within one residue and the allele-specific binding-rank bin and were ordered by amino-acid composition distance, peptide length, and numeric IEDB identifier. Robustness inputs use fixed seeds 104729 and 104759.",
        "",
        "## Current results",
        "",
        f"The current download set contains {complete_count}/150 canonical complete jobs and {len(duplicates)} duplicate sensitivity run. Missing jobs were not imputed. Exact sequence QA was performed for {len(job_rows)} jobs and {len(sample_rows)} samples.",
        "",
        *[f"- {allele}: {register_counts[allele]}/50 modeled peptides had a unique fully contained predicted core." for allele in ALLELES],
        "",
        "Primary anchor geometry from the discovery-seed matrix:",
        "",
        *anchor_lines,
        "",
        "These results describe computational pMHC geometry only. They do not establish natural presentation, TCR binding, T-cell activation, cross-reactivity, molecular mimicry, or an MS disease mechanism. The primary anchor-versus-control result remains incomplete until the fixed-seed robustness jobs are modeled and analyzed.",
        "",
    ]), encoding="utf-8")

    (OUT / "README.md").write_text("\n".join([
        "# Multi-allele EBV–MS analysis package (2026-08-21)",
        "",
        f"- Canonical AF3 jobs complete: **{complete_count}/150**",
        f"- Sample QC rows currently available: **{len(sample_rows)}/750**",
        f"- Allele-specific panel predictions: **{len(panel_registers)}/150**",
        f"- Pair-universe rows: **{len(pair_summaries)}/1,875**",
        f"- Geometry sample combinations: **{len(geometry_rows):,}**",
        f"- Frozen controls: **{len(controls)}**",
        f"- Fixed-seed robustness jobs prepared: **{len(robustness_jobs)}/30 maximum**",
        "",
        "The package is reproducible from saved inputs and raw IEDB responses. Downloads are referenced read-only; all derived files live here. The 20-job retry JSON is copied into this package when available, but a prepared/uploaded JSON is not a completed AlphaFold run.",
        "",
        CLAIM_BOUNDARY,
        "",
        "Reproduce:",
        "",
        "```bash",
        "PYTHONPATH=src python3 src/run_multiallele_manuscript_analysis.py",
        "PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py'",
        "```",
        "",
    ]), encoding="utf-8")
    if RETRY_SOURCE.exists():
        (OUT / "alphafold_retry_20_jobs.json").write_bytes(RETRY_SOURCE.read_bytes())
    (OUT / "ALPHAFOLD_EXTERNAL_STATUS.md").write_text("\n".join([
        "# AlphaFold Server external status",
        "",
        "- The original 20-job retry was partially completed by AlphaFold Server.",
        f"- Thirteen additional retry outputs were downloaded on 2026-08-22; the canonical matrix is now **{complete_count}/150**.",
        f"- The remaining **{len(missing)}** jobs are treated as persistent missingness and are never imputed.",
        "- The 30-job fixed-seed robustness JSON is prepared but intentionally remains unsubmitted until its frozen control manifest is reviewed.",
        "- A prepared JSON is not counted as a completed AlphaFold job; only downloaded folders passing reinventory enter the model matrix.",
        "",
    ]), encoding="utf-8")

    checksum_targets = sorted(path for path in OUT.rglob("*") if path.is_file() and path.name != "SHA256SUMS.txt")
    (OUT / "SHA256SUMS.txt").write_text("\n".join(f"{sha256_file(path)}  {path.relative_to(OUT)}" for path in checksum_targets) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUT),
        "complete_jobs": complete_count,
        "sample_rows": len(sample_rows),
        "panel_predictions": len(panel_registers),
        "pair_rows": len(pair_summaries),
        "geometry_rows": len(geometry_rows),
        "controls": len(controls),
        "robustness_jobs": len(robustness_jobs),
    }, indent=2))


if __name__ == "__main__":
    main()
