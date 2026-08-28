"""Record IEDB MHC-II binding/register hypotheses for the BALF5--MBP control.

This is a reproducible computational context check for the two DR15-haplotype
class-II molecules used by the established BALF5--MBP calibration system.  It
does not establish a shared TCR surface or biological cross-reactivity.
"""

from __future__ import annotations

import csv
import urllib.parse
import urllib.request
from pathlib import Path

from premeeting_rigor import parse_iedb_mhcii_tsv


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "register_sensitivity"
IEDB_ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
METHOD = "recommended_binding"
CONTROLS = [
    {
        "control_id": "BALF5_on_DRB5_0101",
        "peptide": "TGGVYHFVKKHVHES",
        "allele": "HLA-DRB5*01:01",
        "interpretation": "DR2a computational control context; not directly comparable to a DRB1-only screen without mentor-approved allele logic.",
    },
    {
        "control_id": "BALF5_on_DRB1_1501",
        "peptide": "TGGVYHFVKKHVHES",
        "allele": "HLA-DRB1*15:01",
        "interpretation": "DR2b computational context only; not evidence that the experimental BALF5 complex is presented by DRB1*15:01.",
    },
    {
        "control_id": "MBP_on_DRB1_1501",
        "peptide": "ENPVVHFFKNIVTPR",
        "allele": "HLA-DRB1*15:01",
        "interpretation": "DR2b computational control context; not a TCR-recognition result.",
    },
]


def fetch(control: dict[str, str]) -> str:
    fasta = f">{control['control_id']}\n{control['peptide']}"
    body = urllib.parse.urlencode({
        "method": METHOD,
        "sequence_text": fasta,
        "allele": control["allele"],
        "length": "asis",
    }).encode("utf-8")
    request = urllib.request.Request(
        IEDB_ENDPOINT,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    output_rows = []
    raw_sections = []
    for control in CONTROLS:
        raw = fetch(control)
        rows = parse_iedb_mhcii_tsv(raw)
        if len(rows) != 1 or rows[0]["peptide"] != control["peptide"]:
            raise ValueError(f"Unexpected IEDB response for {control['control_id']}")
        result = rows[0]
        output_rows.append({
            **control,
            "method_requested": METHOD,
            "endpoint": IEDB_ENDPOINT,
            "predicted_core_peptide": result["core_peptide"],
            "predicted_ic50_nM": result.get("ic50", ""),
            "predicted_percentile_rank": result["rank"],
            "claim_boundary": "Computational binding/register hypothesis only; does not test cross-allele structural equivalence, TCR binding, activation, or disease mechanism.",
        })
        raw_sections.append(f"# {control['control_id']}\n{raw.rstrip()}\n")
    with (OUT / "iedb_mhcii_positive_control_raw.txt").open("w", encoding="utf-8") as handle:
        handle.write("\n".join(raw_sections))
    with (OUT / "positive_control_allele_context.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {OUT / 'positive_control_allele_context.csv'}")


if __name__ == "__main__":
    main()
