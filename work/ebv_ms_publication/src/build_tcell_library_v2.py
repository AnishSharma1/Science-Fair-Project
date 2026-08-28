"""Build the dated EBV-MS T-cell library and positive-recovery package.

This is an additive V2 build.  It reads existing V1 tables, IEDB, UniProt and
RCSB references, but writes only under processed/tcell_library_v2_2026-08-22.
Prepared AlphaFold inputs are not represented as completed models.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import shutil
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

from premeeting_rigor import binding_rank_bin
from tcell_library_v2 import (
    ALLELES,
    CALIBRATION_SEEDS,
    LIBRARY_VERSION,
    PANEL_VERSION,
    build_af3_jobs,
    build_calibration_comparison_universe,
    build_discovery_pairs,
    build_native_calibration_jobs,
    freeze_panel,
    select_native_controls,
    validate_registry,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "tcell_library_v2_2026-08-22"
AF3_OUT = OUT / "alphafold_jobs"
RAW_OUT = OUT / "raw_responses"
SEARCH_FREEZE_DATE = "2026-08-22"
IEDB_QUERY_ENDPOINT = "https://query-api.iedb.org/tcell_search"
IEDB_PREDICTION_ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb"
CLAIM_BOUNDARY = (
    "Computational pMHC geometry only; not evidence of presentation, TCR binding, "
    "activation, cross-reactivity, molecular mimicry, or MS disease mechanism."
)


DRA_SEQUENCE = (
    "EHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKAN"
    "LEIMTKRSNYTPITNVPPEVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLRNGKPVTTGVSETV"
    "FLPREDHLFRKFHYLPFLPSTEDVYDCRVEHWGLDEPLLKHWEFD"
)
DRB1_1501_SEQUENCE = (
    "DTRPRFLWQPKRECHFFNGTERVRFLDRYFYNQEESVRFDSDVGEFRAVTELGRPDAEYWNSQKDIL"
    "EQARAAVDTYCRHNYGVVESFTVQRRVQPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFLNG"
    "QEEKAGMVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWRA"
)
DRB5_0101_SEQUENCE = (
    "DTRPRFLQQDKYECHFFNGTERVRFLHRDIYNQEEDLRFDSDVGEYRAVTELGRPDAEYWNSQKDFL"
    "EDRRAAVDTYCRHNYGVGESFTVQRRVEPKVTVYPARTQTLQHHNLLVCSVNGFYPGSIEVRWFRNS"
    "QEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWRA"
)


SELF_PROTEINS = {
    "MBP": {"accession": "P02686", "name": "myelin basic protein"},
    "PLP1": {"accession": "P60201", "name": "proteolipid protein 1"},
    "MOG": {"accession": "Q16653", "name": "myelin oligodendrocyte glycoprotein"},
    "ANO2": {"accession": "Q9NQ90", "name": "anoctamin-2"},
    "CNP": {"accession": "P09543", "name": "2',3'-cyclic-nucleotide 3'-phosphodiesterase"},
    "MOBP": {"accession": "Q13875", "name": "myelin-associated oligodendrocyte basic protein"},
    "MAG": {"accession": "P20916", "name": "myelin-associated glycoprotein"},
    "CNTN2": {"accession": "Q02246", "name": "contactin-2"},
    "CLDN11": {"accession": "O75508", "name": "claudin-11"},
    "TALDO1": {"accession": "P37837", "name": "transaldolase 1"},
    "CRYAB": {"accession": "P02511", "name": "alpha-crystallin B chain"},
}


def is_allowed_non_cns_control_source(accession: str, source_name: str) -> bool:
    """Return False for every study/CNS protein, including legacy names."""
    base_accession = re.split(r"[.-]", accession.upper())[0]
    study_accessions = {value["accession"] for value in SELF_PROTEINS.values()} | {"Q14CZ8"}
    if base_accession in study_accessions:
        return False
    normalized_name = source_name.lower()
    excluded_terms = {
        "mbp", "myelin", "proteolipid", "mog", "glialcam", "glial cam",
        "claudin", "mobp", "crystallin", "anoctamin", "contactin",
        "transaldolase", "cyclic-nucleotide phosphodiesterase",
    }
    return not any(term in normalized_name for term in excluded_terms)

EBV_SYMBOL_BY_ACCESSION = {
    "P03211": "EBNA1", "P12978": "EBNA2", "P12977": "EBNA3A",
    "P03203": "EBNA3B", "P03204": "EBNA3C", "P03230": "LMP1",
    "P13285": "LMP2", "P03206": "BZLF1", "P03198": "BALF5",
    "P03191": "BMRF1", "P03182": "BHRF1", "P03188": "BALF4_gB",
    "P03231": "BXLF2_gH", "P03200": "BLLF1_gp350", "P03187": "BFRF3",
    "P03186": "BPLF1", "P03228": "BARF1", "P03179": "BNRF1",
    "P03226": "BcLF1", "P03227": "BALF2", "P03185": "BFRF1",
    "P0CAP6": "BaRF1", "Q8AZK7": "EBNA-LP",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Optional[Sequence[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fields:
        raise ValueError(f"refusing to write empty table without fields: {path}")
    fieldnames = list(fields) if fields else sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_epitope_coordinates(summary: str) -> tuple[Optional[int], Optional[int]]:
    matches = re.findall(r"\((\d+)-(\d+)\)", summary or "")
    if not matches:
        return None, None
    start, end = matches[-1]
    return int(start), int(end)


def _accession_from_iri(value: Any) -> str:
    text = "|".join(map(str, value)) if isinstance(value, list) else str(value or "")
    match = re.search(r"UNIPROT:([A-Z0-9]+)", text)
    return match.group(1) if match else ""


def normalize_iedb_tcell_rows(rows: Sequence[dict[str, Any]], kingdom: str) -> list[dict[str, Any]]:
    """Collapse assay duplicates while retaining distinct exact peptide records."""
    grouped: dict[tuple[str, Optional[int], Optional[int], str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        sequence = str(row.get("linear_sequence") or "").upper()
        if not 11 <= len(sequence) <= 30 or re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
            continue
        if str(row.get("epitope_structure_defined", "")) != "Exact Epitope":
            continue
        accession = _accession_from_iri(row.get("parent_source_antigen_iri"))
        start, end = parse_epitope_coordinates(str(row.get("epitope_summary", "")))
        grouped[(accession, start, end, sequence)].append(row)
    output = []
    ordered_groups = sorted(
        grouped.items(),
        key=lambda item: (
            item[0][0],
            item[0][1] if item[0][1] is not None else 10**9,
            item[0][2] if item[0][2] is not None else 10**9,
            item[0][3],
        ),
    )
    for (accession, start, end, sequence), records in ordered_groups:
        first = records[0]
        symbol = EBV_SYMBOL_BY_ACCESSION.get(accession, accession or "unresolved_source")
        if kingdom == "human_self":
            symbol = next((key for key, value in SELF_PROTEINS.items() if value["accession"] == accession), symbol)
        digest = hashlib.sha256(f"{accession}|{start}|{end}|{sequence}".encode()).hexdigest()[:12]
        human_host = any("Homo sapiens" in str(record.get("host_organism_name", "")) for record in records)
        alleles = sorted({str(record.get("mhc_allele_name", "")) for record in records if record.get("mhc_allele_name")})
        pmids = sorted({str(record.get("pubmed_id", "")) for record in records if record.get("pubmed_id")})
        assays = sorted({str(record.get("assay_names", "")) for record in records if record.get("assay_names")})
        output.append({
            "candidate_id": f"{'EBV' if kingdom == 'EBV' else 'SELF'}_IEDB_{digest}",
            "library_version": LIBRARY_VERSION,
            "kingdom": kingdom,
            "species": "Human herpesvirus 4" if kingdom == "EBV" else "Homo sapiens",
            "reference_strain_or_isoform": "B95-8 reference mapping; study variant retained" if kingdom == "EBV" else "canonical reviewed UniProt entry",
            "protein_symbol": symbol,
            "protein_name": str(first.get("parent_source_antigen_name", "")),
            "source_accession": accession,
            "accession": accession,
            "source_start_1_based": start if start is not None else "",
            "source_end_1_based": end if end is not None else "",
            "start": start if start is not None else 10**9,
            "end": end if end is not None else 10**9,
            "sequence": sequence,
            "sequence_length": len(sequence),
            "modification": "unmodified",
            "evidence_priority": 1 if human_host else 3,
            "native_hla_evidence": any("HLA-DR" in allele for allele in alleles),
            "source_certainty": "exact_iedb_positive",
            "evidence_scope": "positive human T-cell MHC-II record" if human_host else "positive animal T-cell MHC-II record",
            "host_species": "Homo sapiens" if human_host else "nonhuman_or_mixed",
            "hla_restrictions": ";".join(alleles),
            "assay_types": ";".join(assays),
            "pmids": ";".join(pmids),
            "supporting_record_count": len(records),
            "required_for_confirmed_system": False,
            "natural_flanks_verified": len(sequence) >= 11,
            "proposed_core": sequence[(len(sequence) - 9) // 2:(len(sequence) - 9) // 2 + 9],
            "selection_reason": "direct positive IEDB T-cell/MHC-II record; not proof of cross-reactivity",
        })
    return output


def canonical_tiles(
    *, protein_symbol: str, accession: str, sequence: str,
    region_start: int = 1, region_end: Optional[int] = None, count: int = 4,
) -> list[dict[str, Any]]:
    sequence = sequence.strip().upper()
    region_end = region_end or len(sequence)
    if region_start < 1 or region_end > len(sequence) or region_end - region_start + 1 < 15:
        raise ValueError(f"invalid tiling region for {protein_symbol}")
    min_start = region_start
    max_start = region_end - 14
    starts = [round(min_start + i * (max_start - min_start) / max(1, count - 1)) for i in range(count)]
    starts = sorted(dict.fromkeys(starts))
    if len(starts) != count:
        raise ValueError(f"could not create {count} unique tiles for {protein_symbol}")
    return [{
        "candidate_id": f"SELF_CANON_{protein_symbol}_{start:04d}_{start + 14:04d}",
        "library_version": LIBRARY_VERSION,
        "kingdom": "human_self",
        "species": "Homo sapiens",
        "reference_strain_or_isoform": "canonical reviewed UniProt entry",
        "protein_symbol": protein_symbol,
        "protein_name": SELF_PROTEINS.get(protein_symbol, {}).get("name", protein_symbol),
        "source_accession": accession,
        "accession": accession,
        "source_start_1_based": start,
        "source_end_1_based": start + 14,
        "start": start,
        "end": start + 14,
        "sequence": sequence[start - 1:start + 14],
        "sequence_length": 15,
        "modification": "unmodified",
        "evidence_priority": 2 if protein_symbol == "ANO2" else 4,
        "native_hla_evidence": False,
        "source_certainty": "region_mapped_primary_source" if protein_symbol == "ANO2" else "canonical_tiling_for_protein_level_evidence",
        "evidence_scope": "region-mapped human dual-reactivity coverage" if protein_symbol == "ANO2" else "exploratory canonical sentinel for protein-level evidence",
        "host_species": "Homo sapiens",
        "hla_restrictions": "unresolved",
        "assay_types": "not an exact positive epitope; deterministic coverage tile",
        "pmids": "41534529" if protein_symbol == "ANO2" else "38843296",
        "supporting_record_count": 0,
        "required_for_confirmed_system": protein_symbol == "ANO2",
        "natural_flanks_verified": True,
        "proposed_core": sequence[start + 2:start + 11],
        "selection_reason": "canonical sequence tile; never classified as an experimentally positive peptide",
    } for start in starts]


def fetch_json(url: str, params: dict[str, str]) -> tuple[list[dict[str, Any]], str]:
    full_url = url + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(full_url, headers={"User-Agent": "EBV-MS-library/2.0"})
    with urllib.request.urlopen(request, timeout=180) as response:
        text = response.read().decode("utf-8")
    value = json.loads(text)
    if not isinstance(value, list):
        raise ValueError(f"expected list response from {url}")
    return value, full_url


def fetch_uniprot_sequence(accession: str) -> str:
    request = urllib.request.Request(f"{UNIPROT_ENDPOINT}/{accession}.fasta", headers={"User-Agent": "EBV-MS-library/2.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        lines = response.read().decode("utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if not sequence or re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
        raise ValueError(f"invalid UniProt FASTA for {accession}")
    return sequence


def fetch_binding_predictions(allele: str, rows: Sequence[dict[str, Any]]) -> tuple[list[dict[str, Any]], str]:
    fasta = "\n".join(f">{row['candidate_id']}\n{row['sequence']}" for row in rows)
    body = urllib.parse.urlencode({
        "method": "recommended_binding", "sequence_text": fasta,
        "allele": allele, "length": "asis",
    }).encode("utf-8")
    request = urllib.request.Request(
        IEDB_PREDICTION_ENDPOINT, data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "EBV-MS-library/2.0"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=240) as response:
        raw = response.read().decode("utf-8")
    response_rows = list(csv.DictReader(io.StringIO(raw), delimiter="\t"))
    required = {"allele", "seq_num", "core_peptide", "peptide", "ic50", "rank"}
    if not response_rows or not required.issubset(response_rows[0]):
        raise ValueError(f"IEDB prediction response incomplete for {allele}: {raw[:300]}")
    by_seq = {int(row["seq_num"]): row for row in response_rows}
    if len(by_seq) != len(rows):
        raise ValueError(f"IEDB returned {len(by_seq)} rows for {len(rows)} sequences on {allele}")
    output = []
    for seq_num, source in enumerate(rows, start=1):
        response = by_seq[seq_num]
        if response["peptide"] != source["sequence"] or response["allele"] != allele:
            raise ValueError(f"IEDB sequence/allele mismatch at {allele} seq_num {seq_num}")
        core = response["core_peptide"]
        starts = [i + 1 for i in range(len(source["sequence"]) - 8) if source["sequence"][i:i + 9] == core]
        resolution = "resolved_unique_fully_contained" if len(starts) == 1 else "unresolved_tied_or_repeated" if starts else "unresolved_flank_dependent"
        output.append({
            "panel_version": source.get("panel_version", PANEL_VERSION),
            "allele": allele,
            "candidate_id": source["candidate_id"],
            "seq_num": seq_num,
            "sequence": source["sequence"],
            "prediction_method": "IEDB recommended_binding",
            "prediction_status": "predicted",
            "raw_response_file": "",
            "percentile_rank": response["rank"],
            "predicted_ic50_nM": response["ic50"],
            "predicted_core": core,
            "core_start": ";".join(map(str, starts)),
            "binding_rank_bin": binding_rank_bin(float(response["rank"])),
            "register_resolution": resolution,
            "eligible_for_p1_p9_geometry": resolution == "resolved_unique_fully_contained",
            "interpretation": "binding/register prediction only; not experimental presentation evidence",
        })
    return output, raw


def literature_registry() -> list[dict[str, Any]]:
    rows = [
        {
            "biological_system_id": "SYS_BALF5_MBP_HY2E11", "evidence_tier": "E1_exact_pmhc_positive",
            "receptor_modality": "human_T_cell_shared_clone", "tcell_positive_denominator": True,
            "viral_protein": "BALF5", "viral_accession": "P03198", "viral_sequence": "TGGVYHFVKKHVHES",
            "viral_coordinates_1_based": "627-641", "viral_hla": "HLA-DRB5*01:01",
            "self_protein": "MBP", "self_accession": "P02686", "self_sequence": "ENPVVHFFKNIVTPR",
            "self_coordinates_1_based": "217-231 canonical P02686; conventional MBP85-99",
            "self_hla": "HLA-DRB1*15:01", "species": "human TCR; EBV and human antigens",
            "strain_or_isoform": "EBV study sequence; human MBP conventional numbering",
            "assay_type": "same clone/TCR recognition plus pMHC crystal structures", "receptor_or_clone_id": "Hy.2E11",
            "ptms": "none declared", "primary_source": "PMID:12244309", "doi": "10.1038/ni835",
            "unresolved_fields": "complete alpha/beta Hy.2E11 sequences not publicly recovered",
            "evaluability_status": "strict_recovery_evaluable_after_models",
        },
        {
            "biological_system_id": "SYS_EBNA1_ANO2_2026", "evidence_tier": "E2_human_tcell_protein_or_region",
            "receptor_modality": "human_T_cell_clones", "tcell_positive_denominator": False,
            "viral_protein": "EBNA1", "viral_accession": "P03211", "viral_sequence": "unresolved",
            "viral_coordinates_1_based": "unresolved", "viral_hla": "incomplete peptide/register resolution",
            "self_protein": "ANO2", "self_accession": "Q9NQ90", "self_sequence": "region 79-168",
            "self_coordinates_1_based": "79-168", "self_hla": "incomplete peptide/register resolution",
            "species": "human", "strain_or_isoform": "study-specific constructs retained in source",
            "assay_type": "dual-reactive human T-cell clones", "receptor_or_clone_id": "study clone set",
            "ptms": "none declared", "primary_source": "PMID:41534529", "doi": "10.1016/j.cell.2025.12.032",
            "unresolved_fields": "exact viral/self peptide pair, register, and native-HLA pair",
            "evaluability_status": "partially_mapped",
        },
        {
            "biological_system_id": "SYS_EBNA1_MOG_2024", "evidence_tier": "E2_human_tcell_protein_or_region",
            "receptor_modality": "human_T_cell_clones", "tcell_positive_denominator": False,
            "viral_protein": "EBNA1", "viral_accession": "P03211", "viral_sequence": "whole-protein MVA",
            "viral_coordinates_1_based": "unresolved", "viral_hla": "unresolved",
            "self_protein": "MOG", "self_accession": "Q16653", "self_sequence": "whole-protein MVA",
            "self_coordinates_1_based": "unresolved", "self_hla": "unresolved", "species": "human",
            "strain_or_isoform": "study MVA constructs", "assay_type": "human dual-reactive T-cell clones",
            "receptor_or_clone_id": "study clone set", "ptms": "none declared",
            "primary_source": "PMID:38843296", "doi": "10.1371/journal.ppat.1012177",
            "unresolved_fields": "exact peptide pair, registers, and HLA", "evaluability_status": "covered_but_not_pmhc_evaluable",
        },
        {
            "biological_system_id": "SYS_EBNA1_MYELIN_POOL_2008", "evidence_tier": "E3_supportive_tcell",
            "receptor_modality": "human_T_cell_clones_pooled_targets", "tcell_positive_denominator": False,
            "viral_protein": "EBNA1", "viral_accession": "P03211", "viral_sequence": "overlapping peptide pool",
            "viral_coordinates_1_based": "multiple", "viral_hla": "incompletely resolved",
            "self_protein": "MBP;PLP1;MOG;CNP", "self_accession": "P02686;P60201;Q16653;P09543",
            "self_sequence": "pooled peptides", "self_coordinates_1_based": "multiple", "self_hla": "incompletely resolved",
            "species": "human", "strain_or_isoform": "source constructs", "assay_type": "clone recognition of peptide pools",
            "receptor_or_clone_id": "study clone set", "ptms": "none declared",
            "primary_source": "PMID:18663124", "doi": "10.1084/jem.20072397",
            "unresolved_fields": "one exact dual-reactive peptide pair and native pMHC registers",
            "evaluability_status": "partially_mapped",
        },
        {
            "biological_system_id": "SYS_EBNA1_CRYAB_2023", "evidence_tier": "E3_supportive_tcell",
            "receptor_modality": "mixed_antibody_and_supportive_T_cell", "tcell_positive_denominator": False,
            "viral_protein": "EBNA1", "viral_accession": "P03211", "viral_sequence": "study regions",
            "viral_coordinates_1_based": "see primary source", "viral_hla": "not an exact shared-human-TCR pMHC pair",
            "self_protein": "CRYAB", "self_accession": "P02511", "self_sequence": "study regions",
            "self_coordinates_1_based": "see primary source", "self_hla": "not an exact shared-human-TCR pMHC pair",
            "species": "human and animal evidence", "strain_or_isoform": "source constructs",
            "assay_type": "antibody cross-reactivity plus supportive reciprocal T-cell evidence",
            "receptor_or_clone_id": "not one exact shared human clone", "ptms": "none declared",
            "primary_source": "PMID:37196088", "doi": "10.1126/sciadv.adg3032",
            "unresolved_fields": "strict exact human TCR/pMHC system", "evaluability_status": "partially_mapped",
        },
        {
            "biological_system_id": "SYS_GB_GH_DRB1501_2024", "evidence_tier": "context_only",
            "receptor_modality": "human_T_cell_presentation_context", "tcell_positive_denominator": False,
            "viral_protein": "BALF4/gB;BXLF2/gH", "viral_accession": "P03188;P03231",
            "viral_sequence": "exact DRB1*15:01-positive epitopes in source", "viral_coordinates_1_based": "see source",
            "viral_hla": "HLA-DRB1*15:01", "self_protein": "none", "self_accession": "none",
            "self_sequence": "none", "self_coordinates_1_based": "none", "self_hla": "none",
            "species": "human", "strain_or_isoform": "source EBV antigens", "assay_type": "T-cell antigen presentation",
            "receptor_or_clone_id": "not a viral-self dual-reactive clone", "ptms": "none declared",
            "primary_source": "PMID:39432795", "doi": "10.1073/pnas.2416097121",
            "unresolved_fields": "no demonstrated self cross-reactive partner", "evaluability_status": "context_only_not_a_pair",
        },
        {
            "biological_system_id": "SYS_EBNA1_GLIALCAM_ANTIBODY", "evidence_tier": "context_only",
            "receptor_modality": "antibody_only", "tcell_positive_denominator": False,
            "viral_protein": "EBNA1", "viral_accession": "P03211", "viral_sequence": "386-405 region",
            "viral_coordinates_1_based": "386-405", "viral_hla": "not applicable to antibody evidence",
            "self_protein": "GlialCAM", "self_accession": "Q14CZ8", "self_sequence": "370-389 region",
            "self_coordinates_1_based": "370-389", "self_hla": "not applicable to antibody evidence",
            "species": "human antibodies", "strain_or_isoform": "source constructs", "assay_type": "B-cell/antibody cross-reactivity",
            "receptor_or_clone_id": "patient-derived B-cell clones", "ptms": "none declared",
            "primary_source": "PMID:35073561", "doi": "10.1038/s41586-022-04432-7",
            "unresolved_fields": "direct T-cell cross-reactivity absent", "evaluability_status": "excluded_antibody_only",
        },
    ]
    validate_registry(rows)
    return rows


def _iedb_query_params(**extra: str) -> dict[str, str]:
    return {"limit": "10000", "qualitative_measure": "like.Positive*", "mhc_class": "eq.II", **extra}


def load_hla_sequences() -> dict[str, dict[str, str]]:
    package = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Downloads/alphafold_multiallele_5x30_2026-08-20/hla_sequence_manifest.csv"
    rows = read_csv(package)
    result = {"HLA-DRB1*15:01": {"dra_sequence": DRA_SEQUENCE, "drb_sequence": DRB1_1501_SEQUENCE}}
    code_map = {"DRB1*13:03": "HLA-DRB1*13:03", "DRB1*03:01": "HLA-DRB1*03:01", "DRB1*08:01": "HLA-DRB1*08:01"}
    for row in rows:
        if row["chain"] != "HLA-DRB":
            continue
        allele = next((full for short, full in code_map.items() if row["allele_or_name"].startswith(short)), None)
        if allele:
            result[allele] = {"dra_sequence": DRA_SEQUENCE, "drb_sequence": row["sequence"]}
    if set(result) != set(ALLELES):
        raise ValueError(f"missing HLA chains: {set(ALLELES) - set(result)}")
    return result


def _ensure_required_candidates(master: list[dict[str, Any]]) -> None:
    for row in master:
        if row["kingdom"] == "EBV" and row["protein_symbol"] == "BALF5" and row["sequence"] == "TGGVYHFVKKHVHES":
            row["required_for_confirmed_system"] = True
            row["evidence_priority"] = 0
            row["selection_reason"] = "E1 exact BALF5 arm of Hy.2E11 system"
        if row["kingdom"] == "human_self" and row["protein_symbol"] == "MBP" and row["sequence"] == "ENPVVHFFKNIVTPR":
            row["required_for_confirmed_system"] = True
            row["evidence_priority"] = 0
            row["selection_reason"] = "E1 exact MBP arm of Hy.2E11 system"
        if row["kingdom"] == "human_self" and row["protein_symbol"] == "MOG" and row["sequence"] == "MEVGWYRPPFSRVVHLYRNGK":
            row["required_for_confirmed_system"] = True
            row["evidence_priority"] = min(int(row["evidence_priority"]), 1)
            row["selection_reason"] = "human Q16653 residues 64-84; conventional human MOG35-55 region"


def build_master_library() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, str], dict[str, str]]:
    search_log = []
    raw_urls = {}
    ebv_raw, ebv_url = fetch_json(IEDB_QUERY_ENDPOINT, _iedb_query_params(
        host_organism_iri="eq.NCBITaxon:9606",
        source_organism_iri_search="cs.{NCBITaxon:10376}",
    ))
    raw_urls["iedb_ebv_human_positive_mhcii"] = ebv_url
    search_log.append({"search_id": "IEDB_EBV_001", "date": SEARCH_FREEZE_DATE, "source": "IEDB Query API", "query_url": ebv_url, "raw_record_count": len(ebv_raw), "purpose": "positive human-host EBV T-cell/MHC-II epitopes"})
    master = normalize_iedb_tcell_rows(ebv_raw, "EBV")

    uniprot_sequences = {}
    for symbol, metadata in SELF_PROTEINS.items():
        accession = metadata["accession"]
        sequence = fetch_uniprot_sequence(accession)
        uniprot_sequences[accession] = sequence
        params = _iedb_query_params(parent_source_antigen_iri_search=f"cs.{{UNIPROT:{accession}}}")
        source_rows, source_url = fetch_json(IEDB_QUERY_ENDPOINT, params)
        raw_urls[f"iedb_self_{symbol}"] = source_url
        search_log.append({"search_id": f"IEDB_SELF_{symbol}", "date": SEARCH_FREEZE_DATE, "source": "IEDB Query API", "query_url": source_url, "raw_record_count": len(source_rows), "purpose": f"positive MHC-II T-cell records for {symbol}"})
        master.extend(normalize_iedb_tcell_rows(source_rows, "human_self"))
        region = (79, 168) if symbol == "ANO2" else (1, len(sequence))
        master.extend(canonical_tiles(protein_symbol=symbol, accession=accession, sequence=sequence, region_start=region[0], region_end=region[1], count=4))

    # Ensure exact calibration/MOG convention records exist even if a live API
    # query represents them under a historical accession.
    explicit = [
        ("EBV", "BALF5", "P03198", 627, "TGGVYHFVKKHVHES", "PMID:12244309"),
        ("human_self", "MBP", "P02686", 217, "ENPVVHFFKNIVTPR", "PMID:12244309"),
        ("human_self", "MOG", "Q16653", 64, "MEVGWYRPPFSRVVHLYRNGK", "MOG35-55 convention"),
    ]
    existing = {(row["kingdom"], row["accession"], row["sequence"]) for row in master}
    for kingdom, symbol, accession, start, sequence, source in explicit:
        if (kingdom, accession, sequence) in existing:
            continue
        master.append({
            "candidate_id": f"{'EBV' if kingdom == 'EBV' else 'SELF'}_LIT_{symbol}_{start:04d}_{start + len(sequence) - 1:04d}",
            "library_version": LIBRARY_VERSION, "kingdom": kingdom,
            "species": "Human herpesvirus 4" if kingdom == "EBV" else "Homo sapiens",
            "reference_strain_or_isoform": "study sequence" if kingdom == "EBV" else "canonical reviewed UniProt entry",
            "protein_symbol": symbol, "protein_name": symbol, "source_accession": accession, "accession": accession,
            "source_start_1_based": start, "source_end_1_based": start + len(sequence) - 1,
            "start": start, "end": start + len(sequence) - 1, "sequence": sequence, "sequence_length": len(sequence),
            "modification": "unmodified", "evidence_priority": 0, "native_hla_evidence": True,
            "source_certainty": "exact_primary_source", "evidence_scope": "exact literature calibration or named convention",
            "host_species": "Homo sapiens", "hla_restrictions": "native HLA in registry" if symbol != "MOG" else "not implied",
            "assay_types": source, "pmids": source.replace("PMID:", ""), "supporting_record_count": 1,
            "required_for_confirmed_system": symbol in {"BALF5", "MBP", "MOG"}, "natural_flanks_verified": True,
            "proposed_core": sequence[(len(sequence) - 9) // 2:(len(sequence) - 9) // 2 + 9],
            "selection_reason": "exact primary-source calibration target" if symbol != "MOG" else "human Q16653 64-84; conventional human MOG35-55",
        })
    _ensure_required_candidates(master)

    # Collapse only exact duplicate biological epitope keys, not alternative lengths.
    deduped = {}
    for row in master:
        key = (row["kingdom"], row["species"], row["accession"], row["start"], row["end"], row["sequence"], row["modification"])
        old = deduped.get(key)
        if old is None or (int(row["evidence_priority"]), row["candidate_id"]) < (int(old["evidence_priority"]), old["candidate_id"]):
            deduped[key] = row
    return sorted(deduped.values(), key=lambda row: (row["kingdom"], row["protein_symbol"], int(row["start"]), row["sequence"])), search_log, raw_urls, uniprot_sequences


def current_panel_audit() -> list[dict[str, Any]]:
    source = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Downloads/alphafold_multiallele_5x30_2026-08-20/peptide_panel_manifest.csv"
    rows = read_csv(source)
    counts = Counter((row.get("arm", ""), row.get("source_antigen_name", "")) for row in rows)
    return [{"panel": "V1_fixed_50", "arm": arm, "protein_or_source": protein, "peptide_count": count, "audit_status": "preserved_unchanged"} for (arm, protein), count in sorted(counts.items())]


def _raw_filename(allele: str, prefix: str) -> str:
    return f"{prefix}_{re.sub(r'[^a-z0-9]+', '_', allele.lower()).strip('_')}.tsv"


def run_predictions(panel: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    predictions = []
    for allele in ALLELES:
        rows, raw = fetch_binding_predictions(allele, panel)
        raw_name = _raw_filename(allele, "iedb_panel80")
        raw_path = RAW_OUT / raw_name
        raw_path.write_text(raw, encoding="utf-8")
        for row in rows:
            row["raw_response_file"] = str(raw_path.relative_to(OUT))
        predictions.extend(rows)
    if len(predictions) != 320:
        raise AssertionError("V2 register table must contain exactly 320 records")
    return predictions


def _prediction_for(predictions: Sequence[dict[str, Any]], allele: str, candidate_id: str) -> dict[str, Any]:
    matches = [row for row in predictions if row["allele"] == allele and row["candidate_id"] == candidate_id]
    if len(matches) != 1:
        raise ValueError(f"expected one prediction for {allele} {candidate_id}")
    return matches[0]


def build_calibration(
    master: Sequence[dict[str, Any]], panel: Sequence[dict[str, Any]], panel_predictions: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    balf5 = next(row for row in panel if row["protein_symbol"] == "BALF5" and row["sequence"] == "TGGVYHFVKKHVHES")
    mbp = next(row for row in panel if row["protein_symbol"] == "MBP" and row["sequence"] == "ENPVVHFFKNIVTPR")

    ebv_pool_sources = [balf5] + [row for row in master if row["kingdom"] == "EBV" and row["candidate_id"] != balf5["candidate_id"] and abs(len(row["sequence"]) - len(balf5["sequence"])) <= 1]
    # Stable deduplication keeps a single candidate per ID.
    ebv_pool_sources = list({row["candidate_id"]: row for row in ebv_pool_sources}.values())
    ebv_predictions, ebv_raw = fetch_binding_predictions("HLA-DRB5*01:01", ebv_pool_sources)
    ebv_raw_path = RAW_OUT / "iedb_calibration_ebv_drb5_0101.tsv"
    ebv_raw_path.write_text(ebv_raw, encoding="utf-8")
    ebv_by_id = {row["candidate_id"]: row for row in ebv_predictions}
    balf5_prediction = ebv_by_id[balf5["candidate_id"]]
    ebv_control_pool = [{
        **source,
        "binding_rank_bin": ebv_by_id[source["candidate_id"]]["binding_rank_bin"],
        "excluded_source": source["candidate_id"] == balf5["candidate_id"],
    } for source in ebv_pool_sources if source["candidate_id"] in ebv_by_id]
    ebv_controls = select_native_controls(
        {"candidate_id": balf5["candidate_id"], "sequence": balf5["sequence"], "binding_rank_bin": balf5_prediction["binding_rank_bin"]},
        ebv_control_pool, count=5,
    )

    universe = read_csv(ROOT / "processed/structural_control_expansion_2026-08-15/frozen_control_universe.csv")
    human_sources = []
    for row in universe:
        sequence = row["peptide"]
        if abs(len(sequence) - len(mbp["sequence"])) > 1 or not is_allowed_non_cns_control_source(row["source_accession"], row["source_antigen_name"]):
            continue
        human_sources.append({
            "candidate_id": row["candidate_id"], "sequence": sequence,
            "protein_symbol": row["source_antigen_name"], "source_accession": row["source_accession"],
            "excluded_source": False,
        })
    human_sources = list({row["candidate_id"]: row for row in human_sources}.values())
    human_predictions, human_raw = fetch_binding_predictions("HLA-DRB1*15:01", human_sources)
    human_raw_path = RAW_OUT / "iedb_calibration_human_non_cns_drb1_1501.tsv"
    human_raw_path.write_text(human_raw, encoding="utf-8")
    human_by_id = {row["candidate_id"]: row for row in human_predictions}
    mbp_prediction = _prediction_for(panel_predictions, "HLA-DRB1*15:01", mbp["candidate_id"])
    for row in human_sources:
        row["binding_rank_bin"] = human_by_id[row["candidate_id"]]["binding_rank_bin"]
    human_controls = select_native_controls(
        {"candidate_id": mbp["candidate_id"], "sequence": mbp["sequence"], "binding_rank_bin": mbp_prediction["binding_rank_bin"]},
        human_sources, count=5,
    )

    entities = [{
        "entity_id": balf5["candidate_id"], "entity_role": "E1_positive_BALF5", "arm": "viral",
        "allele": "HLA-DRB5*01:01", "sequence": balf5["sequence"], "dra_sequence": DRA_SEQUENCE,
        "drb_sequence": DRB5_0101_SEQUENCE,
    }] + [{
        "entity_id": row["candidate_id"], "entity_role": "score_blind_matched_EBV_control", "arm": "viral",
        "allele": "HLA-DRB5*01:01", "sequence": row["sequence"], "dra_sequence": DRA_SEQUENCE,
        "drb_sequence": DRB5_0101_SEQUENCE,
    } for row in ebv_controls] + [{
        "entity_id": mbp["candidate_id"], "entity_role": "E1_positive_MBP", "arm": "self",
        "allele": "HLA-DRB1*15:01", "sequence": mbp["sequence"], "dra_sequence": DRA_SEQUENCE,
        "drb_sequence": DRB1_1501_SEQUENCE,
    }] + [{
        "entity_id": row["candidate_id"], "entity_role": "score_blind_matched_non_CNS_human_control", "arm": "self",
        "allele": "HLA-DRB1*15:01", "sequence": row["sequence"], "dra_sequence": DRA_SEQUENCE,
        "drb_sequence": DRB1_1501_SEQUENCE,
    } for row in human_controls]
    jobs, manifest = build_native_calibration_jobs(entities)
    controls = []
    for arm, target, target_pred, selected in (
        ("viral", balf5, balf5_prediction, ebv_controls),
        ("self", mbp, mbp_prediction, human_controls),
    ):
        for rank, row in enumerate(selected, start=1):
            controls.append({
                "arm": arm, "target_candidate_id": target["candidate_id"], "target_sequence": target["sequence"],
                "target_binding_rank_bin": target_pred["binding_rank_bin"], "control_rank": rank,
                "control_candidate_id": row["candidate_id"], "control_sequence": row["sequence"],
                "control_source": row.get("protein_symbol", ""), "control_accession": row.get("source_accession", row.get("accession", "")),
                "control_binding_rank_bin": row["binding_rank_bin"], "composition_distance": row["composition_distance"],
                "length_difference": row["length_difference"], "selection_key": row["selection_key"],
                "score_blind_freeze": True,
            })
    return jobs, manifest, controls, ebv_predictions + human_predictions


def write_coverage_svg(panel: Sequence[dict[str, Any]], path: Path) -> None:
    counts = Counter((row["kingdom"], row["protein_symbol"]) for row in panel)
    labels = [("EBV", p, c) for (_, p), c in sorted(counts.items()) if _ == "EBV"] + [("Self", p, c) for (k, p), c in sorted(counts.items()) if k == "human_self"]
    width, row_h = 900, 24
    height = 90 + row_h * len(labels)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', '<style>text{font-family:Arial,sans-serif;fill:#17233b}.title{font-size:22px;font-weight:700}.label{font-size:13px}.count{font-size:12px;fill:#4a5568}</style>', '<text x="30" y="34" class="title">Frozen V2 protein coverage (40 EBV + 40 self)</text>']
    for i, (arm, protein, count) in enumerate(labels):
        y = 68 + i * row_h
        color = "#5B5BD6" if arm == "EBV" else "#0F9D8A"
        parts.append(f'<text x="30" y="{y + 13}" class="label">{arm}: {protein}</text>')
        parts.append(f'<rect x="270" y="{y}" width="{count * 110}" height="16" rx="3" fill="{color}"/>')
        parts.append(f'<text x="{282 + count * 110}" y="{y + 13}" class="count">{count}</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    AF3_OUT.mkdir(parents=True)
    RAW_OUT.mkdir(parents=True)

    registry = literature_registry()
    master, search_log, raw_urls, uniprot_sequences = build_master_library()
    panel = freeze_panel(master)
    predictions = run_predictions(panel)
    hla_sequences = load_hla_sequences()
    jobs, model_inventory, batches = build_af3_jobs(panel, hla_sequences)
    pairs = build_discovery_pairs(panel)
    calibration_jobs, calibration_manifest, controls, calibration_predictions = build_calibration(master, panel, predictions)

    write_csv(OUT / "current_panel_protein_coverage_audit.csv", current_panel_audit())
    write_csv(OUT / "literature_search_log.csv", search_log + [
        {"search_id": "PRIMARY_LANG_2002", "date": SEARCH_FREEZE_DATE, "source": "PubMed/primary study", "query_url": "https://pubmed.ncbi.nlm.nih.gov/12244309/", "raw_record_count": 1, "purpose": "E1 exact pMHC/shared-clone calibration"},
        {"search_id": "PRIMARY_ANO2_2026", "date": SEARCH_FREEZE_DATE, "source": "PubMed/primary study", "query_url": "https://pubmed.ncbi.nlm.nih.gov/41534529/", "raw_record_count": 1, "purpose": "E2 EBNA1-ANO2 human T-cell system"},
        {"search_id": "PRIMARY_MOG_2024", "date": SEARCH_FREEZE_DATE, "source": "PubMed/primary study", "query_url": "https://pubmed.ncbi.nlm.nih.gov/38843296/", "raw_record_count": 1, "purpose": "E2 EBNA1-MOG human T-cell system"},
        {"search_id": "PRIMARY_LUNEMANN_2008", "date": SEARCH_FREEZE_DATE, "source": "PubMed/primary study", "query_url": "https://pubmed.ncbi.nlm.nih.gov/18663124/", "raw_record_count": 1, "purpose": "E3 pooled myelin T-cell system"},
    ])
    write_csv(OUT / "literature_tcell_pair_registry.csv", registry)
    write_csv(OUT / "master_protein_epitope_library.csv", master)
    write_csv(OUT / "frozen_v2_80_peptide_panel.csv", panel)
    write_csv(OUT / "allele_register_predictions_320.csv", predictions)
    write_csv(OUT / "model_inventory_320.csv", model_inventory)
    sample_qc = [{**row, "sample_index": sample, "qc_status": "pending_external_alphafold", "ranking_score": "", "iptm": "", "ptm": "", "clashes": "", "peptide_plddt": "", "contacts": "", "within_job_pose_rmsd": ""} for row in model_inventory for sample in range(5)]
    write_csv(OUT / "model_sample_qc_1600.csv", sample_qc)
    write_csv(OUT / "within_allele_pair_universe_6400.csv", pairs)
    write_json(AF3_OUT / "all_320_jobs.json", jobs)
    for index, batch in enumerate(batches, start=1):
        write_json(AF3_OUT / f"ebvms_v2_batch_{index:02d}_{len(batch)}_jobs.json", batch)
    write_json(AF3_OUT / "native_hla_calibration_24_jobs.json", calibration_jobs)
    write_csv(OUT / "native_hla_calibration_manifest_24.csv", calibration_manifest)
    write_csv(OUT / "calibration_comparison_universe_72.csv", build_calibration_comparison_universe(calibration_manifest))
    write_csv(OUT / "frozen_native_hla_controls.csv", controls)
    write_csv(OUT / "calibration_control_binding_predictions.csv", calibration_predictions)
    write_csv(OUT / "unresolved_evidence_report.csv", [{
        "biological_system_id": row["biological_system_id"], "evidence_tier": row["evidence_tier"],
        "status": row["evaluability_status"], "unresolved_fields": row["unresolved_fields"],
    } for row in registry if row["evidence_tier"] != "E1_exact_pmhc_positive"])
    write_csv(OUT / "positive_recovery_report.csv", [{
        "biological_system_id": "SYS_BALF5_MBP_HY2E11", "recovery_status": "not_evaluable_models_not_yet_run",
        "required_rule": "top 3 of 26 on both seeds and below equal-weight control median",
        "full_decoy_combinations": 25, "single_arm_sensitivity_combinations": 10,
        "controls_frozen": True, "claim_boundary": CLAIM_BOUNDARY,
    }])
    write_csv(OUT / "reference_sequence_manifest.csv", [{
        "entity": accession, "source": f"UniProt {accession}", "sequence_length": len(sequence),
        "sha256": hashlib.sha256(sequence.encode()).hexdigest(), "sequence": sequence,
    } for accession, sequence in sorted(uniprot_sequences.items())] + [{
        "entity": allele, "source": "existing validated project HLA input", "sequence_length": len(value["drb_sequence"]),
        "sha256": hashlib.sha256(value["drb_sequence"].encode()).hexdigest(), "sequence": value["drb_sequence"],
    } for allele, value in sorted(hla_sequences.items())] + [{
        "entity": "HLA-DRB5*01:01", "source": "RCSB PDB 1H15 chain B/E mature extracellular construct", "sequence_length": len(DRB5_0101_SEQUENCE),
        "sha256": hashlib.sha256(DRB5_0101_SEQUENCE.encode()).hexdigest(), "sequence": DRB5_0101_SEQUENCE,
    }])
    write_json(OUT / "source_urls.json", raw_urls)
    write_coverage_svg(panel, OUT / "protein_coverage_figure.svg")

    summary = {
        "library_version": LIBRARY_VERSION, "panel_version": PANEL_VERSION, "search_frozen_through": SEARCH_FREEZE_DATE,
        "registry_systems": len(registry), "strict_e1_systems": 1, "master_epitopes": len(master),
        "panel_ebv": sum(row["kingdom"] == "EBV" for row in panel), "panel_self": sum(row["kingdom"] == "human_self" for row in panel),
        "panel_ebv_proteins": len({row["protein_symbol"] for row in panel if row["kingdom"] == "EBV"}),
        "panel_self_proteins": len({row["protein_symbol"] for row in panel if row["kingdom"] == "human_self"}),
        "prediction_records": len(predictions), "prepared_model_jobs": len(jobs), "pending_sample_qc_rows": len(sample_qc),
        "within_allele_pairs": len(pairs), "calibration_jobs": len(calibration_jobs), "calibration_seeds": list(CALIBRATION_SEEDS),
        "external_state": "AlphaFold jobs prepared but not submitted; recovery remains non-evaluable until downloads exist",
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(OUT / "validation_summary.json", summary)
    (OUT / "README.md").write_text(f"""# EBV-MS T-cell Library V2 ({SEARCH_FREEZE_DATE})

This additive package leaves all V1 analyses unchanged. It freezes a 40-EBV/40-self breadth panel, live IEDB register predictions, 320 AlphaFold Server inputs, and a separate 24-job native-HLA positive-control calibration.

## Current status

- Literature registry: {len(registry)} independent biological systems; one strict E1 system.
- Master library: {len(master)} exact IEDB/literature records plus explicitly labeled canonical coverage tiles.
- Frozen panel: 40 EBV peptides across {summary['panel_ebv_proteins']} proteins and 40 self peptides across {summary['panel_self_proteins']} proteins.
- IEDB: 320 completed prediction records with exact `seq_num` and saved raw responses.
- AlphaFold: 320 discovery jobs and 24 calibration jobs are prepared, **not submitted**.
- Geometry/recovery: pending model downloads; no positive has been declared recovered.

Canonical tiling rows for protein/region-level evidence are not experimentally positive epitopes. Antibody-only EBNA1-GlialCAM is documented but excluded from the T-cell denominator. The EBNA1-MBP structural lead remains a computational hypothesis and is separate from positive recovery.

{CLAIM_BOUNDARY}
""", encoding="utf-8")

    files = sorted(path for path in OUT.rglob("*") if path.is_file() and path.name != "SHA256SUMS.csv")
    write_csv(OUT / "SHA256SUMS.csv", [{"relative_path": str(path.relative_to(OUT)), "sha256": sha256_file(path), "bytes": path.stat().st_size} for path in files])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
