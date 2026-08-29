"""Prepare, fetch, and analyze the additive eight-candidate evidence dossier."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import ssl
import subprocess
import time
from typing import Any, Iterable, Mapping, Optional, Sequence
import urllib.parse
import urllib.request

from high_yield_candidate_evidence import (
    CLAIM_BOUNDARY,
    build_stage2_gate,
    classify_assay_evidence,
    classify_ligand_hit,
    classify_stage1,
    normalize_hla,
    presentation_conditioned_rarity,
    scan_similarity_rarity_fast,
    sequence_relation,
    summarize_conservation,
    summarize_predictor_evidence,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_UPSTREAM = ROOT / "processed/high_yield_control_validation_2026-08-28"
DEFAULT_REGISTER = ROOT / "processed/high_yield_register_resolution_2026-08-28"
DEFAULT_V3 = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
DEFAULT_OUT = ROOT / "processed/high_yield_candidate_evidence_2026-08-28"
HUMAN_SOURCE_FASTA = ROOT / "processed/protein_region_annotations/source_records/human_parent_proteins_uniprot.fasta"
EBV_SOURCE_FASTA = ROOT / "processed/protein_region_annotations/source_records/ebv_parent_proteins.fasta"

HLA_ATLAS_RELEASE = "2020.12"
MIXMHCIIPRED_RELEASE = "v2.1.beta1.2"
MIXMHCIIPRED_ARCHIVE_SHA256 = "fe86e0390c96ca7b4e7b8b68d563c717ce43f4091fe460fc283a33ccd716be74"
HLA_ATLAS_BASE = f"https://hla-ligand-atlas.org/rel/{HLA_ATLAS_RELEASE}"
IEDB_QUERY_BASE = "https://query-api.iedb.org"
IEDB_PREDICT_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
UNIPROT_HUMAN_REVIEWED_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query="
    "%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29"
)
UNIPROT_EBV_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query="
    "%28organism_id%3A10376%29"
)
PXD068488_URL = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD068488"
PXD068488_FTP_BASE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD068488"
GNOMAD_GRAPHQL_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
COMMON_VARIANT_AF_THRESHOLD = 0.01

CLAIM_FIELDS = {
    "specificity_claim_allowed": False,
    "cross_reactivity_claim_allowed": False,
    "molecular_mimicry_claim_allowed": False,
    "discovery_unlock_allowed": False,
    "weights_frozen": False,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] = (),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or sorted({key for row in rows for key in row}))
    if not fieldnames:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_fasta(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    records: list[dict[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append(_fasta_record(header, "".join(sequence)))
            header = line[1:].strip()
            sequence = []
        else:
            sequence.append(line.strip())
    if header:
        records.append(_fasta_record(header, "".join(sequence)))
    return records


def _fasta_record(header: str, sequence: str) -> dict[str, str]:
    token = header.split()[0]
    accession = token.split("|")[1] if token.count("|") >= 2 else token
    return {
        "accession": accession,
        "protein": header,
        "header": header,
        "sequence": re.sub(r"[^A-Za-z]", "", sequence).upper(),
    }


def _source_match(
    peptide: str,
    records: Sequence[Mapping[str, str]],
) -> tuple[str, str, int]:
    matches = []
    for record in records:
        start = str(record["sequence"]).find(peptide)
        if start >= 0:
            matches.append((str(record["accession"]), str(record["sequence"]), start))
    if not matches:
        return "", "", -1
    return sorted(matches, key=lambda item: item[0])[0]


def _mix_context(protein_sequence: str, peptide: str, start: int) -> tuple[str, str]:
    if start < 0 or not protein_sequence:
        return f"XXX{peptide[:3]}{peptide[-3:]}XXX", "unknown_external_flanks"
    end = start + len(peptide)
    upstream = protein_sequence[max(0, start - 3) : start].rjust(3, "-")
    downstream = protein_sequence[end : end + 3].ljust(3, "-")
    return upstream + peptide[:3] + peptide[-3:] + downstream, "exact_parent_context"


def _frozen_targets(upstream_dir: Path) -> list[dict[str, str]]:
    rows = [row for row in read_csv(upstream_dir / "frozen_target_registry.csv") if row.get("lane") == "sequence"]
    rows.sort(key=lambda row: row["target_id"])
    if len(rows) != 8:
        raise ValueError("candidate dossier requires exactly eight sequence-lane targets")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["allele"]] = counts.get(row["allele"], 0) + 1
        for side in ("ebv", "self"):
            sequence = row[f"{side}_sequence"].strip().upper()
            core = row[f"{side}_core"].strip().upper()
            if len(core) != 9 or core not in sequence:
                raise ValueError(f"{row['target_id']} {side} core is not an exact contained nonamer")
    if sorted(counts.values()) != [2, 2, 2, 2]:
        raise ValueError("candidate dossier requires two targets per HLA")
    return rows


def _arm_registry(targets: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    human_records = parse_fasta(HUMAN_SOURCE_FASTA)
    ebv_records = parse_fasta(EBV_SOURCE_FASTA)
    arms = []
    for target in targets:
        for side, kingdom, records in (
            ("ebv", "EBV", ebv_records),
            ("self", "human", human_records),
        ):
            peptide = target[f"{side}_sequence"].strip().upper()
            core = target[f"{side}_core"].strip().upper()
            accession, parent_sequence, start = _source_match(peptide, records)
            context, context_status = _mix_context(parent_sequence, peptide, start)
            arms.append(
                {
                    "arm_id": f"{target['target_id']}__{side}",
                    "target_id": target["target_id"],
                    "pair_id": target["pair_id"],
                    "side": side,
                    "kingdom": kingdom,
                    "allele": target["allele"],
                    "candidate_id": target[f"{side}_candidate_id"],
                    "protein": target[f"{side}_protein"],
                    "sequence": peptide,
                    "core": core,
                    "declared_core_start_1_based": peptide.index(core) + 1,
                    "source_accession": accession,
                    "source_record_status": "exact_parent_sequence_match" if accession else "not_evaluable",
                    "mixmhc_context": context,
                    "mixmhc_context_status": context_status,
                }
            )
    arms.sort(key=lambda row: row["arm_id"])
    if len(arms) != 16:
        raise ValueError("candidate dossier requires exactly sixteen peptide arms")
    return arms


def _upstream_checksums(upstream_dir: Path, register_dir: Path, v3_dir: Path) -> dict[str, str]:
    paths = [
        upstream_dir / "frozen_target_registry.csv",
        upstream_dir / "ranking_context_gate.json",
        register_dir / "register_resolution_gate.json",
        v3_dir / "protocol_lock.json",
        v3_dir / "v3_all_hla_ranked_pairs.csv",
        v3_dir / "combined_drb1501_v3_ranked_pairs.csv",
    ]
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing frozen upstream artifacts: {missing}")
    return {str(path.relative_to(ROOT)): sha256_file(path) for path in paths}


def _source_manifest(output_dir: Path) -> dict[str, Any]:
    sources = [
        ("iedb_assays", "raw_responses/iedb_assay_records.csv", IEDB_QUERY_BASE),
        ("netmhciipan_4_3", "raw_responses/predictor_records.csv", IEDB_PREDICT_URL),
        ("mixmhc2pred_2_1", "raw_responses/predictor_records.csv", "local_version_pinned_binary"),
        ("hla_ligand_atlas_2020_12", "raw_responses/hla_ligand_atlas_hits.csv", HLA_ATLAS_BASE),
        ("pxd068488_published_tables", "raw_responses/pxd068488_published_hits.csv", PXD068488_URL),
        ("human_reviewed_reference_proteome", "raw_responses/human_reviewed_canonical.fasta", UNIPROT_HUMAN_REVIEWED_URL),
        ("ebv_sequence_collection", "raw_responses/ebv_uniprot_sequences.fasta", UNIPROT_EBV_URL),
        ("human_common_variants", "raw_responses/human_common_variant_records.csv", GNOMAD_GRAPHQL_URL),
    ]
    rows = []
    for source_id, relative_path, url in sources:
        path = output_dir / relative_path
        rows.append(
            {
                "source_id": source_id,
                "relative_path": relative_path,
                "source_url": url,
                "status": "cached" if path.exists() else "not_evaluable_not_cached",
                "sha256": sha256_file(path) if path.exists() else "",
            }
        )
    return {"sources": rows, "claim_boundary": CLAIM_BOUNDARY}


def prepare_package(
    *,
    output_dir: Path = DEFAULT_OUT,
    upstream_dir: Path = DEFAULT_UPSTREAM,
    register_dir: Path = DEFAULT_REGISTER,
    v3_dir: Path = DEFAULT_V3,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = _frozen_targets(upstream_dir)
    arms = _arm_registry(targets)
    upstream_checksums = _upstream_checksums(upstream_dir, register_dir, v3_dir)
    write_csv(output_dir / "frozen_candidate_registry.csv", targets)
    write_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv", arms)

    predictor_queries = []
    iedb_queries = []
    for arm in arms:
        predictor_queries.append(
            {
                **arm,
                "netmhcii_sequence_format": "fasta_single_exact_peptide",
                "netmhcii_requested_length": len(arm["sequence"]),
                "mixmhc_context_mode": "natural_parent_context",
            }
        )
        for endpoint in ("mhc_search", "tcell_search"):
            for query_scope, query_sequence in (("exact_peptide", arm["sequence"]), ("declared_core", arm["core"])):
                iedb_queries.append(
                    {
                        "arm_id": arm["arm_id"],
                        "target_id": arm["target_id"],
                        "endpoint": endpoint,
                        "query_scope": query_scope,
                        "query_sequence": query_sequence,
                        "target_sequence": arm["sequence"],
                        "target_hla": arm["allele"],
                    }
                )
    write_csv(output_dir / "prepared_inputs/predictor_queries.csv", predictor_queries)
    write_csv(output_dir / "prepared_inputs/iedb_query_manifest.csv", iedb_queries)

    with (output_dir / "prepared_inputs/mixmhc2pred_context.tsv").open("w", encoding="utf-8") as handle:
        for arm in arms:
            handle.write(f"{arm['sequence']}\t{arm['mixmhc_context']}\n")
    with (output_dir / "prepared_inputs/mixmhc2pred_no_context.tsv").open("w", encoding="utf-8") as handle:
        for arm in arms:
            handle.write(f"{arm['sequence']}\n")

    protocol = {
        "protocol_id": "high_yield_candidate_evidence_2026-08-28",
        "status": "prepared_for_cache_backed_fetch_and_analysis",
        "target_count": len(targets),
        "arm_count": len(arms),
        "alleles": sorted({row["allele"] for row in targets}),
        "frozen_candidate_registry_sha256": sha256_file(output_dir / "frozen_candidate_registry.csv"),
        "upstream_checksums": upstream_checksums,
        "tool_versions": {
            "netmhciipan_el": "4.3",
            "netmhciipan_ba": "4.3",
            "mixmhc2pred": MIXMHCIIPRED_RELEASE,
            "mixmhc2pred_archive_sha256": MIXMHCIIPRED_ARCHIVE_SHA256,
            "hla_ligand_atlas_release": HLA_ATLAS_RELEASE,
            "iedb_query_api": "IQ-API snapshot queried 2026-08-28",
            "gnomad_dataset": GNOMAD_DATASET,
        },
        "database_releases": {
            "hla_ligand_atlas": HLA_ATLAS_RELEASE,
            "human_proteome": "UniProt reviewed reference proteome UP000005640",
            "ebv_sequences": "UniProt organism taxonomy 10376 collection",
            "gnomad": GNOMAD_DATASET,
            "pxd068488": "PRIDE archive 2026-03 release",
        },
        "query_date": "2026-08-28",
        "binding_consensus_rule": "NetMHCIIpan-4.3 EL percentile <=5 and MixMHC2pred-2.1 context percentile <=5",
        "binding_support_rule": "at least one independent percentile <=5 or both independent percentiles <=20",
        "register_consensus_rule": "NetMHCIIpan-4.3 EL and MixMHC2pred-2.1 context select the same P1-P9 core",
        "rarity_rule": "frozen V3 lexicographic sequence order; high-priority threshold <=1 empirical percentile",
        "presentation_conditioned_rarity_rule": "exact-HLA V3 library; 2043-pair combined library for DRB1*15:01 and 1600-pair V3 libraries otherwise; diagnostic only",
        "immunopeptidome_absence_rule": "not_observed_is_missing_evidence_not_a_negative",
        "raw_pxd068488_reprocessing": False,
        "common_variant_rule": "canonical-transcript missense variant with maximum reported exome/genome allele frequency >=0.01 in gnomAD r4",
        "single_composite_score_created": False,
        "discovery_rankings_modified": False,
        **CLAIM_FIELDS,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "protocol_lock.json", protocol)
    write_json(output_dir / "source_manifest.json", _source_manifest(output_dir))
    return {
        "target_count": len(targets),
        "arm_count": len(arms),
        "upstream_checksums": upstream_checksums,
    }


def _ssl_context() -> ssl.SSLContext:
    try:
        import certifi  # type: ignore

        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def _download(url: str, path: Path, timeout: int = 180) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "EBV-MS-evidence-dossier/1.0"})
    with urllib.request.urlopen(request, timeout=timeout, context=_ssl_context()) as response:
        path.write_bytes(response.read())


def _fetch_iedb_assays(output_dir: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    queries = read_csv(output_dir / "prepared_inputs/iedb_query_manifest.csv")
    raw_dir = output_dir / "raw_responses/iedb_query_api"
    raw_dir.mkdir(parents=True, exist_ok=True)
    standardized = []
    failures = []

    def fetch_one(index: int, query: Mapping[str, str]) -> tuple[Mapping[str, str], Path, list[dict[str, Any]], Optional[str]]:
        params = urllib.parse.urlencode(
            {
                "linear_sequence": f"like.*{query['query_sequence']}*",
                "limit": "500",
            }
        )
        url = f"{IEDB_QUERY_BASE}/{query['endpoint']}?{params}"
        raw_path = raw_dir / f"{index:03d}_{query['arm_id']}_{query['endpoint']}_{query['query_scope']}.json"
        try:
            if not raw_path.exists():
                _download(url, raw_path, timeout=45)
            records = json.loads(raw_path.read_text(encoding="utf-8"))
        except Exception as error:
            return query, raw_path, [], f"{type(error).__name__}: {error}"
        return query, raw_path, records, None

    fetched = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(fetch_one, index, query): (index, query)
            for index, query in enumerate(queries, start=1)
        }
        for future in as_completed(futures):
            query, raw_path, records, error = future.result()
            if error:
                failures.append({"query": dict(query), "error": error})
            else:
                fetched.append((dict(query), raw_path, records))

    for query, raw_path, records in sorted(fetched, key=lambda item: str(item[1])):
        for record_index, record in enumerate(records):
            observed = record.get("linear_sequence", "")
            standardized.append(
                {
                    "arm_id": query["arm_id"],
                    "target_id": query["target_id"],
                    "endpoint": query["endpoint"],
                    "query_scope": query["query_scope"],
                    "record_index": record_index,
                    "epitope_sequence": observed,
                    "sequence_relation": sequence_relation(query["target_sequence"], observed),
                    "mhc_allele": record.get("mhc_allele_name", ""),
                    "mhc_class": record.get("mhc_class", ""),
                    "host_organism": record.get("host_organism_name", ""),
                    "qualitative_measure": record.get("qualitative_measure", ""),
                    "assay_name": record.get("assay_description") or record.get("assay_names") or "",
                    "parent_source_antigen": record.get("parent_source_antigen_name", ""),
                    "source_organism": record.get("parent_source_antigen_source_org_name", ""),
                    "pubmed_id": record.get("pubmed_id", ""),
                    "pdb_id": record.get("pdb_id", ""),
                    "evidence_class": classify_assay_evidence(record, query["target_sequence"], query["target_hla"]),
                    "raw_response": str(raw_path.relative_to(output_dir)),
                }
            )
    unique = {}
    for row in standardized:
        key = (
            row["arm_id"], row["endpoint"], row["epitope_sequence"], row["mhc_allele"],
            row["qualitative_measure"], row["pubmed_id"], row["assay_name"],
        )
        unique[key] = row
    rows = sorted(unique.values(), key=lambda row: (row["arm_id"], row["endpoint"], row["epitope_sequence"], str(row["mhc_allele"])))
    fields = [
        "arm_id", "target_id", "endpoint", "query_scope", "record_index", "epitope_sequence",
        "sequence_relation", "mhc_allele", "mhc_class", "host_organism", "qualitative_measure",
        "assay_name", "parent_source_antigen", "source_organism", "pubmed_id", "pdb_id",
        "evidence_class", "raw_response",
    ]
    write_csv(output_dir / "raw_responses/iedb_assay_records.csv", rows, fields)
    write_json(output_dir / "raw_responses/iedb_fetch_status.json", {"query_count": len(queries), "record_count": len(rows), "failures": failures})
    return rows, {"failure_count": len(failures), "record_count": len(rows)}


def _post_form(url: str, payload: Mapping[str, str], timeout: int = 180) -> bytes:
    data = urllib.parse.urlencode(payload).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=data,
        headers={"User-Agent": "EBV-MS-evidence-dossier/1.0", "Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(request, timeout=timeout, context=_ssl_context()) as response:
        return response.read()


def _parse_iedb_predictor_tsv(text: str, predictor: str, arm: Mapping[str, str]) -> dict[str, Any]:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise ValueError(f"empty {predictor} response")
    exact = [row for row in rows if row.get("peptide", "").upper() == arm["sequence"].upper()]
    row = exact[0] if exact else rows[0]
    return {
        "arm_id": arm["arm_id"],
        "target_id": arm["target_id"],
        "allele": arm["allele"],
        "predictor": predictor,
        "percentile_rank": row.get("rank", ""),
        "core": row.get("core_peptide", ""),
        "core_reliability": row.get("core_rel", ""),
        "orientation": "1",
        "score": row.get("score", row.get("ic50", "")),
        "context_status": "exact_peptide_query",
        "raw_response": "",
    }


def parse_netmhcii_batch(
    text: str,
    predictor: str,
    arms: Sequence[Mapping[str, str]],
    raw_response: str,
) -> list[dict[str, Any]]:
    table = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    results = []
    for seq_num, arm in enumerate(arms, start=1):
        matches = [
            row for row in table
            if str(row.get("seq_num", "")) == str(seq_num)
            and str(row.get("peptide", "")).upper() == arm["sequence"].upper()
        ]
        if len(matches) != 1:
            raise ValueError(
                f"{predictor} batch has {len(matches)} exact rows for seq_num {seq_num}"
            )
        row = matches[0]
        results.append(
            {
                "arm_id": arm["arm_id"],
                "target_id": arm["target_id"],
                "allele": arm["allele"],
                "predictor": predictor,
                "percentile_rank": row.get("rank", ""),
                "core": row.get("core_peptide", ""),
                "core_reliability": row.get("core_rel", ""),
                "orientation": "1",
                "score": row.get("score", row.get("ic50", "")),
                "context_status": "exact_peptide_query",
                "raw_response": raw_response,
            }
        )
    return results


def _fetch_netmhcii(output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    arms = read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
    raw_dir = output_dir / "raw_responses/netmhciipan_4_3"
    raw_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    failures = []
    groups: dict[tuple[str, int], list[dict[str, str]]] = {}
    for arm in arms:
        groups.setdefault((arm["allele"], len(arm["sequence"])), []).append(arm)
    for group_arms in groups.values():
        group_arms.sort(key=lambda arm: arm["arm_id"])
    for (allele, length), group_arms in sorted(groups.items()):
        allele_token = allele.replace("HLA-", "").replace("*", "_").replace(":", "_")
        sequence_text = "\n".join(
            f">{arm['arm_id']}\n{arm['sequence']}" for arm in group_arms
        )
        for predictor, method in (
            ("netmhciipan_4_3_el", "netmhciipan_el-4.3"),
            ("netmhciipan_4_3_ba", "netmhciipan_ba-4.3"),
        ):
            raw_path = raw_dir / f"batch__{allele_token}__len{length}__{predictor}.tsv"
            try:
                if raw_path.exists():
                    data = raw_path.read_bytes()
                else:
                    last_error: Optional[Exception] = None
                    data = b""
                    for delay in (0, 10, 30):
                        if delay:
                            time.sleep(delay)
                        try:
                            data = _post_form(
                                IEDB_PREDICT_URL,
                                {
                                    "method": method,
                                    "sequence_text": sequence_text,
                                    "allele": allele,
                                    "length": str(length),
                                },
                                timeout=120,
                            )
                            break
                        except Exception as error:
                            last_error = error
                    if not data:
                        raise last_error or RuntimeError("empty NetMHCIIpan response")
                    raw_path.write_bytes(data)
                rows.extend(
                    parse_netmhcii_batch(
                        data.decode("utf-8", "replace"),
                        predictor,
                        group_arms,
                        str(raw_path.relative_to(output_dir)),
                    )
                )
            except Exception as error:
                failures.append(
                    {
                        "arm_ids": ";".join(arm["arm_id"] for arm in group_arms),
                        "predictor": predictor,
                        "error": f"{type(error).__name__}: {error}",
                    }
                )
    rows.sort(key=lambda row: (row["arm_id"], row["predictor"]))
    failures.sort(key=lambda row: (row.get("arm_ids", ""), row["predictor"]))
    return rows, failures


def _parse_mixmhc_output(
    path: Path,
    arms: Sequence[Mapping[str, str]],
    predictor: str,
    context_status: str,
) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        table = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if len(table) != len(arms):
        raise ValueError(f"MixMHC2pred returned {len(table)} rows for {len(arms)} arms")
    results = []
    for arm, row in zip(arms, table):
        allele_key = arm["allele"].replace("HLA-", "").replace("*", "_").replace(":", "_")
        rank = row.get(f"%Rank_{allele_key}") or row.get(f"%Rank_{arm['allele']}") or row.get("%Rank_best") or ""
        core_p1 = row.get(f"CoreP1_{allele_key}") or row.get(f"CoreP1_{arm['allele']}") or ""
        orientation = row.get(f"subSpec_{allele_key}") or "1"
        try:
            start = int(float(core_p1)) - 1
            core = arm["sequence"][start : start + 9]
        except (TypeError, ValueError):
            core = row.get("Core_best", "")
        results.append(
            {
                "arm_id": arm["arm_id"],
                "target_id": arm["target_id"],
                "allele": arm["allele"],
                "predictor": predictor,
                "percentile_rank": rank,
                "core": core,
                "core_reliability": "",
                "orientation": orientation,
                "score": "",
                "context_status": context_status,
                "raw_response": str(path.relative_to(path.parents[2])),
            }
        )
    return results


def _run_mixmhc(output_dir: Path, mixmhc_binary: Optional[Path]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    arms = read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
    failures = []
    if mixmhc_binary is None or not mixmhc_binary.exists():
        return [], [{"predictor": "mixmhc2pred_2_1", "error": "version-pinned binary not available"}]
    raw_dir = output_dir / "raw_responses/mixmhc2pred_2_1"
    raw_dir.mkdir(parents=True, exist_ok=True)
    allele_args = [arm["allele"] for arm in arms]
    unique_alleles = list(dict.fromkeys(allele_args))
    results = []
    for mode, input_name, predictor, extra in (
        ("context", "mixmhc2pred_context.tsv", "mixmhc2pred_2_1_context", []),
        ("no_context", "mixmhc2pred_no_context.tsv", "mixmhc2pred_2_1_no_context", ["--no_context"]),
    ):
        output_path = raw_dir / f"{mode}.tsv"
        command = [str(mixmhc_binary), "-i", str(output_dir / "prepared_inputs" / input_name), "-o", str(output_path), "-a", *unique_alleles, *extra]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        (raw_dir / f"{mode}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (raw_dir / f"{mode}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            failures.append({"predictor": predictor, "error": f"exit {completed.returncode}: {completed.stderr.strip()}"})
            continue
        try:
            results.extend(_parse_mixmhc_output(output_path, arms, predictor, mode))
        except Exception as error:
            failures.append({"predictor": predictor, "error": f"{type(error).__name__}: {error}"})
    return results, failures


def _fetch_hla_atlas(output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cache_dir = output_dir / "raw_responses/hla_ligand_atlas_2020_12"
    cache_dir.mkdir(parents=True, exist_ok=True)
    failures = []
    for table in ("peptides", "donors", "sample_hits", "protein_map"):
        path = cache_dir / f"{table}.tsv.gz"
        if path.exists():
            continue
        try:
            _download(f"{HLA_ATLAS_BASE}/{table}.tsv.gz", path, timeout=300)
        except Exception as error:
            failures.append({"table": table, "error": f"{type(error).__name__}: {error}"})
    required = [cache_dir / f"{table}.tsv.gz" for table in ("peptides", "donors", "sample_hits")]
    if not all(path.exists() for path in required):
        return [], failures
    def read_gzip_table(path: Path) -> list[dict[str, str]]:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    peptides = read_gzip_table(cache_dir / "peptides.tsv.gz")
    donors = read_gzip_table(cache_dir / "donors.tsv.gz")
    sample_hits = read_gzip_table(cache_dir / "sample_hits.tsv.gz")
    sequence_by_id = {row["peptide_sequence_id"]: row["peptide_sequence"] for row in peptides}
    alleles_by_donor: dict[str, set[str]] = {}
    for row in donors:
        alleles_by_donor.setdefault(row["donor"], set()).add(normalize_hla(row["hla_allele"]))
    arms = read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
    results = []
    for hit in sample_hits:
        if hit.get("hla_class") != "HLA-II":
            continue
        observed = sequence_by_id.get(hit["peptide_sequence_id"], "")
        if not observed:
            continue
        donor_alleles = alleles_by_donor.get(hit["donor"], set())
        for arm in arms:
            relation = sequence_relation(arm["sequence"], observed)
            if relation == "none" and sequence_relation(arm["core"], observed) == "none":
                continue
            exact_hla = normalize_hla(arm["allele"]) in donor_alleles
            class_ii_alleles = [item for item in donor_alleles if item.startswith(("HLA-DR", "HLA-DQ", "HLA-DP"))]
            monoallelic = exact_hla and len(class_ii_alleles) == 1
            results.append(
                {
                    "arm_id": arm["arm_id"],
                    "target_id": arm["target_id"],
                    "source": "HLA_Ligand_Atlas_2020.12",
                    "observed_sequence": observed,
                    "sequence_relation": relation if relation != "none" else "core_overlap",
                    "donor": hit["donor"],
                    "tissue": hit["tissue"],
                    "hla_class": hit["hla_class"],
                    "donor_hla_alleles": ";".join(sorted(donor_alleles)),
                    "exact_hla_compatible": exact_hla,
                    "monoallelic": monoallelic,
                    "hit_class": classify_ligand_hit(
                        arm["sequence"], observed, exact_hla=exact_hla,
                        monoallelic=monoallelic, target_core=arm["core"],
                    ),
                    "absence_is_not_a_negative": True,
                    "raw_response": "raw_responses/hla_ligand_atlas_2020_12/sample_hits.tsv.gz",
                }
            )
    results.sort(key=lambda row: (row["arm_id"], row["observed_sequence"], row["donor"], row["tissue"]))
    return results, failures


def _fetch_pxd068488_processed(output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Search the eight submitter-provided processed peptide tables without raw reprocessing."""
    arms = read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
    raw_dir = output_dir / "raw_responses/pxd068488"
    raw_dir.mkdir(parents=True, exist_ok=True)
    failures = []
    hits = []
    for donor in ("a", "b", "c", "d"):
        for fraction in ("DR2a", "DR2b"):
            filename = f"HD_{donor}_EBV_B_cells_{fraction}_peptide.tsv"
            path = raw_dir / filename
            try:
                if not path.exists():
                    _download(f"{PXD068488_FTP_BASE}/{filename}", path, timeout=180)
                with path.open(newline="", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
            except Exception as error:
                failures.append({"file": filename, "error": f"{type(error).__name__}: {error}"})
                continue
            for observed_row in rows:
                observed = str(observed_row.get("Peptide") or "").strip().upper()
                if not observed:
                    continue
                for arm in arms:
                    peptide_relation = sequence_relation(arm["sequence"], observed)
                    core_relation = sequence_relation(arm["core"], observed)
                    if peptide_relation == "none" and core_relation == "none":
                        continue
                    exact_hla = (
                        fraction == "DR2b"
                        and normalize_hla(arm["allele"]) == "HLA-DRB1*15:01"
                    )
                    hits.append(
                        {
                            "arm_id": arm["arm_id"],
                            "target_id": arm["target_id"],
                            "source": f"PXD068488_processed_{fraction}",
                            "observed_sequence": observed,
                            "sequence_relation": peptide_relation if peptide_relation != "none" else "core_overlap",
                            "donor": f"HD_{donor}",
                            "tissue": "EBV_infected_B_cells",
                            "hla_class": "HLA-II",
                            "donor_hla_alleles": "HLA-DRB1*15:01" if fraction == "DR2b" else "HLA-DRB5*01:01",
                            "exact_hla_compatible": exact_hla,
                            "monoallelic": False,
                            "hit_class": classify_ligand_hit(
                                arm["sequence"], observed, exact_hla=exact_hla,
                                monoallelic=False, target_core=arm["core"],
                            ),
                            "absence_is_not_a_negative": True,
                            "protein": observed_row.get("Protein", ""),
                            "protein_id": observed_row.get("Protein ID", ""),
                            "qvalue": observed_row.get("Qvalue", ""),
                            "spectral_count": observed_row.get("Spectral Count", ""),
                            "raw_response": str(path.relative_to(output_dir)),
                        }
                    )
    hits.sort(key=lambda row: (row["arm_id"], row["source"], row["donor"], row["observed_sequence"]))
    return hits, failures


def _fetch_human_common_variants(output_dir: Path) -> dict[str, Any]:
    """Cache canonical-transcript missense variants overlapping each self peptide."""
    arms = [
        row for row in read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
        if row.get("side") == "self"
    ]
    human_path = output_dir / "raw_responses/human_reviewed_canonical.fasta"
    human_records = parse_fasta(human_path) if human_path.exists() else parse_fasta(HUMAN_SOURCE_FASTA)
    raw_dir = output_dir / "raw_responses/gnomad"
    raw_dir.mkdir(parents=True, exist_ok=True)
    query = """
query GeneVariants($geneSymbol: String!) {
  gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
    symbol
    canonical_transcript_id
    variants(dataset: gnomad_r4) {
      variant_id
      hgvsp
      consequence
      transcript_id
      exome { ac an af }
      genome { ac an af }
    }
  }
}
""".strip()
    gene_payloads: dict[str, Mapping[str, Any]] = {}
    failures = []
    for gene in sorted({row["protein"] for row in arms}):
        raw_path = raw_dir / f"{gene}.json"
        try:
            if raw_path.exists():
                payload = json.loads(raw_path.read_text(encoding="utf-8"))
            else:
                body = json.dumps({"query": query, "variables": {"geneSymbol": gene}}).encode("utf-8")
                request = urllib.request.Request(
                    GNOMAD_GRAPHQL_URL,
                    data=body,
                    headers={"Content-Type": "application/json", "User-Agent": "EBV-MS-evidence-dossier/1.0"},
                )
                with urllib.request.urlopen(request, timeout=120, context=_ssl_context()) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                write_json(raw_path, payload)
            if payload.get("errors") or not payload.get("data", {}).get("gene"):
                raise ValueError(payload.get("errors") or "gene not resolved")
            gene_payloads[gene] = payload["data"]["gene"]
        except Exception as error:
            failures.append({"gene": gene, "error": f"{type(error).__name__}: {error}"})

    rows = []
    arm_status = []
    for arm in arms:
        gene = arm["protein"]
        payload = gene_payloads.get(gene)
        accession, _, start = _source_match(arm["sequence"], human_records)
        if not payload or start < 0:
            arm_status.append({"arm_id": arm["arm_id"], "gene": gene, "status": "not_evaluable"})
            continue
        canonical_transcript = str(payload.get("canonical_transcript_id") or "")
        peptide_start = start + 1
        peptide_end = start + len(arm["sequence"])
        overlap_count = 0
        common_count = 0
        for variant in payload.get("variants", []):
            if variant.get("transcript_id") != canonical_transcript:
                continue
            if variant.get("consequence") != "missense_variant":
                continue
            match = re.search(r"p\.[A-Za-z*]{3}(\d+)", str(variant.get("hgvsp") or ""))
            if not match:
                continue
            protein_position = int(match.group(1))
            if not peptide_start <= protein_position <= peptide_end:
                continue
            frequencies = []
            for sequencing_type in ("exome", "genome"):
                value = (variant.get(sequencing_type) or {}).get("af")
                if value is not None:
                    frequencies.append(float(value))
            max_af = max(frequencies) if frequencies else None
            is_common = max_af is not None and max_af >= COMMON_VARIANT_AF_THRESHOLD
            overlap_count += 1
            common_count += int(is_common)
            rows.append(
                {
                    "arm_id": arm["arm_id"],
                    "target_id": arm["target_id"],
                    "gene": gene,
                    "source_accession": accession,
                    "canonical_transcript_id": canonical_transcript,
                    "variant_id": variant.get("variant_id", ""),
                    "hgvsp": variant.get("hgvsp", ""),
                    "protein_position_1_based": protein_position,
                    "peptide_position_1_based": protein_position - peptide_start + 1,
                    "exome_af": (variant.get("exome") or {}).get("af"),
                    "genome_af": (variant.get("genome") or {}).get("af"),
                    "maximum_af": max_af,
                    "common_af_at_least_0_01": is_common,
                    "dataset": GNOMAD_DATASET,
                    "raw_response": str(raw_dir.relative_to(output_dir) / f"{gene}.json"),
                }
            )
        arm_status.append(
            {
                "arm_id": arm["arm_id"],
                "gene": gene,
                "status": "evaluable",
                "overlapping_missense_count": overlap_count,
                "overlapping_common_missense_count": common_count,
            }
        )
    rows.sort(key=lambda row: (row["arm_id"], int(row["protein_position_1_based"]), row["variant_id"]))
    write_csv(output_dir / "raw_responses/human_common_variant_records.csv", rows, [
        "arm_id", "target_id", "gene", "source_accession", "canonical_transcript_id",
        "variant_id", "hgvsp", "protein_position_1_based", "peptide_position_1_based",
        "exome_af", "genome_af", "maximum_af", "common_af_at_least_0_01", "dataset", "raw_response",
    ])
    write_csv(output_dir / "raw_responses/human_common_variant_status.csv", arm_status)
    return {
        "status": "complete" if not failures else "partial",
        "arm_count": len(arms),
        "variant_record_count": len(rows),
        "failure_count": len(failures),
        "failures": failures,
    }


def fetch_evidence(
    *,
    output_dir: Path = DEFAULT_OUT,
    mixmhc_binary: Optional[Path] = None,
) -> dict[str, Any]:
    if not (output_dir / "protocol_lock.json").exists():
        prepare_package(output_dir=output_dir)
    statuses: dict[str, Any] = {}
    _, statuses["iedb"] = _fetch_iedb_assays(output_dir)
    net_rows, net_failures = _fetch_netmhcii(output_dir)
    mix_rows, mix_failures = _run_mixmhc(output_dir, mixmhc_binary)
    predictor_rows = sorted(net_rows + mix_rows, key=lambda row: (row["arm_id"], row["predictor"]))
    predictor_fields = [
        "arm_id", "target_id", "allele", "predictor", "percentile_rank", "core",
        "core_reliability", "orientation", "score", "context_status", "raw_response",
    ]
    write_csv(output_dir / "raw_responses/predictor_records.csv", predictor_rows, predictor_fields)
    statuses["predictors"] = {"record_count": len(predictor_rows), "failures": net_failures + mix_failures}
    atlas_rows, atlas_failures = _fetch_hla_atlas(output_dir)
    atlas_fields = [
        "arm_id", "target_id", "source", "observed_sequence", "sequence_relation", "donor", "tissue",
        "hla_class", "donor_hla_alleles", "exact_hla_compatible", "monoallelic", "hit_class",
        "absence_is_not_a_negative", "raw_response",
    ]
    write_csv(output_dir / "raw_responses/hla_ligand_atlas_hits.csv", atlas_rows, atlas_fields)
    statuses["hla_ligand_atlas"] = {"hit_count": len(atlas_rows), "failures": atlas_failures}

    pxd_rows, pxd_failures = _fetch_pxd068488_processed(output_dir)
    pxd_fields = atlas_fields + ["protein", "protein_id", "qvalue", "spectral_count"]
    write_csv(output_dir / "raw_responses/pxd068488_published_hits.csv", pxd_rows, pxd_fields)

    for source_id, url, name in (
        ("human_reviewed_reference_proteome", UNIPROT_HUMAN_REVIEWED_URL, "human_reviewed_canonical.fasta"),
        ("ebv_sequence_collection", UNIPROT_EBV_URL, "ebv_uniprot_sequences.fasta"),
    ):
        path = output_dir / "raw_responses" / name
        if not path.exists():
            try:
                _download(url, path, timeout=300)
                statuses[source_id] = {"status": "cached", "sha256": sha256_file(path)}
            except Exception as error:
                statuses[source_id] = {"status": "not_evaluable", "error": f"{type(error).__name__}: {error}"}
        else:
            statuses[source_id] = {"status": "cached", "sha256": sha256_file(path)}

    statuses["human_common_variants"] = _fetch_human_common_variants(output_dir)

    pxd_status = {
        "source": PXD068488_URL,
        "status": "processed_tables_searched" if not pxd_failures else "partial",
        "published_result_scope": "eight submitter-provided HD EBV B-cell DR2a/DR2b peptide TSV tables",
        "processed_hit_count": len(pxd_rows),
        "failure_count": len(pxd_failures),
        "failures": pxd_failures,
        "raw_reprocessing_performed": False,
    }
    write_json(output_dir / "raw_responses/pxd068488_status.json", pxd_status)
    statuses["pxd068488"] = pxd_status
    write_json(output_dir / "fetch_status.json", statuses)
    write_json(output_dir / "source_manifest.json", _source_manifest(output_dir))
    return statuses


def _assay_summary(arm_id: str, rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    selected = [row for row in rows if row.get("arm_id") == arm_id]
    positives = [row for row in selected if "positive" in str(row.get("qualitative_measure", "")).lower()]
    classes = [row.get("evidence_class", "none") for row in positives]
    return {
        "iedb_record_count": len(selected),
        "iedb_positive_record_count": len(positives),
        "iedb_best_evidence_class": next(
            (item for item in ("exact_sequence_exact_hla", "overlap_exact_hla", "other_human_hla", "class_i", "nonhuman", "untyped") if item in classes),
            "none",
        ),
        "iedb_exact_hla_positive": any(item in {"exact_sequence_exact_hla", "overlap_exact_hla"} for item in classes),
    }


def _ligand_summary(arm_id: str, rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    selected = [row for row in rows if row.get("arm_id") == arm_id]
    exact = [row for row in selected if row.get("hit_class") == "exact_sequence_monoallelic_exact_hla"]
    compatible = [row for row in selected if str(row.get("exact_hla_compatible", "")).lower() == "true"]
    return {
        "immunopeptidome_hit_count": len(selected),
        "immunopeptidome_exact_monoallelic_hit_count": len(exact),
        "immunopeptidome_exact_hla_compatible_hit_count": len(compatible),
        "immunopeptidome_status": "observed_support" if compatible else "not_observed_missing_evidence",
        "immunopeptidome_absence_is_negative": False,
    }


def _homologous_sequences(peptide: str, records: Sequence[Mapping[str, str]]) -> list[str]:
    homologs = []
    for record in records:
        sequence = record["sequence"]
        if peptide in sequence:
            homologs.append(sequence)
            continue
        best_identity = 0
        for start in range(max(0, len(sequence) - len(peptide) + 1)):
            window = sequence[start : start + len(peptide)]
            if len(window) == len(peptide):
                best_identity = max(best_identity, sum(a == b for a, b in zip(peptide, window)))
        if best_identity / len(peptide) >= 0.6:
            homologs.append(sequence)
    return homologs


def _read_cached_records(primary: Path, fallback: Path) -> tuple[list[dict[str, str]], str]:
    if primary.exists():
        return parse_fasta(primary), "full_cached_public_source"
    return parse_fasta(fallback), "fallback_candidate_source_records"


def _blank_predictor_summary() -> dict[str, Any]:
    return summarize_predictor_evidence([], "AAAAAAAAA")


def analyze_package(
    *,
    output_dir: Path = DEFAULT_OUT,
    v3_dir: Path = DEFAULT_V3,
) -> dict[str, Any]:
    if not (output_dir / "protocol_lock.json").exists():
        prepare_package(output_dir=output_dir)
    targets = read_csv(output_dir / "frozen_candidate_registry.csv")
    arms = read_csv(output_dir / "prepared_inputs/peptide_arm_registry.csv")
    iedb_rows = read_csv(output_dir / "raw_responses/iedb_assay_records.csv")
    predictor_rows = read_csv(output_dir / "raw_responses/predictor_records.csv")
    ligand_rows = (
        read_csv(output_dir / "raw_responses/hla_ligand_atlas_hits.csv")
        + read_csv(output_dir / "raw_responses/pxd068488_published_hits.csv")
    )
    variant_rows = read_csv(output_dir / "raw_responses/human_common_variant_records.csv")
    variant_status_rows = read_csv(output_dir / "raw_responses/human_common_variant_status.csv")
    human_records, human_scope = _read_cached_records(
        output_dir / "raw_responses/human_reviewed_canonical.fasta", HUMAN_SOURCE_FASTA
    )
    ebv_records, ebv_scope = _read_cached_records(
        output_dir / "raw_responses/ebv_uniprot_sequences.fasta", EBV_SOURCE_FASTA
    )

    arm_evidence = []
    predictor_comparisons = []
    conservation_rows = []
    for arm in arms:
        records = [row for row in predictor_rows if row.get("arm_id") == arm["arm_id"]]
        predictor = summarize_predictor_evidence(records, arm["core"])
        assay = _assay_summary(arm["arm_id"], iedb_rows)
        ligand = _ligand_summary(arm["arm_id"], ligand_rows)
        source_collection = ebv_records if arm["side"] == "ebv" else human_records
        resolved_accession, resolved_parent, resolved_start = _source_match(
            arm["sequence"], source_collection
        )
        source_record_status = (
            "exact_parent_sequence_match" if resolved_accession else "not_evaluable"
        )
        identity_conflict = not bool(resolved_accession)
        arm_row = {
            **arm,
            "source_accession": resolved_accession,
            "source_record_status": source_record_status,
            "source_start_1_based": resolved_start + 1 if resolved_start >= 0 else "",
            **assay,
            **predictor,
            **ligand,
            "identity_or_hla_conflict": identity_conflict,
            "claim_boundary": CLAIM_BOUNDARY,
        }
        arm_evidence.append(arm_row)
        predictor_comparisons.append(
            {
                "arm_id": arm["arm_id"],
                "target_id": arm["target_id"],
                "allele": arm["allele"],
                "declared_core": arm["core"],
                **predictor,
            }
        )
        sequence_collection = ebv_records if arm["side"] == "ebv" else [
            row for row in human_records if arm["sequence"] in row["sequence"] or arm["core"] in row["sequence"]
        ]
        homologs = _homologous_sequences(arm["sequence"], sequence_collection)
        arm_variants = [row for row in variant_rows if row.get("arm_id") == arm["arm_id"]]
        common_variants = [
            row for row in arm_variants
            if str(row.get("common_af_at_least_0_01", "")).lower() == "true"
        ]
        variant_status = next(
            (row.get("status", "not_evaluable") for row in variant_status_rows if row.get("arm_id") == arm["arm_id"]),
            "not_applicable_viral_arm" if arm["side"] == "ebv" else "not_evaluable",
        )
        changed_core_positions = []
        changed_tcr_positions = []
        core_start = int(arm["declared_core_start_1_based"])
        for variant in common_variants:
            peptide_position = int(variant["peptide_position_1_based"])
            core_position = peptide_position - core_start + 1
            if 1 <= core_position <= 9:
                changed_core_positions.append(f"P{core_position}")
                if core_position in {2, 3, 5, 7, 8}:
                    changed_tcr_positions.append(f"P{core_position}")
        conservation_rows.append(
            {
                "arm_id": arm["arm_id"],
                "target_id": arm["target_id"],
                "side": arm["side"],
                "sequence_collection_scope": ebv_scope if arm["side"] == "ebv" else human_scope,
                **summarize_conservation(arm["sequence"], arm["core"], homologs),
                "common_variant_status": variant_status,
                "overlapping_missense_variant_count": len(arm_variants),
                "overlapping_common_missense_variant_count": len(common_variants),
                "changed_register_positions": ";".join(sorted(set(changed_core_positions))),
                "changed_tcr_facing_positions": ";".join(sorted(set(changed_tcr_positions))),
            }
        )
    arm_evidence.sort(key=lambda row: row["arm_id"])
    write_csv(output_dir / "peptide_arm_evidence.csv", arm_evidence)
    write_csv(output_dir / "predictor_register_comparison.csv", predictor_comparisons)
    write_csv(output_dir / "iedb_assay_provenance.csv", iedb_rows, [
        "arm_id", "target_id", "endpoint", "query_scope", "record_index", "epitope_sequence",
        "sequence_relation", "mhc_allele", "mhc_class", "host_organism", "qualitative_measure",
        "assay_name", "parent_source_antigen", "source_organism", "pubmed_id", "pdb_id",
        "evidence_class", "raw_response",
    ])
    write_csv(output_dir / "immunopeptidome_hits.csv", ligand_rows, [
        "arm_id", "target_id", "source", "observed_sequence", "sequence_relation", "donor", "tissue",
        "hla_class", "donor_hla_alleles", "exact_hla_compatible", "monoallelic", "hit_class",
        "absence_is_not_a_negative", "protein", "protein_id", "qvalue", "spectral_count", "raw_response",
    ])
    write_csv(output_dir / "conservation_results.csv", conservation_rows)

    rarity_rows = []
    conditioned_rarity_rows = []
    nearest_rows = []
    candidate_rows = []
    recommendations = []
    arm_by_id = {row["arm_id"]: row for row in arm_evidence}
    v3_rows = read_csv(v3_dir / "v3_all_hla_ranked_pairs.csv")
    combined_drb1501_rows = read_csv(v3_dir / "combined_drb1501_v3_ranked_pairs.csv")
    for target in targets:
        ebv_arm = arm_by_id[f"{target['target_id']}__ebv"]
        self_arm = arm_by_id[f"{target['target_id']}__self"]
        forward = scan_similarity_rarity_fast(
            query_core=target["ebv_core"],
            paired_core=target["self_core"],
            database_records=human_records,
            exclude_accession=self_arm.get("source_accession", ""),
            exclude_core=target["self_core"],
        )
        reciprocal = scan_similarity_rarity_fast(
            query_core=target["self_core"],
            paired_core=target["ebv_core"],
            database_records=ebv_records,
            exclude_accession=ebv_arm.get("source_accession", ""),
            exclude_core=target["ebv_core"],
        )
        percentiles = [value for value in (forward["empirical_percentile"], reciprocal["empirical_percentile"]) if value is not None]
        conservative_percentile = max(percentiles) if len(percentiles) == 2 else None
        rarity_status = "evaluable" if conservative_percentile is not None else "not_evaluable"
        conditioned_library = (
            combined_drb1501_rows
            if normalize_hla(target["allele"]) == "HLA-DRB1*15:01"
            else v3_rows
        )
        conditioned = presentation_conditioned_rarity(
            target_pair_id=target["pair_id"],
            target_allele=target["allele"],
            candidate_rows=conditioned_library,
        )
        conditioned_rarity_rows.append(
            {
                "target_id": target["target_id"],
                "pair_id": target["pair_id"],
                "allele": target["allele"],
                "library_scope": "combined_drb1501_2043" if normalize_hla(target["allele"]) == "HLA-DRB1*15:01" else "v3_exact_hla_1600",
                **conditioned,
                "used_in_stage1_gate": False,
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )
        for direction, result, scope in (
            ("ebv_core_vs_human_proteome", forward, human_scope),
            ("self_core_vs_ebv_sequences", reciprocal, ebv_scope),
        ):
            rarity_rows.append(
                {
                    "target_id": target["target_id"],
                    "allele": target["allele"],
                    "direction": direction,
                    "database_scope": scope,
                    "query_core": result["query_core"],
                    "paired_core": result["paired_core"],
                    "evaluated_window_count": result["evaluated_window_count"],
                    "excluded_target_window_count": result["excluded_target_window_count"],
                    "at_least_as_good_count": result["at_least_as_good_count"],
                    "empirical_percentile": result["empirical_percentile"],
                    "presentation_conditioned": False,
                    "claim_boundary": CLAIM_BOUNDARY,
                }
            )
            for neighbor in result["nearest_neighbors"]:
                nearest_rows.append(
                    {
                        "target_id": target["target_id"],
                        "allele": target["allele"],
                        "direction": direction,
                        **neighbor,
                    }
                )
        stage1 = classify_stage1(ebv_arm, self_arm, rarity_percentile=conservative_percentile)
        candidate = {
            "target_id": target["target_id"],
            "pair_id": target["pair_id"],
            "allele": target["allele"],
            "ebv_protein": target["ebv_protein"],
            "ebv_sequence": target["ebv_sequence"],
            "ebv_core": target["ebv_core"],
            "self_protein": target["self_protein"],
            "self_sequence": target["self_sequence"],
            "self_core": target["self_core"],
            "ebv_binding_consensus": ebv_arm["binding_consensus"],
            "self_binding_consensus": self_arm["binding_consensus"],
            "ebv_register_consensus_matches_declared": ebv_arm["register_consensus_matches_declared"],
            "self_register_consensus_matches_declared": self_arm["register_consensus_matches_declared"],
            "iedb_exact_hla_positive_arm_count": sum(bool(row["iedb_exact_hla_positive"]) for row in (ebv_arm, self_arm)),
            "immunopeptidome_exact_hla_compatible_arm_count": sum(int(row["immunopeptidome_exact_hla_compatible_hit_count"]) > 0 for row in (ebv_arm, self_arm)),
            "forward_rarity_percentile": forward["empirical_percentile"],
            "reciprocal_rarity_percentile": reciprocal["empirical_percentile"],
            "conservative_rarity_percentile": conservative_percentile,
            "rarity_status": rarity_status,
            "presentation_conditioned_rarity_status": conditioned["status"],
            "presentation_conditioned_rarity_percentile": conditioned["empirical_percentile"],
            "stage1_status": stage1,
            "stage2_status": "not_evaluable_pending_experimental_binding_and_register",
            "single_composite_score": "not_created",
            "claim_boundary": CLAIM_BOUNDARY,
        }
        candidate_rows.append(candidate)
        recommendations.append(
            {
                **candidate,
                "recommended_next_experiment": (
                    "exact-HLA peptide binding plus nested-peptide register mapping"
                    if stage1 != "stage1_hold"
                    else "resolve predictor, identity, or HLA conflict before ordering assays"
                ),
                "tcell_assay_recommended_now": False,
            }
        )
    candidate_rows.sort(key=lambda row: row["target_id"])
    recommendations.sort(key=lambda row: row["target_id"])
    write_csv(output_dir / "proteome_rarity_summary.csv", rarity_rows)
    write_csv(output_dir / "presentation_conditioned_rarity.csv", conditioned_rarity_rows)
    write_csv(output_dir / "proteome_nearest_neighbors.csv", nearest_rows)
    write_csv(output_dir / "candidate_evidence_matrix.csv", candidate_rows)
    write_csv(output_dir / "stage1_assay_recommendations.csv", recommendations)
    stage1_gate = {
        "gate_name": "binding_and_register_assay_prioritization",
        "status": "complete" if all(row["rarity_status"] == "evaluable" for row in candidate_rows) else "partial",
        "candidate_count": len(candidate_rows),
        "high_priority_count": sum(row["stage1_status"] == "stage1_high_priority" for row in candidate_rows),
        "medium_priority_count": sum(row["stage1_status"] == "stage1_medium_priority" for row in candidate_rows),
        "hold_count": sum(row["stage1_status"] == "stage1_hold" for row in candidate_rows),
        **CLAIM_FIELDS,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "stage1_assay_gate.json", stage1_gate)
    write_json(output_dir / "stage2_tcell_gate.json", build_stage2_gate(candidate_rows))
    write_json(output_dir / "source_manifest.json", _source_manifest(output_dir))
    _write_raw_response_manifest(output_dir, arms)
    _write_candidate_dossiers(output_dir, candidate_rows, arm_by_id)
    _write_readme(output_dir, stage1_gate, human_scope, ebv_scope)
    checksums = _write_checksums(output_dir)
    return {
        "target_count": len(candidate_rows),
        "arm_count": len(arm_evidence),
        "stage1_gate": stage1_gate,
        "file_checksums": checksums,
    }


def _write_candidate_dossiers(
    output_dir: Path,
    candidates: Sequence[Mapping[str, Any]],
    arm_by_id: Mapping[str, Mapping[str, Any]],
) -> None:
    dossier_dir = output_dir / "candidate_dossiers"
    dossier_dir.mkdir(parents=True, exist_ok=True)
    for candidate in candidates:
        ebv = arm_by_id[f"{candidate['target_id']}__ebv"]
        self_arm = arm_by_id[f"{candidate['target_id']}__self"]
        lines = [
            f"# {candidate['target_id']}: {candidate['ebv_protein']} - {candidate['self_protein']}",
            "",
            f"- HLA: `{candidate['allele']}`",
            f"- EBV peptide/core: `{candidate['ebv_sequence']}` / `{candidate['ebv_core']}`",
            f"- Self peptide/core: `{candidate['self_sequence']}` / `{candidate['self_core']}`",
            f"- Stage-one status: `{candidate['stage1_status']}`",
            f"- Stage-two status: `{candidate['stage2_status']}`",
            "",
            "## Evidence",
            "",
            f"- EBV independent binding consensus: `{ebv['binding_consensus']}`; declared-register consensus: `{ebv['register_consensus_matches_declared']}`.",
            f"- Self independent binding consensus: `{self_arm['binding_consensus']}`; declared-register consensus: `{self_arm['register_consensus_matches_declared']}`.",
            f"- Exact-HLA positive IEDB arms: `{candidate['iedb_exact_hla_positive_arm_count']}/2`.",
            f"- Exact-HLA-compatible immunopeptidome arms: `{candidate['immunopeptidome_exact_hla_compatible_arm_count']}/2`.",
            f"- Conservative sequence-rarity percentile: `{candidate['conservative_rarity_percentile']}`.",
            "",
            "## Next experiment",
            "",
            "Measure binding of both exact peptides to the modeled HLA and map the P1-P9 register with nested peptides. Do not proceed to a cross-reactive T-cell claim from this dossier alone.",
            "",
            CLAIM_BOUNDARY,
            "",
        ]
        (dossier_dir / f"{candidate['target_id']}.md").write_text("\n".join(lines), encoding="utf-8")


def _write_readme(output_dir: Path, stage1_gate: Mapping[str, Any], human_scope: str, ebv_scope: str) -> None:
    text = f"""# Eight-candidate computational evidence dossier

This additive package audits the eight sequence-supported high-yield HLA-II candidates. It does not rerank the discovery universe or create a composite score.

## Result state

- Stage-one gate: `{stage1_gate['status']}`
- High priority: `{stage1_gate['high_priority_count']}`
- Medium priority: `{stage1_gate['medium_priority_count']}`
- Hold: `{stage1_gate['hold_count']}`
- Stage two: `not_evaluable_pending_experimental_binding_and_register`

## Evidence layers

- IEDB assay provenance is classified by exact sequence, exact HLA, MHC class, and host.
- NetMHCIIpan 4.3 EL/BA and MixMHC2pred 2.1 remain separate predictors.
- HLA Ligand Atlas release {HLA_ATLAS_RELEASE} hits distinguish exact, nested, monoallelic, and multiallelic evidence.
- All eight submitter-provided PXD068488 DR2a/DR2b processed peptide tables were searched; raw spectra were not reprocessed.
- gnomAD r4 canonical-transcript missense variants were checked with a locked common-frequency threshold of 1%.
- Human rarity database scope: `{human_scope}`.
- EBV reciprocal rarity scope: `{ebv_scope}`.
- Immunopeptidome absence is missing evidence, never a negative.

## Files

- `protocol_lock.json` and `source_manifest.json`: frozen rules, versions, sources, and checksums.
- `raw_response_manifest.csv`: checksum and exact target linkage for every cached raw file.
- `peptide_arm_evidence.csv`: sixteen peptide-arm evidence records.
- `candidate_evidence_matrix.csv`: eight candidate-level records.
- `iedb_assay_provenance.csv`: exact assay classifications with raw-response linkage.
- `predictor_register_comparison.csv`: independent predictor scores and core agreement.
- `immunopeptidome_hits.csv`: observed exact/nested ligand records.
- `proteome_rarity_summary.csv` and `proteome_nearest_neighbors.csv`: unconditioned sequence rarity diagnostics.
- `presentation_conditioned_rarity.csv`: separate exact-HLA candidate-library rarity diagnostic.
- `conservation_results.csv`: sequence conservation plus mapped common overlapping missense variants.
- `stage1_assay_recommendations.csv`, `stage1_assay_gate.json`, and `stage2_tcell_gate.json`: experimental funnel.
- `candidate_dossiers/`: human-readable evidence sheets.
- `SHA256SUMS.csv`: deterministic checksums.

## Claim boundary

{CLAIM_BOUNDARY}
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def _write_raw_response_manifest(
    output_dir: Path,
    arms: Sequence[Mapping[str, Any]],
) -> None:
    linkage: dict[str, set[str]] = {}
    for table in (
        "raw_responses/iedb_assay_records.csv",
        "raw_responses/predictor_records.csv",
        "raw_responses/hla_ligand_atlas_hits.csv",
        "raw_responses/pxd068488_published_hits.csv",
        "raw_responses/human_common_variant_records.csv",
    ):
        for row in read_csv(output_dir / table):
            raw_path = str(row.get("raw_response") or "")
            arm_id = str(row.get("arm_id") or "")
            if raw_path and arm_id:
                linkage.setdefault(raw_path, set()).add(arm_id)
    arms_by_gene: dict[str, set[str]] = {}
    for arm in arms:
        if arm.get("side") == "self":
            arms_by_gene.setdefault(str(arm["protein"]), set()).add(str(arm["arm_id"]))
    all_arm_ids = {str(arm["arm_id"]) for arm in arms}
    rows = []
    raw_root = output_dir / "raw_responses"
    for path in sorted(item for item in raw_root.rglob("*") if item.is_file()):
        relative = str(path.relative_to(output_dir))
        linked = set(linkage.get(relative, set()))
        if relative.startswith("raw_responses/gnomad/"):
            linked.update(arms_by_gene.get(path.stem, set()))
        if not linked and (
            relative.endswith((".fasta", ".tsv.gz"))
            or "/mixmhc2pred_2_1/" in relative
            or "/pxd068488/" in relative
        ):
            linked.update(all_arm_ids)
        target_ids = sorted({arm_id.split("__", 1)[0] for arm_id in linked})
        source_id = path.relative_to(raw_root).parts[0]
        rows.append(
            {
                "relative_path": relative,
                "source_id": source_id,
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
                "linked_arm_ids": ";".join(sorted(linked)),
                "linked_target_ids": ";".join(target_ids),
                "linkage_status": "exact_or_query_manifest" if linked else "package_level_status_or_standardized_output",
            }
        )
    write_csv(output_dir / "raw_response_manifest.csv", rows, [
        "relative_path", "source_id", "sha256", "size_bytes", "linked_arm_ids",
        "linked_target_ids", "linkage_status",
    ])


def _write_checksums(output_dir: Path) -> dict[str, str]:
    checksum_path = output_dir / "SHA256SUMS.csv"
    files = sorted(
        path for path in output_dir.rglob("*")
        if path.is_file() and path != checksum_path and "raw_responses" not in path.parts
    )
    rows = []
    checksums = {}
    for path in files:
        relative = str(path.relative_to(output_dir))
        digest = sha256_file(path)
        checksums[relative] = digest
        rows.append({"relative_path": relative, "sha256": digest})
    write_csv(checksum_path, rows, ["relative_path", "sha256"])
    return checksums


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "fetch", "analyze", "all"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--mixmhc-binary", type=Path, default=None)
    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare_package(output_dir=args.output_dir)
    elif args.command == "fetch":
        result = fetch_evidence(output_dir=args.output_dir, mixmhc_binary=args.mixmhc_binary)
    elif args.command == "analyze":
        result = analyze_package(output_dir=args.output_dir)
    else:
        prepare_package(output_dir=args.output_dir)
        fetch_evidence(output_dir=args.output_dir, mixmhc_binary=args.mixmhc_binary)
        result = analyze_package(output_dir=args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
