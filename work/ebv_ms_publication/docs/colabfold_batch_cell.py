"""Paste this entire cell into the ColabFold notebook after installation.

It uploads either one multi-entry FASTA or one ``id,sequence`` CSV, runs all
pMHC-II candidates in one ColabFold-Multimer workflow, and copies the results
to Google Drive.
"""

# Upload one prepared batch file:
# - processed/pmhc_colabfold_inputs.fasta
# - processed/pmhc_colabfold_batch.csv
from google.colab import files, drive
import csv
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

uploaded = files.upload()
if len(uploaded) != 1:
    raise ValueError(f"Upload exactly one batch file; found {list(uploaded)}")
uploaded_file = Path(next(iter(uploaded)))


def parse_fasta(path: Path) -> list[dict[str, str]]:
    rows = []
    header = None
    seq_lines = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                rows.append({"id": header, "sequence": "".join(seq_lines)})
            header = line[1:].split("|", 1)[0].strip()
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        rows.append({"id": header, "sequence": "".join(seq_lines)})
    return rows


def parse_batch(path: Path) -> list[dict[str, str]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows or set(rows[0]) != {"id", "sequence"}:
            raise ValueError("The CSV must have exactly two columns: id,sequence")
        return rows
    if suffix in {".fasta", ".fa", ".faa"}:
        rows = parse_fasta(path)
        if not rows:
            raise ValueError("The FASTA file did not contain any sequences")
        return rows
    raise ValueError("Upload a .csv, .fasta, .fa, or .faa batch file")


rows = parse_batch(uploaded_file)
if len({row["id"] for row in rows}) != len(rows):
    raise ValueError("Candidate IDs must be unique")
for row in rows:
    row["id"] = row["id"].strip()
    row["sequence"] = "".join(row["sequence"].split()).upper()
    if not row["id"]:
        raise ValueError("Every candidate must have a non-empty ID")
    if row["sequence"].count(":") != 2:
        raise ValueError(f"{row['id']} does not contain DRA:DRB:peptide chains")
    if any(aa not in "ACDEFGHIKLMNPQRSTVWY:" for aa in row["sequence"]):
        raise ValueError(f"{row['id']} contains a non-standard amino-acid symbol")

BATCH_NAME = "ebv_ms_pmhc_batch"
RUN_DIR = Path("/content") / BATCH_NAME
RUN_DIR.mkdir(parents=True, exist_ok=True)
batch_csv = RUN_DIR / "pmhc_colabfold_batch.csv"
shutil.copy2(uploaded_file, RUN_DIR / uploaded_file.name)
with batch_csv.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["id", "sequence"])
    writer.writeheader()
    writer.writerows(rows)

settings = {
    "model_type": "alphafold2_multimer_v3",
    "msa_mode": "mmseqs2_uniref_env",
    "pair_mode": "unpaired_paired",
    "num_models": 5,
    "num_recycles": 3,
    "num_seeds": 1,
    "num_relax": 0,
    "save_all": True,
    "calc_extra_ptm": True,
}
(RUN_DIR / "run_metadata.json").write_text(
    json.dumps(
        {
            "candidate_count": len(rows),
            "candidate_ids": [row["id"] for row in rows],
            "input_sha256": hashlib.sha256(batch_csv.read_bytes()).hexdigest(),
            "settings": settings,
        },
        indent=2,
    )
)

def ensure_colabfold():
    try:
        from colabfold.batch import get_queries, run, set_model_type
        from colabfold.download import download_alphafold_params
        from colabfold.utils import setup_logging
        return get_queries, run, set_model_type, download_alphafold_params, setup_logging
    except ModuleNotFoundError:
        print("ColabFold is missing in this runtime. Installing it now...")
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "-q",
                "--no-warn-conflicts",
                "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold",
            ],
            check=True,
        )
        for pattern in [
            "/usr/local/lib/python*/dist-packages/tensorflow/core/kernels/libtfkernel_sobol_op.so",
            "/usr/local/lib/python*/dist-packages/tensorflow/lite/python/*/*.so",
        ]:
            for match in glob.glob(pattern):
                try:
                    os.remove(match)
                except OSError:
                    pass
        from colabfold.batch import get_queries, run, set_model_type
        from colabfold.download import download_alphafold_params
        from colabfold.utils import setup_logging
        return get_queries, run, set_model_type, download_alphafold_params, setup_logging


get_queries, run, set_model_type, download_alphafold_params, setup_logging = ensure_colabfold()

queries, is_complex = get_queries(str(batch_csv))
if not is_complex:
    raise ValueError("The uploaded batch was not recognized as a complex")
model_type = set_model_type(is_complex, settings["model_type"])
setup_logging(RUN_DIR / "log.txt")
download_alphafold_params(model_type, Path("/content"))

run(
    queries=queries,
    result_dir=str(RUN_DIR),
    use_templates=False,
    custom_template_path=None,
    num_relax=settings["num_relax"],
    msa_mode=settings["msa_mode"],
    model_type=model_type,
    num_models=settings["num_models"],
    num_recycles=settings["num_recycles"],
    relax_max_iterations=200,
    recycle_early_stop_tolerance=0.5,
    num_seeds=settings["num_seeds"],
    use_dropout=False,
    model_order=[1, 2, 3, 4, 5],
    is_complex=is_complex,
    data_dir=Path("/content"),
    keep_existing_results=False,
    rank_by="auto",
    pair_mode=settings["pair_mode"],
    pairing_strategy="greedy",
    stop_at_score=100.0,
    prediction_callback=None,
    dpi=200,
    zip_results=False,
    save_all=settings["save_all"],
    max_msa=None,
    use_cluster_profile=False,
    input_features_callback=None,
    save_recycles=False,
    user_agent="colabfold/google-colab-batch",
    calc_extra_ptm=settings["calc_extra_ptm"],
)

# Preserve the full run in Drive for later QA and publication analysis.
drive.mount("/content/drive")
drive_dir = Path("/content/drive/MyDrive") / BATCH_NAME
if drive_dir.exists():
    from datetime import datetime
    drive_dir = drive_dir.with_name(
        f"{BATCH_NAME}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
shutil.copytree(RUN_DIR, drive_dir)
archive = shutil.make_archive(str(drive_dir), "zip", root_dir=drive_dir.parent, base_dir=drive_dir.name)
print(f"Completed {len(rows)} candidates")
print(f"Drive folder: {drive_dir}")
print(f"Drive archive: {archive}")
