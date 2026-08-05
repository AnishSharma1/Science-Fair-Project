import fs from "node:fs/promises";

const root = "/Users/anishsharma/Documents/New project/outputs/ebv_ms_model_package";
const source = `${root}/source_fasta`;
const out = `${root}/alphafold_server_batch.json`;
const account2Out = `${root}/alphafold_server_account_2_remaining_30.json`;
const account3Out = `${root}/alphafold_server_account_3_remaining_30.json`;
const manifestPath = `${root}/alphafold_server_batch_manifest.tsv`;
const seeds = [101, 202, 303, 404, 505];

function parseFasta(text) {
  const records = new Map();
  let header = null;
  let sequence = [];
  for (const line of text.trim().split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (header) records.set(header, sequence.join(""));
      header = line.slice(1);
      sequence = [];
    } else sequence.push(line.trim());
  }
  if (header) records.set(header, sequence.join(""));
  return records;
}
function pick(records, prefix) {
  const item = [...records.entries()].find(([header]) => header.startsWith(prefix));
  if (!item) throw new Error(`Missing FASTA record: ${prefix}`);
  return item[1];
}
function protein(sequence) {
  return {proteinChain: {sequence, count: 1, useStructureTemplate: false}};
}
function job(name, sequences) {
  return seeds.map(seed => ({name: `${name}_S${seed}`, modelSeeds: [seed], sequences: sequences.map(protein)}));
}

async function main() {
  const [f1ymm, f1zgl, f2wbj, f1h15] = await Promise.all(
    ["1YMM", "1ZGL", "2WBJ", "1H15"].map(async id => parseFasta(await fs.readFile(`${source}/${id}.fasta`, "utf8")))
  );
  const hyAlpha = "SQQGEEDPQALSIQEGENATMNCSYKTSINNLQWYRQNSGRGLVHLILIRSNEREKHSGRLRVTLDTSKKSSSLLITASRAADTASYFCATDSGGSYIPTFGRGTSLIVHPY";
  const hyBeta = "SAVISQKPSRDICQRGTSLTIQCQVDSQVTMIFWYRQQPGQSLTLIATANQGSEATYESGFVIDKFPISRPNLTFSTLTVSNMSPEDSSIYLCSAWPSGQGTYGYTFGSGTRLTVV";
  const obAlpha = pick(f1ymm, "1YMM_4|");
  const obBeta = pick(f1ymm, "1YMM_5|");
  const a3Alpha = pick(f1zgl, "1ZGL_4|");
  const a3Beta = pick(f1zgl, "1ZGL_5|");
  const ob2wbjAlpha = pick(f2wbj, "2WBJ_3|").slice(1); // Remove initial cloning Met only.
  const ob2wbjBeta = pick(f2wbj, "2WBJ_4|").slice(pick(f2wbj, "2WBJ_4|").indexOf("AVVSQHPS"));
  const mhcDrb1 = [pick(f1ymm, "1YMM_1|"), pick(f1ymm, "1YMM_2|")];
  const mhcDrb5 = [pick(f1h15, "1H15_1|"), pick(f1h15, "1H15_2|")];
  const mhc2wbj = [pick(f2wbj, "2WBJ_1|"), pick(f2wbj, "2WBJ_2|")];
  const mhc1zgl = [pick(f1zgl, "1ZGL_1|"), pick(f1zgl, "1ZGL_2|")];
  const pmhc = {
    MBP_DRB1: [...mhcDrb1, "ENPVVHFFKNIVTPR"],
    BALF5_FULL15: [...mhcDrb5, "TGGVYHFVKKHVHES"],
    BALF5_CORE14: [...mhcDrb5, "GGVYHFVKKHVHES"],
    MBP_DRB5: [...mhc1zgl, "VHFFKNIVTPRTP"],
    ENGA_2WBJ: [...mhc2wbj, "MDFARVHFISALHG"],
  };
  const tcr = {HY: [hyAlpha, hyBeta], OB: [obAlpha, obBeta], A3: [a3Alpha, a3Beta]};
  const jobs = [
    job("PMHC_MBP_DRB1", pmhc.MBP_DRB1),
    job("PMHC_BALF5_FULL15", pmhc.BALF5_FULL15),
    job("PMHC_BALF5_CORE14", pmhc.BALF5_CORE14),
    job("PMHC_MBP_DRB5", pmhc.MBP_DRB5),
    job("PMHC_ENGA_2WBJ", pmhc.ENGA_2WBJ),
    job("CAL_1YMM", [...tcr.OB, ...pmhc.MBP_DRB1]),
    job("CAL_1ZGL", [...tcr.A3, ...pmhc.MBP_DRB5]),
    job("CAL_2WBJ", [ob2wbjAlpha, ob2wbjBeta, ...pmhc.ENGA_2WBJ]),
    job("HY_MBP_DRB1", [...tcr.HY, ...pmhc.MBP_DRB1]),
    job("HY_BALF5_FULL15", [...tcr.HY, ...pmhc.BALF5_FULL15]),
    job("HY_BALF5_CORE14", [...tcr.HY, ...pmhc.BALF5_CORE14]),
    job("DECOY_01_HY_MBP_DRB5", [...tcr.HY, ...pmhc.MBP_DRB5]),
    job("DECOY_02_HY_ENGA_DRB1", [...tcr.HY, ...pmhc.ENGA_2WBJ]),
    job("DECOY_03_OB_BALF5", [...tcr.OB, ...pmhc.BALF5_FULL15]),
    job("DECOY_04_OB_MBP_DRB5", [...tcr.OB, ...pmhc.MBP_DRB5]),
    job("DECOY_05_3A6_MBP_DRB1", [...tcr.A3, ...pmhc.MBP_DRB1]),
    job("DECOY_06_3A6_BALF5", [...tcr.A3, ...pmhc.BALF5_FULL15]),
    job("DECOY_07_3A6_ENGA", [...tcr.A3, ...pmhc.ENGA_2WBJ]),
  ].flat();
  if (jobs.length !== 90 || jobs.some(j => ![3, 5].includes(j.sequences.length)) || jobs.some(j => j.modelSeeds.length !== 1)) {
    throw new Error("Batch validation failed: expected 90 jobs, 3/5 chains, and exactly one seed each.");
  }
  if (JSON.stringify(jobs).includes("SGGGSGGGGG")) throw new Error("Batch validation failed: 2WBJ fusion linker leaked into input.");
  await fs.writeFile(out, `${JSON.stringify(jobs, null, 2)}\n`, "utf8");
  const account2Bases = [
    "CAL_1ZGL", "CAL_2WBJ", "HY_MBP_DRB1", "HY_BALF5_FULL15", "HY_BALF5_CORE14", "DECOY_01_HY_MBP_DRB5",
  ];
  const account3Bases = [
    "DECOY_02_HY_ENGA_DRB1", "DECOY_03_OB_BALF5", "DECOY_04_OB_MBP_DRB5", "DECOY_05_3A6_MBP_DRB1", "DECOY_06_3A6_BALF5", "DECOY_07_3A6_ENGA",
  ];
  const selectBases = bases => jobs.filter(item => bases.some(base => item.name.startsWith(`${base}_S`)));
  const account2Jobs = selectBases(account2Bases);
  const account3Jobs = selectBases(account3Bases);
  if (account2Jobs.length !== 30 || account3Jobs.length !== 30 || new Set([...account2Jobs, ...account3Jobs].map(item => item.name)).size !== 60) {
    throw new Error("Split validation failed: expected two disjoint 30-job files.");
  }
  await fs.writeFile(account2Out, `${JSON.stringify(account2Jobs, null, 2)}\n`, "utf8");
  await fs.writeFile(account3Out, `${JSON.stringify(account3Jobs, null, 2)}\n`, "utf8");
  const rows = [
    ["job_name", "category", "interpretation"],
    ["PMHC_MBP_DRB1", "pMHC control", "pMHC geometry against 1BX2"],
    ["PMHC_BALF5_FULL15", "pMHC control", "Biological 15-mer"],
    ["PMHC_BALF5_CORE14", "pMHC control", "Direct 1H15 deposited-core comparison"],
    ["PMHC_MBP_DRB5", "pMHC control", "pMHC geometry against 1ZGL"],
    ["PMHC_ENGA_2WBJ", "pMHC control", "Strict deposited 2WBJ peptide condition"],
    ["CAL_1YMM", "ternary calibration", "Experimental OB.1A12/MBP/DRB1 reference"],
    ["CAL_1ZGL", "ternary calibration", "Experimental 3A6/MBP/DRB5 reference"],
    ["CAL_2WBJ", "ternary calibration", "Experimental OB.1A12/ENGA/DRB1 reference; fusion linker removed"],
    ["HY_MBP_DRB1", "functional positive", "Hy.2E11 self pMHC"],
    ["HY_BALF5_FULL15", "functional positive", "Hy.2E11 biological EBV 15-mer"],
    ["HY_BALF5_CORE14", "structural sensitivity", "Hy.2E11 1H15 deposited-core condition"],
    ...Array.from({length: 7}, (_, i) => [`DECOY_${String(i + 1).padStart(2, "0")}_`, "decoy / unknown", "Selectivity stress test only; not a functional negative"]),
  ];
  await fs.writeFile(manifestPath, `${rows.map(r => r.join("\t")).join("\n")}\n`, "utf8");
  console.log(JSON.stringify({out, jobs: jobs.length, account2Out, account2Jobs: account2Jobs.length, account3Out, account3Jobs: account3Jobs.length, seedsPerJob: 1}, null, 2));
}
await main();
