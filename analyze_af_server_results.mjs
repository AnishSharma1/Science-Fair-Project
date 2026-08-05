import fs from "node:fs/promises";
import path from "node:path";

const downloads = "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads";
const packageDir = "/Users/anishsharma/Documents/New project/outputs/ebv_ms_model_package";
const batchPath = `${packageDir}/alphafold_server_batch.json`;
const outputDir = `${packageDir}/results_analysis`;

async function walk(dir, found = []) {
  for (const entry of await fs.readdir(dir, {withFileTypes: true})) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) await walk(full, found);
    else if (entry.name.endsWith("_job_request.json")) found.push(full);
  }
  return found;
}
function number(value) { return typeof value === "number" && Number.isFinite(value) ? value : null; }
function mean(values) { const use = values.filter(v => v !== null); return use.length ? use.reduce((a, b) => a + b, 0) / use.length : null; }
function fmt(value) { return typeof value === "number" && Number.isFinite(value) ? value.toFixed(3) : ""; }
function conditionOf(name) { return name.replace(/_S\d+$/i, ""); }
function categoryOf(name) {
  if (name.startsWith("PMHC_")) return "pMHC control";
  if (name.startsWith("CAL_")) return "ternary calibration";
  if (name.startsWith("HY_")) return "Hy.2E11 positive";
  if (name.startsWith("DECOY_")) return "decoy / unknown";
  return "other";
}
function interfaceScores(summary, chainCount) {
  const matrix = summary.chain_pair_iptm;
  if (!Array.isArray(matrix)) return {tcrPmhcIptm: null, tcrPeptideIptm: null, peptideMhcIptm: null};
  const pair = (i, j) => number(matrix?.[i]?.[j]);
  if (chainCount === 5) {
    const tcrPmhc = [pair(0, 2), pair(0, 3), pair(0, 4), pair(1, 2), pair(1, 3), pair(1, 4)];
    return {tcrPmhcIptm: mean(tcrPmhc), tcrPeptideIptm: mean([pair(0, 4), pair(1, 4)]), peptideMhcIptm: mean([pair(2, 4), pair(3, 4)])};
  }
  return {tcrPmhcIptm: null, tcrPeptideIptm: null, peptideMhcIptm: mean([pair(0, 2), pair(1, 2)])};
}
async function requestRecords(files, expected) {
  const records = [];
  for (const file of files) {
    try {
      const raw = JSON.parse(await fs.readFile(file, "utf8"));
      for (const job of (Array.isArray(raw) ? raw : [raw])) {
        const name = String(job.name || "");
        if (expected.has(name.toLowerCase())) records.push({name, dir: path.dirname(file), requestFile: file, chainCount: job.sequences?.length || null});
      }
    } catch { /* Ignore malformed or unrelated download metadata. */ }
  }
  return records;
}
async function summarize(record) {
  const entries = await fs.readdir(record.dir);
  const summaryFiles = entries.filter(name => /_summary_confidences_\d+\.json$/.test(name));
  const summaries = [];
  for (const name of summaryFiles) {
    try {
      const payload = JSON.parse(await fs.readFile(path.join(record.dir, name), "utf8"));
      summaries.push({model: Number(name.match(/_(\d+)\.json$/)?.[1]), ...payload});
    } catch { /* incomplete download */ }
  }
  const ranked = summaries.map(s => ({summary: s, ranking: number(s.ranking_score)})).filter(x => x.ranking !== null).sort((a, b) => b.ranking - a.ranking);
  const best = ranked[0]?.summary || null;
  const metrics = best ? interfaceScores(best, record.chainCount) : interfaceScores({}, record.chainCount);
  return {
    ...record,
    modelsReturned: summaries.length,
    bestModel: best?.model ?? null,
    bestRanking: best ? number(best.ranking_score) : null,
    meanRanking: mean(summaries.map(s => number(s.ranking_score))),
    bestIptm: best ? number(best.iptm) : null,
    meanIptm: mean(summaries.map(s => number(s.iptm))),
    bestPtm: best ? number(best.ptm) : null,
    meanPtm: mean(summaries.map(s => number(s.ptm))),
    hasClash: best ? number(best.has_clash) : null,
    ...metrics,
  };
}
async function main() {
  const batch = JSON.parse(await fs.readFile(batchPath, "utf8"));
  const expected = new Map(batch.map(job => [job.name.toLowerCase(), job]));
  const requestFiles = await walk(downloads);
  const rawRecords = await requestRecords(requestFiles, expected);
  const byName = new Map();
  for (const record of rawRecords) {
    const key = record.name.toLowerCase();
    const current = byName.get(key);
    if (!current || record.dir.length < current.dir.length) byName.set(key, record);
  }
  const rows = [];
  for (const job of batch) {
    const record = byName.get(job.name.toLowerCase());
    if (!record) rows.push({name: job.name, condition: conditionOf(job.name), category: categoryOf(job.name), status: "missing", chainCount: job.sequences.length, modelsReturned: 0});
    else rows.push({...(await summarize(record)), condition: conditionOf(record.name), category: categoryOf(record.name), status: "completed"});
  }
  const headers = ["job_name", "condition", "category", "status", "chain_count", "models_returned", "best_model", "best_ranking_score", "mean_ranking_score", "best_iptm", "mean_iptm", "best_ptm", "mean_ptm", "tcr_pmhc_pair_iptm", "tcr_peptide_pair_iptm", "peptide_mhc_pair_iptm", "best_has_clash", "result_folder"];
  const tsv = [headers, ...rows.map(r => [r.name, r.condition, r.category, r.status, r.chainCount, r.modelsReturned, r.bestModel ?? "", fmt(r.bestRanking), fmt(r.meanRanking), fmt(r.bestIptm), fmt(r.meanIptm), fmt(r.bestPtm), fmt(r.meanPtm), fmt(r.tcrPmhcIptm), fmt(r.tcrPeptideIptm), fmt(r.peptideMhcIptm), r.hasClash ?? "", r.dir || ""])].map(row => row.join("\t")).join("\n") + "\n";
  const groups = new Map();
  for (const row of rows) {
    if (!groups.has(row.condition)) groups.set(row.condition, []);
    groups.get(row.condition).push(row);
  }
  const groupRows = [...groups.entries()].map(([condition, set]) => ({
    condition,
    category: set[0].category,
    completedSeeds: set.filter(r => r.status === "completed").length,
    missingSeeds: set.filter(r => r.status === "missing").length,
    meanBestRanking: mean(set.map(r => r.bestRanking)),
    meanBestIptm: mean(set.map(r => r.bestIptm)),
    meanTcrPmhcIptm: mean(set.map(r => r.tcrPmhcIptm)),
    meanTcrPeptideIptm: mean(set.map(r => r.tcrPeptideIptm)),
    meanPeptideMhcIptm: mean(set.map(r => r.peptideMhcIptm)),
  }));
  const groupTsv = [["condition", "category", "completed_seeds", "missing_seeds", "mean_best_ranking", "mean_best_iptm", "mean_tcr_pmhc_pair_iptm", "mean_tcr_peptide_pair_iptm", "mean_peptide_mhc_pair_iptm"], ...groupRows.map(r => [r.condition, r.category, r.completedSeeds, r.missingSeeds, fmt(r.meanBestRanking), fmt(r.meanBestIptm), fmt(r.meanTcrPmhcIptm), fmt(r.meanTcrPeptideIptm), fmt(r.meanPeptideMhcIptm)])].map(row => row.join("\t")).join("\n") + "\n";
  const missing = rows.filter(r => r.status === "missing").map(r => r.name);
  const retryJobs = batch.filter(job => missing.includes(job.name));
  if (retryJobs.some(job => job.modelSeeds.length !== 1)) throw new Error("Retry validation failed: every rerun must have exactly one seed.");
  const duplicateCount = rawRecords.length - byName.size;
  const report = `# AlphaFold Server download inventory\n\n- Expected seed jobs: ${batch.length}\n- Unique completed seed jobs found: ${rows.filter(r => r.status === "completed").length}\n- Missing seed jobs: ${missing.length}\n- Duplicate downloaded job folders: ${duplicateCount}\n- Completed jobs with five summary-confidence files: ${rows.filter(r => r.status === "completed" && r.modelsReturned === 5).length}\n\n## Missing seed jobs\n\n${missing.map(name => `- ${name}`).join("\n") || "- None"}\n\n## Interpretation boundary\n\nThe TSV files summarize AlphaFold confidence metrics only. They do not establish binding, cross-reactivity, or a functional negative. The ternary-calibration cases still require structural alignment to their experimental PDB references before interpreting Hy.2E11 or decoy rankings.\n`;
  await fs.mkdir(outputDir, {recursive: true});
  await fs.writeFile(`${outputDir}/af3_seed_job_metrics.tsv`, tsv, "utf8");
  await fs.writeFile(`${outputDir}/af3_condition_summary.tsv`, groupTsv, "utf8");
  await fs.writeFile(`${outputDir}/AF3_DOWNLOAD_INVENTORY.md`, report, "utf8");
  await fs.writeFile(`${outputDir}/af3_seed_job_metrics.json`, JSON.stringify(rows, null, 2) + "\n", "utf8");
  await fs.writeFile(`${packageDir}/alphafold_server_rerun_missing_${retryJobs.length}.json`, JSON.stringify(retryJobs, null, 2) + "\n", "utf8");
  console.log(JSON.stringify({completed: rows.filter(r => r.status === "completed").length, missing, retryBatch: `${packageDir}/alphafold_server_rerun_missing_${retryJobs.length}.json`, duplicateDownloadedFolders: duplicateCount, outputDir}, null, 2));
}
await main();
