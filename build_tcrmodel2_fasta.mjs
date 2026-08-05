import fs from "node:fs/promises";

const project = "/Users/anishsharma/Documents/New project";
const modelPackage = `${project}/outputs/ebv_ms_model_package`;
const source = `${modelPackage}/source_fasta`;
const out = `${modelPackage}/tcrmodel2_fasta_submissions`;

function parseFasta(text) {
  const records = new Map();
  let header;
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
  const match = [...records.entries()].find(([header]) => header.startsWith(prefix));
  if (!match) throw new Error(`Missing ${prefix}`);
  return match[1];
}

function between(sequence, start, end) {
  const a = sequence.indexOf(start);
  const b = sequence.indexOf(end, a);
  if (a < 0 || b < 0) throw new Error(`Could not trim sequence using ${start} … ${end}`);
  return sequence.slice(a, b + end.length);
}

function fasta(name, records) {
  return records.map(([header, sequence]) => `>${header}\n${sequence}`).join("\n") + "\n";
}

function mhcBindingDomain(sequence, start, length) {
  const offset = sequence.indexOf(start);
  if (offset < 0) throw new Error(`Could not find MHC domain start ${start}`);
  const domain = sequence.slice(offset, offset + length);
  if (domain.length !== length) throw new Error(`MHC domain ${start} was unexpectedly short`);
  return domain;
}

async function main() {
  const [f1ymm, f1zgl, f2wbj, f1h15] = await Promise.all(
    ["1YMM", "1ZGL", "2WBJ", "1H15"].map(async id => parseFasta(await fs.readFile(`${source}/${id}.fasta`, "utf8")))
  );

  // Variable domains only; tags, constant domains, and the 2WBJ peptide-linker fusion are excluded.
  const obAlpha = between(pick(f1ymm, "1YMM_4|"), "SQQGEED", "IFGTGTRLKVLA");
  const obBeta = between(pick(f1ymm, "1YMM_5|"), "GAVVSQHPS", "GPGTRLTVL");
  const a3Alpha = between(pick(f1zgl, "1ZGL_4|"), "GDSVTQME", "IFGSGTRLLVRP");
  const a3Beta = between(pick(f1zgl, "1ZGL_5|"), "VTQTPRYL", "GQGTRLTVV");
  const ob2wbjAlpha = between(pick(f2wbj, "2WBJ_3|"), "QQGEED", "IFGTGTRLKVLPN");
  const ob2wbjBeta = between(pick(f2wbj, "2WBJ_4|"), "AVVSQHPS", "GPGTRLTVT");
  const hyAlpha = "SQQGEEDPQALSIQEGENATMNCSYKTSINNLQWYRQNSGRGLVHLILIRSNEREKHSGRLRVTLDTSKKSSSLLITASRAADTASYFCATDSGGSYIPTFGRGTSLIVHPY";
  const hyBeta = "SAVISQKPSRDICQRGTSLTIQCQVDSQVTMIFWYRQQPGQSLTLIATANQGSEATYESGFVIDKFPISRPNLTFSTLTVSNMSPEDSSIYLCSAWPSGQGTYGYTFGSGTRLTVV";

  // Match the server's provided class-II FASTA example: alpha1 (83 aa) and beta1 (91 aa) binding domains.
  const mhc2 = (alpha, beta) => [mhcBindingDomain(alpha, "IKEEHV", 83), mhcBindingDomain(beta, "TRPRFL", 91)];
  const drb1 = mhc2(pick(f1ymm, "1YMM_1|"), pick(f1ymm, "1YMM_2|"));
  const drb5 = mhc2(pick(f1h15, "1H15_1|"), pick(f1h15, "1H15_2|"));
  const drb5From1zgl = mhc2(pick(f1zgl, "1ZGL_1|"), pick(f1zgl, "1ZGL_2|"));
  const drb1From2wbj = mhc2(pick(f2wbj, "2WBJ_1|"), pick(f2wbj, "2WBJ_2|"));

  const submissions = [
    ["01_CAL_1YMM", "calibration", "DRB1*15:01", "1YMM", "Known OB.1A12–MBP–DRB1 structure; exclude 1YMM in Advanced options.", [obAlpha, obBeta, "ENPVVHFFKNIVTPR", ...drb1]],
    ["02_CAL_1ZGL", "calibration", "DRB5*01:01", "1ZGL", "Known 3A6–MBP–DRB5 structure; 11-mer V-HFFKNIVTP-R; exclude 1ZGL in Advanced options.", [a3Alpha, a3Beta, "VHFFKNIVTPR", ...drb5From1zgl]],
    ["03_CAL_2WBJ", "calibration", "DRB1*15:01", "2WBJ", "Known OB.1A12–ENGA–DRB1 structure; 11-mer D-FARVHFISA-L centered on the deposited 14-mer; exclude 2WBJ in Advanced options.", [ob2wbjAlpha, ob2wbjBeta, "DFARVHFISAL", ...drb1From2wbj]],
    ["04_HY_MBP_DRB1", "Hy.2E11 test", "DRB1*15:01", "1YMM", "Hy.2E11 functional positive, MBP 11-mer V-HFFKNIVTP-R / DRB1*15:01; exclude 1YMM, 1ZGL, 2WBJ, and 1H15.", [hyAlpha, hyBeta, "VHFFKNIVTPR", ...drb1]],
    ["05_HY_BALF5_FULL15", "Hy.2E11 test", "DRB5*01:01", "1H15", "Hy.2E11 functional positive, BALF5 11-mer G-VYHFVKKHV-H / DRB5*01:01; exclude 1H15, 1YMM, 1ZGL, and 2WBJ.", [hyAlpha, hyBeta, "GVYHFVKKHVH", ...drb5]],
    ["06_HY_BALF5_CORE14", "Hy.2E11 test", "DRB5*01:01", "1H15", "Hy.2E11 structural-sensitivity condition, same BALF5 11-mer G-VYHFVKKHV-H / DRB5*01:01; exclude 1H15, 1YMM, 1ZGL, and 2WBJ.", [hyAlpha, hyBeta, "GVYHFVKKHVH", ...drb5]],
  ];

  const alternative = [
    "03b_CAL_2WBJ_CANONICAL_OB",
    "Canonical-OB calibration",
    "DRB1*15:01",
    "2WBJ",
    "Uses the unmodified OB.1A12 variable domains from 1YMM instead of the 2WBJ peptide-linker fusion construct; retains 2WBJ DRB1 and 11-mer D-FARVHFISA-L. Exclude 2WBJ in Advanced options.",
    [obAlpha, obBeta, "DFARVHFISAL", ...drb1From2wbj],
  ];

  await fs.mkdir(out, {recursive: true});
  for (const [id, , , , , sequences] of submissions) {
    const headers = ["TCRa", "TCRb", "pep", "MHC1", "MHC2"];
    await fs.writeFile(`${out}/${id}.fasta`, fasta(id, headers.map((header, index) => [header, sequences[index]])));
  }
  {
    const [id, , , , , sequences] = alternative;
    const headers = ["TCRa", "TCRb", "pep", "MHC1", "MHC2"];
    await fs.writeFile(`${out}/${id}.fasta`, fasta(id, headers.map((header, index) => [header, sequences[index]])));
  }
  const rows = [
    ["submission", "role", "peptide_submitted", "MHC", "direct_reference_PDB_to_exclude", "interpretation"],
    ...submissions.map(([id, role, mhc, pdb, note, sequences]) => [id, role, sequences[2], mhc, pdb, note]),
    [alternative[0], alternative[1], alternative[5][2], alternative[2], alternative[3], alternative[4]],
  ];
  await fs.writeFile(`${out}/submission_manifest.tsv`, `${rows.map(row => row.join("\t")).join("\n")}\n`);
  await fs.writeFile(`${out}/README.md`, `# TCRmodel2 submission bundle\n\nUse the **TCR-pMHCII** tab at [TCRmodel2](https://tcrmodel.ibbr.umd.edu/). Upload one FASTA file per job. The five records are deliberately ordered: TCR alpha, TCR beta, peptide, MHC alpha, MHC beta. Headers are descriptive only; the server uses record order.\n\nAlthough the current server help describes accepting class-II peptides of at least 9 aa, the failed submissions show that this path is enforcing 11 aa. Every pending FASTA therefore uses exactly 11 residues: a 9-aa core with one N- and one C-terminal flank. The already successful 01_CAL_1YMM file is left unchanged.\n\n## Submit first\n\n1. 02_CAL_1ZGL.fasta — exclude 1ZGL.\n2. 03b_CAL_2WBJ_CANONICAL_OB.fasta — revised 2WBJ calibration using canonical OB.1A12 variable domains; exclude 2WBJ.\n3. 04_HY_MBP_DRB1.fasta — exclude 1YMM, 1ZGL, 2WBJ, 1H15.\n\nOnly interpret the Hy.2E11 jobs if the template-excluded calibrators produce credible structures. No resulting score proves binding or cross-reactivity.\n\nDownload all five models and result JSON for each submitted job. Record model confidence, TCR-pMHC ipTM, I-pLDDT, template list, and the unusual-docking warning, if any.\n`);
  console.log(`Wrote ${submissions.length + 1} FASTA submissions to ${out}`);
}

await main();
