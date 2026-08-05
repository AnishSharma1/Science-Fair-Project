import fs from "node:fs/promises";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const outDir = "/Users/anishsharma/Documents/New project/outputs/ebv_ms_benchmark_manifest";
const outPath = `${outDir}/EBV_MS_TCR_pMHC_Benchmark_Manifest.xlsx`;
const navy = "#17365D", blue = "#2E74B5", paleBlue = "#E8EEF5", paleYellow = "#FFF4CC", paleRed = "#FCE8E6", paleGreen = "#E6F4EA", border = "#C7D5E5";

function fmt(range, {fill, bold, color, size, wrap = true, align = "left"} = {}) {
  range.format.wrapText = wrap;
  range.format.horizontalAlignment = align;
  range.format.verticalAlignment = "center";
  if (fill) range.format.fill = fill;
  if (bold || color || size) range.format.font = {bold, color, size};
}
function box(range, preset = "all") { range.format.borders = {preset, style: "thin", color: border}; }
function widths(sheet, spec) { for (const [c, width] of Object.entries(spec)) sheet.getRange(`${c}:${c}`).format.columnWidth = width; }
function table(sheet, range, headers, rows, colWidths) {
  sheet.getRange(range).values = [headers, ...rows];
  const endRow = rows.length + 1;
  const endCol = String.fromCharCode(64 + headers.length);
  fmt(sheet.getRange(`A1:${endCol}1`), {fill: navy, bold: true, color: "#FFFFFF", align: "center"});
  fmt(sheet.getRange(`A2:${endCol}${endRow}`), {size: 9});
  box(sheet.getRange(`A1:${endCol}${endRow}`));
  sheet.freezePanes.freezeRows(1);
  widths(sheet, colWidths);
  sheet.getRange(`A1:${endCol}${endRow}`).format.autofitRows();
}

async function main() {
  await fs.mkdir(outDir, {recursive: true});
  const wb = Workbook.create();
  const readme = wb.worksheets.add("Read Me");
  const pmhc = wb.worksheets.add("pMHC References");
  const tcr = wb.worksheets.add("TCR References");
  const pairs = wb.worksheets.add("Benchmark Matrix");
  const decoys = wb.worksheets.add("Decoy Panel");
  const inputs = wb.worksheets.add("Frozen Model Inputs");
  for (const s of [readme, pmhc, tcr, pairs, decoys, inputs]) s.showGridLines = false;

  readme.mergeCells("A1:B1");
  readme.getRange("A1").values = [["EBV–MS TCR–pMHC Benchmark Manifest"]];
  fmt(readme.getRange("A1:B1"), {fill: navy, bold: true, color: "#FFFFFF", size: 16});
  readme.getRange("A1:B1").format.rowHeight = 30;
  readme.mergeCells("A2:B2");
  readme.getRange("A2").values = [["Frozen evidence map for staged pMHC and ternary-complex evaluation • updated 2026-08-04"]];
  fmt(readme.getRange("A2:B2"), {fill: paleBlue, color: navy, size: 10});
  readme.getRange("A4:B9").values = [
    ["Status", "Meaning"],
    ["Ready: pMHC", "Exact peptide, allele and experimental pMHC structure are available."],
    ["Ready: ternary calibration", "Published full TCR–pMHC structure permits geometric evaluation."],
    ["Functional only", "T-cell activation is published, but allele/structure is not yet sufficient for an allele-specific model run."],
    ["Decoy / unknown", "Pre-specified noncognate-looking pairing. It is not a demonstrated non-binder."],
    ["True functional negative", "Requires an explicitly absent response in a matched TCR–peptide–HLA assay. None are claimed here."],
  ];
  fmt(readme.getRange("A4:B4"), {fill: paleBlue, bold: true, color: navy});
  fmt(readme.getRange("A5:A9"), {bold: true, color: navy});
  box(readme.getRange("A4:B9"));
  readme.getRange("A8:B8").format.fill = paleYellow;
  readme.getRange("A9:B9").format.fill = paleRed;
  readme.mergeCells("A11:B11");
  readme.getRange("A11").values = [["Resolved inputs and interpretation"]];
  fmt(readme.getRange("A11:B11"), {fill: navy, bold: true, color: "#FFFFFF"});
  readme.mergeCells("A12:B15");
  readme.getRange("A12").values = [["Hy.2E11 alpha/beta variable-region amino-acid sequences are captured from US20220364057A1 Table 1 and its original junctions agree with the 1997 clone report. Use full BALF5(627–641) TGGVYHFVKKHVHES for biological reproduction. PDB 1H15 contains the resolved/deposited 14-residue core GGVYHFVKKHVHES (628–641); when comparing a model to 1H15, score only the shared core plus MHC. MBP(85–99)/DRB1*15:01 and BALF5/DRB5*01:01 are a cross-allele biological positive, not a same-allele structural pair."]];
  fmt(readme.getRange("A12:B15"), {fill: paleYellow, size: 10});
  box(readme.getRange("A11:B15"), "outside");
  readme.mergeCells("A17:B17");
  readme.getRange("A17").values = [["Required before an efficacy claim"]];
  fmt(readme.getRange("A17:B17"), {fill: navy, bold: true, color: "#FFFFFF"});
  readme.getRange("A18:B22").values = [
    ["1", "Submit the frozen input records identically to both methods; log chain boundaries, template use and version."],
    ["2", "Evaluate pMHC geometry separately from TCR–pMHC geometry. Confidence is not binding or cross-reactivity evidence."],
    ["3", "Use 1YMM, 1ZGL and 2WBJ as exact ternary-structure calibration cases before interpreting Hy.2E11 models."],
    ["4", "Use the Decoy Panel only as a selectivity stress test. Do not report its entries as true negatives."],
    ["5", "A 10–20 member matched functional-negative set remains an experimental-data requirement; literature-free pairings do not meet it."],
  ];
  fmt(readme.getRange("A18:A22"), {fill: paleBlue, bold: true, color: navy, align: "center"});
  fmt(readme.getRange("B18:B22"), {size: 10});
  box(readme.getRange("A18:B22"));
  widths(readme, {A: 16, B: 62});
  readme.getRange("A1:B22").format.autofitRows();

  const pmhcHeaders = ["pMHC_ID", "Role", "Peptide source", "Residues / record", "Biological peptide sequence", "Exact structural comparison sequence", "Class II alpha", "Class II beta", "PDB", "Task status", "Evidence / rule", "Source URL"];
  const pmhcRows = [
    ["PMHC_MBP_DRB1", "Hy.2E11 self component; calibration", "Human MBP", "85–99", "ENPVVHFFKNIVTPR", "ENPVVHFFKNIVTPR", "DRA*01:01", "DRB1*15:01", "1BX2; 1YMM", "Ready: pMHC", "1BX2 is pMHC; 1YMM is OB.1A12 ternary with same pMHC.", "https://www.rcsb.org/structure/1BX2"],
    ["PMHC_BALF5_DRB5", "Hy.2E11 viral component", "EBV BALF5 DNA polymerase", "Biological: 627–641; 1H15 deposit: 628–641", "TGGVYHFVKKHVHES", "GGVYHFVKKHVHES", "DRA*01:01", "DRB5*01:01", "1H15", "Ready: pMHC", "Model the full 15-mer; compare the shared 14-mer core to 1H15.", "https://www.rcsb.org/structure/1H15"],
    ["PMHC_FLU_HA_DR2", "Hy.2E11 functional mimic", "Influenza A hemagglutinin", "reported peptide", "YRNLVWFIKKNTRYP", "", "Not allele-resolved", "DR2 family", "None", "Functional only", "Activation reported in the 1995 panel; no allele-specific structural run until restriction is resolved.", "https://pubmed.ncbi.nlm.nih.gov/7534214/"],
    ["PMHC_REO_SIGMA2_DR2", "Hy.2E11 functional mimic", "Reovirus type 3 sigma 2", "reported peptide", "MARAAFLFKTVGFGG", "", "Not allele-resolved", "DR2 family", "None", "Functional only", "Activation reported in the 1995 panel; no allele-specific structural run until restriction is resolved.", "https://pubmed.ncbi.nlm.nih.gov/7534214/"],
    ["PMHC_MBP_DRB5", "3A6 calibration", "Human MBP", "89–101", "VHFFKNIVTPRTP", "VHFFKNIVTPRTP", "DRA*01:01", "DRB5*01:01", "1ZGL", "Ready: ternary calibration", "Experimental 3A6/MBP/DR2a ternary complex.", "https://www.rcsb.org/structure/1ZGL"],
    ["PMHC_ENGA_DRB1", "OB.1A12 calibration", "Bacterial enolase (ENGA)", "core: FARVHFISALHG; deposited construct includes N-terminal MD", "FARVHFISALHG", "MDFARVHFISALHG", "DRA*01:01", "DRB1*15:01", "2WBJ", "Ready: ternary calibration", "Use deposited MDFARVHFISALHG for strict 2WBJ recovery; do not include its artificial TCR-fusion linker.", "https://www.rcsb.org/structure/2WBJ"],
  ];
  table(pmhc, `A1:L${pmhcRows.length + 1}`, pmhcHeaders, pmhcRows, {A: 20, B: 25, C: 28, D: 34, E: 24, F: 28, G: 16, H: 17, I: 15, J: 27, K: 52, L: 42});
  pmhc.getRange("J2:J3").format.fill = paleGreen;
  pmhc.getRange("J4:J5").format.fill = paleYellow;
  pmhc.getRange("J6:J7").format.fill = paleGreen;

  const alpha = "SQQGEEDPQALSIQEGENATMNCSYKTSINNLQWYRQNSGRGLVHLILIRSNEREKHSGRLRVTLDTSKKSSSLLITASRAADTASYFCATDSGGSYIPTFGRGTSLIVHPY";
  const beta = "SAVISQKPSRDICQRGTSLTIQCQVDSQVTMIFWYRQQPGQSLTLIATANQGSEATYESGFVIDKFPISRPNLTFSTLTVSNMSPEDSSIYLCSAWPSGQGTYGYTFGSGTRLTVV";
  const tcrHeaders = ["TCR_ID", "Clone", "Use", "Alpha variable AA sequence", "Beta variable AA sequence", "Sequence scope", "Sequence evidence", "Ternary reference", "Run status", "Source URL"];
  const tcrRows = [
    ["TCR_HY_2E11", "Hy.2E11", "Functional cross-reactivity reproduction", alpha, beta, "Variable domains only; add canonical human constants only if server input requires full chains.", "Exact V-region rows: US20220364057A1 Table 1; original 1997 report corroborates junctions.", "None", "Ready for variable-domain modeling", "https://patents.google.com/patent/US20220364057A1/en"],
    ["TCR_OB_1A12", "OB.1A12", "Structural calibration", "Extract deposited chain from 1YMM / 2WBJ", "Extract deposited chain from 1YMM / 2WBJ", "Use deposited construct for RMSD calibration.", "Experimental ternary entries.", "1YMM; 2WBJ", "Ready: reference structure", "https://www.rcsb.org/structure/1YMM"],
    ["TCR_3A6", "3A6", "Structural calibration", "Extract deposited chain from 1ZGL", "Extract deposited chain from 1ZGL", "Use deposited construct for RMSD calibration.", "Experimental ternary entry.", "1ZGL", "Ready: reference structure", "https://www.rcsb.org/structure/1ZGL"],
  ];
  table(tcr, `A1:J${tcrRows.length + 1}`, tcrHeaders, tcrRows, {A: 18, B: 16, C: 27, D: 54, E: 54, F: 44, G: 43, H: 20, I: 29, J: 46});
  tcr.getRange("I2:I4").format.fill = paleGreen;

  const pairHeaders = ["Benchmark_ID", "Benchmark class", "TCR", "pMHC A", "pMHC B / complex", "Expected result", "Evidence tier", "Primary evaluation", "Decision rule", "Source URL"];
  const pairRows = [
    ["HY_BIO_001", "Functional positive", "TCR_HY_2E11", "PMHC_MBP_DRB1", "PMHC_BALF5_DRB5", "Same TCR recognizes self and EBV pMHC", "Direct functional + two pMHC structures", "Separate pMHC quality; ternary interface is prediction-only", "Do not compute ternary RMSD; assess concordance across methods and decoys.", "https://pubmed.ncbi.nlm.nih.gov/12244309/"],
    ["HY_BIO_002", "Functional positive", "TCR_HY_2E11", "PMHC_MBP_DRB1", "PMHC_FLU_HA_DR2", "Same TCR recognizes MBP and influenza peptide", "Direct functional panel", "Functional reproduction only", "Hold structural run until HLA restriction is resolved.", "https://pubmed.ncbi.nlm.nih.gov/7534214/"],
    ["HY_BIO_003", "Functional positive", "TCR_HY_2E11", "PMHC_MBP_DRB1", "PMHC_REO_SIGMA2_DR2", "Same TCR recognizes MBP and reovirus peptide", "Direct functional panel", "Functional reproduction only", "Hold structural run until HLA restriction is resolved.", "https://pubmed.ncbi.nlm.nih.gov/7534214/"],
    ["CAL_1YMM", "Exact ternary calibration", "TCR_OB_1A12", "PMHC_MBP_DRB1", "1YMM", "Recover experimental ternary geometry", "Experimental X-ray, 3.50 Å", "TCR–pMHC RMSD / interface", "Calibration only; document chain mapping and comparison region.", "https://www.rcsb.org/structure/1YMM"],
    ["CAL_1ZGL", "Exact ternary calibration", "TCR_3A6", "PMHC_MBP_DRB5", "1ZGL", "Recover experimental ternary geometry", "Experimental X-ray, 2.80 Å", "TCR–pMHC RMSD / interface", "Calibration only; document chain mapping and comparison region.", "https://www.rcsb.org/structure/1ZGL"],
    ["CAL_2WBJ", "Exact ternary calibration", "TCR_OB_1A12", "PMHC_ENGA_DRB1", "2WBJ", "Recover microbial-mimic ternary geometry", "Experimental X-ray, 3.00 Å", "TCR–pMHC RMSD / interface", "Calibration only; document chain mapping and comparison region.", "https://www.rcsb.org/structure/2WBJ"],
  ];
  table(pairs, `A1:J${pairRows.length + 1}`, pairHeaders, pairRows, {A: 18, B: 25, C: 18, D: 22, E: 24, F: 38, G: 31, H: 34, I: 46, J: 44});
  pairs.getRange("B2:B4").format.fill = paleGreen;
  pairs.getRange("B5:B7").format.fill = paleBlue;

  const decoyHeaders = ["Decoy_ID", "TCR", "pMHC", "Why included", "Status", "Allowed analysis", "Not allowed", "Evidence requirement before calling negative"];
  const decoySeeds = [
    ["TCR_HY_2E11", "PMHC_MBP_DRB5"], ["TCR_HY_2E11", "PMHC_ENGA_DRB1"],
    ["TCR_OB_1A12", "PMHC_BALF5_DRB5"], ["TCR_OB_1A12", "PMHC_MBP_DRB5"],
    ["TCR_3A6", "PMHC_MBP_DRB1"], ["TCR_3A6", "PMHC_BALF5_DRB5"], ["TCR_3A6", "PMHC_ENGA_DRB1"],
  ];
  const decoyRows = decoySeeds.map(([clone, target], i) => [
    `DECOY_${String(i + 1).padStart(2, "0")}`, clone, target,
    "Cross-pairing of a published complex with a TCR that is not its documented cognate.", "Decoy / unknown",
    "Selectivity stress test; compare ranking and interface confidence only.", "Do not count as a non-binder, false positive, or specificity estimate.",
    "Matched assay showing no response for this exact TCR, peptide and HLA combination."
  ]);
  table(decoys, `A1:H${decoyRows.length + 1}`, decoyHeaders, decoyRows, {A: 16, B: 20, C: 23, D: 48, E: 20, F: 43, G: 49, H: 53});
  decoys.getRange(`E2:H${decoyRows.length + 1}`).format.fill = paleYellow;

  const inputHeaders = ["Run_ID", "Method", "Goal", "TCR alpha input", "TCR beta input", "pMHC alpha", "pMHC beta", "Peptide input", "Comparison region", "Interpretation gate"];
  const inputRows = [
    ["HY_EBV_FULL15", "AlphaFold 3 / TCRmodel2", "Hy.2E11 biological reproduction", alpha, beta, "DRA*01:01", "DRB5*01:01", "TGGVYHFVKKHVHES", "For 1H15: MHC + peptide residues 2–15 only", "No ternary RMSD claim; compare methods and decoy ranks."],
    ["HY_MBP", "AlphaFold 3 / TCRmodel2", "Hy.2E11 biological reproduction", alpha, beta, "DRA*01:01", "DRB1*15:01", "ENPVVHFFKNIVTPR", "pMHC quality against 1BX2", "No ternary RMSD claim; compare methods and decoy ranks."],
    ["CAL_1YMM", "AlphaFold 3 / TCRmodel2", "Ternary geometry calibration", "Deposited 1YMM TCR alpha", "Deposited 1YMM TCR beta", "Deposited 1YMM DRA", "Deposited 1YMM DRB1", "ENPVVHFFKNIVTPR", "Defined shared atoms after chain mapping", "Report RMSD/interface recovery with chain mapping."],
    ["CAL_1ZGL", "AlphaFold 3 / TCRmodel2", "Ternary geometry calibration", "Deposited 1ZGL TCR alpha", "Deposited 1ZGL TCR beta", "Deposited 1ZGL DRA", "Deposited 1ZGL DRB5", "VHFFKNIVTPRTP", "Defined shared atoms after chain mapping", "Report RMSD/interface recovery with chain mapping."],
    ["CAL_2WBJ", "AlphaFold 3 / TCRmodel2", "Ternary geometry calibration", "Deposited 2WBJ TCR alpha", "Deposited 2WBJ TCR beta, linker removed", "Deposited 2WBJ DRA", "Deposited 2WBJ DRB1", "MDFARVHFISALHG", "Defined shared atoms after chain mapping", "Report RMSD/interface recovery; do not submit the peptide–TCR fusion linker."],
  ];
  table(inputs, `A1:J${inputRows.length + 1}`, inputHeaders, inputRows, {A: 17, B: 25, C: 31, D: 52, E: 52, F: 21, G: 22, H: 26, I: 39, J: 43});
  inputs.getRange("A2:J3").format.fill = paleGreen;
  inputs.getRange("A4:J6").format.fill = paleBlue;

  const check = await wb.inspect({kind: "table", range: "Frozen Model Inputs!A1:J6", include: "values,formulas", tableMaxRows: 8, tableMaxCols: 10});
  console.log(check.ndjson);
  const file = await SpreadsheetFile.exportXlsx(wb);
  await file.save(outPath);
  console.log(outPath);
}

await main();
