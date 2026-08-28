# Brian Pierce call prep

Meeting: Tuesday, August 4, 2026, 2:00 PM CT  
Zoom: https://umd.zoom.us/my/pierce  
Thread: Brian Pierce, University of Maryland, TCR-pMHC modeling

## Why this call matters

Dr. Pierce is directly in the TCR-pMHC structural modeling lane. This is the best current chance to decide whether the project should:

1. keep TCR-level modeling as rejected/future work, or
2. add a calibrated, receptor-level follow-up using appropriate tools and controls.

The scientific posture should stay humble and precise: the current project is a calibrated pMHC prioritization workflow, not proof of TCR binding.

## Brian's explicit question

He asked:

> Just to check, you've been checking out the relevant literature as background or reference to help guide your work, and if so, which studies have you been using for reference on the TCR complex structural modeling so far?

## Best short answer to give him

So far, I have mainly used structural positive controls rather than trusting receptor prediction scores directly. The key TCR-positive control in my project is Ob.1A12 bound to HLA-DRB1*15:01/MBP in PDB 1YMM. I also used the experimental pMHC-only structures 1H15 and 1BX2 as the BALF5/MBP molecular-mimicry anchor.

For TCR-complex modeling, I initially tried exploratory five-chain ColabFold-Multimer, but I checked it against the known 1YMM positive-control orientation and it failed. The predicted TCR chains looked folded, but the TCR pose did not recover the experimental orientation, so I excluded those results as evidence of binding.

Since then, I have been shifting toward the actual TCR-pMHC modeling literature: TCRmodel2 from your lab, Bradley's AlphaFold/TCRdock work, and docking/benchmark papers that emphasize controls, docking geometry, and limits of AlphaFold-style methods. My main goal for the call is to learn what benchmark would be defensible before making any receptor-level claim.

## Literature / methods currently relevant

### 1. Ob.1A12 experimental positive control

Use in this project:
- Biological receptor-level anchor.
- PDB 1YMM: Ob.1A12 TCR bound to HLA-DR2b / HLA-DRB1*15:01 with MBP peptide.
- Current conclusion: our five-chain ColabFold model failed to recover the 1YMM TCR orientation, so ternary predictions are not evidence.

Why it matters:
- This is the exact positive-control orientation any receptor-level method should recover before being trusted on EBV candidates.

Relevant source:
- RCSB 1YMM: https://www.rcsb.org/structure/1YMM

### 2. TCRmodel2

Use in next step:
- Candidate replacement for generic ColabFold-Multimer.
- It is specifically designed for TCR and TCR-pMHC modeling.
- It gives confidence scores and is benchmarked against other methods.

Questions for Brian:
- Is TCRmodel2 appropriate for HLA class II / DRB1*15:01 cases?
- For class II, TCRmodel2 expects an 11-aa peptide input: a 9-aa core plus one residue on each side. How should we define the 11-aa core for MBP and EBV candidates?
- Can TCRmodel2 handle the unusual Ob.1A12 docking mode, or does that make it a stress test?

Relevant sources:
- TCRmodel2 paper / PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10320165/
- TCRmodel2 PubMed: https://pubmed.ncbi.nlm.nih.gov/37140040/
- TCRmodel2 web server: https://tcrmodel.ibbr.umd.edu/

### 3. Bradley / TCRdock / specialized AlphaFold-style TCR-pMHC modeling

Use in next step:
- Shows that specialized TCR-pMHC structure prediction can sometimes discriminate correct from incorrect pMHC targets.
- More appropriate than treating generic ColabFold-Multimer ipTM as a binding score.

Questions for Brian:
- Would TCRdock-style docking geometry be a better benchmark than raw contact counts?
- Should we compare predicted docking geometries to the 1YMM experimental geometry?
- Is this approach reliable enough for a student computational paper, or only hypothesis-generating?

Relevant sources:
- eLife paper: https://elifesciences.org/articles/82813
- PMC full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC9859041/
- TCRdock GitHub: https://github.com/phbradley/TCRdock

### 4. Information-driven / restrained docking

Use in next step:
- Alternative if template-guided or restrained docking is more defensible than unconstrained five-chain prediction.
- Literature suggests docking improves when using interface information/constraints, with HADDOCK often performing well in benchmark settings.

Questions for Brian:
- If using the 1YMM geometry as a template, what counts as legitimate restraint vs circular overfitting?
- Would restrained docking be defensible if evaluated against decoys and negative controls?
- What independent score should be used besides contact counts?

Relevant sources:
- Information-driven docking paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8219952/
- UCL summary: https://discovery.ucl.ac.uk/id/eprint/10130964/

### 5. AlphaFold-Multimer limitations

Use in current claim boundary:
- Supports why the project correctly rejected generic five-chain ColabFold as binding evidence.
- In this project, ColabFold's global scores looked good but failed receptor-pose calibration.

Questions for Brian:
- Is the failure mode we saw typical for TCR-pMHC?
- Are pLDDT/pTM/ipTM useful at all in this context, or should they be treated only as weak secondary QC?
- Is the safest publication framing to leave receptor modeling out unless a specialized benchmark passes?

Relevant source:
- AlphaFold-Multimer benchmark PDF noting low near-native success on a small TCR-pMHC benchmark: https://drum.lib.umd.edu/bitstreams/48b7f9d4-d0c9-4d5a-b760-ff33a613af11/download

## Current project evidence to show Brian

### What worked

- Experimental pMHC anchor:
  - 1H15: EBV BALF5 peptide on HLA-DRB5*01:01 / DR2a.
  - 1BX2: MBP peptide on HLA-DRB1*15:01 / DR2b.
- 86 modeled rank-1 pMHC structures passed peptide/chain-layout QA.
- 32 EBV-CNS candidate pairs ranked by local pMHC geometry and physicochemical similarity.
- Classic BALF5-MBP positive controls were strongly recovered in external validation.

### What failed

- Five-chain ColabFold ternary test using Ob.1A12 failed calibration:
  - HLA fit was good.
  - TCR chain pLDDT was decent.
  - But the TCR pose did not recover the experimental 1YMM orientation.
  - Median HLA-aligned TCR alpha RMSD: 10.227 Å.
  - Median HLA-aligned TCR beta RMSD: 18.040 Å.
- Therefore, no TCR-recognition claim is currently made.

## Questions to ask Brian

1. Is it scientifically reasonable to keep the paper as pMHC prioritization only, with receptor modeling as future work?
2. If we do add receptor modeling, should the starting point be TCRmodel2, TCRdock, template-guided docking, or something else?
3. What exact positive and negative controls would make the receptor-level branch defensible?
4. For HLA class II, how should the 9-aa binding core plus flanking residues be selected for TCRmodel2-style input?
5. What output metrics should be reported?
   - TCR pose RMSD to 1YMM?
   - docking geometry?
   - CDR3/interface pLDDT?
   - peptide-contact patterns?
   - decoy enrichment?
6. Does Ob.1A12's unusual autoimmune docking topology make it a bad test case for generic tools, or an ideal stress test?
7. Would he recommend any specific benchmark datasets or negative controls?
8. What claim boundary would he be comfortable seeing in a student paper?

## Draft reply to Brian

Dear Dr. Pierce,

Thank you — Tuesday, August 4 at 2 PM CT works for me, and I’ll use the Zoom link you sent.

On the literature side, I’ve mainly been using the experimental structures as controls so far: PDB 1YMM for Ob.1A12 bound to HLA-DRB1*15:01/MBP, plus 1H15 and 1BX2 for the BALF5/MBP pMHC molecular-mimicry anchor. I initially tried exploratory five-chain ColabFold-Multimer, but after benchmarking it against the known 1YMM orientation, it failed to recover the experimental TCR pose, so I excluded it as evidence.

Before tomorrow, I’m reviewing TCRmodel2, Bradley’s TCRdock/AlphaFold-based TCR-pMHC work, and benchmark/docking papers on TCR-pMHC complex prediction. My main goal is to understand what controls and metrics would make any receptor-level follow-up defensible, or whether the paper should stop at pMHC prioritization.

Thank you again,
Anish

