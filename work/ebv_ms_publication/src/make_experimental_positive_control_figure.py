"""Create a compact SVG figure from experimental PDB 1H15 and 1BX2."""
from __future__ import annotations

import html
import sys
from pathlib import Path

import numpy as np

from triage_colabfold_pmhc import ca_coordinates, kabsch, parse_pdb, sequence

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "experimental_positive_control"


def circle(x: float, y: float, color: str, label: str) -> str:
    return (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{color}"/>'
            f'<text x="{x + 8:.1f}" y="{y + 4:.1f}" class="res">{html.escape(label)}</text>')


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: make_experimental_positive_control_figure.py 1H15.pdb 1BX2.pdb")
    ebv, mbp = parse_pdb(Path(sys.argv[1])), parse_pdb(Path(sys.argv[2]))
    ebv_frame = np.vstack((ca_coordinates(ebv['A'][:85]), ca_coordinates(ebv['B'][:85])))
    mbp_frame = np.vstack((ca_coordinates(mbp['A'][:85]), ca_coordinates(mbp['B'][:85])))
    rotation, translation, _ = kabsch(ebv_frame, mbp_frame)
    e_seq, m_seq = sequence(ebv['C']), sequence(mbp['C'])
    e_all = ca_coordinates(ebv['C']) @ rotation + translation
    m_all = ca_coordinates(mbp['C'])
    e0, m0 = e_seq.index('YHFVKKH'), m_seq.index('VHFFKNI')
    e, m = e_all[e0:e0 + 7], m_all[m0:m0 + 7]
    # Principal-component projection makes a deterministic two-dimensional
    # depiction; it is a figure aid, not an additional structural metric.
    combined = np.vstack((e, m)); centered = combined - combined.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    p = centered @ vt[:2].T
    lo, hi = p.min(axis=0), p.max(axis=0)
    xy = 315 + 260 * (p - lo) / np.maximum(hi - lo, 1e-6)
    e_xy, m_xy = xy[:7], xy[7:]
    rmsd = float(np.sqrt(np.mean(np.sum((e - m) ** 2, axis=1))))
    svg = ['''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="660" viewBox="0 0 1100 660">
<style>.title{font:700 27px Arial,sans-serif;fill:#172033}.sub{font:16px Arial,sans-serif;fill:#3e4b63}.head{font:700 19px Arial,sans-serif;fill:#172033}.body{font:16px Arial,sans-serif;fill:#25324a}.small{font:14px Arial,sans-serif;fill:#3e4b63}.res{font:13px Arial,sans-serif;fill:#172033}</style>
<rect width="1100" height="660" fill="#ffffff"/><text x="45" y="52" class="title">Experimental positive control: EBV BALF5 and MBP pMHC surfaces</text>
<text x="45" y="80" class="sub">Crystal structures, not predicted complexes. Geometric similarity does not itself demonstrate TCR binding.</text>
<line x1="45" y1="110" x2="1055" y2="110" stroke="#cfd7e6"/>
<text x="45" y="150" class="head">Established cross-reactive system</text>
<rect x="45" y="175" width="430" height="185" rx="12" fill="#eef6ff" stroke="#8ebfe8"/>
<text x="68" y="210" class="body">MBP(85–99) → HLA-DRB1*15:01 (DR2b)</text>
<text x="68" y="242" class="body">EBV BALF5(627–641) → HLA-DRB5*01:01 (DR2a)</text>
<text x="68" y="290" class="body">Both ligands are recognized by Hy.2E11</text>
<text x="68" y="318" class="small">Allele-pair switching is essential biology.</text>
<text x="45" y="405" class="head">Experimental coordinates used</text>
<text x="45" y="437" class="body">1H15: DR2a–BALF5</text><text x="45" y="465" class="body">1BX2: DR2b–MBP</text>
<text x="45" y="515" class="small">Core metric: peptide Cα RMSD after HLA-groove fit</text>
<text x="45" y="557" class="title">''' + f'{rmsd:.3f} Å' + '''</text>
<text x="550" y="150" class="head">Aligned seven-position peptide core</text>
<text x="550" y="178" class="small">2D projection after experimental HLA peptide-binding-platform alignment</text>
<rect x="535" y="200" width="520" height="385" rx="12" fill="#fbfcff" stroke="#cfd7e6"/>
<line x1="570" y1="230" x2="590" y2="230" stroke="#1177cc" stroke-width="3"/><text x="598" y="235" class="small">BALF5, 1H15</text>
<line x1="735" y1="230" x2="755" y2="230" stroke="#e66a2c" stroke-width="3"/><text x="763" y="235" class="small">MBP, 1BX2</text>''']
    for i, ((ex, ey), (mx, my)) in enumerate(zip(e_xy, m_xy), 1):
        svg.append(f'<line x1="{ex:.1f}" y1="{ey:.1f}" x2="{mx:.1f}" y2="{my:.1f}" stroke="#b8c4d8" stroke-width="1.4"/>')
        svg.append(circle(ex, ey, '#1177cc', f'{i}:{e_seq[e0+i-1]}'))
        svg.append(circle(mx, my, '#e66a2c', f'{i}:{m_seq[m0+i-1]}'))
    svg.append('''<text x="560" y="560" class="small">Each grey connector joins structurally aligned core positions.</text>
<text x="45" y="625" class="small">Sources: Lang et al., Nat. Immunol. 2002; RCSB PDB 1H15 and 1BX2.</text></svg>''')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'figure_1_experimental_positive_control.svg').write_text('\n'.join(svg))


if __name__ == '__main__':
    main()
