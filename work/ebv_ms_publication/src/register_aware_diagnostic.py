"""Register-position diagnostics for the predeclared modeled-pMHC shortlist.

These functions report whether residues that were already locally aligned in
the original screen occupy the same P1--P9 position under a selected or
possible peptide core.  This is a sensitivity/feasibility artifact, not a new
rank, predictor, binding result, or cross-reactivity claim.
"""

from __future__ import annotations

from itertools import product


def parse_local_alignment_positions(value: str) -> list[tuple[int, str, int, str]]:
    """Parse the geometry-matrix ``positionResidue:positionResidue`` entries."""
    if not value:
        return []
    rows = []
    for item in value.split(";"):
        left, right = item.split(":", maxsplit=1)
        rows.append((int(left[:-1]), left[-1], int(right[:-1]), right[-1]))
    return rows


def same_register_alignment_count(
    alignment: list[tuple[int, str, int, str]], ebv_core_start: int, human_core_start: int
) -> int:
    """Count existing local alignments whose residues share a P1--P9 index."""
    return sum(
        ebv_position - ebv_core_start == human_position - human_core_start
        and 0 <= ebv_position - ebv_core_start < 9
        and 0 <= human_position - human_core_start < 9
        for ebv_position, _ebv_residue, human_position, _human_residue in alignment
    )


def window_pair_sensitivity(
    ebv_peptide: str, human_peptide: str, alignment: list[tuple[int, str, int, str]]
) -> list[dict[str, int]]:
    """Retain every pair of manifest-contained 9-mer registers for sensitivity."""
    ebv_starts = range(1, max(0, len(ebv_peptide) - 8) + 1)
    human_starts = range(1, max(0, len(human_peptide) - 8) + 1)
    return [
        {
            "ebv_core_start_1_based": ebv_start,
            "human_core_start_1_based": human_start,
            "same_register_alignment_count": same_register_alignment_count(
                alignment, ebv_start, human_start
            ),
        }
        for ebv_start, human_start in product(ebv_starts, human_starts)
    ]
