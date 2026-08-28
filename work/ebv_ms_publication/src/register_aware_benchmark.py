"""Pure helpers for a register-aware, score-blind pMHC benchmark.

The functions in this module organize computational hypotheses and matched
decoy covariates. They do not infer peptide presentation, TCR recognition,
cross-reactivity, T-cell activation, or multiple-sclerosis mechanism.
"""

from __future__ import annotations

from typing import Any

from premeeting_rigor import eligible_decoys
from register_aware_diagnostic import same_register_alignment_count


def is_assessable_same_register_pair(
    alignment: list[tuple[int, str, int, str]],
    ebv_core_start: int,
    human_core_start: int,
) -> bool:
    """Return whether an existing local alignment includes a shared P1--P9 index."""
    return bool(alignment) and same_register_alignment_count(
        alignment, ebv_core_start, human_core_start
    ) > 0


def strict_eligible_decoys(
    target: dict[str, Any], candidates: list[dict[str, Any]], limit: int = 5
) -> tuple[list[dict[str, object]], int]:
    """Return strict score-blind decoys under the project length/bin rule."""
    return eligible_decoys(target, candidates, limit=limit)
