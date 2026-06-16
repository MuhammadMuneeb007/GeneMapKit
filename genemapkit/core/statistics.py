"""Statistical helpers for benchmark reports."""

from __future__ import annotations

import math
from statistics import NormalDist
from typing import Dict, Optional


def wilson_interval(
    successes: int,
    total: int,
    *,
    confidence_level: float = 0.95,
) -> tuple[Optional[float], Optional[float]]:
    """Return a Wilson score confidence interval as percentages."""
    if total <= 0:
        return None, None
    if successes < 0 or successes > total:
        raise ValueError("successes must be between zero and total")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between zero and one")
    z = NormalDist().inv_cdf(0.5 + confidence_level / 2.0)
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    margin = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / total
            + z * z / (4.0 * total * total)
        )
        / denominator
    )
    return (
        round(100.0 * max(0.0, center - margin), 4),
        round(100.0 * min(1.0, center + margin), 4),
    )


def interval_fields(
    prefix: str,
    successes: int,
    total: int,
    *,
    confidence_level: float = 0.95,
) -> Dict[str, object]:
    """Return explicit numerator, denominator, and Wilson interval report fields."""
    low, high = wilson_interval(
        successes, total, confidence_level=confidence_level
    )
    return {
        f"{prefix}_ci_method": "Wilson score",
        f"{prefix}_confidence_level": confidence_level,
        f"{prefix}_numerator": successes,
        f"{prefix}_denominator": total,
        f"{prefix}_ci_low_pct": low,
        f"{prefix}_ci_high_pct": high,
    }
