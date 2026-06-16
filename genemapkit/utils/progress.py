"""Shared progress-bar helper for long-running GeneMapKit operations."""

from __future__ import annotations

from typing import Iterable, Optional, TypeVar

from tqdm.auto import tqdm


T = TypeVar("T")


def track(
    iterable: Iterable[T],
    *,
    enabled: bool = False,
    description: str,
    total: Optional[int] = None,
    unit: str = "items",
) -> Iterable[T]:
    """Wrap an iterable in a stable stderr progress bar when requested."""
    if not enabled:
        return iterable
    return tqdm(
        iterable,
        desc=description,
        total=total,
        unit=unit,
        dynamic_ncols=True,
        mininterval=0.5,
        leave=True,
    )
