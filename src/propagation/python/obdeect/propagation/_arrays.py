"""Argument helpers shared by the propagation Python API (private).

They convert Python arguments to the types of the C++ core without silent changes: integers
are taken with ``operator.index`` (1.7 is rejected, not truncated), photon columns must already
be float32 (no silent cast) and pairs must have exactly two elements.
"""

import operator
from collections.abc import Callable, Mapping

import numpy as np


def broadcast_call(function: Callable, *arguments):
    """Broadcast float64 arguments to one shape and apply function to the flat arrays.

    A function returning one array gives one array of the broadcast shape; a function returning
    several arrays (a tuple or list) gives a tuple of such arrays.
    """
    arrays = np.broadcast_arrays(*(np.asarray(a, dtype=np.float64) for a in arguments))
    shape = arrays[0].shape
    flat = [np.ascontiguousarray(a).reshape(-1) for a in arrays]
    result = function(*flat)
    if isinstance(result, tuple | list):
        return tuple(np.asarray(r).reshape(shape) for r in result)
    return np.asarray(result).reshape(shape)


def integer(value, name: str) -> int:
    """The integer value, rejecting non-integers (floats, strings) instead of truncating them."""
    try:
        return operator.index(value)
    except TypeError:
        raise TypeError(f"{name} must be an integer, not {value!r}") from None


def float_pair(value, name: str) -> tuple[float, float]:
    """Exactly two numbers as (float, float)."""
    if isinstance(value, str | bytes):
        raise TypeError(f"{name} must be a pair of numbers, got {value!r}")
    items = tuple(value)
    if len(items) != 2:
        raise ValueError(f"{name} must have exactly two elements, got {len(items)}")
    return float(items[0]), float(items[1])


def float32_columns(columns: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Contiguous float32 photon columns; a column of another dtype is rejected (no cast)."""
    out = {}
    for key, column in columns.items():
        array = np.asarray(column)
        if array.dtype != np.float32:
            raise TypeError(f"photon column {key!r} must be float32, got {array.dtype}")
        out[key] = np.ascontiguousarray(array)
    return out
