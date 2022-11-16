"""
This type stub file was generated by pyright.
"""

import numpy as np
from astropy.utils.compat import NUMPY_LT_1_19, NUMPY_LT_1_20, NUMPY_LT_1_23

"""Helpers for letting numpy functions interact with Masked arrays.

The module supplies helper routines for numpy functions that propagate
masks appropriately., for use in the ``__array_function__``
implementation of `~astropy.utils.masked.MaskedNDArray`.  They are not
very useful on their own, but the ones with docstrings are included in
the documentation so that there is a place to find out how the mask is
interpreted.

"""
__all__ = ['MASKED_SAFE_FUNCTIONS', 'APPLY_TO_BOTH_FUNCTIONS', 'DISPATCHED_FUNCTIONS', 'UNSUPPORTED_FUNCTIONS']
MASKED_SAFE_FUNCTIONS = ...
APPLY_TO_BOTH_FUNCTIONS = ...
DISPATCHED_FUNCTIONS = ...
UNSUPPORTED_FUNCTIONS = ...
IGNORED_FUNCTIONS = ...
if NUMPY_LT_1_20:
    ...
if NUMPY_LT_1_23:
    ...
apply_to_both = ...
dispatched_function = ...
@dispatched_function
def datetime_as_string(arr, *args, **kwargs): # -> tuple[str_, Unknown, None]:
    ...

@dispatched_function
def sinc(x): # -> tuple[floating[Any], Unknown, None]:
    ...

@dispatched_function
def iscomplex(x): # -> tuple[bool_, Unknown, None]:
    ...

@dispatched_function
def unwrap(p, *args, **kwargs): # -> tuple[NDArray[floating[Any]], Unknown, None]:
    ...

@dispatched_function
def nan_to_num(x, copy=..., nan=..., posinf=..., neginf=...): # -> tuple[Unknown, Unknown, None]:
    ...

@apply_to_both(np.copy, np.asfarray, np.resize, np.moveaxis, np.rollaxis, np.roll)
def masked_a_helper(a, *args, **kwargs): # -> tuple[tuple[Unknown | Any | None, ...], tuple[Any | NDArray[Any], ...], dict[str, Unknown], None]:
    ...

@apply_to_both(np.flip, np.flipud, np.fliplr, np.rot90, np.triu, np.tril)
def masked_m_helper(m, *args, **kwargs): # -> tuple[tuple[Unknown | Any | None, ...], tuple[Any | NDArray[Any], ...], dict[str, Unknown], None]:
    ...

@apply_to_both(np.diag, np.diagflat)
def masked_v_helper(v, *args, **kwargs): # -> tuple[tuple[Unknown | Any | None, ...], tuple[Any | NDArray[Any], ...], dict[str, Unknown], None]:
    ...

@apply_to_both(np.delete)
def masked_arr_helper(array, *args, **kwargs): # -> tuple[tuple[Unknown | Any | None, ...], tuple[Any | NDArray[Any], ...], dict[str, Unknown], None]:
    ...

@apply_to_both
def broadcast_to(array, shape, subok=...): # -> tuple[tuple[Unknown | Any | None, ...], tuple[Any | NDArray[Any], ...], dict[str, bool], None]:
    """Broadcast array to the given shape.

    Like `numpy.broadcast_to`, and applied to both unmasked data and mask.
    Note that ``subok`` is taken to mean whether or not subclasses of
    the unmasked data and mask are allowed, i.e., for ``subok=False``,
    a `~astropy.utils.masked.MaskedNDArray` will be returned.
    """
    ...

@dispatched_function
def outer(a, b, out=...): # -> NDArray[Any]:
    ...

@dispatched_function
def empty_like(prototype, dtype=..., order=..., subok=..., shape=...): # -> tuple[Unknown, Unknown, None]:
    """Return a new array with the same shape and type as a given array.

    Like `numpy.empty_like`, but will add an empty mask.
    """
    ...

@dispatched_function
def zeros_like(a, dtype=..., order=..., subok=..., shape=...): # -> tuple[Unknown, Literal[False], None]:
    """Return an array of zeros with the same shape and type as a given array.

    Like `numpy.zeros_like`, but will add an all-false mask.
    """
    ...

@dispatched_function
def ones_like(a, dtype=..., order=..., subok=..., shape=...): # -> tuple[Unknown, Literal[False], None]:
    """Return an array of ones with the same shape and type as a given array.

    Like `numpy.ones_like`, but will add an all-false mask.
    """
    ...

@dispatched_function
def full_like(a, fill_value, dtype=..., order=..., subok=..., shape=...):
    """Return a full array with the same shape and type as a given array.

    Like `numpy.full_like`, but with a mask that is also set.
    If ``fill_value`` is `numpy.ma.masked`, the data will be left unset
    (i.e., as created by `numpy.empty_like`).
    """
    ...

@dispatched_function
def put(a, ind, v, mode=...): # -> None:
    """Replaces specified elements of an array with given values.

    Like `numpy.put`, but for masked array ``a`` and possibly masked
    value ``v``.  Masked indices ``ind`` are not supported.
    """
    ...

@dispatched_function
def putmask(a, mask, values): # -> None:
    """Changes elements of an array based on conditional and input values.

    Like `numpy.putmask`, but for masked array ``a`` and possibly masked
    ``values``.  Masked ``mask`` is not supported.
    """
    ...

@dispatched_function
def place(arr, mask, vals): # -> None:
    """Change elements of an array based on conditional and input values.

    Like `numpy.place`, but for masked array ``a`` and possibly masked
    ``values``.  Masked ``mask`` is not supported.
    """
    ...

@dispatched_function
def copyto(dst, src, casting=..., where=...): # -> None:
    """Copies values from one array to another, broadcasting as necessary.

    Like `numpy.copyto`, but for masked destination ``dst`` and possibly
    masked source ``src``.
    """
    ...

@dispatched_function
def packbits(a, *args, **kwargs): # -> tuple[NDArray[uint8], NDArray[Any], None]:
    ...

@dispatched_function
def unpackbits(a, *args, **kwargs): # -> tuple[NDArray[uint8], NDArray[Any], None]:
    ...

@dispatched_function
def bincount(x, weights=..., minlength=...): # -> tuple[NDArray[intp], NDArray[Any] | None, None]:
    """Count number of occurrences of each value in array of non-negative ints.

    Like `numpy.bincount`, but masked entries in ``x`` will be skipped.
    Any masked entries in ``weights`` will lead the corresponding bin to
    be masked.
    """
    ...

@dispatched_function
def msort(a):
    ...

@dispatched_function
def sort_complex(a):
    ...

if NUMPY_LT_1_20:
    @apply_to_both
    def concatenate(arrays, axis=..., out=...): # -> tuple[tuple[tuple[Unknown | Any | None, ...]], tuple[tuple[Any | NDArray[Any], ...]], dict[str, int], Unknown | None]:
        ...
    
else:
    @dispatched_function
    def concatenate(arrays, axis=..., out=..., dtype=..., casting=...): # -> tuple[Unknown, NDArray[Any], None] | Masked:
        ...
    
@apply_to_both
def append(arr, values, axis=...): # -> tuple[tuple[Unknown | Any | None, ...], tuple[Any | NDArray[Any], ...], dict[str, Unknown | None], None]:
    ...

@dispatched_function
def block(arrays): # -> Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any:
    ...

@dispatched_function
def broadcast_arrays(*args, subok=...): # -> list[Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Unknown | Masked | Any | NDArray[Any]] | Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any | NDArray[Any]:
    """Broadcast arrays to a common shape.

    Like `numpy.broadcast_arrays`, applied to both unmasked data and masks.
    Note that ``subok`` is taken to mean whether or not subclasses of
    the unmasked data and masks are allowed, i.e., for ``subok=False``,
    `~astropy.utils.masked.MaskedNDArray` instances will be returned.
    """
    ...

@apply_to_both
def insert(arr, obj, values, axis=...): # -> tuple[tuple[Unknown | Any | None, Unknown, Unknown | Any | None, Unknown | None], tuple[Any | NDArray[Any], Unknown, Any | NDArray[Any], Unknown | None], dict[Unknown, Unknown], None]:
    """Insert values along the given axis before the given indices.

    Like `numpy.insert` but for possibly masked ``arr`` and ``values``.
    Masked ``obj`` is not supported.
    """
    ...

if NUMPY_LT_1_19:
    @dispatched_function
    def count_nonzero(a, axis=...): # -> int:
        """Counts the number of non-zero values in the array ``a``.

        Like `numpy.count_nonzero`, with masked values counted as 0 or `False`.
        """
        ...
    
else:
    @dispatched_function
    def count_nonzero(a, axis=..., *, keepdims=...):
        """Counts the number of non-zero values in the array ``a``.

        Like `numpy.count_nonzero`, with masked values counted as 0 or `False`.
        """
        ...
    
if NUMPY_LT_1_19:
    ...
else:
    _zeros_like = ...
@dispatched_function
def median(a, axis=..., out=..., overwrite_input=..., keepdims=...): # -> Masked:
    ...

@dispatched_function
def quantile(a, q, axis=..., out=..., **kwargs): # -> Masked:
    ...

@dispatched_function
def percentile(a, q, *args, **kwargs):
    ...

@dispatched_function
def array_equal(a1, a2, equal_nan=...): # -> bool:
    ...

@dispatched_function
def array_equiv(a1, a2): # -> bool:
    ...

@dispatched_function
def where(condition, *args): # -> tuple[Unknown, None, None] | Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any:
    ...

@dispatched_function
def choose(a, choices, out=..., mode=...): # -> Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any:
    """Construct an array from an index array and a set of arrays to choose from.

    Like `numpy.choose`.  Masked indices in ``a`` will lead to masked output
    values and underlying data values are ignored if out of bounds (for
    ``mode='raise'``).  Any values masked in ``choices`` will be propagated
    if chosen.

    """
    ...

@apply_to_both
def select(condlist, choicelist, default=...): # -> tuple[tuple[list[Unknown], tuple[Unknown | Any | None, ...], property | Unknown | Any], tuple[list[Unknown], tuple[Any | NDArray[Any], ...], Any | property | Unknown], dict[Unknown, Unknown], None]:
    """Return an array drawn from elements in choicelist, depending on conditions.

    Like `numpy.select`, with masks in ``choicelist`` are propagated.
    Any masks in ``condlist`` are ignored.

    """
    ...

@dispatched_function
def piecewise(x, condlist, funclist, *args, **kw):
    """Evaluate a piecewise-defined function.

    Like `numpy.piecewise` but for masked input array ``x``.
    Any masks in ``condlist`` are ignored.

    """
    ...

@dispatched_function
def interp(x, xp, fp, *args, **kwargs): # -> Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any:
    """One-dimensional linear interpolation.

    Like `numpy.interp`, but any masked points in ``xp`` and ``fp``
    are ignored.  Any masked values in ``x`` will still be evaluated,
    but masked on output.
    """
    ...

@dispatched_function
def lexsort(keys, axis=...): # -> Any:
    """Perform an indirect stable sort using a sequence of keys.

    Like `numpy.lexsort` but for possibly masked ``keys``.  Masked
    values are sorted towards the end for each key.
    """
    ...

@dispatched_function
def apply_over_axes(func, a, axes): # -> NDArray[Unknown]:
    ...

class MaskedFormat:
    """Formatter for masked array scalars.

    For use in `numpy.array2string`, wrapping the regular formatters such
    that if a value is masked, its formatted string is replaced.

    Typically initialized using the ``from_data`` class method.
    """
    def __init__(self, format_function) -> None:
        ...
    
    def __call__(self, x): # -> LiteralString:
        ...
    
    @classmethod
    def from_data(cls, data, **options): # -> Self@MaskedFormat:
        ...
    


@dispatched_function
def array2string(a, max_line_width=..., precision=..., suppress_small=..., separator=..., prefix=..., style=..., formatter=..., threshold=..., edgeitems=..., sign=..., floatmode=..., suffix=...): # -> Literal['[]']:
    ...

@dispatched_function
def array_str(a, max_line_width=..., precision=..., suppress_small=...): # -> Literal['[]']:
    ...

_nanfunc_fill_values = ...
def masked_nanfunc(nanfuncname): # -> (a: Unknown, *args: Unknown, **kwargs: Unknown) -> Any:
    ...

__all__ += sorted(helper.__name__ for helper in (set(APPLY_TO_BOTH_FUNCTIONS.values()) | set(DISPATCHED_FUNCTIONS.values())) if helper.__doc__)