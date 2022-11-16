"""
This type stub file was generated by pyright.
"""

import numpy as np
from astropy.utils.shapes import NDArrayShapeMethods
from astropy.utils.data_info import ParentDtypeInfo

"""
Built-in mask mixin class.

The design uses `Masked` as a factory class which automatically
generates new subclasses for any data class that is itself a
subclass of a predefined masked class, with `MaskedNDArray`
providing such a predefined class for `~numpy.ndarray`.

Generally, any new predefined class should override the
``from_unmasked(data, mask, copy=False)`` class method that
creates an instance from unmasked data and a mask, as well as
the ``unmasked`` property that returns just the data.
The `Masked` class itself provides a base ``mask`` property,
which can also be overridden if needed.

"""
__all__ = ['Masked', 'MaskedNDArray']
get__doc__ = ...
class Masked(NDArrayShapeMethods):
    """A scalar value or array of values with associated mask.

    The resulting instance will take its exact type from whatever the
    contents are, with the type generated on the fly as needed.

    Parameters
    ----------
    data : array-like
        The data for which a mask is to be added.  The result will be a
        a subclass of the type of ``data``.
    mask : array-like of bool, optional
        The initial mask to assign.  If not given, taken from the data.
    copy : bool
        Whether the data and mask should be copied. Default: `False`.

    """
    _base_classes = ...
    _masked_classes = ...
    def __new__(cls, *args, **kwargs): # -> Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any | Self@Masked:
        ...
    
    def __init_subclass__(cls, base_cls=..., data_cls=..., **kwargs): # -> None:
        """Register a Masked subclass.

        Parameters
        ----------
        base_cls : type, optional
            If given, it is taken to mean that ``cls`` can be used as
            a base for masked versions of all subclasses of ``base_cls``,
            so it is registered as such in ``_base_classes``.
        data_cls : type, optional
            If given, ``cls`` should will be registered as the masked version of
            ``data_cls``.  Will set the private ``cls._data_cls`` attribute,
            and auto-generate a docstring if not present already.
        **kwargs
            Passed on for possible further initialization by superclasses.

        """
        ...
    
    @classmethod
    def from_unmasked(cls, data, mask=..., copy=...): # -> Self@Masked:
        """Create an instance from unmasked data and a mask."""
        ...
    
    mask = ...
    @property
    def unmasked(self):
        """The unmasked values.

        See Also
        --------
        astropy.utils.masked.Masked.filled
        """
        ...
    
    def filled(self, fill_value):
        """Get a copy of the underlying data, with masked values filled in.

        Parameters
        ----------
        fill_value : object
            Value to replace masked values with.

        See Also
        --------
        astropy.utils.masked.Masked.unmasked
        """
        ...
    
    def __setitem__(self, item, value): # -> None:
        ...
    


class MaskedInfoBase:
    mask_val = ...
    def __init__(self, bound=...) -> None:
        ...
    


class MaskedNDArrayInfo(MaskedInfoBase, ParentDtypeInfo):
    """
    Container for meta information like name, description, format.
    """
    attr_names = ...
    _represent_as_dict_primary_data = ...


class MaskedArraySubclassInfo(MaskedInfoBase):
    """Mixin class to create a subclasses such as MaskedQuantityInfo."""
    ...


class MaskedIterator:
    """
    Flat iterator object to iterate over Masked Arrays.

    A `~astropy.utils.masked.MaskedIterator` iterator is returned by ``m.flat``
    for any masked array ``m``.  It allows iterating over the array as if it
    were a 1-D array, either in a for-loop or by calling its `next` method.

    Iteration is done in C-contiguous style, with the last index varying the
    fastest. The iterator can also be indexed using basic slicing or
    advanced indexing.

    Notes
    -----
    The design of `~astropy.utils.masked.MaskedIterator` follows that of
    `~numpy.ma.core.MaskedIterator`.  It is not exported by the
    `~astropy.utils.masked` module.  Instead of instantiating directly,
    use the ``flat`` method in the masked array instance.
    """
    def __init__(self, m) -> None:
        ...
    
    def __iter__(self): # -> Self@MaskedIterator:
        ...
    
    def __getitem__(self, indx):
        ...
    
    def __setitem__(self, index, value): # -> None:
        ...
    
    def __next__(self):
        """
        Return the next value, or raise StopIteration.
        """
        ...
    
    next = ...


class MaskedNDArray(Masked, np.ndarray, base_cls=np.ndarray, data_cls=np.ndarray):
    _mask = ...
    info = ...
    def __new__(cls, *args, mask=..., **kwargs): # -> Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_] | Masked | Any | Self@MaskedNDArray:
        """Get data class instance from arguments and then set mask."""
        ...
    
    def __init_subclass__(cls, **kwargs): # -> None:
        ...
    
    @classmethod
    def from_unmasked(cls, data, mask=..., copy=...):
        ...
    
    @property
    def unmasked(self):
        ...
    
    @property
    def flat(self): # -> MaskedIterator:
        """A 1-D iterator over the Masked array.

        This returns a ``MaskedIterator`` instance, which behaves the same
        as the `~numpy.flatiter` instance returned by `~numpy.ndarray.flat`,
        and is similar to Python's built-in iterator, except that it also
        allows assignment.
        """
        ...
    
    def view(self, dtype=..., type=...): # -> NDArray[Any]:
        """New view of the masked array.

        Like `numpy.ndarray.view`, but always returning a masked array subclass.
        """
        ...
    
    def __array_finalize__(self, obj): # -> None:
        ...
    
    @property
    def shape(self): # -> _Shape:
        """The shape of the data and the mask.

        Usually used to get the current shape of an array, but may also be
        used to reshape the array in-place by assigning a tuple of array
        dimensions to it.  As with `numpy.reshape`, one of the new shape
        dimensions can be -1, in which case its value is inferred from the
        size of the array and the remaining dimensions.

        Raises
        ------
        AttributeError
            If a copy is required, of either the data or the mask.

        """
        ...
    
    @shape.setter
    def shape(self, shape): # -> None:
        ...
    
    _eq_simple = ...
    _ne_simple = ...
    __lt__ = ...
    __le__ = ...
    __gt__ = ...
    __ge__ = ...
    def __eq__(self, other) -> bool:
        ...
    
    def __ne__(self, other) -> bool:
        ...
    
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        ...
    
    def __array_function__(self, function, types, args, kwargs): # -> Any | _NotImplementedType | Masked | tuple[Unknown, ...] | Type[Masked] | Type[MaskedArray[Unknown, Unknown]] | Type[MaskedNDArray] | Type[_]:
        ...
    
    def trace(self, offset=..., axis1=..., axis2=..., dtype=..., out=...):
        ...
    
    def min(self, axis=..., out=..., **kwargs): # -> Any:
        ...
    
    def max(self, axis=..., out=..., **kwargs): # -> Any:
        ...
    
    def nonzero(self): # -> tuple[Unknown, ...] | tuple[NDArray[intp], ...]:
        ...
    
    def compress(self, condition, axis=..., out=...): # -> Self@MaskedNDArray:
        ...
    
    def repeat(self, repeats, axis=...): # -> Self@MaskedNDArray:
        ...
    
    def choose(self, choices, out=..., mode=...):
        ...
    
    def argmin(self, axis=..., out=...): # -> Any:
        ...
    
    def argmax(self, axis=..., out=...): # -> Any:
        ...
    
    def argsort(self, axis=..., kind=..., order=...): # -> Any:
        """Returns the indices that would sort an array.

        Perform an indirect sort along the given axis on both the array
        and the mask, with masked items being sorted to the end.

        Parameters
        ----------
        axis : int or None, optional
            Axis along which to sort.  The default is -1 (the last axis).
            If None, the flattened array is used.
        kind : str or None, ignored.
            The kind of sort.  Present only to allow subclasses to work.
        order : str or list of str.
            For an array with fields defined, the fields to compare first,
            second, etc.  A single field can be specified as a string, and not
            all fields need be specified, but unspecified fields will still be
            used, in dtype order, to break ties.

        Returns
        -------
        index_array : ndarray, int
            Array of indices that sorts along the specified ``axis``.  Use
            ``np.take_along_axis(self, index_array, axis=axis)`` to obtain
            the sorted array.

        """
        ...
    
    def sort(self, axis=..., kind=..., order=...): # -> None:
        """Sort an array in-place. Refer to `numpy.sort` for full documentation."""
        ...
    
    def argpartition(self, kth, axis=..., kind=..., order=...): # -> Any:
        ...
    
    def partition(self, kth, axis=..., kind=..., order=...): # -> None:
        ...
    
    def cumsum(self, axis=..., dtype=..., out=...): # -> NDArray[Any]:
        ...
    
    def cumprod(self, axis=..., dtype=..., out=...): # -> NDArray[Any]:
        ...
    
    def clip(self, min=..., max=..., out=..., **kwargs): # -> NDArray[Any]:
        """Return an array whose values are limited to ``[min, max]``.

        Like `~numpy.clip`, but any masked values in ``min`` and ``max``
        are ignored for clipping.  The mask of the input array is propagated.
        """
        ...
    
    def mean(self, axis=..., dtype=..., out=..., keepdims=..., *, where=...): # -> Any:
        ...
    
    def var(self, axis=..., dtype=..., out=..., ddof=..., keepdims=..., *, where=...): # -> Any:
        ...
    
    def std(self, axis=..., dtype=..., out=..., ddof=..., keepdims=..., *, where=...): # -> Any:
        ...
    
    def __bool__(self): # -> bool:
        ...
    
    def any(self, axis=..., out=..., keepdims=..., *, where=...): # -> Any:
        ...
    
    def all(self, axis=..., out=..., keepdims=..., *, where=...): # -> Any:
        ...
    
    def __str__(self) -> str:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __format__(self, format_spec): # -> LiteralString | str:
        ...
    


class MaskedRecarray(np.recarray, MaskedNDArray, data_cls=np.recarray):
    def __array_finalize__(self, obj): # -> None:
        ...
    
    def getfield(self, dtype, offset=...): # -> Any:
        ...
    
    def setfield(self, val, dtype, offset=...): # -> None:
        ...
    


