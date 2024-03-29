"""
This type stub file was generated by pyright.
"""

import numpy as np

class FITS_record:
    """
    FITS record class.

    `FITS_record` is used to access records of the `FITS_rec` object.
    This will allow us to deal with scaled columns.  It also handles
    conversion/scaling of columns in ASCII tables.  The `FITS_record`
    class expects a `FITS_rec` object as input.
    """
    def __init__(self, input, row=..., start=..., end=..., step=..., base=..., **kwargs) -> None:
        """
        Parameters
        ----------
        input : array
            The array to wrap.
        row : int, optional
            The starting logical row of the array.
        start : int, optional
            The starting column in the row associated with this object.
            Used for subsetting the columns of the `FITS_rec` object.
        end : int, optional
            The ending column in the row associated with this object.
            Used for subsetting the columns of the `FITS_rec` object.
        """
        ...
    
    def __getitem__(self, key): # -> Self@FITS_record:
        ...
    
    def __setitem__(self, key, value): # -> None:
        ...
    
    def __len__(self): # -> int:
        ...
    
    def __repr__(self): # -> LiteralString:
        """
        Display a single row.
        """
        ...
    
    def field(self, field): # -> Self@FITS_record:
        """
        Get the field data of the record.
        """
        ...
    
    def setfield(self, field, value): # -> None:
        """
        Set the field data of the record.
        """
        ...
    


class FITS_rec(np.recarray):
    """
    FITS record array class.

    `FITS_rec` is the data part of a table HDU's data part.  This is a layer
    over the `~numpy.recarray`, so we can deal with scaled columns.

    It inherits all of the standard methods from `numpy.ndarray`.
    """
    _record_type = FITS_record
    _character_as_bytes = ...
    def __new__(subtype, input): # -> recarray[Any, dtype[Any]]:
        """
        Construct a FITS record array from a recarray.
        """
        ...
    
    def __setstate__(self, state): # -> None:
        ...
    
    def __reduce__(self): # -> tuple[str | Any, str | Any, Unknown]:
        """
        Return a 3-tuple for pickling a FITS_rec. Use the super-class
        functionality but then add in a tuple of FITS_rec-specific
        values that get used in __setstate__.
        """
        ...
    
    def __array_finalize__(self, obj): # -> None:
        ...
    
    @classmethod
    def from_columns(cls, columns, nrows=..., fill=..., character_as_bytes=...): # -> Self@FITS_rec:
        """
        Given a `ColDefs` object of unknown origin, initialize a new `FITS_rec`
        object.

        .. note::

            This was originally part of the ``new_table`` function in the table
            module but was moved into a class method since most of its
            functionality always had more to do with initializing a `FITS_rec`
            object than anything else, and much of it also overlapped with
            ``FITS_rec._scale_back``.

        Parameters
        ----------
        columns : sequence of `Column` or a `ColDefs`
            The columns from which to create the table data.  If these
            columns have data arrays attached that data may be used in
            initializing the new table.  Otherwise the input columns
            will be used as a template for a new table with the requested
            number of rows.

        nrows : int
            Number of rows in the new table.  If the input columns have data
            associated with them, the size of the largest input column is used.
            Otherwise the default is 0.

        fill : bool
            If `True`, will fill all cells with zeros or blanks.  If
            `False`, copy the data from input, undefined cells will still
            be filled with zeros/blanks.
        """
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __getattribute__(self, attr): # -> Any | NDArray[Any] | _VLF:
        ...
    
    def __getitem__(self, key): # -> Any | NDArray[Any] | _VLF | _record_type | Self@FITS_rec:
        ...
    
    def __setitem__(self, key, value): # -> None:
        ...
    
    def copy(self, order=...): # -> FITS_rec:
        """
        The Numpy documentation lies; `numpy.ndarray.copy` is not equivalent to
        `numpy.copy`.  Differences include that it re-views the copied array as
        self's ndarray subclass, as though it were taking a slice; this means
        ``__array_finalize__`` is called and the copy shares all the array
        attributes (including ``._converted``!).  So we need to make a deep
        copy of all those attributes so that the two arrays truly do not share
        any data.
        """
        ...
    
    @property
    def columns(self): # -> Any | None:
        """A user-visible accessor for the coldefs."""
        ...
    
    def __del__(self): # -> None:
        ...
    
    @property
    def names(self): # -> list[Unknown] | Any | None:
        """List of column names."""
        ...
    
    @property
    def formats(self): # -> Any | None:
        """List of column FITS formats."""
        ...
    
    def field(self, key): # -> NDArray[Any] | _VLF:
        """
        A view of a `Column`'s data as an array.
        """
        ...
    
    def tolist(self): # -> list[list[list[Any]]]:
        ...
    


class _UnicodeArrayEncodeError(UnicodeEncodeError):
    def __init__(self, encoding, object_, start, end, reason, index) -> None:
        ...
    


