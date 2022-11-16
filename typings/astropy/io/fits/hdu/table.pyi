"""
This type stub file was generated by pyright.
"""

import csv
from .base import ExtensionHDU, _ValidHDU
from astropy.io.fits.column import ColDefs, _AsciiColDefs
from astropy.io.fits.fitsrec import FITS_rec
from astropy.utils import lazyproperty

class FITSTableDumpDialect(csv.excel):
    """
    A CSV dialect for the Astropy format of ASCII dumps of FITS tables.
    """
    delimiter = ...
    lineterminator = ...
    quotechar = ...
    quoting = ...
    skipinitialspace = ...


class _TableLikeHDU(_ValidHDU):
    """
    A class for HDUs that have table-like data.  This is used for both
    Binary/ASCII tables as well as Random Access Group HDUs (which are
    otherwise too dissimilar for tables to use _TableBaseHDU directly).
    """
    _data_type = FITS_rec
    _columns_type = ColDefs
    _uint = ...
    @classmethod
    def match_header(cls, header):
        """
        This is an abstract HDU type for HDUs that contain table-like data.
        This is even more abstract than _TableBaseHDU which is specifically for
        the standard ASCII and Binary Table types.
        """
        ...
    
    @classmethod
    def from_columns(cls, columns, header=..., nrows=..., fill=..., character_as_bytes=..., **kwargs): # -> Self@_TableLikeHDU:
        """
        Given either a `ColDefs` object, a sequence of `Column` objects,
        or another table HDU or table data (a `FITS_rec` or multi-field
        `numpy.ndarray` or `numpy.recarray` object, return a new table HDU of
        the class this method was called on using the column definition from
        the input.

        See also `FITS_rec.from_columns`.

        Parameters
        ----------
        columns : sequence of `Column`, `ColDefs` -like
            The columns from which to create the table data, or an object with
            a column-like structure from which a `ColDefs` can be instantiated.
            This includes an existing `BinTableHDU` or `TableHDU`, or a
            `numpy.recarray` to give some examples.

            If these columns have data arrays attached that data may be used in
            initializing the new table.  Otherwise the input columns will be
            used as a template for a new table with the requested number of
            rows.

        header : `Header`
            An optional `Header` object to instantiate the new HDU yet.  Header
            keywords specifically related to defining the table structure (such
            as the "TXXXn" keywords like TTYPEn) will be overridden by the
            supplied column definitions, but all other informational and data
            model-specific keywords are kept.

        nrows : int
            Number of rows in the new table.  If the input columns have data
            associated with them, the size of the largest input column is used.
            Otherwise the default is 0.

        fill : bool
            If `True`, will fill all cells with zeros or blanks.  If `False`,
            copy the data from input, undefined cells will still be filled with
            zeros/blanks.

        character_as_bytes : bool
            Whether to return bytes for string columns when accessed from the
            HDU. By default this is `False` and (unicode) strings are returned,
            but for large tables this may use up a lot of memory.

        Notes
        -----
        Any additional keyword arguments accepted by the HDU class's
        ``__init__`` may also be passed in as keyword arguments.
        """
        ...
    
    @lazyproperty
    def columns(self): # -> ColDefs:
        """
        The :class:`ColDefs` objects describing the columns in this table.
        """
        ...
    


class _TableBaseHDU(ExtensionHDU, _TableLikeHDU):
    """
    FITS table extension base HDU class.

    Parameters
    ----------
    data : array
        Data to be used.
    header : `Header` instance
        Header to be used. If the ``data`` is also specified, header keywords
        specifically related to defining the table structure (such as the
        "TXXXn" keywords like TTYPEn) will be overridden by the supplied column
        definitions, but all other informational and data model-specific
        keywords are kept.
    name : str
        Name to be populated in ``EXTNAME`` keyword.
    uint : bool, optional
        Set to `True` if the table contains unsigned integer columns.
    ver : int > 0 or None, optional
        The ver of the HDU, will be the value of the keyword ``EXTVER``.
        If not given or None, it defaults to the value of the ``EXTVER``
        card of the ``header`` or 1.
        (default: None)
    character_as_bytes : bool
        Whether to return bytes for string columns. By default this is `False`
        and (unicode) strings are returned, but this does not respect memory
        mapping and loads the whole column in memory when accessed.
    """
    _manages_own_heap = ...
    def __init__(self, data=..., header=..., name=..., uint=..., ver=..., character_as_bytes=...) -> None:
        ...
    
    @classmethod
    def match_header(cls, header):
        """
        This is an abstract type that implements the shared functionality of
        the ASCII and Binary Table HDU types, which should be used instead of
        this.
        """
        ...
    
    @lazyproperty
    def columns(self): # -> Any | _columns_type | None:
        """
        The :class:`ColDefs` objects describing the columns in this table.
        """
        ...
    
    @lazyproperty
    def data(self): # -> _data_type:
        ...
    
    @data.setter
    def data(self, data): # -> FITS_rec | None:
        ...
    
    def update(self): # -> None:
        """
        Update header keywords to reflect recent changes of columns.
        """
        ...
    
    def copy(self): # -> Self@_TableBaseHDU:
        """
        Make a copy of the table HDU, both header and data are copied.
        """
        ...
    


class TableHDU(_TableBaseHDU):
    """
    FITS ASCII table extension HDU class.

    Parameters
    ----------
    data : array or `FITS_rec`
        Data to be used.
    header : `Header`
        Header to be used.
    name : str
        Name to be populated in ``EXTNAME`` keyword.
    ver : int > 0 or None, optional
        The ver of the HDU, will be the value of the keyword ``EXTVER``.
        If not given or None, it defaults to the value of the ``EXTVER``
        card of the ``header`` or 1.
        (default: None)
    character_as_bytes : bool
        Whether to return bytes for string columns. By default this is `False`
        and (unicode) strings are returned, but this does not respect memory
        mapping and loads the whole column in memory when accessed.

    """
    _extension = ...
    _ext_comment = ...
    _padding_byte = ...
    _columns_type = _AsciiColDefs
    __format_RE = ...
    def __init__(self, data=..., header=..., name=..., ver=..., character_as_bytes=...) -> None:
        ...
    
    @classmethod
    def match_header(cls, header): # -> bool:
        ...
    


class BinTableHDU(_TableBaseHDU):
    """
    Binary table HDU class.

    Parameters
    ----------
    data : array, `FITS_rec`, or `~astropy.table.Table`
        Data to be used.
    header : `Header`
        Header to be used.
    name : str
        Name to be populated in ``EXTNAME`` keyword.
    uint : bool, optional
        Set to `True` if the table contains unsigned integer columns.
    ver : int > 0 or None, optional
        The ver of the HDU, will be the value of the keyword ``EXTVER``.
        If not given or None, it defaults to the value of the ``EXTVER``
        card of the ``header`` or 1.
        (default: None)
    character_as_bytes : bool
        Whether to return bytes for string columns. By default this is `False`
        and (unicode) strings are returned, but this does not respect memory
        mapping and loads the whole column in memory when accessed.

    """
    _extension = ...
    _ext_comment = ...
    def __init__(self, data=..., header=..., name=..., uint=..., ver=..., character_as_bytes=...) -> None:
        ...
    
    @classmethod
    def match_header(cls, header): # -> bool:
        ...
    
    _tdump_file_format = ...
    def dump(self, datafile=..., cdfile=..., hfile=..., overwrite=...): # -> None:
        """
        Dump the table HDU to a file in ASCII format.  The table may be dumped
        in three separate files, one containing column definitions, one
        containing header parameters, and one for table data.

        Parameters
        ----------
        datafile : path-like or file-like, optional
            Output data file.  The default is the root name of the
            fits file associated with this HDU appended with the
            extension ``.txt``.

        cdfile : path-like or file-like, optional
            Output column definitions file.  The default is `None`, no
            column definitions output is produced.

        hfile : path-like or file-like, optional
            Output header parameters file.  The default is `None`,
            no header parameters output is produced.

        overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` if ``False`` and the output file exists. Default is
            ``False``.

        Notes
        -----
        The primary use for the `dump` method is to allow viewing and editing
        the table data and parameters in a standard text editor.
        The `load` method can be used to create a new table from the three
        plain text (ASCII) files.
        """
        ...
    
    if isinstance(dump.__doc__, str):
        ...
    def load(cls, datafile, cdfile=..., hfile=..., replace=..., header=...):
        """
        Create a table from the input ASCII files.  The input is from up to
        three separate files, one containing column definitions, one containing
        header parameters, and one containing column data.

        The column definition and header parameters files are not required.
        When absent the column definitions and/or header parameters are taken
        from the header object given in the header argument; otherwise sensible
        defaults are inferred (though this mode is not recommended).

        Parameters
        ----------
        datafile : path-like or file-like
            Input data file containing the table data in ASCII format.

        cdfile : path-like or file-like, optional
            Input column definition file containing the names,
            formats, display formats, physical units, multidimensional
            array dimensions, undefined values, scale factors, and
            offsets associated with the columns in the table.  If
            `None`, the column definitions are taken from the current
            values in this object.

        hfile : path-like or file-like, optional
            Input parameter definition file containing the header
            parameter definitions to be associated with the table.  If
            `None`, the header parameter definitions are taken from
            the current values in this objects header.

        replace : bool, optional
            When `True`, indicates that the entire header should be
            replaced with the contents of the ASCII file instead of
            just updating the current header.

        header : `~astropy.io.fits.Header`, optional
            When the cdfile and hfile are missing, use this Header object in
            the creation of the new table and HDU.  Otherwise this Header
            supersedes the keywords from hfile, which is only used to update
            values not present in this Header, unless ``replace=True`` in which
            this Header's values are completely replaced with the values from
            hfile.

        Notes
        -----
        The primary use for the `load` method is to allow the input of ASCII
        data that was edited in a standard text editor of the table data and
        parameters.  The `dump` method can be used to create the initial ASCII
        files.
        """
        ...
    
    if isinstance(load.__doc__, str):
        ...
    load = ...


