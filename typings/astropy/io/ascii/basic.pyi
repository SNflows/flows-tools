"""
This type stub file was generated by pyright.
"""

from . import core

"""An extensible ASCII table reader and writer.

basic.py:
  Basic table read / write functionality for simple character
  delimited files with various options for column header definition.

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""
class BasicHeader(core.BaseHeader):
    """
    Basic table Header Reader

    Set a few defaults for common ascii table formats
    (start at line 0, comments begin with ``#`` and possibly white space)
    """
    start_line = ...
    comment = ...
    write_comment = ...


class BasicData(core.BaseData):
    """
    Basic table Data Reader

    Set a few defaults for common ascii table formats
    (start at line 1, comments begin with ``#`` and possibly white space)
    """
    start_line = ...
    comment = ...
    write_comment = ...


class Basic(core.BaseReader):
    r"""Character-delimited table with a single header line at the top.

    Lines beginning with a comment character (default='#') as the first
    non-whitespace character are comments.

    Example table::

      # Column definition is the first uncommented line
      # Default delimiter is the space character.
      apples oranges pears

      # Data starts after the header column definition, blank lines ignored
      1 2 3
      4 5 6
    """
    _format_name = ...
    _description = ...
    _io_registry_format_aliases = ...
    header_class = BasicHeader
    data_class = BasicData


class NoHeaderHeader(BasicHeader):
    """
    Reader for table header without a header

    Set the start of header line number to `None`, which tells the basic
    reader there is no header line.
    """
    start_line = ...


class NoHeaderData(BasicData):
    """
    Reader for table data without a header

    Data starts at first uncommented line since there is no header line.
    """
    start_line = ...


class NoHeader(Basic):
    """Character-delimited table with no header line.

    When reading, columns are autonamed using header.auto_format which defaults
    to "col%d".  Otherwise this reader the same as the :class:`Basic` class
    from which it is derived.  Example::

      # Table data
      1 2 "hello there"
      3 4 world

    """
    _format_name = ...
    _description = ...
    header_class = NoHeaderHeader
    data_class = NoHeaderData


class CommentedHeaderHeader(BasicHeader):
    """
    Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """
    def process_lines(self, lines): # -> Generator[Unknown, None, None]:
        """
        Return only lines that start with the comment regexp.  For these
        lines strip out the matching characters.
        """
        ...
    
    def write(self, lines): # -> None:
        ...
    


class CommentedHeader(Basic):
    """Character-delimited table with column names in a comment line.

    When reading, ``header_start`` can be used to specify the
    line index of column names, and it can be a negative index (for example -1
    for the last commented line).  The default delimiter is the <space>
    character.

    This matches the format produced by ``np.savetxt()``, with ``delimiter=','``,
    and ``header='<comma-delimited-column-names-list>'``.

    Example::

      # col1 col2 col3
      # Comment line
      1 2 3
      4 5 6

    """
    _format_name = ...
    _description = ...
    header_class = CommentedHeaderHeader
    data_class = NoHeaderData
    def read(self, table): # -> Table:
        """
        Read input data (file-like object, filename, list of strings, or
        single string) into a Table and return the result.
        """
        ...
    
    def write_header(self, lines, meta): # -> None:
        """
        Write comment lines after, rather than before, the header.
        """
        ...
    


class TabHeaderSplitter(core.DefaultSplitter):
    """Split lines on tab and do not remove whitespace"""
    delimiter = ...
    def process_line(self, line):
        ...
    


class TabDataSplitter(TabHeaderSplitter):
    """
    Don't strip data value whitespace since that is significant in TSV tables
    """
    process_val = ...
    skipinitialspace = ...


class TabHeader(BasicHeader):
    """
    Reader for header of tables with tab separated header
    """
    splitter_class = TabHeaderSplitter


class TabData(BasicData):
    """
    Reader for data of tables with tab separated data
    """
    splitter_class = TabDataSplitter


class Tab(Basic):
    """Tab-separated table.

    Unlike the :class:`Basic` reader, whitespace is not stripped from the
    beginning and end of either lines or individual column values.

    Example::

      col1 <tab> col2 <tab> col3
      # Comment line
      1 <tab> 2 <tab> 5

    """
    _format_name = ...
    _description = ...
    header_class = TabHeader
    data_class = TabData


class CsvSplitter(core.DefaultSplitter):
    """
    Split on comma for CSV (comma-separated-value) tables
    """
    delimiter = ...


class CsvHeader(BasicHeader):
    """
    Header that uses the :class:`astropy.io.ascii.basic.CsvSplitter`
    """
    splitter_class = CsvSplitter
    comment = ...
    write_comment = ...


class CsvData(BasicData):
    """
    Data that uses the :class:`astropy.io.ascii.basic.CsvSplitter`
    """
    splitter_class = CsvSplitter
    fill_values = ...
    comment = ...
    write_comment = ...


class Csv(Basic):
    """CSV (comma-separated-values) table.

    This file format may contain rows with fewer entries than the number of
    columns, a situation that occurs in output from some spreadsheet editors.
    The missing entries are marked as masked in the output table.

    Masked values (indicated by an empty '' field value when reading) are
    written out in the same way with an empty ('') field.  This is different
    from the typical default for `astropy.io.ascii` in which missing values are
    indicated by ``--``.

    Since the `CSV format <https://tools.ietf.org/html/rfc4180>`_ does not
    formally support comments, any comments defined for the table via
    ``tbl.meta['comments']`` are ignored by default. If you would still like to
    write those comments then include a keyword ``comment='#'`` to the
    ``write()`` call.

    Example::

      num,ra,dec,radius,mag
      1,32.23222,10.1211
      2,38.12321,-88.1321,2.2,17.0

    """
    _format_name = ...
    _io_registry_format_aliases = ...
    _io_registry_can_write = ...
    _io_registry_suffix = ...
    _description = ...
    header_class = CsvHeader
    data_class = CsvData
    def inconsistent_handler(self, str_vals, ncols):
        """
        Adjust row if it is too short.

        If a data row is shorter than the header, add empty values to make it the
        right length.
        Note that this will *not* be called if the row already matches the header.

        Parameters
        ----------
        str_vals : list
            A list of value strings from the current row of the table.
        ncols : int
            The expected number of entries from the table header.

        Returns
        -------
        str_vals : list
            List of strings to be parsed into data entries in the output table.
        """
        ...
    


class RdbHeader(TabHeader):
    """
    Header for RDB tables
    """
    col_type_map = ...
    def get_type_map_key(self, col):
        ...
    
    def get_cols(self, lines): # -> None:
        """
        Initialize the header Column objects from the table ``lines``.

        This is a specialized get_cols for the RDB type:
        Line 0: RDB col names
        Line 1: RDB col definitions
        Line 2+: RDB data rows

        Parameters
        ----------
        lines : list
            List of table lines

        Returns
        -------
        None

        """
        ...
    
    def write(self, lines): # -> None:
        ...
    


class RdbData(TabData):
    """
    Data reader for RDB data. Starts reading at line 2.
    """
    start_line = ...


class Rdb(Tab):
    """Tab-separated file with an extra line after the column definition line that
    specifies either numeric (N) or string (S) data.

    See: https://www.drdobbs.com/rdb-a-unix-command-line-database/199101326

    Example::

      col1 <tab> col2 <tab> col3
      N <tab> S <tab> N
      1 <tab> 2 <tab> 5

    """
    _format_name = ...
    _io_registry_format_aliases = ...
    _io_registry_suffix = ...
    _description = ...
    header_class = RdbHeader
    data_class = RdbData

