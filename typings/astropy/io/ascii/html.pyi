"""
This type stub file was generated by pyright.
"""

from . import core

"""An extensible HTML table reader and writer.

html.py:
  Classes to read and write HTML tables

`BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_
must be installed to read HTML tables.
"""
class SoupString(str):
    """
    Allows for strings to hold BeautifulSoup data.
    """
    def __new__(cls, *args, **kwargs): # -> Self@SoupString:
        ...
    
    def __init__(self, val) -> None:
        ...
    


class ListWriter:
    """
    Allows for XMLWriter to write to a list instead of a file.
    """
    def __init__(self, out) -> None:
        ...
    
    def write(self, data): # -> None:
        ...
    


def identify_table(soup, htmldict, numtable): # -> Literal[False]:
    """
    Checks whether the given BeautifulSoup tag is the table
    the user intends to process.
    """
    ...

class HTMLInputter(core.BaseInputter):
    """
    Input lines of HTML in a valid form.

    This requires `BeautifulSoup
    <http://www.crummy.com/software/BeautifulSoup/>`_ to be installed.
    """
    def process_lines(self, lines): # -> list[SoupString]:
        """
        Convert the given input into a list of SoupString rows
        for further processing.
        """
        ...
    


class HTMLSplitter(core.BaseSplitter):
    """
    Split HTML table data.
    """
    def __call__(self, lines): # -> Generator[list[tuple[Unknown, Unknown] | Unknown] | list[Unknown], None, None]:
        """
        Return HTML data from lines as a generator.
        """
        ...
    


class HTMLOutputter(core.TableOutputter):
    """
    Output the HTML data as an ``astropy.table.Table`` object.

    This subclass allows for the final table to contain
    multidimensional columns (defined using the colspan attribute
    of <th>).
    """
    default_converters = ...
    def __call__(self, cols, meta): # -> Table:
        """
        Process the data in multidimensional columns.
        """
        ...
    


class HTMLHeader(core.BaseHeader):
    splitter_class = HTMLSplitter
    def start_line(self, lines): # -> int | None:
        """
        Return the line number at which header data begins.
        """
        ...
    


class HTMLData(core.BaseData):
    splitter_class = HTMLSplitter
    def start_line(self, lines): # -> int:
        """
        Return the line number at which table data begins.
        """
        ...
    
    def end_line(self, lines): # -> int | None:
        """
        Return the line number at which table data ends.
        """
        ...
    


class HTML(core.BaseReader):
    """HTML format table.

    In order to customize input and output, a dict of parameters may
    be passed to this class holding specific customizations.

    **htmldict** : Dictionary of parameters for HTML input/output.

        * css : Customized styling
            If present, this parameter will be included in a <style>
            tag and will define stylistic attributes of the output.

        * table_id : ID for the input table
            If a string, this defines the HTML id of the table to be processed.
            If an integer, this specifies the index of the input table in the
            available tables. Unless this parameter is given, the reader will
            use the first table found in the input file.

        * multicol : Use multi-dimensional columns for output
            The writer will output tuples as elements of multi-dimensional
            columns if this parameter is true, and if not then it will
            use the syntax 1.36583e-13 .. 1.36583e-13 for output. If not
            present, this parameter will be true by default.

        * raw_html_cols : column name or list of names with raw HTML content
            This allows one to include raw HTML content in the column output,
            for instance to include link references in a table.  This option
            requires that the bleach package be installed.  Only whitelisted
            tags are allowed through for security reasons (see the
            raw_html_clean_kwargs arg).

        * raw_html_clean_kwargs : dict of keyword args controlling HTML cleaning
            Raw HTML will be cleaned to prevent unsafe HTML from ending up in
            the table output.  This is done by calling ``bleach.clean(data,
            **raw_html_clean_kwargs)``.  For details on the available options
            (e.g. tag whitelist) see:
            https://bleach.readthedocs.io/en/latest/clean.html

        * parser : Specific HTML parsing library to use
            If specified, this specifies which HTML parsing library
            BeautifulSoup should use as a backend. The options to choose
            from are 'html.parser' (the standard library parser), 'lxml'
            (the recommended parser), 'xml' (lxml's XML parser), and
            'html5lib'. html5lib is a highly lenient parser and therefore
            might work correctly for unusual input if a different parser
            fails.

        * jsfiles : list of js files to include when writing table.

        * cssfiles : list of css files to include when writing table.

        * js : js script to include in the body when writing table.

        * table_class : css class for the table

    """
    _format_name = ...
    _io_registry_format_aliases = ...
    _io_registry_suffix = ...
    _description = ...
    header_class = HTMLHeader
    data_class = HTMLData
    inputter_class = HTMLInputter
    max_ndim = ...
    def __init__(self, htmldict=...) -> None:
        """
        Initialize classes for HTML reading and writing.
        """
        ...
    
    def read(self, table): # -> Table:
        """
        Read the ``table`` in HTML format and return a resulting ``Table``.
        """
        ...
    
    def write(self, table): # -> list[str]:
        """
        Return data in ``table`` converted to HTML as a list of strings.
        """
        ...
    
    def fill_values(self, col, col_str_iters): # -> Generator[Unknown, None, None]:
        """
        Return an iterator of the values with replacements based on fill_values
        """
        ...
    

