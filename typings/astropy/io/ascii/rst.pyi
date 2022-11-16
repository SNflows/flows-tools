"""
This type stub file was generated by pyright.
"""

from .core import DefaultSplitter
from .fixedwidth import FixedWidth, FixedWidthData, FixedWidthHeader, FixedWidthTwoLineDataSplitter

"""
:Author: Simon Gibbons (simongibbons@gmail.com)
"""
class SimpleRSTHeader(FixedWidthHeader):
    position_line = ...
    start_line = ...
    splitter_class = DefaultSplitter
    position_char = ...
    def get_fixedwidth_params(self, line): # -> tuple[list[Unknown], list[Unknown] | list[int], list[Unknown | None] | list[Unknown]]:
        ...
    


class SimpleRSTData(FixedWidthData):
    start_line = ...
    end_line = ...
    splitter_class = FixedWidthTwoLineDataSplitter


class RST(FixedWidth):
    """reStructuredText simple format table.

    See: https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#simple-tables

    Example::

        ==== ===== ======
        Col1  Col2  Col3
        ==== ===== ======
          1    2.3  Hello
          2    4.5  Worlds
        ==== ===== ======

    Currently there is no support for reading tables which utilize continuation lines,
    or for ones which define column spans through the use of an additional
    line of dashes in the header.

    """
    _format_name = ...
    _description = ...
    data_class = SimpleRSTData
    header_class = SimpleRSTHeader
    def __init__(self) -> None:
        ...
    
    def write(self, lines): # -> list[Unknown]:
        ...
    


