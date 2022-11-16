"""
This type stub file was generated by pyright.
"""

from .core import AllType, BaseData, BaseHeader, BaseInputter, BaseOutputter, BaseReader, BaseSplitter, Column, ContinuationLinesInputter, DefaultSplitter, FloatType, InconsistentTableError, IntType, NoType, NumType, ParameterError, StrType, TableOutputter, WhitespaceSplitter, convert_numpy, masked
from .basic import Basic, BasicData, BasicHeader, CommentedHeader, Csv, NoHeader, Rdb, Tab
from .fastbasic import FastBasic, FastCommentedHeader, FastCsv, FastNoHeader, FastRdb, FastTab
from .cds import Cds
from .mrt import Mrt
from .ecsv import Ecsv
from .latex import AASTex, Latex, latexdicts
from .html import HTML
from .ipac import Ipac
from .daophot import Daophot
from .qdp import QDP
from .sextractor import SExtractor
from .fixedwidth import FixedWidth, FixedWidthData, FixedWidthHeader, FixedWidthNoHeader, FixedWidthSplitter, FixedWidthTwoLine
from .rst import RST
from .ui import get_read_trace, get_reader, get_writer, read, set_guess, write
from . import connect

""" An extensible ASCII table reader and writer.

"""
