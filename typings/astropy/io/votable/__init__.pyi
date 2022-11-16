"""
This type stub file was generated by pyright.
"""

from .table import from_table, is_votable, parse, parse_single_table, validate, writeto
from .exceptions import IOWarning, UnimplementedWarning, VOTableChangeWarning, VOTableSpecError, VOTableSpecWarning, VOWarning
from astropy import config as _config

"""
This package reads and writes data formats used by the Virtual
Observatory (VO) initiative, particularly the VOTable XML format.
"""
__all__ = ['Conf', 'conf', 'parse', 'parse_single_table', 'validate', 'from_table', 'is_votable', 'writeto', 'VOWarning', 'VOTableChangeWarning', 'VOTableSpecWarning', 'UnimplementedWarning', 'IOWarning', 'VOTableSpecError']
class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.io.votable`.
    """
    verify = ...


conf = ...