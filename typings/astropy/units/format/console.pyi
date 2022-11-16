"""
This type stub file was generated by pyright.
"""

from . import base

"""
Handles the "Console" unit format.
"""
class Console(base.Base):
    """
    Output-only format for to display pretty formatting at the
    console.

    For example::

      >>> import astropy.units as u
      >>> print(u.Ry.decompose().to_string('console'))  # doctest: +FLOAT_CMP
                       m^2 kg
      2.1798721*10^-18 ------
                        s^2
    """
    _times = ...
    _line = ...
    @classmethod
    def format_exponential_notation(cls, val): # -> str:
        ...
    
    @classmethod
    def to_string(cls, unit): # -> str:
        ...
    


