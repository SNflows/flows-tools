"""
This type stub file was generated by pyright.
"""

from .base import Base

"""
Handles a "generic" string format for units
"""
class Generic(Base):
    """
    A "generic" format.

    The syntax of the format is based directly on the FITS standard,
    but instead of only supporting the units that FITS knows about, it
    supports any unit available in the `astropy.units` namespace.
    """
    _show_scale = ...
    _tokens = ...
    _unit_symbols = ...
    _prefixable_unit_symbols = ...
    _unit_suffix_symbols = ...
    _translations = ...
    _superscripts = ...
    _superscript_translations = ...
    _regex_superscript = ...
    _regex_deg = ...
    @classmethod
    def parse(cls, s, debug=...): # -> Any:
        ...
    
    @classmethod
    def to_string(cls, unit): # -> LiteralString | None:
        ...
    


class Unscaled(Generic):
    """
    A format that doesn't display the scale part of the unit, other
    than that, it is identical to the `Generic` format.

    This is used in some error messages where the scale is irrelevant.
    """
    _show_scale = ...


