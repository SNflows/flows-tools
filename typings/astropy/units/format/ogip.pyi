"""
This type stub file was generated by pyright.
"""

from . import generic

"""
Handles units in `Office of Guest Investigator Programs (OGIP)
FITS files
<https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
"""
class OGIP(generic.Generic):
    """
    Support the units in `Office of Guest Investigator Programs (OGIP)
    FITS files
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
    """
    _tokens = ...
    @classmethod
    def parse(cls, s, debug=...): # -> Any:
        ...
    
    @classmethod
    def to_string(cls, unit): # -> LiteralString | None:
        ...
    


