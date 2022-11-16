"""
This type stub file was generated by pyright.
"""

from . import generic

"""
Handles the "VOUnit" unit format.
"""
class VOUnit(generic.Generic):
    """
    The IVOA standard for units used by the VO.

    This is an implementation of `Units in the VO 1.0
    <http://www.ivoa.net/documents/VOUnits/>`_.
    """
    _explicit_custom_unit_regex = ...
    _custom_unit_regex = ...
    _custom_units = ...
    @classmethod
    def parse(cls, s, debug=...): # -> Any | None:
        ...
    
    @classmethod
    def to_string(cls, unit): # -> str | bytes:
        ...
    

