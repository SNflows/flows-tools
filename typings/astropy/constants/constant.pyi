"""
This type stub file was generated by pyright.
"""

from astropy.units.quantity import Quantity

__all__ = ['Constant', 'EMConstant']
class ConstantMeta(type):
    """Metaclass for `~astropy.constants.Constant`. The primary purpose of this
    is to wrap the double-underscore methods of `~astropy.units.Quantity`
    which is the superclass of `~astropy.constants.Constant`.

    In particular this wraps the operator overloads such as `__add__` to
    prevent their use with constants such as ``e`` from being used in
    expressions without specifying a system.  The wrapper checks to see if the
    constant is listed (by name) in ``Constant._has_incompatible_units``, a set
    of those constants that are defined in different systems of units are
    physically incompatible.  It also performs this check on each `Constant` if
    it hasn't already been performed (the check is deferred until the
    `Constant` is actually used in an expression to speed up import times,
    among other reasons).
    """
    def __new__(mcls, name, bases, d): # -> Self@ConstantMeta:
        ...
    


class Constant(Quantity, metaclass=ConstantMeta):
    """A physical or astronomical constant.

    These objects are quantities that are meant to represent physical
    constants.
    """
    _registry = ...
    _has_incompatible_units = ...
    def __new__(cls, abbrev, name, value, unit, uncertainty, reference=..., system=...): # -> Self@Constant:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __str__(self) -> str:
        ...
    
    def __quantity_subclass__(self, unit): # -> tuple[Type[Quantity], Literal[False]]:
        ...
    
    def copy(self): # -> Self@Constant:
        """
        Return a copy of this `Constant` instance.  Since they are by
        definition immutable, this merely returns another reference to
        ``self``.
        """
        ...
    
    __deepcopy__ = ...
    @property
    def abbrev(self): # -> float | NDArray[void] | Any:
        """A typical ASCII text abbreviation of the constant, also generally
        the same as the Python variable used for this constant.
        """
        ...
    
    @property
    def name(self): # -> float | NDArray[void] | Any:
        """The full name of the constant."""
        ...
    
    @property
    def uncertainty(self): # -> float | NDArray[void] | Any:
        """The known absolute uncertainty in this constant's value."""
        ...
    
    @property
    def reference(self): # -> float | NDArray[void] | Any:
        """The source used for the value of this constant."""
        ...
    
    @property
    def system(self): # -> float | NDArray[void] | Any:
        """The system of units in which this constant is defined (typically
        `None` so long as the constant's units can be directly converted
        between systems).
        """
        ...
    
    @property
    def si(self): # -> Any:
        """If the Constant is defined in the SI system return that instance of
        the constant, else convert to a Quantity in the appropriate SI units.
        """
        ...
    
    @property
    def cgs(self): # -> Any:
        """If the Constant is defined in the CGS system return that instance of
        the constant, else convert to a Quantity in the appropriate CGS units.
        """
        ...
    
    def __array_finalize__(self, obj): # -> None:
        ...
    


class EMConstant(Constant):
    """An electromagnetic constant."""
    @property
    def cgs(self):
        """Overridden for EMConstant to raise a `TypeError`
        emphasizing that there are multiple EM extensions to CGS.
        """
        ...
    


