"""
This type stub file was generated by pyright.
"""

from .core import FittableModel, Model

"""
Special models useful for complex compound models where control is needed over
which outputs from a source model are mapped to which inputs of a target model.
"""
__all__ = ['Mapping', 'Identity', 'UnitsMapping']
class Mapping(FittableModel):
    """
    Allows inputs to be reordered, duplicated or dropped.

    Parameters
    ----------
    mapping : tuple
        A tuple of integers representing indices of the inputs to this model
        to return and in what order to return them.  See
        :ref:`astropy:compound-model-mappings` for more details.
    n_inputs : int
        Number of inputs; if `None` (default) then ``max(mapping) + 1`` is
        used (i.e. the highest input index used in the mapping).
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).
    meta : dict-like
        Free-form metadata to associate with this model.

    Raises
    ------
    TypeError
        Raised when number of inputs is less that ``max(mapping)``.

    Examples
    --------

    >>> from astropy.modeling.models import Polynomial2D, Shift, Mapping
    >>> poly1 = Polynomial2D(1, c0_0=1, c1_0=2, c0_1=3)
    >>> poly2 = Polynomial2D(1, c0_0=1, c1_0=2.4, c0_1=2.1)
    >>> model = (Shift(1) & Shift(2)) | Mapping((0, 1, 0, 1)) | (poly1 & poly2)
    >>> model(1, 2)  # doctest: +FLOAT_CMP
    (17.0, 14.2)
    """
    linear = ...
    def __init__(self, mapping, n_inputs=..., name=..., meta=...) -> None:
        ...
    
    @property
    def n_inputs(self):
        ...
    
    @property
    def n_outputs(self): # -> int:
        ...
    
    @property
    def mapping(self): # -> Unknown:
        """Integers representing indices of the inputs."""
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def evaluate(self, *args): # -> tuple[Unknown, ...]:
        ...
    
    @property
    def inverse(self): # -> Self@Mapping:
        """
        A `Mapping` representing the inverse of the current mapping.

        Raises
        ------
        `NotImplementedError`
            An inverse does no exist on mappings that drop some of its inputs
            (there is then no way to reconstruct the inputs that were dropped).
        """
        ...
    


class Identity(Mapping):
    """
    Returns inputs unchanged.

    This class is useful in compound models when some of the inputs must be
    passed unchanged to the next model.

    Parameters
    ----------
    n_inputs : int
        Specifies the number of inputs this identity model accepts.
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).
    meta : dict-like
        Free-form metadata to associate with this model.

    Examples
    --------

    Transform ``(x, y)`` by a shift in x, followed by scaling the two inputs::

        >>> from astropy.modeling.models import (Polynomial1D, Shift, Scale,
        ...                                      Identity)
        >>> model = (Shift(1) & Identity(1)) | Scale(1.2) & Scale(2)
        >>> model(1,1)  # doctest: +FLOAT_CMP
        (2.4, 2.0)
        >>> model.inverse(2.4, 2) # doctest: +FLOAT_CMP
        (1.0, 1.0)
    """
    linear = ...
    def __init__(self, n_inputs, name=..., meta=...) -> None:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    @property
    def inverse(self): # -> Self@Identity:
        """
        The inverse transformation.

        In this case of `Identity`, ``self.inverse is self``.
        """
        ...
    


class UnitsMapping(Model):
    """
    Mapper that operates on the units of the input, first converting to
    canonical units, then assigning new units without further conversion.
    Used by Model.coerce_units to support units on otherwise unitless models
    such as Polynomial1D.

    Parameters
    ----------
    mapping : tuple
        A tuple of (input_unit, output_unit) pairs, one per input, matched to the
        inputs by position.  The first element of the each pair is the unit that
        the model will accept (specify ``dimensionless_unscaled``
        to accept dimensionless input).  The second element is the unit that the
        model will return.  Specify ``dimensionless_unscaled``
        to return dimensionless Quantity, and `None` to return raw values without
        Quantity.
    input_units_equivalencies : dict, optional
        Default equivalencies to apply to input values.  If set, this should be a
        dictionary where each key is a string that corresponds to one of the
        model inputs.
    input_units_allow_dimensionless : dict or bool, optional
        Allow dimensionless input. If this is True, input values to evaluate will
        gain the units specified in input_units. If this is a dictionary then it
        should map input name to a bool to allow dimensionless numbers for that
        input.
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).
    meta : dict-like, optional
        Free-form metadata to associate with this model.

    Examples
    --------

    Wrapping a unitless model to require and convert units:

    >>> from astropy.modeling.models import Polynomial1D, UnitsMapping
    >>> from astropy import units as u
    >>> poly = Polynomial1D(1, c0=1, c1=2)
    >>> model = UnitsMapping(((u.m, None),)) | poly
    >>> model = model | UnitsMapping(((None, u.s),))
    >>> model(u.Quantity(10, u.m))  # doctest: +FLOAT_CMP
    <Quantity 21. s>
    >>> model(u.Quantity(1000, u.cm)) # doctest: +FLOAT_CMP
    <Quantity 21. s>
    >>> model(u.Quantity(10, u.cm)) # doctest: +FLOAT_CMP
    <Quantity 1.2 s>

    Wrapping a unitless model but still permitting unitless input:

    >>> from astropy.modeling.models import Polynomial1D, UnitsMapping
    >>> from astropy import units as u
    >>> poly = Polynomial1D(1, c0=1, c1=2)
    >>> model = UnitsMapping(((u.m, None),), input_units_allow_dimensionless=True) | poly
    >>> model = model | UnitsMapping(((None, u.s),))
    >>> model(u.Quantity(10, u.m))  # doctest: +FLOAT_CMP
    <Quantity 21. s>
    >>> model(10)  # doctest: +FLOAT_CMP
    <Quantity 21. s>
    """
    def __init__(self, mapping, input_units_equivalencies=..., input_units_allow_dimensionless=..., name=..., meta=...) -> None:
        ...
    
    @property
    def n_inputs(self): # -> int:
        ...
    
    @property
    def n_outputs(self): # -> int:
        ...
    
    @property
    def inputs(self): # -> tuple[Literal['x']] | tuple[Literal['x'], Literal['y']] | tuple[str, ...] | tuple[()]:
        ...
    
    @inputs.setter
    def inputs(self, value): # -> None:
        ...
    
    @property
    def outputs(self): # -> tuple[Literal['y']] | tuple[Literal['z']] | tuple[str, ...] | tuple[()]:
        ...
    
    @outputs.setter
    def outputs(self, value): # -> None:
        ...
    
    @property
    def input_units(self): # -> dict[str, Unknown]:
        ...
    
    @property
    def mapping(self): # -> Unknown:
        ...
    
    def evaluate(self, *args): # -> tuple[Unknown, ...]:
        ...
    
    def __repr__(self): # -> str:
        ...
    


