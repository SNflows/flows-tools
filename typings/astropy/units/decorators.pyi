"""
This type stub file was generated by pyright.
"""

__all__ = ['quantity_input']
NoneType = ...
class QuantityInput:
    @classmethod
    def as_decorator(cls, func=..., **kwargs): # -> ((*func_args: Unknown, **func_kwargs: Unknown) -> Unknown) | Self@QuantityInput:
        r"""
        A decorator for validating the units of arguments to functions.

        Unit specifications can be provided as keyword arguments to the
        decorator, or by using function annotation syntax. Arguments to the
        decorator take precedence over any function annotations present.

        A `~astropy.units.UnitsError` will be raised if the unit attribute of
        the argument is not equivalent to the unit specified to the decorator or
        in the annotation. If the argument has no unit attribute, i.e. it is not
        a Quantity object, a `ValueError` will be raised unless the argument is
        an annotation. This is to allow non Quantity annotations to pass
        through.

        Where an equivalency is specified in the decorator, the function will be
        executed with that equivalency in force.

        Notes
        -----

        The checking of arguments inside variable arguments to a function is not
        supported (i.e. \*arg or \**kwargs).

        The original function is accessible by the attributed ``__wrapped__``.
        See :func:`functools.wraps` for details.

        Examples
        --------

        .. code-block:: python

            import astropy.units as u
            @u.quantity_input(myangle=u.arcsec)
            def myfunction(myangle):
                return myangle**2


        .. code-block:: python

            import astropy.units as u
            @u.quantity_input
            def myfunction(myangle: u.arcsec):
                return myangle**2

        Or using a unit-aware Quantity annotation.

        .. code-block:: python

            @u.quantity_input
            def myfunction(myangle: u.Quantity[u.arcsec]):
                return myangle**2

        Also you can specify a return value annotation, which will
        cause the function to always return a `~astropy.units.Quantity` in that
        unit.

        .. code-block:: python

            import astropy.units as u
            @u.quantity_input
            def myfunction(myangle: u.arcsec) -> u.deg**2:
                return myangle**2

        Using equivalencies::

            import astropy.units as u
            @u.quantity_input(myenergy=u.eV, equivalencies=u.mass_energy())
            def myfunction(myenergy):
                return myenergy**2

        """
        ...
    
    def __init__(self, func=..., strict_dimensionless=..., **kwargs) -> None:
        ...
    
    def __call__(self, wrapped_function): # -> (*func_args: Unknown, **func_kwargs: Unknown) -> Unknown:
        ...
    


quantity_input = ...
