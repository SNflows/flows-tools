"""
This type stub file was generated by pyright.
"""

"""
This package defines colloquially used Imperial units.  They are
available in the `astropy.units.imperial` namespace, but not in the
top-level `astropy.units` namespace, e.g.::

    >>> import astropy.units as u
    >>> mph = u.imperial.mile / u.hour
    >>> mph
    Unit("mi / h")

To include them in `~astropy.units.UnitBase.compose` and the results of
`~astropy.units.UnitBase.find_equivalent_units`, do::

    >>> import astropy.units as u
    >>> u.imperial.enable()  # doctest: +SKIP
"""
_ns = ...
if __doc__ is not None:
    ...
def enable(): # -> _UnitContext:
    """
    Enable Imperial units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable Imperial
    units only temporarily.
    """
    ...

