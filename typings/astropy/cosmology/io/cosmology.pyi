"""
This type stub file was generated by pyright.
"""

"""
The following are private functions. These functions are registered into
:meth:`~astropy.cosmology.Cosmology.to_format` and
:meth:`~astropy.cosmology.Cosmology.from_format` and should only be accessed
via these methods.
"""
__all__ = []
def from_cosmology(cosmo, /, **kwargs):
    """Return the |Cosmology| unchanged.

    Parameters
    ----------
    cosmo : `~astropy.cosmology.Cosmology`
        The cosmology to return.
    **kwargs
        This argument is required for compatibility with the standard set of
        keyword arguments in format `~astropy.cosmology.Cosmology.from_format`,
        e.g. "cosmology". If "cosmology" is included and is not `None`,
        ``cosmo`` is checked for correctness.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance
        Just ``cosmo`` passed through.

    Raises
    ------
    TypeError
        If the |Cosmology| object is not an instance of ``cosmo`` (and
        ``cosmology`` is not `None`).
    """
    ...

def to_cosmology(cosmo, *args):
    """Return the |Cosmology| unchanged.

    Parameters
    ----------
    cosmo : `~astropy.cosmology.Cosmology`
        The cosmology to return.
    *args
        Not used.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance
        Just ``cosmo`` passed through.
    """
    ...

def cosmology_identify(origin, format, *args, **kwargs): # -> bool:
    """Identify if object is a `~astropy.cosmology.Cosmology`.

    Returns
    -------
    bool
    """
    ...
