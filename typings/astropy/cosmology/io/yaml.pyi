"""
This type stub file was generated by pyright.
"""

"""
The following are private functions, included here **FOR REFERENCE ONLY** since
the io registry cannot be displayed. These functions are registered into
:meth:`~astropy.cosmology.Cosmology.to_format` and
:meth:`~astropy.cosmology.Cosmology.from_format` and should only be accessed
via these methods.
"""
__all__ = []
def yaml_representer(tag): # -> (dumper: Unknown, obj: Unknown) -> Unknown:
    """:mod:`yaml` representation of |Cosmology| object.

    Parameters
    ----------
    tag : str
        The class tag, e.g. '!astropy.cosmology.LambdaCDM'

    Returns
    -------
    representer : callable[[`~astropy.io.misc.yaml.AstropyDumper`, |Cosmology|], str]
        Function to construct :mod:`yaml` representation of |Cosmology| object.
    """
    ...

def yaml_constructor(cls): # -> (loader: Unknown, node: Unknown) -> Unknown:
    """Cosmology| object from :mod:`yaml` representation.

    Parameters
    ----------
    cls : type
        The class type, e.g. `~astropy.cosmology.LambdaCDM`.

    Returns
    -------
    constructor : callable
        Function to construct |Cosmology| object from :mod:`yaml` representation.
    """
    ...

def register_cosmology_yaml(cosmo_cls): # -> None:
    """Register :mod:`yaml` for Cosmology class.

    Parameters
    ----------
    cosmo_cls : `~astropy.cosmology.Cosmology` class
    """
    ...

def from_yaml(yml, *, cosmology=...): # -> Any:
    """Load `~astropy.cosmology.Cosmology` from :mod:`yaml` object.

    Parameters
    ----------
    yml : str
        :mod:`yaml` representation of |Cosmology| object
    cosmology : str, `~astropy.cosmology.Cosmology` class, or None (optional, keyword-only)
        The expected cosmology class (or string name thereof). This argument is
        is only checked for correctness if not `None`.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Raises
    ------
    TypeError
        If the |Cosmology| object loaded from ``yml`` is not an instance of
        the ``cosmology`` (and ``cosmology`` is not `None`).
    """
    ...

def to_yaml(cosmology, *args): # -> _Yaml:
    """Return the cosmology class, parameters, and metadata as a :mod:`yaml` object.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`

    Returns
    -------
    str
        :mod:`yaml` representation of |Cosmology| object
    """
    ...
