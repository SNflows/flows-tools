"""
This type stub file was generated by pyright.
"""

import abc
from typing import Mapping, Optional, Type, TypeVar

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin"]
__doctest_requires__ = ...
_COSMOLOGY_CLASSES = ...
_CosmoT = TypeVar("_CosmoT", bound="Cosmology")
_FlatCosmoT = TypeVar("_FlatCosmoT", bound="FlatCosmologyMixin")
class CosmologyError(Exception):
    ...


class Cosmology(metaclass=abc.ABCMeta):
    """Base-class for all Cosmologies.

    Parameters
    ----------
    *args
        Arguments into the cosmology; used by subclasses, not this base class.
    name : str or None (optional, keyword-only)
        The name of the cosmology.
    meta : dict or None (optional, keyword-only)
        Metadata for the cosmology, e.g., a reference.
    **kwargs
        Arguments into the cosmology; used by subclasses, not this base class.

    Notes
    -----
    Class instances are static -- you cannot (and should not) change the values
    of the parameters.  That is, all of the above attributes (except meta) are
    read only.

    For details on how to create performant custom subclasses, see the
    documentation on :ref:`astropy-cosmology-fast-integrals`.
    """
    meta = ...
    from_format = ...
    to_format = ...
    read = ...
    write = ...
    __parameters__ = ...
    __all_parameters__ = ...
    def __init_subclass__(cls): # -> None:
        ...
    
    def __init__(self, name=..., meta=...) -> None:
        ...
    
    @property
    def name(self): # -> str | None:
        """The name of the Cosmology instance."""
        ...
    
    @property
    @abc.abstractmethod
    def is_flat(self):
        """
        Return bool; `True` if the cosmology is flat.
        This is abstract and must be defined in subclasses.
        """
        ...
    
    def clone(self, *, meta=..., **kwargs): # -> Self@Cosmology:
        """Returns a copy of this object with updated parameters, as specified.

        This cannot be used to change the type of the cosmology, so ``clone()``
        cannot be used to change between flat and non-flat cosmologies.

        Parameters
        ----------
        meta : mapping or None (optional, keyword-only)
            Metadata that will update the current metadata.
        **kwargs
            Cosmology parameter (and name) modifications. If any parameter is
            changed and a new name is not given, the name will be set to "[old
            name] (modified)".

        Returns
        -------
        newcosmo : `~astropy.cosmology.Cosmology` subclass instance
            A new instance of this class with updated parameters as specified.
            If no arguments are given, then a reference to this object is
            returned instead of copy.

        Examples
        --------
        To make a copy of the ``Planck13`` cosmology with a different matter
        density (``Om0``), and a new name:

            >>> from astropy.cosmology import Planck13
            >>> Planck13.clone(name="Modified Planck 2013", Om0=0.35)
            FlatLambdaCDM(name="Modified Planck 2013", H0=67.77 km / (Mpc s),
                          Om0=0.35, ...

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'
        """
        ...
    
    def is_equivalent(self, other, *, format=...): # -> _NotImplementedType | bool:
        r"""Check equivalence between Cosmologies.

        Two cosmologies may be equivalent even if not the same class.
        For example, an instance of ``LambdaCDM`` might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like ``FlatLambdaCDM``.

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object in which to compare.
        format : bool or None or str, optional keyword-only
            Whether to allow, before equivalence is checked, the object to be
            converted to a |Cosmology|. This allows, e.g. a |Table| to be
            equivalent to a Cosmology.
            `False` (default) will not allow conversion. `True` or `None` will,
            and will use the auto-identification to try to infer the correct
            format. A `str` is assumed to be the correct format to use when
            converting.

        Returns
        -------
        bool
            True if cosmologies are equivalent, False otherwise.

        Examples
        --------
        Two cosmologies may be equivalent even if not of the same class.
        In this examples the ``LambdaCDM`` has ``Ode0`` set to the same value
        calculated in ``FlatLambdaCDM``.

            >>> import astropy.units as u
            >>> from astropy.cosmology import LambdaCDM, FlatLambdaCDM
            >>> cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
            >>> cosmo2 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3)
            >>> cosmo1.is_equivalent(cosmo2)
            True

        While in this example, the cosmologies are not equivalent.

            >>> cosmo3 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, Tcmb0=3 * u.K)
            >>> cosmo3.is_equivalent(cosmo2)
            False

        Also, using the keyword argument, the notion of equivalence is extended
        to any Python object that can be converted to a |Cosmology|.

            >>> from astropy.cosmology import Planck18
            >>> tbl = Planck18.to_format("astropy.table")
            >>> Planck18.is_equivalent(tbl, format=True)
            True

        The list of valid formats, e.g. the |Table| in this example, may be
        checked with ``Cosmology.from_format.list_formats()``.

        As can be seen in the list of formats, not all formats can be
        auto-identified by ``Cosmology.from_format.registry``. Objects of
        these kinds can still be checked for equivalence, but the correct
        format string must be used.

            >>> tbl = Planck18.to_format("yaml")
            >>> Planck18.is_equivalent(tbl, format="yaml")
            True
        """
        ...
    
    def __equiv__(self, other): # -> _NotImplementedType | bool:
        """Cosmology equivalence. Use ``.is_equivalent()`` for actual check!

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object in which to compare.

        Returns
        -------
        bool or `NotImplemented`
            `NotImplemented` if 'other' is from a different class.
            `True` if 'other' is of the same class and has matching parameters
            and parameter values. `False` otherwise.
        """
        ...
    
    def __eq__(self, other) -> bool:
        """Check equality between Cosmologies.

        Checks the Parameters and immutable fields (i.e. not "meta").

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object in which to compare.

        Returns
        -------
        bool
            `True` if Parameters and names are the same, `False` otherwise.
        """
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __astropy_table__(self, cls, copy, **kwargs):
        """Return a `~astropy.table.Table` of type ``cls``.

        Parameters
        ----------
        cls : type
            Astropy ``Table`` class or subclass.
        copy : bool
            Ignored.
        **kwargs : dict, optional
            Additional keyword arguments. Passed to ``self.to_format()``.
            See ``Cosmology.to_format.help("astropy.table")`` for allowed kwargs.

        Returns
        -------
        `astropy.table.Table` or subclass instance
            Instance of type ``cls``.
        """
        ...
    


class FlatCosmologyMixin(metaclass=abc.ABCMeta):
    """
    Mixin class for flat cosmologies. Do NOT instantiate directly.
    Note that all instances of ``FlatCosmologyMixin`` are flat, but not all
    flat cosmologies are instances of ``FlatCosmologyMixin``. As example,
    ``LambdaCDM`` **may** be flat (for the a specific set of parameter values),
    but ``FlatLambdaCDM`` **will** be flat.
    """
    def __init_subclass__(cls: Type[_FlatCosmoT]) -> None:
        ...
    
    _nonflat_cls_ = ...
    @property
    def is_flat(self): # -> Literal[True]:
        """Return `True`, the cosmology is flat."""
        ...
    
    @abc.abstractmethod
    def nonflat(self: _FlatCosmoT) -> _CosmoT:
        """Return the equivalent non-flat-class instance of this cosmology."""
        ...
    
    def clone(self, *, meta: Optional[Mapping] = ..., to_nonflat: bool = ..., **kwargs): # -> Any:
        """Returns a copy of this object with updated parameters, as specified.

        This cannot be used to change the type of the cosmology, except for
        changing to the non-flat version of this cosmology.

        Parameters
        ----------
        meta : mapping or None (optional, keyword-only)
            Metadata that will update the current metadata.
        to_nonflat : bool, optional keyword-only
            Whether to change to the non-flat version of this cosmology.
        **kwargs
            Cosmology parameter (and name) modifications. If any parameter is
            changed and a new name is not given, the name will be set to "[old
            name] (modified)".

        Returns
        -------
        newcosmo : `~astropy.cosmology.Cosmology` subclass instance
            A new instance of this class with updated parameters as specified.
            If no arguments are given, then a reference to this object is
            returned instead of copy.

        Examples
        --------
        To make a copy of the ``Planck13`` cosmology with a different matter
        density (``Om0``), and a new name:

            >>> from astropy.cosmology import Planck13
            >>> Planck13.clone(name="Modified Planck 2013", Om0=0.35)
            FlatLambdaCDM(name="Modified Planck 2013", H0=67.77 km / (Mpc s),
                          Om0=0.35, ...

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'

        The keyword 'to_nonflat' can be used to clone on the non-flat equivalent
        cosmology.

            >>> Planck13.clone(to_nonflat=True)
            LambdaCDM(name="Planck13", ...

            >>> Planck13.clone(H0=70, to_nonflat=True)
            LambdaCDM(name="Planck13 (modified)", H0=70.0 km / (Mpc s), ...
        """
        ...
    


def __getattr__(attr): # -> Any:
    ...
