"""
This type stub file was generated by pyright.
"""

from astropy import units as u

"""
This module contains the classes and utility functions for distance and
cartesian coordinates.
"""
__all__ = ['Distance']
__doctest_requires__ = ...
class Distance(u.SpecificTypeQuantity):
    """
    A one-dimensional distance.

    This can be initialized by providing one of the following:

    * Distance ``value`` (array or float) and a ``unit``
    * |Quantity| object with dimensionality of length
    * Redshift and (optionally) a `~astropy.cosmology.Cosmology`
    * Distance modulus
    * Parallax

    Parameters
    ----------
    value : scalar or `~astropy.units.Quantity` ['length']
        The value of this distance.
    unit : `~astropy.units.UnitBase` ['length']
        The unit for this distance.
    z : float
        A redshift for this distance.  It will be converted to a distance
        by computing the luminosity distance for this redshift given the
        cosmology specified by ``cosmology``. Must be given as a keyword
        argument.
    cosmology : `~astropy.cosmology.Cosmology` or None
        A cosmology that will be used to compute the distance from ``z``.
        If `None`, the current cosmology will be used (see
        `astropy.cosmology` for details).
    distmod : float or `~astropy.units.Quantity`
        The distance modulus for this distance. Note that if ``unit`` is not
        provided, a guess will be made at the unit between AU, pc, kpc, and Mpc.
    parallax : `~astropy.units.Quantity` or `~astropy.coordinates.Angle`
        The parallax in angular units.
    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.
    copy : bool, optional
        See `~astropy.units.Quantity`.
    order : {'C', 'F', 'A'}, optional
        See `~astropy.units.Quantity`.
    subok : bool, optional
        See `~astropy.units.Quantity`.
    ndmin : int, optional
        See `~astropy.units.Quantity`.
    allow_negative : bool, optional
        Whether to allow negative distances (which are possible in some
        cosmologies).  Default: `False`.

    Raises
    ------
    `~astropy.units.UnitsError`
        If the ``unit`` is not a length unit.
    ValueError
        If value specified is less than 0 and ``allow_negative=False``.

        If ``cosmology`` is provided when ``z`` is *not* given.

        If either none or more than one of ``value``, ``z``, ``distmod``,
        or ``parallax`` were given.


    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.cosmology import WMAP5
    >>> Distance(10, u.Mpc)
    <Distance 10. Mpc>
    >>> Distance(40*u.pc, unit=u.kpc)
    <Distance 0.04 kpc>
    >>> Distance(z=0.23)                      # doctest: +FLOAT_CMP
    <Distance 1184.01657566 Mpc>
    >>> Distance(z=0.23, cosmology=WMAP5)     # doctest: +FLOAT_CMP
    <Distance 1147.78831918 Mpc>
    >>> Distance(distmod=24.47*u.mag)         # doctest: +FLOAT_CMP
    <Distance 783.42964277 kpc>
    >>> Distance(parallax=21.34*u.mas)        # doctest: +FLOAT_CMP
    <Distance 46.86035614 pc>
    """
    _equivalent_unit = ...
    _include_easy_conversion_members = ...
    def __new__(cls, value=..., unit=..., z=..., cosmology=..., distmod=..., parallax=..., dtype=..., copy=..., order=..., subok=..., ndmin=..., allow_negative=...):
        ...
    
    @property
    def z(self): # -> Any:
        """Short for ``self.compute_z()``"""
        ...
    
    def compute_z(self, cosmology=..., **atzkw): # -> Any:
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` or None
            The cosmology to assume for this calculation, or `None` to use the
            current cosmology (see `astropy.cosmology` for details).
        **atzkw
            keyword arguments for :func:`~astropy.cosmology.z_at_value`

        Returns
        -------
        z : `~astropy.units.Quantity`
            The redshift of this distance given the provided ``cosmology``.

        Warnings
        --------
        This method can be slow for large arrays.
        The redshift is determined using :func:`astropy.cosmology.z_at_value`,
        which handles vector inputs (e.g. an array of distances) by
        element-wise calling of :func:`scipy.optimize.minimize_scalar`.
        For faster results consider using an interpolation table;
        :func:`astropy.cosmology.z_at_value` provides details.

        See Also
        --------
        :func:`astropy.cosmology.z_at_value` : Find the redshift corresponding to a
            :meth:`astropy.cosmology.FLRW.luminosity_distance`.
        """
        ...
    
    @property
    def distmod(self): # -> Quantity:
        """The distance modulus as a `~astropy.units.Quantity`"""
        ...
    
    @property
    def parallax(self): # -> Angle:
        """The parallax angle as an `~astropy.coordinates.Angle` object"""
        ...
    


