"""
This type stub file was generated by pyright.
"""

"""
This module contains convenience functions for coordinate-related functionality.

This is generally just wrapping around the object-oriented coordinates
framework, but it is useful for some users who are used to more functional
interfaces.
"""
__all__ = ['cartesian_to_spherical', 'spherical_to_cartesian', 'get_sun', 'get_constellation', 'concatenate_representations', 'concatenate']
def cartesian_to_spherical(x, y, z): # -> tuple[Unknown | Distance, Unknown, Unknown]:
    """
    Converts 3D rectangular cartesian coordinates to spherical polar
    coordinates.

    Note that the resulting angles are latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This function simply wraps functionality provided by the
        `~astropy.coordinates.CartesianRepresentation` and
        `~astropy.coordinates.SphericalRepresentation` classes.  In general,
        for both performance and readability, we suggest using these classes
        directly.  But for situations where a quick one-off conversion makes
        sense, this function is provided.

    Parameters
    ----------
    x : scalar, array-like, or `~astropy.units.Quantity`
        The first Cartesian coordinate.
    y : scalar, array-like, or `~astropy.units.Quantity`
        The second Cartesian coordinate.
    z : scalar, array-like, or `~astropy.units.Quantity`
        The third Cartesian coordinate.

    Returns
    -------
    r : `~astropy.units.Quantity`
        The radial coordinate (in the same units as the inputs).
    lat : `~astropy.units.Quantity` ['angle']
        The latitude in radians
    lon : `~astropy.units.Quantity` ['angle']
        The longitude in radians
    """
    ...

def spherical_to_cartesian(r, lat, lon): # -> tuple[Unknown, Unknown, Unknown]:
    """
    Converts spherical polar coordinates to rectangular cartesian
    coordinates.

    Note that the input angles should be in latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This is a low-level function used internally in
        `astropy.coordinates`.  It is provided for users if they really
        want to use it, but it is recommended that you use the
        `astropy.coordinates` coordinate systems.

    Parameters
    ----------
    r : scalar, array-like, or `~astropy.units.Quantity`
        The radial coordinate (in the same units as the inputs).
    lat : scalar, array-like, or `~astropy.units.Quantity` ['angle']
        The latitude (in radians if array or scalar)
    lon : scalar, array-like, or `~astropy.units.Quantity` ['angle']
        The longitude (in radians if array or scalar)

    Returns
    -------
    x : float or array
        The first cartesian coordinate.
    y : float or array
        The second cartesian coordinate.
    z : float or array
        The third cartesian coordinate.

    """
    ...

def get_sun(time): # -> SkyCoord:
    """
    Determines the location of the sun at a given time (or times, if the input
    is an array `~astropy.time.Time` object), in geocentric coordinates.

    Parameters
    ----------
    time : `~astropy.time.Time`
        The time(s) at which to compute the location of the sun.

    Returns
    -------
    newsc : `~astropy.coordinates.SkyCoord`
        The location of the sun as a `~astropy.coordinates.SkyCoord` in the
        `~astropy.coordinates.GCRS` frame.

    Notes
    -----
    The algorithm for determining the sun/earth relative position is based
    on the simplified version of VSOP2000 that is part of ERFA. Compared to
    JPL's ephemeris, it should be good to about 4 km (in the Sun-Earth
    vector) from 1900-2100 C.E., 8 km for the 1800-2200 span, and perhaps
    250 km over the 1000-3000.

    """
    ...

_constellation_data = ...
def get_constellation(coord, short_name=..., constellation_list=...): # -> Any:
    """
    Determines the constellation(s) a given coordinate object contains.

    Parameters
    ----------
    coord : coordinate-like
        The object to determine the constellation of.
    short_name : bool
        If True, the returned names are the IAU-sanctioned abbreviated
        names.  Otherwise, full names for the constellations are used.
    constellation_list : str
        The set of constellations to use.  Currently only ``'iau'`` is
        supported, meaning the 88 "modern" constellations endorsed by the IAU.

    Returns
    -------
    constellation : str or string array
        If ``coords`` contains a scalar coordinate, returns the name of the
        constellation.  If it is an array coordinate object, it returns an array
        of names.

    Notes
    -----
    To determine which constellation a point on the sky is in, this precesses
    to B1875, and then uses the Delporte boundaries of the 88 modern
    constellations, as tabulated by
    `Roman 1987 <http://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/42>`_.
    """
    ...

def concatenate_representations(reps): # -> Any:
    """
    Combine multiple representation objects into a single instance by
    concatenating the data in each component.

    Currently, all of the input representations have to be the same type. This
    properly handles differential or velocity data, but all input objects must
    have the same differential object type as well.

    Parameters
    ----------
    reps : sequence of `~astropy.coordinates.BaseRepresentation`
        The objects to concatenate

    Returns
    -------
    rep : `~astropy.coordinates.BaseRepresentation` subclass instance
        A single representation object with its data set to the concatenation of
        all the elements of the input sequence of representations.

    """
    ...

def concatenate(coords): # -> SkyCoord:
    """
    Combine multiple coordinate objects into a single
    `~astropy.coordinates.SkyCoord`.

    "Coordinate objects" here mean frame objects with data,
    `~astropy.coordinates.SkyCoord`, or representation objects.  Currently,
    they must all be in the same frame, but in a future version this may be
    relaxed to allow inhomogeneous sequences of objects.

    Parameters
    ----------
    coords : sequence of coordinate-like
        The objects to concatenate

    Returns
    -------
    cskycoord : SkyCoord
        A single sky coordinate with its data set to the concatenation of all
        the elements in ``coords``
    """
    ...
