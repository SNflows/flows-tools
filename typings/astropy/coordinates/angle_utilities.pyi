"""
This type stub file was generated by pyright.
"""

"""
This module contains utility functions for working with angles. These are both
used internally in astropy.coordinates.angles, and of possible
"""
__all__ = ['angular_separation', 'position_angle', 'offset_by', 'golden_spiral_grid', 'uniform_spherical_random_surface', 'uniform_spherical_random_volume']
def angular_separation(lon1, lat1, lon2, lat2): # -> Any:
    """
    Angular separation between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    angular separation : `~astropy.units.Quantity` ['angle'] or float
        Type depends on input; ``Quantity`` in angular units, or float in
        radians.

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1]_,
    which is slightly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    .. [1] https://en.wikipedia.org/wiki/Great-circle_distance
    """
    ...

def position_angle(lon1, lat1, lon2, lat2): # -> Angle | None:
    """
    Position Angle (East of North) between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    pa : `~astropy.coordinates.Angle`
        The (positive) position angle of the vector pointing from position 1 to
        position 2.  If any of the angles are arrays, this will contain an array
        following the appropriate `numpy` broadcasting rules.

    """
    ...

def offset_by(lon, lat, posang, distance): # -> tuple[Unknown, Self@Quantity | Quantity | Any | Unknown]:
    """
    Point with the given offset from the given point.

    Parameters
    ----------
    lon, lat, posang, distance : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the starting point,
        position angle and distance to the final point.
        Quantities should be in angular units; floats in radians.
        Polar points at lat= +/-90 are treated as limit of +/-(90-epsilon) and same lon.

    Returns
    -------
    lon, lat : `~astropy.coordinates.Angle`
        The position of the final point.  If any of the angles are arrays,
        these will contain arrays following the appropriate `numpy` broadcasting rules.
        0 <= lon < 2pi.
    """
    ...

def golden_spiral_grid(size): # -> UnitSphericalRepresentation:
    """Generate a grid of points on the surface of the unit sphere using the
    Fibonacci or Golden Spiral method.

    .. seealso::

        `Evenly distributing points on a sphere <https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere>`_

    Parameters
    ----------
    size : int
        The number of points to generate.

    Returns
    -------
    rep : `~astropy.coordinates.UnitSphericalRepresentation`
        The grid of points.
    """
    ...

def uniform_spherical_random_surface(size=...): # -> UnitSphericalRepresentation:
    """Generate a random sampling of points on the surface of the unit sphere.

    Parameters
    ----------
    size : int
        The number of points to generate.

    Returns
    -------
    rep : `~astropy.coordinates.UnitSphericalRepresentation`
        The random points.
    """
    ...

def uniform_spherical_random_volume(size=..., max_radius=...): # -> SphericalRepresentation:
    """Generate a random sampling of points that follow a uniform volume
    density distribution within a sphere.

    Parameters
    ----------
    size : int
        The number of points to generate.
    max_radius : number, quantity-like, optional
        A dimensionless or unit-ful factor to scale the random distances.

    Returns
    -------
    rep : `~astropy.coordinates.SphericalRepresentation`
        The random points.
    """
    ...

__old_angle_utilities_funcs = ...
