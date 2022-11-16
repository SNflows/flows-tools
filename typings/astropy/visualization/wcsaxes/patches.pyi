"""
This type stub file was generated by pyright.
"""

from matplotlib.patches import Polygon

__all__ = ['Quadrangle', 'SphericalCircle']
class SphericalCircle(Polygon):
    """
    Create a patch representing a spherical circle - that is, a circle that is
    formed of all the points that are within a certain angle of the central
    coordinates on a sphere. Here we assume that latitude goes from -90 to +90

    This class is needed in cases where the user wants to add a circular patch
    to a celestial image, since otherwise the circle will be distorted, because
    a fixed interval in longitude corresponds to a different angle on the sky
    depending on the latitude.

    Parameters
    ----------
    center : tuple or `~astropy.units.Quantity` ['angle']
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements
        or a `~astropy.coordinates.SkyCoord` object.
    radius : `~astropy.units.Quantity` ['angle']
        The radius of the circle
    resolution : int, optional
        The number of points that make up the circle - increase this to get a
        smoother circle.
    vertex_unit : `~astropy.units.Unit`
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.

    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """
    def __init__(self, center, radius, resolution=..., vertex_unit=..., **kwargs) -> None:
        ...
    


class Quadrangle(Polygon):
    """
    Create a patch representing a latitude-longitude quadrangle.

    The edges of the quadrangle lie on two lines of constant longitude and two
    lines of constant latitude (or the equivalent component names in the
    coordinate frame of interest, such as right ascension and declination).
    Note that lines of constant latitude are not great circles.

    Unlike `matplotlib.patches.Rectangle`, the edges of this patch will render
    as curved lines if appropriate for the WCS transformation.

    Parameters
    ----------
    anchor : tuple or `~astropy.units.Quantity` ['angle']
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements.
    width : `~astropy.units.Quantity` ['angle']
        The width of the quadrangle in longitude (or, e.g., right ascension)
    height : `~astropy.units.Quantity` ['angle']
        The height of the quadrangle in latitude (or, e.g., declination)
    resolution : int, optional
        The number of points that make up each side of the quadrangle -
        increase this to get a smoother quadrangle.
    vertex_unit : `~astropy.units.Unit` ['angle']
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.

    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """
    def __init__(self, anchor, width, height, resolution=..., vertex_unit=..., **kwargs) -> None:
        ...
    


