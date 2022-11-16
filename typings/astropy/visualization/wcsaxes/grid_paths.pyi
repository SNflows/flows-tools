"""
This type stub file was generated by pyright.
"""

ROUND_TRIP_RTOL = ...
DISCONT_FACTOR = ...
def get_lon_lat_path(lon_lat, pixel, lon_lat_check):
    """
    Draw a curve, taking into account discontinuities.

    Parameters
    ----------
    lon_lat : ndarray
        The longitude and latitude values along the curve, given as a (n,2)
        array.
    pixel : ndarray
        The pixel coordinates corresponding to ``lon_lat``
    lon_lat_check : ndarray
        The world coordinates derived from converting from ``pixel``, which is
        used to ensure round-tripping.
    """
    ...

def get_gridline_path(world, pixel):
    """
    Draw a grid line

    Parameters
    ----------
    world : ndarray
        The longitude and latitude values along the curve, given as a (n,2)
        array.
    pixel : ndarray
        The pixel coordinates corresponding to ``lon_lat``
    """
    ...

