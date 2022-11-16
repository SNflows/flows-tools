"""
This type stub file was generated by pyright.
"""

from .core import *
from .coordinate_helpers import CoordinateHelper
from .coordinates_map import CoordinatesMap
from .patches import *
from astropy import config as _config

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.visualization.wcsaxes`.
    """
    coordinate_range_samples = ...
    frame_boundary_samples = ...
    grid_samples = ...
    contour_grid_samples = ...


conf = ...
