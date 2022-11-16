"""
This type stub file was generated by pyright.
"""

from astropy.utils.decorators import format_doc
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc

__all__ = ['Supergalactic']
doc_components = ...
@format_doc(base_doc, components=doc_components, footer="")
class Supergalactic(BaseCoordinateFrame):
    """
    Supergalactic Coordinates
    (see Lahav et al. 2000, <https://ui.adsabs.harvard.edu/abs/2000MNRAS.312..166L>,
    and references therein).
    """
    frame_specific_representation_info = ...
    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential
    _nsgp_gal = ...


