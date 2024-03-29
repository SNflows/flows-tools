"""
This type stub file was generated by pyright.
"""

from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import DynamicMatrixTransform
from .fk5 import FK5
from .fk4 import FK4NoETerms
from .galactic import Galactic

@frame_transform_graph.transform(DynamicMatrixTransform, FK5, Galactic)
def fk5_to_gal(fk5coord, galframe): # -> Any:
    ...

@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, Galactic)
def fk4_to_gal(fk4coords, galframe): # -> Any:
    ...

@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK4NoETerms)
def gal_to_fk4(galcoords, fk4frame): # -> Any:
    ...

