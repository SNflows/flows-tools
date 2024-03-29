"""
This type stub file was generated by pyright.
"""

from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from .gcrs import GCRS, PrecessedGeocentric
from .cirs import CIRS
from .itrs import ITRS
from .equatorial import TEME, TETE

"""
Contains the transformation functions for getting to/from ITRS, TEME, GCRS, and CIRS.
These are distinct from the ICRS and AltAz functions because they are just
rotations without aberration corrections or offsets.
"""
def teme_to_itrs_mat(time):
    ...

def gcrs_to_cirs_mat(time):
    ...

def cirs_to_itrs_mat(time):
    ...

def tete_to_itrs_mat(time, rbpn=...):
    """Compute the polar motion p-matrix at the given time.

    If the nutation-precession matrix is already known, it should be passed in,
    as this is by far the most expensive calculation.
    """
    ...

def gcrs_precession_mat(equinox):
    ...

def get_location_gcrs(location, obstime, ref_to_itrs, gcrs_to_ref): # -> GCRS:
    """Create a GCRS frame at the location and obstime.

    The reference frame z axis must point to the Celestial Intermediate Pole
    (as is the case for CIRS and TETE).

    This function is here to avoid location.get_gcrs(obstime), which would
    recalculate matrices that are already available below (and return a GCRS
    coordinate, rather than a frame with obsgeoloc and obsgeovel).  Instead,
    it uses the private method that allows passing in the matrices.

    """
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, TETE)
def gcrs_to_tete(gcrs_coo, tete_frame):
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, GCRS)
def tete_to_gcrs(tete_coo, gcrs_frame): # -> Any:
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, ITRS)
def tete_to_itrs(tete_coo, itrs_frame):
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TETE)
def itrs_to_tete(itrs_coo, tete_frame): # -> Any:
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, CIRS)
def gcrs_to_cirs(gcrs_coo, cirs_frame):
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, GCRS)
def cirs_to_gcrs(cirs_coo, gcrs_frame): # -> Any:
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, ITRS)
def cirs_to_itrs(cirs_coo, itrs_frame):
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, CIRS)
def itrs_to_cirs(itrs_coo, cirs_frame): # -> Any:
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, PrecessedGeocentric)
def gcrs_to_precessedgeo(from_coo, to_frame):
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, PrecessedGeocentric, GCRS)
def precessedgeo_to_gcrs(from_coo, to_frame): # -> Any:
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TEME, ITRS)
def teme_to_itrs(teme_coo, itrs_frame): # -> Any:
    ...

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TEME)
def itrs_to_teme(itrs_coo, teme_frame):
    ...

