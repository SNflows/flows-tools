"""
This type stub file was generated by pyright.
"""

from astropy.utils.decorators import format_doc
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc, frame_transform_graph
from astropy.coordinates.transformations import AffineTransform
from .baseradec import BaseRADecFrame, doc_components as doc_components_radec
from .icrs import ICRS
from .galactic import Galactic

J2000 = ...
v_bary_Schoenrich2010 = ...
__all__ = ['LSR', 'GalacticLSR', 'LSRK', 'LSRD']
doc_footer_lsr = ...
@format_doc(base_doc, components=doc_components_radec, footer=doc_footer_lsr)
class LSR(BaseRADecFrame):
    r"""A coordinate or frame in the Local Standard of Rest (LSR).

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSR. Roughly, the LSR is the mean
    velocity of the stars in the solar neighborhood, but the precise definition
    of which depends on the study. As defined in Schönrich et al. (2010):
    "The LSR is the rest frame at the location of the Sun of a star that would
    be on a circular orbit in the gravitational potential one would obtain by
    azimuthally averaging away non-axisymmetric features in the actual Galactic
    potential." No such orbit truly exists, but it is still a commonly used
    velocity frame.

    We use default values from Schönrich et al. (2010) for the barycentric
    velocity relative to the LSR, which is defined in Galactic (right-handed)
    cartesian velocity components
    :math:`(U, V, W) = (11.1, 12.24, 7.25)~{{\rm km}}~{{\rm s}}^{{-1}}`. These
    values are customizable via the ``v_bary`` argument which specifies the
    velocity of the solar system barycenter with respect to the LSR.

    The frame attributes are listed under **Other Parameters**.
    """
    v_bary = ...


@frame_transform_graph.transform(AffineTransform, ICRS, LSR)
def icrs_to_lsr(icrs_coord, lsr_frame): # -> tuple[None, CartesianRepresentation]:
    ...

@frame_transform_graph.transform(AffineTransform, LSR, ICRS)
def lsr_to_icrs(lsr_coord, icrs_frame): # -> tuple[None, CartesianRepresentation]:
    ...

doc_components_gal = ...
@format_doc(base_doc, components=doc_components_gal, footer=doc_footer_lsr)
class GalacticLSR(BaseCoordinateFrame):
    r"""A coordinate or frame in the Local Standard of Rest (LSR), axis-aligned
    to the `Galactic` frame.

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSR. Roughly, the LSR is the mean
    velocity of the stars in the solar neighborhood, but the precise definition
    of which depends on the study. As defined in Schönrich et al. (2010):
    "The LSR is the rest frame at the location of the Sun of a star that would
    be on a circular orbit in the gravitational potential one would obtain by
    azimuthally averaging away non-axisymmetric features in the actual Galactic
    potential." No such orbit truly exists, but it is still a commonly used
    velocity frame.

    We use default values from Schönrich et al. (2010) for the barycentric
    velocity relative to the LSR, which is defined in Galactic (right-handed)
    cartesian velocity components
    :math:`(U, V, W) = (11.1, 12.24, 7.25)~{{\rm km}}~{{\rm s}}^{{-1}}`. These
    values are customizable via the ``v_bary`` argument which specifies the
    velocity of the solar system barycenter with respect to the LSR.

    The frame attributes are listed under **Other Parameters**.
    """
    frame_specific_representation_info = ...
    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential
    v_bary = ...


@frame_transform_graph.transform(AffineTransform, Galactic, GalacticLSR)
def galactic_to_galacticlsr(galactic_coord, lsr_frame): # -> tuple[None, CartesianRepresentation]:
    ...

@frame_transform_graph.transform(AffineTransform, GalacticLSR, Galactic)
def galacticlsr_to_galactic(lsr_coord, galactic_frame): # -> tuple[None, CartesianRepresentation]:
    ...

class LSRK(BaseRADecFrame):
    r"""
    A coordinate or frame in the Kinematic Local Standard of Rest (LSR).

    This frame is defined as having a velocity of 20 km/s towards RA=270 Dec=30
    (B1900) relative to the solar system Barycenter. This is defined in:

        Gordon 1975, Methods of Experimental Physics: Volume 12:
        Astrophysics, Part C: Radio Observations - Section 6.1.5.

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSRK.
    """
    ...


V_OFFSET_LSRK = ...
ICRS_LSRK_OFFSET = ...
LSRK_ICRS_OFFSET = ...
@frame_transform_graph.transform(AffineTransform, ICRS, LSRK)
def icrs_to_lsrk(icrs_coord, lsr_frame): # -> tuple[None, CartesianRepresentation]:
    ...

@frame_transform_graph.transform(AffineTransform, LSRK, ICRS)
def lsrk_to_icrs(lsr_coord, icrs_frame): # -> tuple[None, CartesianRepresentation]:
    ...

class LSRD(BaseRADecFrame):
    r"""
    A coordinate or frame in the Dynamical Local Standard of Rest (LSRD)

    This frame is defined as a velocity of U=9 km/s, V=12 km/s,
    and W=7 km/s in Galactic coordinates or 16.552945 km/s
    towards l=53.13 b=25.02. This is defined in:

       Delhaye 1965, Solar Motion and Velocity Distribution of
       Common Stars.

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSRD.
    """
    ...


V_OFFSET_LSRD = ...
ICRS_LSRD_OFFSET = ...
LSRD_ICRS_OFFSET = ...
@frame_transform_graph.transform(AffineTransform, ICRS, LSRD)
def icrs_to_lsrd(icrs_coord, lsr_frame): # -> tuple[None, CartesianRepresentation]:
    ...

@frame_transform_graph.transform(AffineTransform, LSRD, ICRS)
def lsrd_to_icrs(lsr_coord, icrs_frame): # -> tuple[None, CartesianRepresentation]:
    ...

