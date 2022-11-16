import abc
from ..core.core import PixelRegion, SkyRegion
from _typeshed import Incomplete

class AnnulusPixelRegion(PixelRegion, abc.ABC, metaclass=abc.ABCMeta):
    @property
    def area(self): ...
    @property
    def bounding_box(self): ...
    def contains(self, pixcoord): ...
    def as_artist(self, origin=..., **kwargs): ...
    def to_mask(self, mode: str = ..., subpixels: int = ...): ...
    def rotate(self, center, angle): ...

class CircleAnnulusPixelRegion(AnnulusPixelRegion):
    center: Incomplete
    inner_radius: Incomplete
    outer_radius: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, inner_radius, outer_radius, meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_sky(self, wcs): ...

class CircleAnnulusSkyRegion(SkyRegion):
    center: Incomplete
    inner_radius: Incomplete
    outer_radius: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, inner_radius, outer_radius, meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_pixel(self, wcs): ...

class AsymmetricAnnulusPixelRegion(AnnulusPixelRegion, metaclass=abc.ABCMeta):
    center: Incomplete
    inner_width: Incomplete
    outer_width: Incomplete
    inner_height: Incomplete
    outer_height: Incomplete
    angle: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, inner_width, outer_width, inner_height, outer_height, angle=..., meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_sky_args(self, wcs): ...

class AsymmetricAnnulusSkyRegion(SkyRegion, metaclass=abc.ABCMeta):
    center: Incomplete
    inner_width: Incomplete
    outer_width: Incomplete
    inner_height: Incomplete
    outer_height: Incomplete
    angle: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, inner_width, outer_width, inner_height, outer_height, angle=..., meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_pixel_args(self, wcs): ...

class EllipseAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    center: Incomplete
    inner_width: Incomplete
    outer_width: Incomplete
    inner_height: Incomplete
    outer_height: Incomplete
    angle: Incomplete
    def to_sky(self, wcs): ...

class EllipseAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    center: Incomplete
    inner_width: Incomplete
    outer_width: Incomplete
    inner_height: Incomplete
    outer_height: Incomplete
    angle: Incomplete
    def to_pixel(self, wcs): ...

class RectangleAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    center: Incomplete
    inner_width: Incomplete
    outer_width: Incomplete
    inner_height: Incomplete
    outer_height: Incomplete
    angle: Incomplete
    def to_sky(self, wcs): ...

class RectangleAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    center: Incomplete
    inner_width: Incomplete
    outer_width: Incomplete
    inner_height: Incomplete
    outer_height: Incomplete
    angle: Incomplete
    def to_pixel(self, wcs): ...
