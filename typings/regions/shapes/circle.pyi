from ..core.core import PixelRegion, SkyRegion
from _typeshed import Incomplete

class CirclePixelRegion(PixelRegion):
    center: Incomplete
    radius: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, radius, meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    @property
    def area(self): ...
    def contains(self, pixcoord): ...
    def to_sky(self, wcs): ...
    @property
    def bounding_box(self): ...
    def to_mask(self, mode: str = ..., subpixels: int = ...): ...
    def as_artist(self, origin=..., **kwargs): ...
    def rotate(self, center, angle): ...

class CircleSkyRegion(SkyRegion):
    center: Incomplete
    radius: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, radius, meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_pixel(self, wcs): ...
