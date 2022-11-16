from ..core.core import PixelRegion, SkyRegion
from _typeshed import Incomplete

class EllipsePixelRegion(PixelRegion):
    center: Incomplete
    width: Incomplete
    height: Incomplete
    angle: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, width, height, angle=..., meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    @property
    def area(self): ...
    def contains(self, pixcoord): ...
    def to_sky(self, wcs): ...
    @property
    def bounding_box(self): ...
    def to_mask(self, mode: str = ..., subpixels: int = ...): ...
    def as_artist(self, origin=..., **kwargs): ...
    def as_mpl_selector(self, ax, active: bool = ..., sync: bool = ..., callback: Incomplete | None = ..., **kwargs): ...
    def rotate(self, center, angle): ...

class EllipseSkyRegion(SkyRegion):
    center: Incomplete
    width: Incomplete
    height: Incomplete
    angle: Incomplete
    meta: Incomplete
    visual: Incomplete
    def __init__(self, center, width, height, angle=..., meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_pixel(self, wcs): ...
