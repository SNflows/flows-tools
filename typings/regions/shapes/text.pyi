from .point import PointPixelRegion, PointSkyRegion
from _typeshed import Incomplete

class TextPixelRegion(PointPixelRegion):
    center: Incomplete
    meta: Incomplete
    visual: Incomplete
    text: Incomplete
    def __init__(self, center, text, meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_sky(self, wcs): ...
    def as_artist(self, origin=..., **kwargs): ...

class TextSkyRegion(PointSkyRegion):
    center: Incomplete
    meta: Incomplete
    visual: Incomplete
    text: Incomplete
    def __init__(self, center, text, meta: Incomplete | None = ..., visual: Incomplete | None = ...) -> None: ...
    def to_pixel(self, wcs): ...
