from ...core.core import PixCoord as PixCoord, SkyRegion as SkyRegion
from ...core.metadata import RegionMeta as RegionMeta, RegionVisual as RegionVisual
from ...shapes import CircleAnnulusPixelRegion as CircleAnnulusPixelRegion, CircleAnnulusSkyRegion as CircleAnnulusSkyRegion, CirclePixelRegion as CirclePixelRegion, CircleSkyRegion as CircleSkyRegion, EllipseAnnulusPixelRegion as EllipseAnnulusPixelRegion, EllipseAnnulusSkyRegion as EllipseAnnulusSkyRegion, EllipsePixelRegion as EllipsePixelRegion, EllipseSkyRegion as EllipseSkyRegion, LinePixelRegion as LinePixelRegion, LineSkyRegion as LineSkyRegion, PointPixelRegion as PointPixelRegion, PointSkyRegion as PointSkyRegion, PolygonPixelRegion as PolygonPixelRegion, PolygonSkyRegion as PolygonSkyRegion, RectangleAnnulusPixelRegion as RectangleAnnulusPixelRegion, RectangleAnnulusSkyRegion as RectangleAnnulusSkyRegion, RectanglePixelRegion as RectanglePixelRegion, RectangleSkyRegion as RectangleSkyRegion, RegularPolygonPixelRegion as RegularPolygonPixelRegion, TextPixelRegion as TextPixelRegion, TextSkyRegion as TextSkyRegion
from ..crtf.core import CRTFRegionParserWarning as CRTFRegionParserWarning
from _typeshed import Incomplete

regions_attributes: Incomplete
reg_mapping: Incomplete
valid_coordsys: Incomplete
coordsys_mapping: Incomplete

class RegionConversionError(ValueError): ...

class _ShapeList(list):
    def to_regions(self): ...
    def to_crtf(self, coordsys: str = ..., fmt: str = ..., radunit: str = ...): ...

class _Shape:
    shape_to_sky_region: Incomplete
    shape_to_pixel_region: Incomplete
    error: Incomplete
    coord: Incomplete
    meta: Incomplete
    composite: Incomplete
    include: Incomplete
    def __init__(self, coordsys, region_type, coord, meta, composite, include) -> None: ...
    @property
    def coordsys(self): ...
    @coordsys.setter
    def coordsys(self, value) -> None: ...
    @property
    def region_type(self): ...
    @region_type.setter
    def region_type(self, value) -> None: ...
    def convert_coords(self): ...
    def to_region(self): ...
    def check_crtf(self) -> None: ...
