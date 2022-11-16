from ...core import PixCoord as PixCoord, RegionMeta as RegionMeta, RegionVisual as RegionVisual, Regions as Regions
from ...core.registry import RegionsRegistry as RegionsRegistry
from .core import DS9ParserError as DS9ParserError, ds9_frame_map as ds9_frame_map, ds9_params_template as ds9_params_template, ds9_shape_to_region as ds9_shape_to_region

class _RegionData:
    frame: str
    region_type: str
    shape: str
    shape_params: str
    raw_meta: dict
    region_str: str
    def __init__(self, frame, region_type, shape, shape_params, raw_meta, region_str) -> None: ...
