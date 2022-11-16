import numpy as np
from ...core import Region as Region, Regions as Regions, SkyRegion as SkyRegion
from ...core.registry import RegionsRegistry as RegionsRegistry
from ...shapes import RegularPolygonPixelRegion as RegularPolygonPixelRegion

class _RegionData:
    shape: str
    x: np.ndarray
    y: np.ndarray
    r: np.ndarray
    rotang: np.ndarray
    component: int
    def __init__(self, shape, x, y, r, rotang, component) -> None: ...
