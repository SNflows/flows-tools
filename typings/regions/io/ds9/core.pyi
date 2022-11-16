from ..._utils.optional_deps import HAS_MATPLOTLIB as HAS_MATPLOTLIB
from ...shapes import CircleAnnulusPixelRegion as CircleAnnulusPixelRegion, CircleAnnulusSkyRegion as CircleAnnulusSkyRegion, CirclePixelRegion as CirclePixelRegion, CircleSkyRegion as CircleSkyRegion, EllipseAnnulusPixelRegion as EllipseAnnulusPixelRegion, EllipseAnnulusSkyRegion as EllipseAnnulusSkyRegion, EllipsePixelRegion as EllipsePixelRegion, EllipseSkyRegion as EllipseSkyRegion, LinePixelRegion as LinePixelRegion, LineSkyRegion as LineSkyRegion, PointPixelRegion as PointPixelRegion, PointSkyRegion as PointSkyRegion, PolygonPixelRegion as PolygonPixelRegion, PolygonSkyRegion as PolygonSkyRegion, RectangleAnnulusPixelRegion as RectangleAnnulusPixelRegion, RectangleAnnulusSkyRegion as RectangleAnnulusSkyRegion, RectanglePixelRegion as RectanglePixelRegion, RectangleSkyRegion as RectangleSkyRegion, TextPixelRegion as TextPixelRegion, TextSkyRegion as TextSkyRegion
from _typeshed import Incomplete

ds9_frame_map: Incomplete
pixel_map: Incomplete
sky_map: Incomplete
ds9_shape_to_region: Incomplete
ds9_params_template: Incomplete
ds9_shape_templates: Incomplete
boxcircle: str
arrow: str
vertices: Incomplete
codes: Incomplete
arrow_verts: Incomplete
ds9_valid_symbols: Incomplete

class DS9ParserError(Exception): ...
