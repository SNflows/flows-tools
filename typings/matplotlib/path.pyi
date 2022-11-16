from .bezier import BezierSegment as BezierSegment
from .cbook import simple_linear_interpolation as simple_linear_interpolation
from _typeshed import Incomplete
from collections.abc import Generator

class Path:
    code_type: Incomplete
    STOP: Incomplete
    MOVETO: Incomplete
    LINETO: Incomplete
    CURVE3: Incomplete
    CURVE4: Incomplete
    CLOSEPOLY: Incomplete
    NUM_VERTICES_FOR_CODE: Incomplete
    def __init__(self, vertices, codes: Incomplete | None = ..., _interpolation_steps: int = ..., closed: bool = ..., readonly: bool = ...) -> None: ...
    @property
    def vertices(self): ...
    @vertices.setter
    def vertices(self, vertices) -> None: ...
    @property
    def codes(self): ...
    @codes.setter
    def codes(self, codes) -> None: ...
    @property
    def simplify_threshold(self): ...
    @simplify_threshold.setter
    def simplify_threshold(self, threshold) -> None: ...
    @property
    def should_simplify(self): ...
    @should_simplify.setter
    def should_simplify(self, should_simplify) -> None: ...
    @property
    def readonly(self): ...
    def copy(self): ...
    def __deepcopy__(self, memo: Incomplete | None = ...): ...
    deepcopy: Incomplete
    @classmethod
    def make_compound_path_from_polys(cls, XY): ...
    @classmethod
    def make_compound_path(cls, *args): ...
    def __len__(self) -> int: ...
    def iter_segments(self, transform: Incomplete | None = ..., remove_nans: bool = ..., clip: Incomplete | None = ..., snap: bool = ..., stroke_width: float = ..., simplify: Incomplete | None = ..., curves: bool = ..., sketch: Incomplete | None = ...) -> Generator[Incomplete, None, None]: ...
    def iter_bezier(self, **kwargs) -> Generator[Incomplete, None, None]: ...
    def cleaned(self, transform: Incomplete | None = ..., remove_nans: bool = ..., clip: Incomplete | None = ..., *, simplify: bool = ..., curves: bool = ..., stroke_width: float = ..., snap: bool = ..., sketch: Incomplete | None = ...): ...
    def transformed(self, transform): ...
    def contains_point(self, point, transform: Incomplete | None = ..., radius: float = ...): ...
    def contains_points(self, points, transform: Incomplete | None = ..., radius: float = ...): ...
    def contains_path(self, path, transform: Incomplete | None = ...): ...
    def get_extents(self, transform: Incomplete | None = ..., **kwargs): ...
    def intersects_path(self, other, filled: bool = ...): ...
    def intersects_bbox(self, bbox, filled: bool = ...): ...
    def interpolated(self, steps): ...
    def to_polygons(self, transform: Incomplete | None = ..., width: int = ..., height: int = ..., closed_only: bool = ...): ...
    @classmethod
    def unit_rectangle(cls): ...
    @classmethod
    def unit_regular_polygon(cls, numVertices): ...
    @classmethod
    def unit_regular_star(cls, numVertices, innerCircle: float = ...): ...
    @classmethod
    def unit_regular_asterisk(cls, numVertices): ...
    @classmethod
    def unit_circle(cls): ...
    @classmethod
    def circle(cls, center=..., radius: float = ..., readonly: bool = ...): ...
    @classmethod
    def unit_circle_righthalf(cls): ...
    @classmethod
    def arc(cls, theta1, theta2, n: Incomplete | None = ..., is_wedge: bool = ...): ...
    @classmethod
    def wedge(cls, theta1, theta2, n: Incomplete | None = ...): ...
    @staticmethod
    def hatch(hatchpattern, density: int = ...): ...
    def clip_to_bbox(self, bbox, inside: bool = ...): ...

def get_path_collection_extents(master_transform, paths, transforms, offsets, offset_transform): ...