from ...core import Regions as Regions
from ...core.registry import RegionsRegistry as RegionsRegistry
from .core import CRTFRegionParserError as CRTFRegionParserError, CRTFRegionParserWarning as CRTFRegionParserWarning, valid_symbols as valid_symbols
from .io_core import reg_mapping as reg_mapping
from _typeshed import Incomplete

regex_begin: Incomplete
regex_comment: Incomplete
regex_global: Incomplete
regex_coordinate: Incomplete
regex_length: Incomplete
regex_meta: Incomplete
regex_region: Incomplete
regex_line: Incomplete

class _CRTFParser:
    valid_definition: Incomplete
    valid_global_keys: Incomplete
    region_string: Incomplete
    errors: Incomplete
    global_meta: Incomplete
    shapes: Incomplete
    def __init__(self, region_string, errors: str = ...) -> None: ...
    def parse_line(self, line) -> None: ...
    def run(self) -> None: ...
    def parse_global_meta(self, global_meta_str) -> None: ...

class _CRTFRegionParser:
    coordinate_systems: Incomplete
    coordsys_mapping: Incomplete
    language_spec: Incomplete
    global_meta: Incomplete
    reg_str: Incomplete
    meta_str: Incomplete
    errors: Incomplete
    coord: Incomplete
    coordsys: Incomplete
    coord_str: Incomplete
    type_: Incomplete
    region_type: Incomplete
    meta: Incomplete
    shape: Incomplete
    include: Incomplete
    def __init__(self, global_meta, include, type_, region_type, reg_str, meta_str, errors: str = ...) -> None: ...
    def parse(self) -> None: ...
    def set_coordsys(self) -> None: ...
    def convert_coordinates(self) -> None: ...
    def convert_meta(self) -> None: ...
    def make_shape(self) -> None: ...

class _CRTFCoordinateParser:
    @staticmethod
    def parse_coordinate(string_rep): ...
    @staticmethod
    def parse_angular_length_quantity(string_rep): ...
