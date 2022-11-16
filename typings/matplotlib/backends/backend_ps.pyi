from . import _backend_pdf_ps
from _typeshed import Incomplete
from enum import Enum
from matplotlib import cbook as cbook
from matplotlib._afm import AFM as AFM
from matplotlib._mathtext_data import uni2type1 as uni2type1
from matplotlib._ttconv import convert_ttf_to_ps as convert_ttf_to_ps
from matplotlib.backend_bases import FigureCanvasBase as FigureCanvasBase, FigureManagerBase as FigureManagerBase, RendererBase as RendererBase, _Backend
from matplotlib.backends.backend_mixed import MixedModeRenderer as MixedModeRenderer
from matplotlib.cbook import file_requires_unicode as file_requires_unicode, is_writable_file_like as is_writable_file_like
from matplotlib.font_manager import get_font as get_font
from matplotlib.ft2font import FT2Font as FT2Font, LOAD_NO_SCALE as LOAD_NO_SCALE
from matplotlib.path import Path as Path
from matplotlib.texmanager import TexManager as TexManager
from matplotlib.transforms import Affine2D as Affine2D

debugPS: bool

class PsBackendHelper:
    def __init__(self) -> None: ...

ps_backend_helper: Incomplete
papersize: Incomplete

def quote_ps_string(s): ...

class RendererPS(_backend_pdf_ps.RendererPDFPSBase):
    textcnt: int
    psfrag: Incomplete
    imagedpi: Incomplete
    color: Incomplete
    linewidth: Incomplete
    linejoin: Incomplete
    linecap: Incomplete
    linedash: Incomplete
    fontname: Incomplete
    fontsize: Incomplete
    image_magnification: Incomplete
    def __init__(self, width, height, pswriter, imagedpi: int = ...) -> None: ...
    def set_color(self, r, g, b, store: bool = ...) -> None: ...
    def set_linewidth(self, linewidth, store: bool = ...) -> None: ...
    def set_linejoin(self, linejoin, store: bool = ...) -> None: ...
    def set_linecap(self, linecap, store: bool = ...) -> None: ...
    def set_linedash(self, offset, seq, store: bool = ...) -> None: ...
    def set_font(self, fontname, fontsize, store: bool = ...) -> None: ...
    def create_hatch(self, hatch): ...
    def get_image_magnification(self): ...
    def draw_image(self, gc, x, y, im, transform: Incomplete | None = ...) -> None: ...
    def draw_path(self, gc, path, transform, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_markers(self, gc, marker_path, marker_trans, path, trans, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_path_collection(self, gc, master_transform, paths, all_transforms, offsets, offset_trans, facecolors, edgecolors, linewidths, linestyles, antialiaseds, urls, offset_position): ...
    def draw_tex(self, gc, x, y, s, prop, angle, *, mtext: Incomplete | None = ...) -> None: ...
    def draw_text(self, gc, x, y, s, prop, angle, ismath: bool = ..., mtext: Incomplete | None = ...): ...
    def draw_mathtext(self, gc, x, y, s, prop, angle) -> None: ...
    def draw_gouraud_triangle(self, gc, points, colors, trans) -> None: ...
    def draw_gouraud_triangles(self, gc, points, colors, trans) -> None: ...

class _Orientation(Enum):
    portrait: Incomplete
    landscape: Incomplete
    def swap_if_landscape(self, shape): ...

class FigureCanvasPS(FigureCanvasBase):
    fixed_dpi: int
    filetypes: Incomplete
    def get_default_filetype(self): ...
    print_ps: Incomplete
    print_eps: Incomplete
    def draw(self): ...

def convert_psfrags(tmpfile, psfrags, font_preamble, custom_preamble, paper_width, paper_height, orientation): ...
def gs_distill(tmpfile, eps: bool = ..., ptype: str = ..., bbox: Incomplete | None = ..., rotated: bool = ...) -> None: ...
def xpdf_distill(tmpfile, eps: bool = ..., ptype: str = ..., bbox: Incomplete | None = ..., rotated: bool = ...) -> None: ...
def get_bbox_header(lbrt, rotated: bool = ...): ...
def pstoeps(tmpfile, bbox: Incomplete | None = ..., rotated: bool = ...) -> None: ...
FigureManagerPS = FigureManagerBase
psDefs: Incomplete

class _BackendPS(_Backend):
    backend_version: str
    FigureCanvas: Incomplete
