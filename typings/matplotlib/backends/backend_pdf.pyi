from . import _backend_pdf_ps
from _typeshed import Incomplete
from enum import Enum
from matplotlib import cbook as cbook, dviread as dviread
from matplotlib._afm import AFM as AFM
from matplotlib._pylab_helpers import Gcf as Gcf
from matplotlib.backend_bases import FigureCanvasBase as FigureCanvasBase, FigureManagerBase as FigureManagerBase, GraphicsContextBase as GraphicsContextBase, RendererBase as RendererBase, _Backend
from matplotlib.backends.backend_mixed import MixedModeRenderer as MixedModeRenderer
from matplotlib.dates import UTC as UTC
from matplotlib.figure import Figure as Figure
from matplotlib.font_manager import get_font as get_font
from matplotlib.ft2font import FIXED_WIDTH as FIXED_WIDTH, FT2Font as FT2Font, ITALIC as ITALIC, KERNING_UNFITTED as KERNING_UNFITTED, LOAD_NO_HINTING as LOAD_NO_HINTING, LOAD_NO_SCALE as LOAD_NO_SCALE
from matplotlib.path import Path as Path
from matplotlib.transforms import Affine2D as Affine2D, BboxBase as BboxBase

def fill(strings, linelen: int = ...): ...
def pdfRepr(obj): ...

class Reference:
    id: Incomplete
    def __init__(self, id) -> None: ...
    def pdfRepr(self): ...
    def write(self, contents, file) -> None: ...

class Name:
    name: Incomplete
    def __init__(self, name) -> None: ...
    def __eq__(self, other): ...
    def __lt__(self, other): ...
    def __hash__(self): ...
    @staticmethod
    def hexify(match): ...
    def pdfRepr(self): ...

class Operator:
    op: Incomplete
    def __init__(self, op) -> None: ...
    def pdfRepr(self): ...

class Verbatim:
    def __init__(self, x) -> None: ...
    def pdfRepr(self): ...

class Op(Enum):
    close_fill_stroke: bytes
    fill_stroke: bytes
    fill: bytes
    closepath: bytes
    close_stroke: bytes
    stroke: bytes
    endpath: bytes
    begin_text: bytes
    end_text: bytes
    curveto: bytes
    rectangle: bytes
    lineto: bytes
    moveto: bytes
    concat_matrix: bytes
    use_xobject: bytes
    setgray_stroke: bytes
    setgray_nonstroke: bytes
    setrgb_stroke: bytes
    setrgb_nonstroke: bytes
    setcolorspace_stroke: bytes
    setcolorspace_nonstroke: bytes
    setcolor_stroke: bytes
    setcolor_nonstroke: bytes
    setdash: bytes
    setlinejoin: bytes
    setlinecap: bytes
    setgstate: bytes
    gsave: bytes
    grestore: bytes
    textpos: bytes
    selectfont: bytes
    textmatrix: bytes
    show: bytes
    showkern: bytes
    setlinewidth: bytes
    clip: bytes
    shading: bytes
    op: Incomplete
    def pdfRepr(self): ...
    @classmethod
    def paint_path(cls, fill, stroke): ...

class Stream:
    id: Incomplete
    len: Incomplete
    pdfFile: Incomplete
    file: Incomplete
    compressobj: Incomplete
    extra: Incomplete
    pos: Incomplete
    def __init__(self, id, len, file, extra: Incomplete | None = ..., png: Incomplete | None = ...) -> None: ...
    def end(self) -> None: ...
    def write(self, data) -> None: ...

class PdfFile:
    xrefTable: Incomplete
    passed_in_file_object: bool
    original_file_like: Incomplete
    tell_base: int
    fh: Incomplete
    currentstream: Incomplete
    rootObject: Incomplete
    pagesObject: Incomplete
    pageList: Incomplete
    fontObject: Incomplete
    hatchObject: Incomplete
    gouraudObject: Incomplete
    XObjectObject: Incomplete
    resourceObject: Incomplete
    infoDict: Incomplete
    fontNames: Incomplete
    dviFontInfo: Incomplete
    type1Descriptors: Incomplete
    alphaStates: Incomplete
    hatchPatterns: Incomplete
    gouraudTriangles: Incomplete
    markers: Incomplete
    multi_byte_charprocs: Incomplete
    paths: Incomplete
    pageAnnotations: Incomplete
    def __init__(self, filename, metadata: Incomplete | None = ...) -> None: ...
    def newPage(self, width, height) -> None: ...
    def newTextnote(self, text, positionRect=...) -> None: ...
    def finalize(self) -> None: ...
    def close(self) -> None: ...
    def write(self, data) -> None: ...
    def output(self, *data) -> None: ...
    def beginStream(self, id, len, extra: Incomplete | None = ..., png: Incomplete | None = ...) -> None: ...
    def endStream(self) -> None: ...
    def outputStream(self, ref, data, *, extra: Incomplete | None = ...) -> None: ...
    def fontName(self, fontprop): ...
    def dviFontName(self, dvifont): ...
    def writeFonts(self) -> None: ...
    def createType1Descriptor(self, t1font, fontfile): ...
    def embedTTF(self, filename, characters): ...
    def alphaState(self, alpha): ...
    def writeExtGSTates(self) -> None: ...
    def hatchPattern(self, hatch_style): ...
    def writeHatches(self) -> None: ...
    def addGouraudTriangles(self, points, colors): ...
    def writeGouraudTriangles(self) -> None: ...
    def imageObject(self, image): ...
    def writeImages(self) -> None: ...
    def markerObject(self, path, trans, fill, stroke, lw, joinstyle, capstyle): ...
    def writeMarkers(self) -> None: ...
    def pathCollectionObject(self, gc, path, trans, padding, filled, stroked): ...
    def writePathCollectionTemplates(self) -> None: ...
    @staticmethod
    def pathOperations(path, transform, clip: Incomplete | None = ..., simplify: Incomplete | None = ..., sketch: Incomplete | None = ...): ...
    def writePath(self, path, transform, clip: bool = ..., sketch: Incomplete | None = ...) -> None: ...
    def reserveObject(self, name: str = ...): ...
    def recordXref(self, id) -> None: ...
    def writeObject(self, object, contents) -> None: ...
    startxref: Incomplete
    def writeXref(self) -> None: ...
    infoObject: Incomplete
    def writeInfoDict(self) -> None: ...
    def writeTrailer(self) -> None: ...

class RendererPdf(_backend_pdf_ps.RendererPDFPSBase):
    file: Incomplete
    gc: Incomplete
    image_dpi: Incomplete
    def __init__(self, file, image_dpi, height, width) -> None: ...
    def finalize(self) -> None: ...
    def check_gc(self, gc, fillcolor: Incomplete | None = ...) -> None: ...
    def get_image_magnification(self): ...
    def draw_image(self, gc, x, y, im, transform: Incomplete | None = ...) -> None: ...
    def draw_path(self, gc, path, transform, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_path_collection(self, gc, master_transform, paths, all_transforms, offsets, offset_trans, facecolors, edgecolors, linewidths, linestyles, antialiaseds, urls, offset_position): ...
    def draw_markers(self, gc, marker_path, marker_trans, path, trans, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_gouraud_triangle(self, gc, points, colors, trans) -> None: ...
    def draw_gouraud_triangles(self, gc, points, colors, trans) -> None: ...
    def draw_mathtext(self, gc, x, y, s, prop, angle) -> None: ...
    def draw_tex(self, gc, x, y, s, prop, angle, *, mtext: Incomplete | None = ...) -> None: ...
    def encode_string(self, s, fonttype): ...
    def draw_text(self, gc, x, y, s, prop, angle, ismath: bool = ..., mtext: Incomplete | None = ...): ...
    def new_gc(self): ...

class GraphicsContextPdf(GraphicsContextBase):
    file: Incomplete
    parent: Incomplete
    def __init__(self, file) -> None: ...
    def stroke(self): ...
    def fill(self, *args): ...
    def paint(self): ...
    capstyles: Incomplete
    joinstyles: Incomplete
    def capstyle_cmd(self, style): ...
    def joinstyle_cmd(self, style): ...
    def linewidth_cmd(self, width): ...
    def dash_cmd(self, dashes): ...
    def alpha_cmd(self, alpha, forced, effective_alphas): ...
    def hatch_cmd(self, hatch, hatch_color): ...
    def rgb_cmd(self, rgb): ...
    def fillcolor_cmd(self, rgb): ...
    def push(self): ...
    def pop(self): ...
    def clip_cmd(self, cliprect, clippath): ...
    commands: Incomplete
    def delta(self, other): ...
    def copy_properties(self, other) -> None: ...
    def finalize(self): ...

class PdfPages:
    keep_empty: Incomplete
    def __init__(self, filename, keep_empty: bool = ..., metadata: Incomplete | None = ...) -> None: ...
    def __enter__(self): ...
    def __exit__(self, exc_type, exc_val, exc_tb) -> None: ...
    def close(self) -> None: ...
    def infodict(self): ...
    def savefig(self, figure: Incomplete | None = ..., **kwargs) -> None: ...
    def get_pagecount(self): ...
    def attach_note(self, text, positionRect=...) -> None: ...

class FigureCanvasPdf(FigureCanvasBase):
    fixed_dpi: int
    filetypes: Incomplete
    def get_default_filetype(self): ...
    def print_pdf(self, filename, *, bbox_inches_restore: Incomplete | None = ..., metadata: Incomplete | None = ...) -> None: ...
    def draw(self): ...
FigureManagerPdf = FigureManagerBase

class _BackendPdf(_Backend):
    FigureCanvas: Incomplete
