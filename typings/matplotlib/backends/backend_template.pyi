from _typeshed import Incomplete
from matplotlib._pylab_helpers import Gcf as Gcf
from matplotlib.backend_bases import FigureCanvasBase as FigureCanvasBase, FigureManagerBase as FigureManagerBase, GraphicsContextBase as GraphicsContextBase, RendererBase as RendererBase
from matplotlib.figure import Figure as Figure

class RendererTemplate(RendererBase):
    dpi: Incomplete
    def __init__(self, dpi) -> None: ...
    def draw_path(self, gc, path, transform, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_image(self, gc, x, y, im) -> None: ...
    def draw_text(self, gc, x, y, s, prop, angle, ismath: bool = ..., mtext: Incomplete | None = ...) -> None: ...
    def flipy(self): ...
    def get_canvas_width_height(self): ...
    def get_text_width_height_descent(self, s, prop, ismath): ...
    def new_gc(self): ...
    def points_to_pixels(self, points): ...

class GraphicsContextTemplate(GraphicsContextBase): ...

def show(*, block: Incomplete | None = ...) -> None: ...

class FigureManagerTemplate(FigureManagerBase): ...

class FigureCanvasTemplate(FigureCanvasBase):
    manager_class: Incomplete
    def draw(self) -> None: ...
    filetypes: Incomplete
    def print_foo(self, filename, *args, **kwargs) -> None: ...
    def get_default_filetype(self): ...
FigureCanvas = FigureCanvasTemplate
FigureManager = FigureManagerTemplate
