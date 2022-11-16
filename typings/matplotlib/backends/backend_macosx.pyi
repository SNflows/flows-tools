from . import _macosx
from .backend_agg import FigureCanvasAgg as FigureCanvasAgg
from _typeshed import Incomplete
from matplotlib import cbook as cbook
from matplotlib._pylab_helpers import Gcf as Gcf
from matplotlib.backend_bases import FigureCanvasBase as FigureCanvasBase, FigureManagerBase as FigureManagerBase, NavigationToolbar2 as NavigationToolbar2, ResizeEvent as ResizeEvent, TimerBase as TimerBase, _Backend
from matplotlib.figure import Figure as Figure
from matplotlib.widgets import SubplotTool as SubplotTool

class TimerMac(_macosx.Timer, TimerBase): ...

class FigureCanvasMac(FigureCanvasAgg, _macosx.FigureCanvas, FigureCanvasBase):
    required_interactive_framework: str
    manager_class: Incomplete
    def __init__(self, figure) -> None: ...
    def draw(self) -> None: ...
    def draw_idle(self) -> None: ...
    def blit(self, bbox: Incomplete | None = ...) -> None: ...
    def resize(self, width, height) -> None: ...

class NavigationToolbar2Mac(_macosx.NavigationToolbar2, NavigationToolbar2):
    def __init__(self, canvas) -> None: ...
    def draw_rubberband(self, event, x0, y0, x1, y1) -> None: ...
    def remove_rubberband(self) -> None: ...
    def save_figure(self, *args) -> None: ...
    def prepare_configure_subplots(self): ...

class FigureManagerMac(_macosx.FigureManager, FigureManagerBase):
    def __init__(self, canvas, num) -> None: ...
    def close(self): ...
    def show(self) -> None: ...

class _BackendMac(_Backend):
    FigureCanvas: Incomplete
    FigureManager: Incomplete
    @staticmethod
    def mainloop() -> None: ...
