from ._backend_tk import FigureCanvasTk as FigureCanvasTk, _BackendTk
from .backend_cairo import FigureCanvasCairo as FigureCanvasCairo, cairo as cairo
from _typeshed import Incomplete

class FigureCanvasTkCairo(FigureCanvasCairo, FigureCanvasTk):
    def draw(self) -> None: ...

class _BackendTkCairo(_BackendTk):
    FigureCanvas: Incomplete
