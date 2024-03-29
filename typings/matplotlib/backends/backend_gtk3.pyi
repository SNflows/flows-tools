from ._backend_gtk import _BackendGTK, _FigureManagerGTK, _NavigationToolbar2GTK
from _typeshed import Incomplete
from gi.repository import Gtk
from matplotlib import backend_tools as backend_tools, cbook as cbook
from matplotlib.backend_bases import CloseEvent as CloseEvent, FigureCanvasBase as FigureCanvasBase, KeyEvent as KeyEvent, LocationEvent as LocationEvent, MouseEvent as MouseEvent, ResizeEvent as ResizeEvent, ToolContainerBase as ToolContainerBase
from matplotlib.backend_tools import Cursors as Cursors

class __getattr__:
    @property
    def cursord(self): ...
    icon_filename: Incomplete
    window_icon: Incomplete

class FigureCanvasGTK3(FigureCanvasBase, Gtk.DrawingArea):
    required_interactive_framework: str
    manager_class: Incomplete
    event_mask: Incomplete
    def __init__(self, figure: Incomplete | None = ...) -> None: ...
    def destroy(self) -> None: ...
    def set_cursor(self, cursor) -> None: ...
    def scroll_event(self, widget, event): ...
    def button_press_event(self, widget, event): ...
    def button_release_event(self, widget, event): ...
    def key_press_event(self, widget, event): ...
    def key_release_event(self, widget, event): ...
    def motion_notify_event(self, widget, event): ...
    def enter_notify_event(self, widget, event) -> None: ...
    def leave_notify_event(self, widget, event) -> None: ...
    def size_allocate(self, widget, allocation) -> None: ...
    def configure_event(self, widget, event): ...
    def on_draw_event(self, widget, ctx) -> None: ...
    def draw(self) -> None: ...
    def draw_idle(self): ...
    def flush_events(self) -> None: ...

class NavigationToolbar2GTK3(_NavigationToolbar2GTK, Gtk.Toolbar):
    message: Incomplete
    def __init__(self, canvas, window: Incomplete | None = ...) -> None: ...
    win: Incomplete
    def save_figure(self, *args) -> None: ...

class ToolbarGTK3(ToolContainerBase, Gtk.Box):
    def __init__(self, toolmanager) -> None: ...
    def add_toolitem(self, name, group, position, image_file, description, toggle) -> None: ...
    def toggle_toolitem(self, name, toggled) -> None: ...
    def remove_toolitem(self, name) -> None: ...
    def set_message(self, s) -> None: ...

class SaveFigureGTK3(backend_tools.SaveFigureBase):
    def trigger(self, *args, **kwargs) -> None: ...

class SetCursorGTK3(backend_tools.SetCursorBase):
    def set_cursor(self, cursor) -> None: ...

class HelpGTK3(backend_tools.ToolHelpBase):
    def trigger(self, *args) -> None: ...

class ToolCopyToClipboardGTK3(backend_tools.ToolCopyToClipboardBase):
    def trigger(self, *args, **kwargs) -> None: ...

def error_msg_gtk(msg, parent: Incomplete | None = ...) -> None: ...
Toolbar = ToolbarGTK3

class FigureManagerGTK3(_FigureManagerGTK): ...

class _BackendGTK3(_BackendGTK):
    FigureCanvas: Incomplete
    FigureManager: Incomplete
