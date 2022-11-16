from _typeshed import Incomplete
from enum import Enum, IntEnum
from matplotlib import cbook as cbook, colors as colors, get_backend as get_backend, is_interactive as is_interactive, rcParams as rcParams, textpath as textpath, transforms as transforms, widgets as widgets
from matplotlib._enums import CapStyle as CapStyle, JoinStyle as JoinStyle
from matplotlib._pylab_helpers import Gcf as Gcf
from matplotlib.backend_managers import ToolManager as ToolManager
from matplotlib.path import Path as Path
from matplotlib.texmanager import TexManager as TexManager
from matplotlib.transforms import Affine2D as Affine2D
from typing import NamedTuple

def register_backend(format, backend, description: Incomplete | None = ...) -> None: ...
def get_registered_canvas_class(format): ...

class RendererBase:
    def __init__(self) -> None: ...
    def open_group(self, s, gid: Incomplete | None = ...) -> None: ...
    def close_group(self, s) -> None: ...
    def draw_path(self, gc, path, transform, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_markers(self, gc, marker_path, marker_trans, path, trans, rgbFace: Incomplete | None = ...) -> None: ...
    def draw_path_collection(self, gc, master_transform, paths, all_transforms, offsets, offset_trans, facecolors, edgecolors, linewidths, linestyles, antialiaseds, urls, offset_position) -> None: ...
    def draw_quad_mesh(self, gc, master_transform, meshWidth, meshHeight, coordinates, offsets, offsetTrans, facecolors, antialiased, edgecolors): ...
    def draw_gouraud_triangle(self, gc, points, colors, transform) -> None: ...
    def draw_gouraud_triangles(self, gc, triangles_array, colors_array, transform) -> None: ...
    def get_image_magnification(self): ...
    def draw_image(self, gc, x, y, im, transform: Incomplete | None = ...) -> None: ...
    def option_image_nocomposite(self): ...
    def option_scale_image(self): ...
    def draw_tex(self, gc, x, y, s, prop, angle, *, mtext: Incomplete | None = ...) -> None: ...
    def draw_text(self, gc, x, y, s, prop, angle, ismath: bool = ..., mtext: Incomplete | None = ...) -> None: ...
    def get_text_width_height_descent(self, s, prop, ismath): ...
    def flipy(self): ...
    def get_canvas_width_height(self): ...
    def get_texmanager(self): ...
    def new_gc(self): ...
    def points_to_pixels(self, points): ...
    def start_rasterizing(self) -> None: ...
    def stop_rasterizing(self) -> None: ...
    def start_filter(self) -> None: ...
    def stop_filter(self, filter_func) -> None: ...

class GraphicsContextBase:
    def __init__(self) -> None: ...
    def copy_properties(self, gc) -> None: ...
    def restore(self) -> None: ...
    def get_alpha(self): ...
    def get_antialiased(self): ...
    def get_capstyle(self): ...
    def get_clip_rectangle(self): ...
    def get_clip_path(self): ...
    def get_dashes(self): ...
    def get_forced_alpha(self): ...
    def get_joinstyle(self): ...
    def get_linewidth(self): ...
    def get_rgb(self): ...
    def get_url(self): ...
    def get_gid(self): ...
    def get_snap(self): ...
    def set_alpha(self, alpha) -> None: ...
    def set_antialiased(self, b) -> None: ...
    def set_capstyle(self, cs) -> None: ...
    def set_clip_rectangle(self, rectangle) -> None: ...
    def set_clip_path(self, path) -> None: ...
    def set_dashes(self, dash_offset, dash_list) -> None: ...
    def set_foreground(self, fg, isRGBA: bool = ...) -> None: ...
    def set_joinstyle(self, js) -> None: ...
    def set_linewidth(self, w) -> None: ...
    def set_url(self, url) -> None: ...
    def set_gid(self, id) -> None: ...
    def set_snap(self, snap) -> None: ...
    def set_hatch(self, hatch) -> None: ...
    def get_hatch(self): ...
    def get_hatch_path(self, density: float = ...): ...
    def get_hatch_color(self): ...
    def set_hatch_color(self, hatch_color) -> None: ...
    def get_hatch_linewidth(self): ...
    def get_sketch_params(self): ...
    def set_sketch_params(self, scale: Incomplete | None = ..., length: Incomplete | None = ..., randomness: Incomplete | None = ...) -> None: ...

class TimerBase:
    callbacks: Incomplete
    def __init__(self, interval: Incomplete | None = ..., callbacks: Incomplete | None = ...) -> None: ...
    def __del__(self) -> None: ...
    def start(self, interval: Incomplete | None = ...) -> None: ...
    def stop(self) -> None: ...
    @property
    def interval(self): ...
    @interval.setter
    def interval(self, interval) -> None: ...
    @property
    def single_shot(self): ...
    @single_shot.setter
    def single_shot(self, ss) -> None: ...
    def add_callback(self, func, *args, **kwargs): ...
    def remove_callback(self, func, *args, **kwargs) -> None: ...

class Event:
    name: Incomplete
    canvas: Incomplete
    guiEvent: Incomplete
    def __init__(self, name, canvas, guiEvent: Incomplete | None = ...) -> None: ...

class DrawEvent(Event):
    renderer: Incomplete
    def __init__(self, name, canvas, renderer) -> None: ...

class ResizeEvent(Event):
    def __init__(self, name, canvas) -> None: ...

class CloseEvent(Event): ...

class LocationEvent(Event):
    lastevent: Incomplete
    x: Incomplete
    y: Incomplete
    inaxes: Incomplete
    xdata: Incomplete
    ydata: Incomplete
    def __init__(self, name, canvas, x, y, guiEvent: Incomplete | None = ...) -> None: ...

class MouseButton(IntEnum):
    LEFT: int
    MIDDLE: int
    RIGHT: int
    BACK: int
    FORWARD: int

class MouseEvent(LocationEvent):
    button: Incomplete
    key: Incomplete
    step: Incomplete
    dblclick: Incomplete
    def __init__(self, name, canvas, x, y, button: Incomplete | None = ..., key: Incomplete | None = ..., step: int = ..., dblclick: bool = ..., guiEvent: Incomplete | None = ...) -> None: ...

class PickEvent(Event):
    mouseevent: Incomplete
    artist: Incomplete
    def __init__(self, name, canvas, mouseevent, artist, guiEvent: Incomplete | None = ..., **kwargs) -> None: ...

class KeyEvent(LocationEvent):
    key: Incomplete
    def __init__(self, name, canvas, key, x: int = ..., y: int = ..., guiEvent: Incomplete | None = ...) -> None: ...

class FigureCanvasBase:
    required_interactive_framework: Incomplete
    manager_class: Incomplete
    events: Incomplete
    fixed_dpi: Incomplete
    filetypes: Incomplete
    def supports_blit(cls): ...
    figure: Incomplete
    manager: Incomplete
    widgetlock: Incomplete
    mouse_grabber: Incomplete
    toolbar: Incomplete
    def __init__(self, figure: Incomplete | None = ...) -> None: ...
    callbacks: Incomplete
    button_pick_id: Incomplete
    scroll_pick_id: Incomplete
    @classmethod
    def new_manager(cls, figure, num): ...
    def is_saving(self): ...
    def pick(self, mouseevent) -> None: ...
    def blit(self, bbox: Incomplete | None = ...) -> None: ...
    def resize(self, w, h): ...
    def draw_event(self, renderer) -> None: ...
    def resize_event(self) -> None: ...
    def close_event(self, guiEvent: Incomplete | None = ...) -> None: ...
    def key_press_event(self, key, guiEvent: Incomplete | None = ...) -> None: ...
    def key_release_event(self, key, guiEvent: Incomplete | None = ...) -> None: ...
    def pick_event(self, mouseevent, artist, **kwargs) -> None: ...
    def scroll_event(self, x, y, step, guiEvent: Incomplete | None = ...) -> None: ...
    def button_press_event(self, x, y, button, dblclick: bool = ..., guiEvent: Incomplete | None = ...) -> None: ...
    def button_release_event(self, x, y, button, guiEvent: Incomplete | None = ...) -> None: ...
    def motion_notify_event(self, x, y, guiEvent: Incomplete | None = ...) -> None: ...
    def leave_notify_event(self, guiEvent: Incomplete | None = ...) -> None: ...
    def enter_notify_event(self, guiEvent: Incomplete | None = ..., xy: Incomplete | None = ...) -> None: ...
    def inaxes(self, xy): ...
    def grab_mouse(self, ax) -> None: ...
    def release_mouse(self, ax) -> None: ...
    def set_cursor(self, cursor) -> None: ...
    def draw(self, *args, **kwargs) -> None: ...
    def draw_idle(self, *args, **kwargs) -> None: ...
    @property
    def device_pixel_ratio(self): ...
    def get_width_height(self, *, physical: bool = ...): ...
    @classmethod
    def get_supported_filetypes(cls): ...
    @classmethod
    def get_supported_filetypes_grouped(cls): ...
    def print_figure(self, filename, dpi: Incomplete | None = ..., facecolor: Incomplete | None = ..., edgecolor: Incomplete | None = ..., orientation: str = ..., format: Incomplete | None = ..., *, bbox_inches: Incomplete | None = ..., pad_inches: Incomplete | None = ..., bbox_extra_artists: Incomplete | None = ..., backend: Incomplete | None = ..., **kwargs): ...
    @classmethod
    def get_default_filetype(cls): ...
    def get_default_filename(self): ...
    def switch_backends(self, FigureCanvasClass): ...
    def mpl_connect(self, s, func): ...
    def mpl_disconnect(self, cid): ...
    def new_timer(self, interval: Incomplete | None = ..., callbacks: Incomplete | None = ...): ...
    def flush_events(self) -> None: ...
    def start_event_loop(self, timeout: int = ...) -> None: ...
    def stop_event_loop(self) -> None: ...

def key_press_handler(event, canvas: Incomplete | None = ..., toolbar: Incomplete | None = ...): ...
def button_press_handler(event, canvas: Incomplete | None = ..., toolbar: Incomplete | None = ...) -> None: ...

class NonGuiException(Exception): ...

class FigureManagerBase:
    canvas: Incomplete
    num: Incomplete
    key_press_handler_id: Incomplete
    button_press_handler_id: Incomplete
    toolmanager: Incomplete
    toolbar: Incomplete
    def __init__(self, canvas, num) -> None: ...
    @classmethod
    def create_with_canvas(cls, canvas_class, figure, num): ...
    def show(self) -> None: ...
    def destroy(self) -> None: ...
    def full_screen_toggle(self) -> None: ...
    def resize(self, w, h) -> None: ...
    def get_window_title(self): ...
    def set_window_title(self, title) -> None: ...

cursors: Incomplete

class _Mode(str, Enum):
    NONE: str
    PAN: str
    ZOOM: str

class NavigationToolbar2:
    toolitems: Incomplete
    canvas: Incomplete
    mode: Incomplete
    def __init__(self, canvas) -> None: ...
    def set_message(self, s) -> None: ...
    def draw_rubberband(self, event, x0, y0, x1, y1) -> None: ...
    def remove_rubberband(self) -> None: ...
    def home(self, *args) -> None: ...
    def back(self, *args) -> None: ...
    def forward(self, *args) -> None: ...
    def mouse_move(self, event) -> None: ...
    def pan(self, *args) -> None: ...

    class _PanInfo(NamedTuple):
        button: Incomplete
        axes: Incomplete
        cid: Incomplete
    def press_pan(self, event) -> None: ...
    def drag_pan(self, event) -> None: ...
    def release_pan(self, event) -> None: ...
    def zoom(self, *args) -> None: ...

    class _ZoomInfo(NamedTuple):
        direction: Incomplete
        start_xy: Incomplete
        axes: Incomplete
        cid: Incomplete
        cbar: Incomplete
    def press_zoom(self, event) -> None: ...
    def drag_zoom(self, event) -> None: ...
    def release_zoom(self, event) -> None: ...
    def push_current(self) -> None: ...
    subplot_tool: Incomplete
    def configure_subplots(self, *args): ...
    def save_figure(self, *args) -> None: ...
    def set_cursor(self, cursor) -> None: ...
    def update(self) -> None: ...
    def set_history_buttons(self) -> None: ...

class ToolContainerBase:
    toolmanager: Incomplete
    def __init__(self, toolmanager) -> None: ...
    def add_tool(self, tool, group, position: int = ...) -> None: ...
    def trigger_tool(self, name) -> None: ...
    def add_toolitem(self, name, group, position, image, description, toggle) -> None: ...
    def toggle_toolitem(self, name, toggled) -> None: ...
    def remove_toolitem(self, name) -> None: ...
    def set_message(self, s) -> None: ...

class _Backend:
    backend_version: str
    FigureCanvas: Incomplete
    FigureManager: Incomplete
    mainloop: Incomplete
    @classmethod
    def new_figure_manager(cls, num, *args, **kwargs): ...
    @classmethod
    def new_figure_manager_given_figure(cls, num, figure): ...
    @classmethod
    def draw_if_interactive(cls) -> None: ...
    @classmethod
    def show(cls, *, block: Incomplete | None = ...) -> None: ...
    @staticmethod
    def export(cls): ...

class ShowBase(_Backend):
    def __call__(self, block: Incomplete | None = ...): ...
