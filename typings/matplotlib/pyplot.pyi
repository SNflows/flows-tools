from .ticker import AutoLocator as AutoLocator, FixedFormatter as FixedFormatter, FixedLocator as FixedLocator, FormatStrFormatter as FormatStrFormatter, Formatter as Formatter, FuncFormatter as FuncFormatter, IndexLocator as IndexLocator, LinearLocator as LinearLocator, Locator as Locator, LogFormatter as LogFormatter, LogFormatterExponent as LogFormatterExponent, LogFormatterMathtext as LogFormatterMathtext, LogLocator as LogLocator, MaxNLocator as MaxNLocator, MultipleLocator as MultipleLocator, NullFormatter as NullFormatter, NullLocator as NullLocator, ScalarFormatter as ScalarFormatter, TickHelper as TickHelper
from _typeshed import Incomplete
#from cycler import cycler as cycler
from matplotlib import cbook as cbook, cm as cm, get_backend as get_backend, interactive as interactive, mlab as mlab, rcParams as rcParams, rcParamsDefault as rcParamsDefault, rcParamsOrig as rcParamsOrig, rcsetup as rcsetup, style as style
from matplotlib.artist import Artist as Artist
from matplotlib.axes import Axes as Axes, Subplot as Subplot
from matplotlib.backend_bases import FigureCanvasBase as FigureCanvasBase, MouseButton as MouseButton
from matplotlib.cm import register_cmap as register_cmap
from matplotlib.colors import Normalize as Normalize
from matplotlib.figure import Figure as Figure, FigureBase as FigureBase, figaspect as figaspect
from matplotlib.gridspec import GridSpec as GridSpec, SubplotSpec as SubplotSpec
from matplotlib.lines import Line2D as Line2D
from matplotlib.patches import Arrow as Arrow, Circle as Circle, Polygon as Polygon, Rectangle as Rectangle
from matplotlib.projections import PolarAxes as PolarAxes
from matplotlib.scale import get_scale_names as get_scale_names
from matplotlib.text import Annotation as Annotation, Text as Text
from matplotlib.widgets import Button as Button, Slider as Slider, Widget as Widget
from numbers import Number as Number

def install_repl_displayhook() -> None: ...
def uninstall_repl_displayhook() -> None: ...

draw_all: Incomplete

def set_loglevel(*args, **kwargs): ...
def findobj(o: Incomplete | None = ..., match: Incomplete | None = ..., include_self: bool = ...): ...
def switch_backend(newbackend): ...
def new_figure_manager(*args, **kwargs): ...
def draw_if_interactive(*args, **kwargs): ...
def show(*args, **kwargs): ...
def isinteractive(): ...
def ioff(): ...
def ion(): ...
def pause(interval) -> None: ...
def rc(group, **kwargs) -> None: ...
def rc_context(rc: Incomplete | None = ..., fname: Incomplete | None = ...): ...
def rcdefaults() -> None: ...
def getp(obj, *args, **kwargs): ...
def get(obj, *args, **kwargs): ...
def setp(obj, *args, **kwargs): ...
def xkcd(scale: int = ..., length: int = ..., randomness: int = ...): ...
def figure(num: Incomplete | None = ..., figsize: Incomplete | None = ..., dpi: Incomplete | None = ..., facecolor: Incomplete | None = ..., edgecolor: Incomplete | None = ..., frameon: bool = ..., FigureClass=..., clear: bool = ..., **kwargs): ...
def gcf(): ...
def fignum_exists(num): ...
def get_fignums(): ...
def get_figlabels(): ...
def get_current_fig_manager(): ...
def connect(s, func): ...
def disconnect(cid): ...
def close(fig: Incomplete | None = ...) -> None: ...
def clf() -> None: ...
def draw() -> None: ...
def savefig(*args, **kwargs): ...
def figlegend(*args, **kwargs): ...
def axes(arg: Incomplete | None = ..., **kwargs): ...
def delaxes(ax: Incomplete | None = ...) -> None: ...
def sca(ax) -> None: ...
def cla(): ...
def subplot(*args, **kwargs): ...
def subplots(nrows: int = ..., ncols: int = ..., *, sharex: bool = ..., sharey: bool = ..., squeeze: bool = ..., width_ratios: Incomplete | None = ..., height_ratios: Incomplete | None = ..., subplot_kw: Incomplete | None = ..., gridspec_kw: Incomplete | None = ..., **fig_kw): ...
def subplot_mosaic(mosaic, *, sharex: bool = ..., sharey: bool = ..., width_ratios: Incomplete | None = ..., height_ratios: Incomplete | None = ..., empty_sentinel: str = ..., subplot_kw: Incomplete | None = ..., gridspec_kw: Incomplete | None = ..., **fig_kw): ...
def subplot2grid(shape, loc, rowspan: int = ..., colspan: int = ..., fig: Incomplete | None = ..., **kwargs): ...
def twinx(ax: Incomplete | None = ...): ...
def twiny(ax: Incomplete | None = ...): ...
def subplot_tool(targetfig: Incomplete | None = ...): ...
def box(on: Incomplete | None = ...) -> None: ...
def xlim(*args, **kwargs): ...
def ylim(*args, **kwargs): ...
def xticks(ticks: Incomplete | None = ..., labels: Incomplete | None = ..., *, minor: bool = ..., **kwargs): ...
def yticks(ticks: Incomplete | None = ..., labels: Incomplete | None = ..., *, minor: bool = ..., **kwargs): ...
def rgrids(radii: Incomplete | None = ..., labels: Incomplete | None = ..., angle: Incomplete | None = ..., fmt: Incomplete | None = ..., **kwargs): ...
def thetagrids(angles: Incomplete | None = ..., labels: Incomplete | None = ..., fmt: Incomplete | None = ..., **kwargs): ...
def get_plot_commands(): ...
def colorbar(mappable: Incomplete | None = ..., cax: Incomplete | None = ..., ax: Incomplete | None = ..., **kwargs): ...
def clim(vmin: Incomplete | None = ..., vmax: Incomplete | None = ...) -> None: ...
def get_cmap(name: Incomplete | None = ..., lut: Incomplete | None = ...): ...
def set_cmap(cmap) -> None: ...
def imread(fname, format: Incomplete | None = ...): ...
def imsave(fname, arr, **kwargs): ...
def matshow(A, fignum: Incomplete | None = ..., **kwargs): ...
def polar(*args, **kwargs): ...
def figimage(X, xo: int = ..., yo: int = ..., alpha: Incomplete | None = ..., norm: Incomplete | None = ..., cmap: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., origin: Incomplete | None = ..., resize: bool = ..., **kwargs): ...
def figtext(x, y, s, fontdict: Incomplete | None = ..., **kwargs): ...
def gca(): ...
def gci(): ...
def ginput(n: int = ..., timeout: int = ..., show_clicks: bool = ..., mouse_add=..., mouse_pop=..., mouse_stop=...): ...
def subplots_adjust(left: Incomplete | None = ..., bottom: Incomplete | None = ..., right: Incomplete | None = ..., top: Incomplete | None = ..., wspace: Incomplete | None = ..., hspace: Incomplete | None = ...): ...
def suptitle(t, **kwargs): ...
def tight_layout(*, pad: float = ..., h_pad: Incomplete | None = ..., w_pad: Incomplete | None = ..., rect: Incomplete | None = ...): ...
def waitforbuttonpress(timeout: int = ...): ...
def acorr(x, *, data: Incomplete | None = ..., **kwargs): ...
def angle_spectrum(x, Fs: Incomplete | None = ..., Fc: Incomplete | None = ..., window: Incomplete | None = ..., pad_to: Incomplete | None = ..., sides: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def annotate(text, xy, xytext: Incomplete | None = ..., xycoords: str = ..., textcoords: Incomplete | None = ..., arrowprops: Incomplete | None = ..., annotation_clip: Incomplete | None = ..., **kwargs): ...
def arrow(x, y, dx, dy, **kwargs): ...
def autoscale(enable: bool = ..., axis: str = ..., tight: Incomplete | None = ...): ...
def axhline(y: int = ..., xmin: int = ..., xmax: int = ..., **kwargs): ...
def axhspan(ymin, ymax, xmin: int = ..., xmax: int = ..., **kwargs): ...
def axis(*args, emit: bool = ..., **kwargs): ...
def axline(xy1, xy2: Incomplete | None = ..., *, slope: Incomplete | None = ..., **kwargs): ...
def axvline(x: int = ..., ymin: int = ..., ymax: int = ..., **kwargs): ...
def axvspan(xmin, xmax, ymin: int = ..., ymax: int = ..., **kwargs): ...
def bar(x, height, width: float = ..., bottom: Incomplete | None = ..., *, align: str = ..., data: Incomplete | None = ..., **kwargs): ...
def barbs(*args, data: Incomplete | None = ..., **kwargs): ...
def barh(y, width, height: float = ..., left: Incomplete | None = ..., *, align: str = ..., data: Incomplete | None = ..., **kwargs): ...
def bar_label(container, labels: Incomplete | None = ..., *, fmt: str = ..., label_type: str = ..., padding: int = ..., **kwargs): ...
def boxplot(x, notch: Incomplete | None = ..., sym: Incomplete | None = ..., vert: Incomplete | None = ..., whis: Incomplete | None = ..., positions: Incomplete | None = ..., widths: Incomplete | None = ..., patch_artist: Incomplete | None = ..., bootstrap: Incomplete | None = ..., usermedians: Incomplete | None = ..., conf_intervals: Incomplete | None = ..., meanline: Incomplete | None = ..., showmeans: Incomplete | None = ..., showcaps: Incomplete | None = ..., showbox: Incomplete | None = ..., showfliers: Incomplete | None = ..., boxprops: Incomplete | None = ..., labels: Incomplete | None = ..., flierprops: Incomplete | None = ..., medianprops: Incomplete | None = ..., meanprops: Incomplete | None = ..., capprops: Incomplete | None = ..., whiskerprops: Incomplete | None = ..., manage_ticks: bool = ..., autorange: bool = ..., zorder: Incomplete | None = ..., capwidths: Incomplete | None = ..., *, data: Incomplete | None = ...): ...
def broken_barh(xranges, yrange, *, data: Incomplete | None = ..., **kwargs): ...
def clabel(CS, levels: Incomplete | None = ..., **kwargs): ...
def cohere(x, y, NFFT: int = ..., Fs: int = ..., Fc: int = ..., detrend=..., window=..., noverlap: int = ..., pad_to: Incomplete | None = ..., sides: str = ..., scale_by_freq: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def contour(*args, data: Incomplete | None = ..., **kwargs): ...
def contourf(*args, data: Incomplete | None = ..., **kwargs): ...
def csd(x, y, NFFT: Incomplete | None = ..., Fs: Incomplete | None = ..., Fc: Incomplete | None = ..., detrend: Incomplete | None = ..., window: Incomplete | None = ..., noverlap: Incomplete | None = ..., pad_to: Incomplete | None = ..., sides: Incomplete | None = ..., scale_by_freq: Incomplete | None = ..., return_line: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def errorbar(x, y, yerr: Incomplete | None = ..., xerr: Incomplete | None = ..., fmt: str = ..., ecolor: Incomplete | None = ..., elinewidth: Incomplete | None = ..., capsize: Incomplete | None = ..., barsabove: bool = ..., lolims: bool = ..., uplims: bool = ..., xlolims: bool = ..., xuplims: bool = ..., errorevery: int = ..., capthick: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def eventplot(positions, orientation: str = ..., lineoffsets: int = ..., linelengths: int = ..., linewidths: Incomplete | None = ..., colors: Incomplete | None = ..., linestyles: str = ..., *, data: Incomplete | None = ..., **kwargs): ...
def fill(*args, data: Incomplete | None = ..., **kwargs): ...
def fill_between(x, y1, y2: int = ..., where: Incomplete | None = ..., interpolate: bool = ..., step: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def fill_betweenx(y, x1, x2: int = ..., where: Incomplete | None = ..., step: Incomplete | None = ..., interpolate: bool = ..., *, data: Incomplete | None = ..., **kwargs): ...
def grid(visible: Incomplete | None = ..., which: str = ..., axis: str = ..., **kwargs): ...
def hexbin(x, y, C: Incomplete | None = ..., gridsize: int = ..., bins: Incomplete | None = ..., xscale: str = ..., yscale: str = ..., extent: Incomplete | None = ..., cmap: Incomplete | None = ..., norm: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., alpha: Incomplete | None = ..., linewidths: Incomplete | None = ..., edgecolors: str = ..., reduce_C_function=..., mincnt: Incomplete | None = ..., marginals: bool = ..., *, data: Incomplete | None = ..., **kwargs): ...
def hist(x, bins: Incomplete | None = ..., range: Incomplete | None = ..., density: bool = ..., weights: Incomplete | None = ..., cumulative: bool = ..., bottom: Incomplete | None = ..., histtype: str = ..., align: str = ..., orientation: str = ..., rwidth: Incomplete | None = ..., log: bool = ..., color: Incomplete | None = ..., label: Incomplete | None = ..., stacked: bool = ..., *, data: Incomplete | None = ..., **kwargs): ...
def stairs(values, edges: Incomplete | None = ..., *, orientation: str = ..., baseline: int = ..., fill: bool = ..., data: Incomplete | None = ..., **kwargs): ...
def hist2d(x, y, bins: int = ..., range: Incomplete | None = ..., density: bool = ..., weights: Incomplete | None = ..., cmin: Incomplete | None = ..., cmax: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def hlines(y, xmin, xmax, colors: Incomplete | None = ..., linestyles: str = ..., label: str = ..., *, data: Incomplete | None = ..., **kwargs): ...
def imshow(X, cmap: Incomplete | None = ..., norm: Incomplete | None = ..., aspect: Incomplete | None = ..., interpolation: Incomplete | None = ..., alpha: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., origin: Incomplete | None = ..., extent: Incomplete | None = ..., *, interpolation_stage: Incomplete | None = ..., filternorm: bool = ..., filterrad: float = ..., resample: Incomplete | None = ..., url: Incomplete | None = ..., data: Incomplete | None = ..., **kwargs): ...
def legend(*args, **kwargs): ...
def locator_params(axis: str = ..., tight: Incomplete | None = ..., **kwargs): ...
def loglog(*args, **kwargs): ...
def magnitude_spectrum(x, Fs: Incomplete | None = ..., Fc: Incomplete | None = ..., window: Incomplete | None = ..., pad_to: Incomplete | None = ..., sides: Incomplete | None = ..., scale: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def margins(*margins, x: Incomplete | None = ..., y: Incomplete | None = ..., tight: bool = ...): ...
def minorticks_off(): ...
def minorticks_on(): ...
def pcolor(*args, shading: Incomplete | None = ..., alpha: Incomplete | None = ..., norm: Incomplete | None = ..., cmap: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., data: Incomplete | None = ..., **kwargs): ...
def pcolormesh(*args, alpha: Incomplete | None = ..., norm: Incomplete | None = ..., cmap: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., shading: Incomplete | None = ..., antialiased: bool = ..., data: Incomplete | None = ..., **kwargs): ...
def phase_spectrum(x, Fs: Incomplete | None = ..., Fc: Incomplete | None = ..., window: Incomplete | None = ..., pad_to: Incomplete | None = ..., sides: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def pie(x, explode: Incomplete | None = ..., labels: Incomplete | None = ..., colors: Incomplete | None = ..., autopct: Incomplete | None = ..., pctdistance: float = ..., shadow: bool = ..., labeldistance: float = ..., startangle: int = ..., radius: int = ..., counterclock: bool = ..., wedgeprops: Incomplete | None = ..., textprops: Incomplete | None = ..., center=..., frame: bool = ..., rotatelabels: bool = ..., *, normalize: bool = ..., data: Incomplete | None = ...): ...
def plot(*args, scalex: bool = ..., scaley: bool = ..., data: Incomplete | None = ..., **kwargs): ...
def plot_date(x, y, fmt: str = ..., tz: Incomplete | None = ..., xdate: bool = ..., ydate: bool = ..., *, data: Incomplete | None = ..., **kwargs): ...
def psd(x, NFFT: Incomplete | None = ..., Fs: Incomplete | None = ..., Fc: Incomplete | None = ..., detrend: Incomplete | None = ..., window: Incomplete | None = ..., noverlap: Incomplete | None = ..., pad_to: Incomplete | None = ..., sides: Incomplete | None = ..., scale_by_freq: Incomplete | None = ..., return_line: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def quiver(*args, data: Incomplete | None = ..., **kwargs): ...
def quiverkey(Q, X, Y, U, label, **kwargs): ...
def scatter(x, y, s: Incomplete | None = ..., c: Incomplete | None = ..., marker: Incomplete | None = ..., cmap: Incomplete | None = ..., norm: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., alpha: Incomplete | None = ..., linewidths: Incomplete | None = ..., *, edgecolors: Incomplete | None = ..., plotnonfinite: bool = ..., data: Incomplete | None = ..., **kwargs): ...
def semilogx(*args, **kwargs): ...
def semilogy(*args, **kwargs): ...
def specgram(x, NFFT: Incomplete | None = ..., Fs: Incomplete | None = ..., Fc: Incomplete | None = ..., detrend: Incomplete | None = ..., window: Incomplete | None = ..., noverlap: Incomplete | None = ..., cmap: Incomplete | None = ..., xextent: Incomplete | None = ..., pad_to: Incomplete | None = ..., sides: Incomplete | None = ..., scale_by_freq: Incomplete | None = ..., mode: Incomplete | None = ..., scale: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., *, data: Incomplete | None = ..., **kwargs): ...
def spy(Z, precision: int = ..., marker: Incomplete | None = ..., markersize: Incomplete | None = ..., aspect: str = ..., origin: str = ..., **kwargs): ...
def stackplot(x, *args, labels=..., colors: Incomplete | None = ..., baseline: str = ..., data: Incomplete | None = ..., **kwargs): ...
def stem(*args, linefmt: Incomplete | None = ..., markerfmt: Incomplete | None = ..., basefmt: Incomplete | None = ..., bottom: int = ..., label: Incomplete | None = ..., use_line_collection=..., orientation: str = ..., data: Incomplete | None = ...): ...
def step(x, y, *args, where: str = ..., data: Incomplete | None = ..., **kwargs): ...
def streamplot(x, y, u, v, density: int = ..., linewidth: Incomplete | None = ..., color: Incomplete | None = ..., cmap: Incomplete | None = ..., norm: Incomplete | None = ..., arrowsize: int = ..., arrowstyle: str = ..., minlength: float = ..., transform: Incomplete | None = ..., zorder: Incomplete | None = ..., start_points: Incomplete | None = ..., maxlength: float = ..., integration_direction: str = ..., broken_streamlines: bool = ..., *, data: Incomplete | None = ...): ...
def table(cellText: Incomplete | None = ..., cellColours: Incomplete | None = ..., cellLoc: str = ..., colWidths: Incomplete | None = ..., rowLabels: Incomplete | None = ..., rowColours: Incomplete | None = ..., rowLoc: str = ..., colLabels: Incomplete | None = ..., colColours: Incomplete | None = ..., colLoc: str = ..., loc: str = ..., bbox: Incomplete | None = ..., edges: str = ..., **kwargs): ...
def text(x, y, s, fontdict: Incomplete | None = ..., **kwargs): ...
def tick_params(axis: str = ..., **kwargs): ...
def ticklabel_format(*, axis: str = ..., style: str = ..., scilimits: Incomplete | None = ..., useOffset: Incomplete | None = ..., useLocale: Incomplete | None = ..., useMathText: Incomplete | None = ...): ...
def tricontour(*args, **kwargs): ...
def tricontourf(*args, **kwargs): ...
def tripcolor(*args, alpha: float = ..., norm: Incomplete | None = ..., cmap: Incomplete | None = ..., vmin: Incomplete | None = ..., vmax: Incomplete | None = ..., shading: str = ..., facecolors: Incomplete | None = ..., **kwargs): ...
def triplot(*args, **kwargs): ...
def violinplot(dataset, positions: Incomplete | None = ..., vert: bool = ..., widths: float = ..., showmeans: bool = ..., showextrema: bool = ..., showmedians: bool = ..., quantiles: Incomplete | None = ..., points: int = ..., bw_method: Incomplete | None = ..., *, data: Incomplete | None = ...): ...
def vlines(x, ymin, ymax, colors: Incomplete | None = ..., linestyles: str = ..., label: str = ..., *, data: Incomplete | None = ..., **kwargs): ...
def xcorr(x, y, normed: bool = ..., detrend=..., usevlines: bool = ..., maxlags: int = ..., *, data: Incomplete | None = ..., **kwargs): ...
def sci(im): ...
def title(label, fontdict: Incomplete | None = ..., loc: Incomplete | None = ..., pad: Incomplete | None = ..., *, y: Incomplete | None = ..., **kwargs): ...
def xlabel(xlabel, fontdict: Incomplete | None = ..., labelpad: Incomplete | None = ..., *, loc: Incomplete | None = ..., **kwargs): ...
def ylabel(ylabel, fontdict: Incomplete | None = ..., labelpad: Incomplete | None = ..., *, loc: Incomplete | None = ..., **kwargs): ...
def xscale(value, **kwargs): ...
def yscale(value, **kwargs): ...
def autumn() -> None: ...
def bone() -> None: ...
def cool() -> None: ...
def copper() -> None: ...
def flag() -> None: ...
def gray() -> None: ...
def hot() -> None: ...
def hsv() -> None: ...
def jet() -> None: ...
def pink() -> None: ...
def prism() -> None: ...
def spring() -> None: ...
def summer() -> None: ...
def winter() -> None: ...
def magma() -> None: ...
def inferno() -> None: ...
def plasma() -> None: ...
def viridis() -> None: ...
def nipy_spectral() -> None: ...
