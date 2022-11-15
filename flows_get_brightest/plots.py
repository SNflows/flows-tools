import copy
import logging
import astropy.visualization as viz
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ZScaleInterval
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .observer import Observer
from .utils import numeric


class Plotter:
    """
    Takes an observer with WCS, Target, Plan, Instrument, Regions, and Corners and makes a finding chart.
    """

    def __init__(self, obs: Observer):
        self.obs = obs

    def plot(self):
        """not implemented yet"""
        pass

    def make_finding_chart(self, plot_refcat: bool = True, plot_simbad: bool = True,
                           savefig: bool = True, radius: numeric = 14) -> matplotlib.axes.Axes:
        """Make finding chart for a given observation."""
        obs = self.obs
        image = obs.get_image(radius=radius)
        zscale = ZScaleInterval()
        vmin, vmax = zscale.get_limits(image.data.flat)
        obs.wcs = obs.get_wcs(header=image.header)
        regions = obs.regions_to_physical()

        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': obs.wcs})

        for i, region in enumerate(regions):
            region.plot(ax=ax, edgecolor='cyan', linestyle='-.', label=obs.ins.region_names[i])

        self.plot_image(image.data, ax=ax, cmap='viridis', scale='linear', vmin=vmin, vmax=vmax)
        tar_pix = obs.wcs.all_world2pix(obs.target.ra, obs.target.dec, 0)
        ax.scatter(tar_pix[0], tar_pix[1], marker='*', s=250, label='SN', color='orange')

        if plot_refcat:
            refcat_stars = obs.refcat_coords.to_pixel(obs.wcs)
            ax.scatter(refcat_stars[0], refcat_stars[1], facecolors='none', edgecolors='red', zorder=4, alpha=0.8,
                       marker='o', s=100, label='Refcat2')
        if plot_simbad:
            simbad_stars = obs.simbad_coords.to_pixel(obs.wcs)
            ax.scatter(simbad_stars[0], simbad_stars[1], facecolors='none', edgecolors='orange', zorder=5, alpha=0.8,
                       marker='o', s=100, label='Simbad')

        ax.legend(fontsize=15)
        ax.set_title(f"{obs.target.info['target_name']} {obs.ins.__class__.__name__} FC")
        if savefig:
            fig.savefig(f"{obs.target.info['target_name']}_{obs.ins.__class__.__name__}_finding_chart.png")
        plt.show()
        return ax

    # Copied from flows image plotter.
    @staticmethod
    def plot_image(image, ax=None, scale='log', cmap=None, origin='lower', xlabel=None, ylabel=None, cbar=None,
                   clabel='Flux ($e^{-}s^{-1}$)', cbar_ticks=None, cbar_ticklabels=None, cbar_pad=None, cbar_size='4%',
                   title=None, percentile=95.0, vmin=None, vmax=None, offset_axes=None, color_bad='k', **kwargs):
        """
        Utility function to plot a 2D image.
        Parameters:
            image (2d array): Image data.
            ax (matplotlib.pyplot.axes, optional): Axes in which to plot.
                Default (None) is to use current active axes.
            scale (str or :py:class:`astropy.visualization.ImageNormalize` object, optional):
                Normalization used to stretch the colormap.
                Options: ``'linear'``, ``'sqrt'``, ``'log'``, ``'asinh'``, ``'histeq'``, ``'sinh'``
                and ``'squared'``.
                Can also be a :py:class:`astropy.visualization.ImageNormalize` object.
                Default is ``'log'``.
            origin (str, optional): The origin of the coordinate system.
            xlabel (str, optional): Label for the x-axis.
            ylabel (str, optional): Label for the y-axis.
            cbar (string, optional): Location of color bar.
                Choises are ``'right'``, ``'left'``, ``'top'``, ``'bottom'``.
                Default is not to create colorbar.
            clabel (str, optional): Label for the color bar.
            cbar_size (float, optional): Fractional size of colorbar compared to axes. Default='4%'.
            cbar_pad (float, optional): Padding between axes and colorbar.
            title (str or None, optional): Title for the plot.
            percentile (float, optional): The fraction of pixels to keep in color-trim.
                If single float given, the same fraction of pixels is eliminated from both ends.
                If tuple of two floats is given, the two are used as the percentiles.
                Default=95.
            cmap (matplotlib colormap, optional): Colormap to use. Default is the ``Blues`` colormap.
            vmin (float, optional): Lower limit to use for colormap.
            vmax (float, optional): Upper limit to use for colormap.
            color_bad (str, optional): Color to apply to bad pixels (NaN). Default is black.
            kwargs (dict, optional): Keyword arguments to be passed to :py:func:`matplotlib.pyplot.imshow`.
        Returns:
            :py:class:`matplotlib.image.AxesImage`: Image from returned
                by :py:func:`matplotlib.pyplot.imshow`.
        .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
        """

        logger = logging.getLogger(__name__)

        # Backward compatible settings:
        make_cbar = kwargs.pop('make_cbar', None)
        if make_cbar:
            raise FutureWarning("'make_cbar' is deprecated. Use 'cbar' instead.")
            if not cbar:
                cbar = make_cbar

        # Special treatment for boolean arrays:
        if isinstance(image, np.ndarray) and image.dtype == 'bool':
            if vmin is None: vmin = 0
            if vmax is None: vmax = 1
            if cbar_ticks is None: cbar_ticks = [0, 1]
            if cbar_ticklabels is None: cbar_ticklabels = ['False', 'True']

        # Calculate limits of color scaling:
        interval = None
        if vmin is None or vmax is None:
            if np.isnan(image).all():
                logger.warning("Image is all NaN")
                vmin = 0
                vmax = 1
                if cbar_ticks is None:
                    cbar_ticks = []
                if cbar_ticklabels is None:
                    cbar_ticklabels = []
            elif isinstance(percentile, (list, tuple, np.ndarray)):
                interval = viz.AsymmetricPercentileInterval(percentile[0], percentile[1])
            else:
                interval = viz.PercentileInterval(percentile)

        # Create ImageNormalize object with extracted limits:
        if scale in ('log', 'linear', 'sqrt', 'asinh', 'histeq', 'sinh', 'squared'):
            if scale == 'log':
                stretch = viz.LogStretch()
            elif scale == 'linear':
                stretch = viz.LinearStretch()
            elif scale == 'sqrt':
                stretch = viz.SqrtStretch()
            elif scale == 'asinh':
                stretch = viz.AsinhStretch()
            elif scale == 'histeq':
                stretch = viz.HistEqStretch(image[np.isfinite(image)])
            elif scale == 'sinh':
                stretch = viz.SinhStretch()
            elif scale == 'squared':
                stretch = viz.SquaredStretch()

            # Create ImageNormalize object. Very important to use clip=False here, otherwise
            # NaN points will not be plotted correctly.
            norm = viz.ImageNormalize(data=image, interval=interval, vmin=vmin, vmax=vmax, stretch=stretch, clip=False)

        elif isinstance(scale, (viz.ImageNormalize, matplotlib.colors.Normalize)):
            norm = scale
        else:
            raise ValueError("scale {} is not available.".format(scale))

        if offset_axes:
            extent = (offset_axes[0] - 0.5, offset_axes[0] + image.shape[1] - 0.5, offset_axes[1] - 0.5,
                      offset_axes[1] + image.shape[0] - 0.5)
        else:
            extent = (-0.5, image.shape[1] - 0.5, -0.5, image.shape[0] - 0.5)

        if ax is None:
            ax = plt.gca()

        # Set up the colormap to use. If a bad color is defined,
        # add it to the colormap:
        if cmap is None:
            cmap = copy.copy(plt.get_cmap('Blues'))
        elif isinstance(cmap, str):
            cmap = copy.copy(plt.get_cmap(cmap))

        if color_bad:
            cmap.set_bad(color_bad, 1.0)

        im = ax.imshow(image, cmap=cmap, norm=norm, origin=origin, extent=extent, interpolation='nearest', **kwargs)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if title is not None:
            ax.set_title(title)
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])

        if cbar:
            fig = ax.figure
            divider = make_axes_locatable(ax)
            if cbar == 'top':
                cbar_pad = 0.05 if cbar_pad is None else cbar_pad
                cax = divider.append_axes('top', size=cbar_size, pad=cbar_pad)
                orientation = 'horizontal'
            elif cbar == 'bottom':
                cbar_pad = 0.35 if cbar_pad is None else cbar_pad
                cax = divider.append_axes('bottom', size=cbar_size, pad=cbar_pad)
                orientation = 'horizontal'
            elif cbar == 'left':
                cbar_pad = 0.35 if cbar_pad is None else cbar_pad
                cax = divider.append_axes('left', size=cbar_size, pad=cbar_pad)
                orientation = 'vertical'
            else:
                cbar_pad = 0.05 if cbar_pad is None else cbar_pad
                cax = divider.append_axes('right', size=cbar_size, pad=cbar_pad)
                orientation = 'vertical'

            cb = fig.colorbar(im, cax=cax, orientation=orientation)

            if cbar == 'top':
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
            elif cbar == 'left':
                cax.yaxis.set_ticks_position('left')
                cax.yaxis.set_label_position('left')

            if clabel is not None:
                cb.set_label(clabel)
            if cbar_ticks is not None:
                cb.set_ticks(cbar_ticks)
            if cbar_ticklabels is not None:
                cb.set_ticklabels(cbar_ticklabels)

            # cax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            # cax.yaxis.set_minor_locator(matplotlib.ticker.AutoLocator())
            cax.tick_params(which='both', direction='out', pad=5)

        # Settings for ticks:
        integer_locator = MaxNLocator(nbins=10, integer=True)
        ax.xaxis.set_major_locator(integer_locator)
        ax.xaxis.set_minor_locator(integer_locator)
        ax.yaxis.set_major_locator(integer_locator)
        ax.yaxis.set_minor_locator(integer_locator)
        ax.tick_params(which='both', direction='out', pad=5)
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()

        return im
