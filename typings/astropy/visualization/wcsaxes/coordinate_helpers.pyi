"""
This type stub file was generated by pyright.
"""

"""
This file defines the classes used to represent a 'coordinate', which includes
axes, ticks, tick labels, and grid lines.
"""
__all__ = ['CoordinateHelper']
LINES_TO_PATCHES_LINESTYLE = ...
def wrap_angle_at(values, coord_wrap):
    ...

class CoordinateHelper:
    """
    Helper class to control one of the coordinates in the
    :class:`~astropy.visualization.wcsaxes.WCSAxes`.

    Parameters
    ----------
    parent_axes : :class:`~astropy.visualization.wcsaxes.WCSAxes`
        The axes the coordinate helper belongs to.
    parent_map : :class:`~astropy.visualization.wcsaxes.CoordinatesMap`
        The :class:`~astropy.visualization.wcsaxes.CoordinatesMap` object this
        coordinate belongs to.
    transform : `~matplotlib.transforms.Transform`
        The transform corresponding to this coordinate system.
    coord_index : int
        The index of this coordinate in the
        :class:`~astropy.visualization.wcsaxes.CoordinatesMap`.
    coord_type : {'longitude', 'latitude', 'scalar'}
        The type of this coordinate, which is used to determine the wrapping and
        boundary behavior of coordinates. Longitudes wrap at ``coord_wrap``,
        latitudes have to be in the range -90 to 90, and scalars are unbounded
        and do not wrap.
    coord_unit : `~astropy.units.Unit`
        The unit that this coordinate is in given the output of transform.
    format_unit : `~astropy.units.Unit`, optional
        The unit to use to display the coordinates.
    coord_wrap : float
        The angle at which the longitude wraps (defaults to 360)
    frame : `~astropy.visualization.wcsaxes.frame.BaseFrame`
        The frame of the :class:`~astropy.visualization.wcsaxes.WCSAxes`.
    """
    def __init__(self, parent_axes=..., parent_map=..., transform=..., coord_index=..., coord_type=..., coord_unit=..., coord_wrap=..., frame=..., format_unit=..., default_label=...) -> None:
        ...
    
    def grid(self, draw_grid=..., grid_type=..., **kwargs): # -> None:
        """
        Plot grid lines for this coordinate.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        grid_type : {'lines', 'contours'}
            Whether to plot the contours by determining the grid lines in
            world coordinates and then plotting them in world coordinates
            (``'lines'``) or by determining the world coordinates at many
            positions in the image and then drawing contours
            (``'contours'``). The first is recommended for 2-d images, while
            for 3-d (or higher dimensional) cubes, the ``'contours'`` option
            is recommended. By default, 'lines' is used if the transform has
            an inverse, otherwise 'contours' is used.
        """
        ...
    
    def set_coord_type(self, coord_type, coord_wrap=...): # -> None:
        """
        Set the coordinate type for the axis.

        Parameters
        ----------
        coord_type : str
            One of 'longitude', 'latitude' or 'scalar'
        coord_wrap : float, optional
            The value to wrap at for angular coordinates
        """
        ...
    
    def set_major_formatter(self, formatter): # -> None:
        """
        Set the formatter to use for the major tick labels.

        Parameters
        ----------
        formatter : str or `~matplotlib.ticker.Formatter`
            The format or formatter to use.
        """
        ...
    
    def format_coord(self, value, format=...): # -> Any | str:
        """
        Given the value of a coordinate, will format it according to the
        format of the formatter_locator.

        Parameters
        ----------
        value : float
            The value to format
        format : {'auto', 'ascii', 'latex'}, optional
            The format to use - by default the formatting will be adjusted
            depending on whether Matplotlib is using LaTeX or MathTex. To
            get plain ASCII strings, use format='ascii'.
        """
        ...
    
    def set_separator(self, separator): # -> None:
        """
        Set the separator to use for the angle major tick labels.

        Parameters
        ----------
        separator : str or tuple or None
            The separator between numbers in sexagesimal representation. Can be
            either a string or a tuple (or `None` for default).
        """
        ...
    
    def set_format_unit(self, unit, decimal=..., show_decimal_unit=...): # -> None:
        """
        Set the unit for the major tick labels.

        Parameters
        ----------
        unit : class:`~astropy.units.Unit`
            The unit to which the tick labels should be converted to.
        decimal : bool, optional
            Whether to use decimal formatting. By default this is `False`
            for degrees or hours (which therefore use sexagesimal formatting)
            and `True` for all other units.
        show_decimal_unit : bool, optional
            Whether to include units when in decimal mode.
        """
        ...
    
    def get_format_unit(self): # -> Unit | None:
        """
        Get the unit for the major tick labels.
        """
        ...
    
    def set_ticks(self, values=..., spacing=..., number=..., size=..., width=..., color=..., alpha=..., direction=..., exclude_overlapping=...): # -> None:
        """
        Set the location and properties of the ticks.

        At most one of the options from ``values``, ``spacing``, or
        ``number`` can be specified.

        Parameters
        ----------
        values : iterable, optional
            The coordinate values at which to show the ticks.
        spacing : float, optional
            The spacing between ticks.
        number : float, optional
            The approximate number of ticks shown.
        size : float, optional
            The length of the ticks in points
        color : str or tuple, optional
            A valid Matplotlib color for the ticks
        alpha : float, optional
            The alpha value (transparency) for the ticks.
        direction : {'in','out'}, optional
            Whether the ticks should point inwards or outwards.
        """
        ...
    
    def set_ticks_position(self, position): # -> None:
        """
        Set where ticks should appear

        Parameters
        ----------
        position : str
            The axes on which the ticks for this coordinate should appear.
            Should be a string containing zero or more of ``'b'``, ``'t'``,
            ``'l'``, ``'r'``. For example, ``'lb'`` will lead the ticks to be
            shown on the left and bottom axis.
        """
        ...
    
    def set_ticks_visible(self, visible): # -> None:
        """
        Set whether ticks are visible or not.

        Parameters
        ----------
        visible : bool
            The visibility of ticks. Setting as ``False`` will hide ticks
            along this coordinate.
        """
        ...
    
    def set_ticklabel(self, color=..., size=..., pad=..., exclude_overlapping=..., **kwargs): # -> None:
        """
        Set the visual properties for the tick labels.

        Parameters
        ----------
        size : float, optional
            The size of the ticks labels in points
        color : str or tuple, optional
            A valid Matplotlib color for the tick labels
        pad : float, optional
            Distance in points between tick and label.
        exclude_overlapping : bool, optional
            Whether to exclude tick labels that overlap over each other.
        **kwargs
            Other keyword arguments are passed to :class:`matplotlib.text.Text`.
        """
        ...
    
    def set_ticklabel_position(self, position): # -> None:
        """
        Set where tick labels should appear

        Parameters
        ----------
        position : str
            The axes on which the tick labels for this coordinate should
            appear. Should be a string containing zero or more of ``'b'``,
            ``'t'``, ``'l'``, ``'r'``. For example, ``'lb'`` will lead the
            tick labels to be shown on the left and bottom axis.
        """
        ...
    
    def set_ticklabel_visible(self, visible): # -> None:
        """
        Set whether the tick labels are visible or not.

        Parameters
        ----------
        visible : bool
            The visibility of ticks. Setting as ``False`` will hide this
            coordinate's tick labels.
        """
        ...
    
    def set_axislabel(self, text, minpad=..., **kwargs): # -> None:
        """
        Set the text and optionally visual properties for the axis label.

        Parameters
        ----------
        text : str
            The axis label text.
        minpad : float, optional
            The padding for the label in terms of axis label font size.
        **kwargs
            Keywords are passed to :class:`matplotlib.text.Text`. These
            can include keywords to set the ``color``, ``size``, ``weight``, and
            other text properties.
        """
        ...
    
    def get_axislabel(self): # -> str:
        """
        Get the text for the axis label

        Returns
        -------
        label : str
            The axis label
        """
        ...
    
    def set_auto_axislabel(self, auto_label): # -> None:
        """
        Render default axis labels if no explicit label is provided.

        Parameters
        ----------
        auto_label : `bool`
            `True` if default labels will be rendered.
        """
        ...
    
    def get_auto_axislabel(self): # -> bool:
        """
        Render default axis labels if no explicit label is provided.

        Returns
        -------
        auto_axislabel : `bool`
            `True` if default labels will be rendered.
        """
        ...
    
    def set_axislabel_position(self, position): # -> None:
        """
        Set where axis labels should appear

        Parameters
        ----------
        position : str
            The axes on which the axis label for this coordinate should
            appear. Should be a string containing zero or more of ``'b'``,
            ``'t'``, ``'l'``, ``'r'``. For example, ``'lb'`` will lead the
            axis label to be shown on the left and bottom axis.
        """
        ...
    
    def set_axislabel_visibility_rule(self, rule): # -> None:
        """
        Set the rule used to determine when the axis label is drawn.

        Parameters
        ----------
        rule : str
            If the rule is 'always' axis labels will always be drawn on the
            axis. If the rule is 'ticks' the label will only be drawn if ticks
            were drawn on that axis. If the rule is 'labels' the axis label
            will only be drawn if tick labels were drawn on that axis.
        """
        ...
    
    def get_axislabel_visibility_rule(self, rule): # -> str:
        """
        Get the rule used to determine when the axis label is drawn.
        """
        ...
    
    @property
    def locator(self): # -> (value_min: Unknown, value_max: Unknown) -> (tuple[Quantity, Unknown] | tuple[Unknown, Unknown]):
        ...
    
    @property
    def formatter(self): # -> ((values: Unknown, spacing: Unknown, format: str = 'auto') -> (list[Unknown] | list[str])) | ((values: Unknown, spacing: Unknown, format: str = 'auto') -> (Any | list[Unknown])):
        ...
    
    def display_minor_ticks(self, display_minor_ticks): # -> None:
        """
        Display minor ticks for this coordinate.

        Parameters
        ----------
        display_minor_ticks : bool
            Whether or not to display minor ticks.
        """
        ...
    
    def get_minor_frequency(self): # -> int:
        ...
    
    def set_minor_frequency(self, frequency): # -> None:
        """
        Set the frequency of minor ticks per major ticks.

        Parameters
        ----------
        frequency : int
            The number of minor ticks per major ticks.
        """
        ...
    
    def tick_params(self, which=..., **kwargs): # -> None:
        """
        Method to set the tick and tick label parameters in the same way as the
        :meth:`~matplotlib.axes.Axes.tick_params` method in Matplotlib.

        This is provided for convenience, but the recommended API is to use
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticks`,
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticklabel`,
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticks_position`,
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticklabel_position`,
        and :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.grid`.

        Parameters
        ----------
        which : {'both', 'major', 'minor'}, optional
            Which ticks to apply the settings to. By default, setting are
            applied to both major and minor ticks. Note that if ``'minor'`` is
            specified, only the length of the ticks can be set currently.
        direction : {'in', 'out'}, optional
            Puts ticks inside the axes, or outside the axes.
        length : float, optional
            Tick length in points.
        width : float, optional
            Tick width in points.
        color : color, optional
            Tick color (accepts any valid Matplotlib color)
        pad : float, optional
            Distance in points between tick and label.
        labelsize : float or str, optional
            Tick label font size in points or as a string (e.g., 'large').
        labelcolor : color, optional
            Tick label color (accepts any valid Matplotlib color)
        colors : color, optional
            Changes the tick color and the label color to the same value
             (accepts any valid Matplotlib color).
        bottom, top, left, right : bool, optional
            Where to draw the ticks. Note that this will not work correctly if
            the frame is not rectangular.
        labelbottom, labeltop, labelleft, labelright : bool, optional
            Where to draw the tick labels. Note that this will not work
            correctly if the frame is not rectangular.
        grid_color : color, optional
            The color of the grid lines (accepts any valid Matplotlib color).
        grid_alpha : float, optional
            Transparency of grid lines: 0 (transparent) to 1 (opaque).
        grid_linewidth : float, optional
            Width of grid lines in points.
        grid_linestyle : str, optional
            The style of the grid lines (accepts any valid Matplotlib line
            style).
        """
        ...
    


