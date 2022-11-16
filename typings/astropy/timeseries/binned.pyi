"""
This type stub file was generated by pyright.
"""

from astropy.timeseries.core import BaseTimeSeries, autocheck_required_columns

__all__ = ['BinnedTimeSeries']
@autocheck_required_columns
class BinnedTimeSeries(BaseTimeSeries):
    """
    A class to represent binned time series data in tabular form.

    `~astropy.timeseries.BinnedTimeSeries` provides a class for
    representing time series as a collection of values of different
    quantities measured in time bins (for time series with values
    sampled at specific times, see the `~astropy.timeseries.TimeSeries`
    class). `~astropy.timeseries.BinnedTimeSeries` is a sub-class of
    `~astropy.table.QTable` and thus provides all the standard table
    maniplation methods available to tables, but it also provides
    additional conveniences for dealing with time series, such as a
    flexible initializer for setting up the times, and attributes to
    access the start/center/end time of bins.

    See also: https://docs.astropy.org/en/stable/timeseries/

    Parameters
    ----------
    data : numpy ndarray, dict, list, table-like object, optional
        Data to initialize time series. This does not need to contain the
        times, which can be provided separately, but if it does contain the
        times they should be in columns called ``'time_bin_start'`` and
        ``'time_bin_size'`` to be automatically recognized.
    time_bin_start : `~astropy.time.Time` or iterable
        The times of the start of each bin - this can be either given
        directly as a `~astropy.time.Time` array or as any iterable that
        initializes the `~astropy.time.Time` class. If this is given, then
        the remaining time-related arguments should not be used. This can also
        be a scalar value if ``time_bin_size`` is provided.
    time_bin_end : `~astropy.time.Time` or iterable
        The times of the end of each bin - this can be either given directly
        as a `~astropy.time.Time` array or as any value or iterable that
        initializes the `~astropy.time.Time` class. If this is given, then the
        remaining time-related arguments should not be used. This can only be
        given if ``time_bin_start`` is an array of values. If ``time_bin_end``
        is a scalar, time bins are assumed to be contiguous, such that the end
        of each bin is the start of the next one, and ``time_bin_end`` gives
        the end time for the last bin. If ``time_bin_end`` is an array, the
        time bins do not need to be contiguous. If this argument is provided,
        ``time_bin_size`` should not be provided.
    time_bin_size : `~astropy.time.TimeDelta` or `~astropy.units.Quantity`
        The size of the time bins, either as a scalar value (in which case all
        time bins will be assumed to have the same duration) or as an array of
        values (in which case each time bin can have a different duration).
        If this argument is provided, ``time_bin_end`` should not be provided.
    n_bins : int
        The number of time bins for the series. This is only used if both
        ``time_bin_start`` and ``time_bin_size`` are provided and are scalar
        values.
    **kwargs : dict, optional
        Additional keyword arguments are passed to `~astropy.table.QTable`.
    """
    _required_columns = ...
    def __init__(self, data=..., *, time_bin_start=..., time_bin_end=..., time_bin_size=..., n_bins=..., **kwargs) -> None:
        ...
    
    @property
    def time_bin_start(self): # -> QTable | _VT@OrderedDict | TableColumns | Row | BinnedTimeSeries:
        """
        The start times of all the time bins.
        """
        ...
    
    @property
    def time_bin_center(self):
        """
        The center times of all the time bins.
        """
        ...
    
    @property
    def time_bin_end(self):
        """
        The end times of all the time bins.
        """
        ...
    
    @property
    def time_bin_size(self): # -> QTable | _VT@OrderedDict | TableColumns | Row | BinnedTimeSeries:
        """
        The sizes of all the time bins.
        """
        ...
    
    def __getitem__(self, item): # -> QTable | _VT@OrderedDict | TableColumns | Row | BinnedTimeSeries:
        ...
    
    @classmethod
    def read(self, filename, time_bin_start_column=..., time_bin_end_column=..., time_bin_size_column=..., time_bin_size_unit=..., time_format=..., time_scale=..., format=..., *args, **kwargs): # -> BinnedTimeSeries:
        """
        Read and parse a file and returns a `astropy.timeseries.BinnedTimeSeries`.

        This method uses the unified I/O infrastructure in Astropy which makes
        it easy to define readers/writers for various classes
        (https://docs.astropy.org/en/stable/io/unified.html). By default, this
        method will try and use readers defined specifically for the
        `astropy.timeseries.BinnedTimeSeries` class - however, it is also
        possible to use the ``format`` keyword to specify formats defined for
        the `astropy.table.Table` class - in this case, you will need to also
        provide the column names for column containing the start times for the
        bins, as well as other column names (see the Parameters section below
        for details)::

            >>> from astropy.timeseries.binned import BinnedTimeSeries
            >>> ts = BinnedTimeSeries.read('binned.dat', format='ascii.ecsv',
            ...                            time_bin_start_column='date_start',
            ...                            time_bin_end_column='date_end')  # doctest: +SKIP

        Parameters
        ----------
        filename : str
            File to parse.
        format : str
            File format specifier.
        time_bin_start_column : str
            The name of the column with the start time for each bin.
        time_bin_end_column : str, optional
            The name of the column with the end time for each bin. Either this
            option or ``time_bin_size_column`` should be specified.
        time_bin_size_column : str, optional
            The name of the column with the size for each bin. Either this
            option or ``time_bin_end_column`` should be specified.
        time_bin_size_unit : `astropy.units.Unit`, optional
            If ``time_bin_size_column`` is specified but does not have a unit
            set in the table, you can specify the unit manually.
        time_format : str, optional
            The time format for the start and end columns.
        time_scale : str, optional
            The time scale for the start and end columns.
        *args : tuple, optional
            Positional arguments passed through to the data reader.
        **kwargs : dict, optional
            Keyword arguments passed through to the data reader.

        Returns
        -------
        out : `astropy.timeseries.binned.BinnedTimeSeries`
            BinnedTimeSeries corresponding to the file.

        """
        ...
    

