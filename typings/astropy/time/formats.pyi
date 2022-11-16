"""
This type stub file was generated by pyright.
"""

import datetime
import astropy.units as u
from astropy.utils.decorators import lazyproperty

__all__ = ['TimeFormat', 'TimeJD', 'TimeMJD', 'TimeFromEpoch', 'TimeUnix', 'TimeUnixTai', 'TimeCxcSec', 'TimeGPS', 'TimeDecimalYear', 'TimePlotDate', 'TimeUnique', 'TimeDatetime', 'TimeString', 'TimeISO', 'TimeISOT', 'TimeFITS', 'TimeYearDayTime', 'TimeEpochDate', 'TimeBesselianEpoch', 'TimeJulianEpoch', 'TimeDeltaFormat', 'TimeDeltaSec', 'TimeDeltaJD', 'TimeEpochDateString', 'TimeBesselianEpochString', 'TimeJulianEpochString', 'TIME_FORMATS', 'TIME_DELTA_FORMATS', 'TimezoneInfo', 'TimeDeltaDatetime', 'TimeDatetime64', 'TimeYMDHMS', 'TimeNumeric', 'TimeDeltaNumeric']
__doctest_skip__ = ...
TIME_FORMATS = ...
TIME_DELTA_FORMATS = ...
FITS_DEPRECATED_SCALES = ...
class TimeFormat:
    """
    Base class for time representations.

    Parameters
    ----------
    val1 : numpy ndarray, list, number, str, or bytes
        Values to initialize the time or times.  Bytes are decoded as ascii.
    val2 : numpy ndarray, list, or number; optional
        Value(s) to initialize the time or times.  Only used for numerical
        input, to help preserve precision.
    scale : str
        Time scale of input value(s)
    precision : int
        Precision for seconds as floating point
    in_subfmt : str
        Select subformat for inputting string times
    out_subfmt : str
        Select subformat for outputting string times
    from_jd : bool
        If true then val1, val2 are jd1, jd2
    """
    _default_scale = ...
    subfmts = ...
    _registry = ...
    def __init__(self, val1, val2, scale, precision, in_subfmt, out_subfmt, from_jd=...) -> None:
        ...
    
    def __init_subclass__(cls, **kwargs): # -> None:
        ...
    
    @property
    def in_subfmt(self):
        ...
    
    @in_subfmt.setter
    def in_subfmt(self, subfmt): # -> None:
        ...
    
    @property
    def out_subfmt(self):
        ...
    
    @out_subfmt.setter
    def out_subfmt(self, subfmt): # -> None:
        ...
    
    @property
    def jd1(self): # -> NDArray[float_] | ndarray[Unknown, Unknown] | NDArray[Unknown] | None:
        ...
    
    @jd1.setter
    def jd1(self, jd1): # -> None:
        ...
    
    @property
    def jd2(self): # -> NDArray[float_] | ndarray[Unknown, Unknown] | NDArray[Unknown] | None:
        ...
    
    @jd2.setter
    def jd2(self, jd2): # -> None:
        ...
    
    def __len__(self): # -> int:
        ...
    
    @property
    def scale(self): # -> str:
        """Time scale"""
        ...
    
    @scale.setter
    def scale(self, val): # -> None:
        ...
    
    def mask_if_needed(self, value):
        ...
    
    @property
    def mask(self): # -> dict[Unknown, Unknown]:
        ...
    
    @property
    def masked(self): # -> bool | dict[Unknown, Unknown]:
        ...
    
    @property
    def jd2_filled(self): # -> NDArray[float_] | ndarray[Unknown, Unknown] | NDArray[Unknown] | None:
        ...
    
    @property
    def precision(self): # -> int:
        ...
    
    @precision.setter
    def precision(self, val): # -> None:
        ...
    
    @lazyproperty
    def cache(self): # -> defaultdict[Unknown, dict[Unknown, Unknown]]:
        """
        Return the cache associated with this instance.
        """
        ...
    
    def set_jds(self, val1, val2):
        """
        Set internal jd1 and jd2 from val1 and val2.  Must be provided
        by derived classes.
        """
        ...
    
    def to_value(self, parent=..., out_subfmt=...):
        """
        Return time representation from internal jd1 and jd2 in specified
        ``out_subfmt``.

        This is the base method that ignores ``parent`` and uses the ``value``
        property to compute the output. This is done by temporarily setting
        ``self.out_subfmt`` and calling ``self.value``. This is required for
        legacy Format subclasses prior to astropy 4.0  New code should instead
        implement the value functionality in ``to_value()`` and then make the
        ``value`` property be a simple call to ``self.to_value()``.

        Parameters
        ----------
        parent : object
            Parent `~astropy.time.Time` object associated with this
            `~astropy.time.TimeFormat` object
        out_subfmt : str or None
            Output subformt (use existing self.out_subfmt if `None`)

        Returns
        -------
        value : numpy.array, numpy.ma.array
            Array or masked array of formatted time representation values
        """
        ...
    
    @property
    def value(self):
        ...
    


class TimeNumeric(TimeFormat):
    subfmts = ...
    def to_value(self, jd1=..., jd2=..., parent=..., out_subfmt=...): # -> Any:
        """
        Return time representation from internal jd1 and jd2.
        Subclasses that require ``parent`` or to adjust the jds should
        override this method.
        """
        ...
    
    value = ...


class TimeJD(TimeNumeric):
    """
    Julian Date time format.
    This represents the number of days since the beginning of
    the Julian Period.
    For example, 2451544.5 in JD is midnight on January 1, 2000.
    """
    name = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    


class TimeMJD(TimeNumeric):
    """
    Modified Julian Date time format.
    This represents the number of days since midnight on November 17, 1858.
    For example, 51544.0 in MJD is midnight on January 1, 2000.
    """
    name = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    def to_value(self, **kwargs): # -> Any:
        ...
    
    value = ...


class TimeDecimalYear(TimeNumeric):
    """
    Time as a decimal year, with integer values corresponding to midnight
    of the first day of each year.  For example 2000.5 corresponds to the
    ISO time '2000-07-02 00:00:00'.
    """
    name = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    def to_value(self, **kwargs): # -> Any:
        ...
    
    value = ...


class TimeFromEpoch(TimeNumeric):
    """
    Base class for times that represent the interval from a particular
    epoch as a floating point multiple of a unit time interval (e.g. seconds
    or days).
    """
    @property
    def epoch(self): # -> () -> Time:
        """Reference epoch time from which the time interval is measured"""
        ...
    
    def set_jds(self, val1, val2): # -> None:
        """
        Initialize the internal jd1 and jd2 attributes given val1 and val2.
        For an TimeFromEpoch subclass like TimeUnix these will be floats giving
        the effective seconds since an epoch time (e.g. 1970-01-01 00:00:00).
        """
        ...
    
    def to_value(self, parent=..., **kwargs): # -> Any:
        ...
    
    value = ...


class TimeUnix(TimeFromEpoch):
    """
    Unix time (UTC): seconds from 1970-01-01 00:00:00 UTC, ignoring leap seconds.

    For example, 946684800.0 in Unix time is midnight on January 1, 2000.

    NOTE: this quantity is not exactly unix time and differs from the strict
    POSIX definition by up to 1 second on days with a leap second.  POSIX
    unix time actually jumps backward by 1 second at midnight on leap second
    days while this class value is monotonically increasing at 86400 seconds
    per UTC day.
    """
    name = ...
    unit = ...
    epoch_val = ...
    epoch_val2 = ...
    epoch_scale = ...
    epoch_format = ...


class TimeUnixTai(TimeUnix):
    """
    Unix time (TAI): SI seconds elapsed since 1970-01-01 00:00:00 TAI (see caveats).

    This will generally differ from standard (UTC) Unix time by the cumulative
    integral number of leap seconds introduced into UTC since 1972-01-01 UTC
    plus the initial offset of 10 seconds at that date.

    This convention matches the definition of linux CLOCK_TAI
    (https://www.cl.cam.ac.uk/~mgk25/posix-clocks.html),
    and the Precision Time Protocol
    (https://en.wikipedia.org/wiki/Precision_Time_Protocol), which
    is also used by the White Rabbit protocol in High Energy Physics:
    https://white-rabbit.web.cern.ch.

    Caveats:

    - Before 1972, fractional adjustments to UTC were made, so the difference
      between ``unix`` and ``unix_tai`` time is no longer an integer.
    - Because of the fractional adjustments, to be very precise, ``unix_tai``
      is the number of seconds since ``1970-01-01 00:00:00 TAI`` or equivalently
      ``1969-12-31 23:59:51.999918 UTC``.  The difference between TAI and UTC
      at that epoch was 8.000082 sec.
    - On the day of a positive leap second the difference between ``unix`` and
      ``unix_tai`` times increases linearly through the day by 1.0. See also the
      documentation for the `~astropy.time.TimeUnix` class.
    - Negative leap seconds are possible, though none have been needed to date.

    Examples
    --------

      >>> # get the current offset between TAI and UTC
      >>> from astropy.time import Time
      >>> t = Time('2020-01-01', scale='utc')
      >>> t.unix_tai - t.unix
      37.0

      >>> # Before 1972, the offset between TAI and UTC was not integer
      >>> t = Time('1970-01-01', scale='utc')
      >>> t.unix_tai - t.unix  # doctest: +FLOAT_CMP
      8.000082

      >>> # Initial offset of 10 seconds in 1972
      >>> t = Time('1972-01-01', scale='utc')
      >>> t.unix_tai - t.unix
      10.0
    """
    name = ...
    epoch_val = ...
    epoch_scale = ...


class TimeCxcSec(TimeFromEpoch):
    """
    Chandra X-ray Center seconds from 1998-01-01 00:00:00 TT.
    For example, 63072064.184 is midnight on January 1, 2000.
    """
    name = ...
    unit = ...
    epoch_val = ...
    epoch_val2 = ...
    epoch_scale = ...
    epoch_format = ...


class TimeGPS(TimeFromEpoch):
    """GPS time: seconds from 1980-01-06 00:00:00 UTC
    For example, 630720013.0 is midnight on January 1, 2000.

    Notes
    =====
    This implementation is strictly a representation of the number of seconds
    (including leap seconds) since midnight UTC on 1980-01-06.  GPS can also be
    considered as a time scale which is ahead of TAI by a fixed offset
    (to within about 100 nanoseconds).

    For details, see https://www.usno.navy.mil/USNO/time/gps/usno-gps-time-transfer
    """
    name = ...
    unit = ...
    epoch_val = ...
    epoch_val2 = ...
    epoch_scale = ...
    epoch_format = ...


class TimePlotDate(TimeFromEpoch):
    """
    Matplotlib `~matplotlib.pyplot.plot_date` input:
    1 + number of days from 0001-01-01 00:00:00 UTC

    This can be used directly in the matplotlib `~matplotlib.pyplot.plot_date`
    function::

      >>> import matplotlib.pyplot as plt
      >>> jyear = np.linspace(2000, 2001, 20)
      >>> t = Time(jyear, format='jyear', scale='utc')
      >>> plt.plot_date(t.plot_date, jyear)
      >>> plt.gcf().autofmt_xdate()  # orient date labels at a slant
      >>> plt.draw()

    For example, 730120.0003703703 is midnight on January 1, 2000.
    """
    name = ...
    unit = ...
    epoch_val = ...
    epoch_val2 = ...
    epoch_scale = ...
    epoch_format = ...
    @lazyproperty
    def epoch(self): # -> Time | (() -> Time):
        """Reference epoch time from which the time interval is measured"""
        ...
    


class TimeStardate(TimeFromEpoch):
    """
    Stardate: date units from 2318-07-05 12:00:00 UTC.
    For example, stardate 41153.7 is 00:52 on April 30, 2363.
    See http://trekguide.com/Stardates.htm#TNG for calculations and reference points
    """
    name = ...
    unit = ...
    epoch_val = ...
    epoch_val2 = ...
    epoch_scale = ...
    epoch_format = ...


class TimeUnique(TimeFormat):
    """
    Base class for time formats that can uniquely create a time object
    without requiring an explicit format specifier.  This class does
    nothing but provide inheritance to identify a class as unique.
    """
    ...


class TimeAstropyTime(TimeUnique):
    """
    Instantiate date from an Astropy Time object (or list thereof).

    This is purely for instantiating from a Time object.  The output
    format is the same as the first time instance.
    """
    name = ...
    def __new__(cls, val1, val2, scale, precision, in_subfmt, out_subfmt, from_jd=...): # -> TimeUnique:
        """
        Use __new__ instead of __init__ to output a class instance that
        is the same as the class of the first Time object in the list.
        """
        ...
    


class TimeDatetime(TimeUnique):
    """
    Represent date as Python standard library `~datetime.datetime` object

    Example::

      >>> from astropy.time import Time
      >>> from datetime import datetime
      >>> t = Time(datetime(2000, 1, 2, 12, 0, 0), scale='utc')
      >>> t.iso
      '2000-01-02 12:00:00.000'
      >>> t.tt.datetime
      datetime.datetime(2000, 1, 2, 12, 1, 4, 184000)
    """
    name = ...
    def set_jds(self, val1, val2): # -> None:
        """Convert datetime object contained in val1 to jd1, jd2"""
        ...
    
    def to_value(self, timezone=..., parent=..., out_subfmt=...): # -> NDArray[Any]:
        """
        Convert to (potentially timezone-aware) `~datetime.datetime` object.

        If ``timezone`` is not ``None``, return a timezone-aware datetime
        object.

        Parameters
        ----------
        timezone : {`~datetime.tzinfo`, None}, optional
            If not `None`, return timezone-aware datetime.

        Returns
        -------
        `~datetime.datetime`
            If ``timezone`` is not ``None``, output will be timezone-aware.
        """
        ...
    
    value = ...


class TimeYMDHMS(TimeUnique):
    """
    ymdhms: A Time format to represent Time as year, month, day, hour,
    minute, second (thus the name ymdhms).

    Acceptable inputs must have keys or column names in the "YMDHMS" set of
    ``year``, ``month``, ``day`` ``hour``, ``minute``, ``second``:

    - Dict with keys in the YMDHMS set
    - NumPy structured array, record array or astropy Table, or single row
      of those types, with column names in the YMDHMS set

    One can supply a subset of the YMDHMS values, for instance only 'year',
    'month', and 'day'.  Inputs have the following defaults::

      'month': 1, 'day': 1, 'hour': 0, 'minute': 0, 'second': 0

    When the input is supplied as a ``dict`` then each value can be either a
    scalar value or an array.  The values will be broadcast to a common shape.

    Example::

      >>> from astropy.time import Time
      >>> t = Time({'year': 2015, 'month': 2, 'day': 3,
      ...           'hour': 12, 'minute': 13, 'second': 14.567},
      ...           scale='utc')
      >>> t.iso
      '2015-02-03 12:13:14.567'
      >>> t.ymdhms.year
      2015
    """
    name = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    @property
    def value(self): # -> recarray[Unknown, Unknown]:
        ...
    


class TimezoneInfo(datetime.tzinfo):
    """
    Subclass of the `~datetime.tzinfo` object, used in the
    to_datetime method to specify timezones.

    It may be safer in most cases to use a timezone database package like
    pytz rather than defining your own timezones - this class is mainly
    a workaround for users without pytz.
    """
    @u.quantity_input(utc_offset=u.day, dst=u.day)
    def __init__(self, utc_offset=..., dst=..., tzname=...) -> None:
        """
        Parameters
        ----------
        utc_offset : `~astropy.units.Quantity`, optional
            Offset from UTC in days. Defaults to zero.
        dst : `~astropy.units.Quantity`, optional
            Daylight Savings Time offset in days. Defaults to zero
            (no daylight savings).
        tzname : str or None, optional
            Name of timezone

        Examples
        --------
        >>> from datetime import datetime
        >>> from astropy.time import TimezoneInfo  # Specifies a timezone
        >>> import astropy.units as u
        >>> utc = TimezoneInfo()    # Defaults to UTC
        >>> utc_plus_one_hour = TimezoneInfo(utc_offset=1*u.hour)  # UTC+1
        >>> dt_aware = datetime(2000, 1, 1, 0, 0, 0, tzinfo=utc_plus_one_hour)
        >>> print(dt_aware)
        2000-01-01 00:00:00+01:00
        >>> print(dt_aware.astimezone(utc))
        1999-12-31 23:00:00+00:00
        """
        ...
    
    def utcoffset(self, dt): # -> timedelta:
        ...
    
    def tzname(self, dt): # -> str:
        ...
    
    def dst(self, dt): # -> timedelta:
        ...
    


class TimeString(TimeUnique):
    """
    Base class for string-like time representations.

    This class assumes that anything following the last decimal point to the
    right is a fraction of a second.

    **Fast C-based parser**

    Time format classes can take advantage of a fast C-based parser if the times
    are represented as fixed-format strings with year, month, day-of-month,
    hour, minute, second, OR year, day-of-year, hour, minute, second. This can
    be a factor of 20 or more faster than the pure Python parser.

    Fixed format means that the components always have the same number of
    characters. The Python parser will accept ``2001-9-2`` as a date, but the C
    parser would require ``2001-09-02``.

    A subclass in this case must define a class attribute ``fast_parser_pars``
    which is a `dict` with all of the keys below. An inherited attribute is not
    checked, only an attribute in the class ``__dict__``.

    - ``delims`` (tuple of int): ASCII code for character at corresponding
      ``starts`` position (0 => no character)

    - ``starts`` (tuple of int): position where component starts (including
      delimiter if present). Use -1 for the month component for format that use
      day of year.

    - ``stops`` (tuple of int): position where component ends. Use -1 to
      continue to end of string, or for the month component for formats that use
      day of year.

    - ``break_allowed`` (tuple of int): if true (1) then the time string can
          legally end just before the corresponding component (e.g. "2000-01-01"
          is a valid time but "2000-01-01 12" is not).

    - ``has_day_of_year`` (int): 0 if dates have year, month, day; 1 if year,
      day-of-year
    """
    def __init_subclass__(cls, **kwargs): # -> None:
        ...
    
    def parse_string(self, timestr, subfmts): # -> list[Any] | list[int]:
        """Read time from a single string, using a set of possible formats."""
        ...
    
    def set_jds(self, val1, val2): # -> None:
        """Parse the time strings contained in val1 and set jd1, jd2"""
        ...
    
    def get_jds_python(self, val1, val2): # -> tuple[NDArray[floating[Any]], Unknown]:
        """Parse the time strings contained in val1 and get jd1, jd2"""
        ...
    
    def get_jds_fast(self, val1, val2): # -> tuple[NDArray[floating[Any]], Unknown]:
        """Use fast C parser to parse time strings in val1 and get jd1, jd2"""
        ...
    
    def str_kwargs(self): # -> Generator[dict[str, int | None], None, None]:
        """
        Generator that yields a dict of values corresponding to the
        calendar date and time for the internal JD values.
        """
        ...
    
    def format_string(self, str_fmt, **kwargs):
        """Write time to a string using a given format.

        By default, just interprets str_fmt as a format string,
        but subclasses can add to this.
        """
        ...
    
    @property
    def value(self): # -> ndarray[Any, dtype[Unknown]]:
        ...
    


class TimeISO(TimeString):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DD HH:MM:SS.sss...".
    For example, 2000-01-01 00:00:00.000 is midnight on January 1, 2000.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """
    name = ...
    subfmts = ...
    fast_parser_pars = ...
    def parse_string(self, timestr, subfmts): # -> list[Any] | list[int]:
        ...
    


class TimeISOT(TimeISO):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DDTHH:MM:SS.sss...".
    This is the same as TimeISO except for a "T" instead of space between
    the date and time.
    For example, 2000-01-01T00:00:00.000 is midnight on January 1, 2000.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """
    name = ...
    subfmts = ...
    fast_parser_pars = ...


class TimeYearDayTime(TimeISO):
    """
    Year, day-of-year and time as "YYYY:DOY:HH:MM:SS.sss...".
    The day-of-year (DOY) goes from 001 to 365 (366 in leap years).
    For example, 2000:001:00:00:00.000 is midnight on January 1, 2000.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """
    name = ...
    subfmts = ...
    fast_parser_pars = ...


class TimeDatetime64(TimeISOT):
    name = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    @property
    def value(self): # -> NDArray[Any]:
        ...
    


class TimeFITS(TimeString):
    """
    FITS format: "[±Y]YYYY-MM-DD[THH:MM:SS[.sss]]".

    ISOT but can give signed five-digit year (mostly for negative years);

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date': date
    - 'longdate_hms': as 'date_hms', but with signed 5-digit year
    - 'longdate': as 'date', but with signed 5-digit year

    See Rots et al., 2015, A&A 574:A36 (arXiv:1409.7583).
    """
    name = ...
    subfmts = ...
    subfmts = ...
    def parse_string(self, timestr, subfmts): # -> list[int | float]:
        """Read time and deprecated scale if present"""
        ...
    
    @property
    def value(self): # -> ndarray[Any, dtype[Unknown]]:
        """Convert times to strings, using signed 5 digit if necessary."""
        ...
    


class TimeEpochDate(TimeNumeric):
    """
    Base class for support floating point Besselian and Julian epoch dates
    """
    _default_scale = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    def to_value(self, **kwargs): # -> Any:
        ...
    
    value = ...


class TimeBesselianEpoch(TimeEpochDate):
    """Besselian Epoch year as floating point value(s) like 1950.0"""
    name = ...
    epoch_to_jd = ...
    jd_to_epoch = ...


class TimeJulianEpoch(TimeEpochDate):
    """Julian Epoch year as floating point value(s) like 2000.0"""
    name = ...
    unit = ...
    epoch_to_jd = ...
    jd_to_epoch = ...


class TimeEpochDateString(TimeString):
    """
    Base class to support string Besselian and Julian epoch dates
    such as 'B1950.0' or 'J2000.0' respectively.
    """
    _default_scale = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    @property
    def value(self): # -> ndarray[Any, dtype[Unknown]]:
        ...
    


class TimeBesselianEpochString(TimeEpochDateString):
    """Besselian Epoch year as string value(s) like 'B1950.0'"""
    name = ...
    epoch_to_jd = ...
    jd_to_epoch = ...
    epoch_prefix = ...


class TimeJulianEpochString(TimeEpochDateString):
    """Julian Epoch year as string value(s) like 'J2000.0'"""
    name = ...
    epoch_to_jd = ...
    jd_to_epoch = ...
    epoch_prefix = ...


class TimeDeltaFormat(TimeFormat):
    """Base class for time delta representations"""
    _registry = ...


class TimeDeltaNumeric(TimeDeltaFormat, TimeNumeric):
    def set_jds(self, val1, val2): # -> None:
        ...
    
    def to_value(self, **kwargs): # -> Any:
        ...
    
    value = ...


class TimeDeltaSec(TimeDeltaNumeric):
    """Time delta in SI seconds"""
    name = ...
    unit = ...


class TimeDeltaJD(TimeDeltaNumeric):
    """Time delta in Julian days (86400 SI seconds)"""
    name = ...
    unit = ...


class TimeDeltaDatetime(TimeDeltaFormat, TimeUnique):
    """Time delta in datetime.timedelta"""
    name = ...
    def set_jds(self, val1, val2): # -> None:
        ...
    
    @property
    def value(self): # -> NDArray[Any]:
        ...
    

