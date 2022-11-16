"""
This type stub file was generated by pyright.
"""

from astropy.utils.state import ScienceState

"""
This module contains convenience functions for getting a coordinate object
for a named object by querying SESAME and getting the first returned result.
Note that this is intended to be a convenience, and is very simple. If you
need precise coordinates for an object you should find the appropriate
reference for that measurement and input the coordinates manually.
"""
__all__ = ["get_icrs_coordinates"]
class sesame_url(ScienceState):
    """
    The URL(s) to Sesame's web-queryable database.
    """
    _value = ...
    @classmethod
    def validate(cls, value):
        ...
    


class sesame_database(ScienceState):
    """
    This specifies the default database that SESAME will query when
    using the name resolve mechanism in the coordinates
    subpackage. Default is to search all databases, but this can be
    'all', 'simbad', 'ned', or 'vizier'.
    """
    _value = ...
    @classmethod
    def validate(cls, value):
        ...
    


class NameResolveError(Exception):
    ...


def get_icrs_coordinates(name, parse=..., cache=...): # -> SkyCoord:
    """
    Retrieve an ICRS object by using an online name resolving service to
    retrieve coordinates for the specified name. By default, this will
    search all available databases until a match is found. If you would like
    to specify the database, use the science state
    ``astropy.coordinates.name_resolve.sesame_database``. You can also
    specify a list of servers to use for querying Sesame using the science
    state ``astropy.coordinates.name_resolve.sesame_url``. This will try
    each one in order until a valid response is returned. By default, this
    list includes the main Sesame host and a mirror at vizier.  The
    configuration item `astropy.utils.data.Conf.remote_timeout` controls the
    number of seconds to wait for a response from the server before giving
    up.

    Parameters
    ----------
    name : str
        The name of the object to get coordinates for, e.g. ``'M42'``.
    parse : bool
        Whether to attempt extracting the coordinates from the name by
        parsing with a regex. For objects catalog names that have
        J-coordinates embedded in their names eg:
        'CRTS SSS100805 J194428-420209', this may be much faster than a
        sesame query for the same object name. The coordinates extracted
        in this way may differ from the database coordinates by a few
        deci-arcseconds, so only use this option if you do not need
        sub-arcsecond accuracy for coordinates.
    cache : bool, str, optional
        Determines whether to cache the results or not. Passed through to
        `~astropy.utils.data.download_file`, so pass "update" to update the
        cached value.

    Returns
    -------
    coord : `astropy.coordinates.ICRS` object
        The object's coordinates in the ICRS frame.

    """
    ...

