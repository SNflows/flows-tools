"""
This type stub file was generated by pyright.
"""

from datetime import datetime
from functools import lru_cache
from typing import Optional, Union, Any
from astropy.time import Time
from astropy.coordinates import SkyCoord

"""
Get information about targets in Flows.
Add target to Flows.
"""

def get_target(target: Union[int, str]) -> dict[str, Any]:
    """
    Get target as json
    Args:
        target: Optional[int, str] targetid

    Returns: json of target info

    """
    ...

def get_targets() -> dict[str, Any]:
    """
    Get json list of all targets and some basic info about them
    Returns: json list of all targets

    """
    ...

def add_target(
    name: str,
    coord: SkyCoord,
    redshift: Optional[float] = ...,
    redshift_error: Optional[float] = ...,
    discovery_date: Union[Optional[Time], Optional[datetime], Optional[str]] = ...,
    discovery_mag: Union[Optional[float], Optional[int]] = ...,
    host_galaxy: Optional[str] = ...,
    ztf: Optional[str] = ...,
    sntype: Optional[str] = ...,
    status: str = ...,
    project: str = ...,
) -> int:
    """
    Add new candidate or target.

    Coordinates are specified using an Astropy SkyCoord object, which can be
    created in the following way:

    coord = SkyCoord(ra=19.1, dec=89.00001, unit='deg', frame='icrs')

    The easiest way is to specify ``discovery_date`` as an Astropy Time object:

    discovery_date = Time('2020-01-02 00:00:00', format='iso', scale='utc')

    Alternatively, you can also specify it as a :class:`datetime.datetime` object,
    but some care has to be taken with specifying the correct timezone:

    discovery_date = datetime.strptime('2020-01-02 00:00:00', '%Y-%m-%d %H:%M:%S%f')
    discovery_date = pytz.timezone('America/New_York').localize(discovery_date)

    Lastly, it can be given as a simple date-string of the following form,
    but here the data has to be given in UTC:

    discovery_date = '2020-01-02 23:56:02.123'

    Parameters:
        name (str): Name of target. Must be of the form "YYYYxyz", where YYYY is the year,
            and xyz are letters.
        coord (:class:ʼastropy.coordinates.SkyCoordʼ): Sky coordinates of target.
        redshift (float, optional): Redshift.
        redshift_error (float, optional): Uncertainty on redshift.
        discovery_date (:class:`astropy.time.Time`, :class:`datetime.datetime` or str, optional):
        discovery_mag (float, int, optional): Magnitude at time of discovery.
        host_galaxy (str, optional): Host galaxy name.
        sntype (str, optional): Supernovae type (e.g. Ia, Ib, II).
        ztf (str, optional): ZTF identifier.
        status (str, optional):
        project (str, optional):

    Returns:
        int: New target identifier in Flows system.
    """
    ...