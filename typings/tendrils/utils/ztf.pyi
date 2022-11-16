"""
This type stub file was generated by pyright.
"""

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from typing import Optional

"""
Query ZTF target information using ALeRCE API.
https://alerceapi.readthedocs.io/
@TODO: REFACTOR TO FIT API! Make Class.
"""

def query_ztf_id(coo_centre: SkyCoord, radius: u.Quantity = ..., discovery_date: Optional[Time] = ...) -> Optional[str]:
    """
    Query ALeRCE ZTF api to lookup ZTF identifier.

    In case multiple identifiers are found within the search cone, the one
    closest to the centre is returned.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (Angle, optional): Search radius. Default 3 arcsec.
        discovery_date (:class:`astropy.time.Time`, optional): Discovery date of target to
            match against ZTF. The date is compared to the ZTF first timestamp and ZTF targets
            are rejected if they are not within 15 days prior to the discovery date
            and 90 days after.

    Returns:
        str: ZTF identifier.
    """
    ...

def download_ztf_photometry(targetid: int) -> Table:
    """
    Download ZTF photometry from ALERCE API.

    Parameters:
        targetid (int): Target identifier.

    Returns:
        :class:`astropy.table.Table`: ZTF photometry table.
    """
    ...