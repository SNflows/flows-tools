"""
This type stub file was generated by pyright.
"""

from astropy.table import Table
from typing import Union

"""
Fetch current lightcurve from Flows API.
"""

def get_lightcurve(target: Union[int, str]) -> Table:
    """
    Retrieve lightcurve from Flows server.

    Parameters:
        target (int): Target to download lightcurve for.

    Returns:
        :class:`astropy.table.Table`: Table containing lightcurve.

    TODO:
        - Enable caching of files.
    """
    ...
