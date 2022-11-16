import os
import argparse
from typing import Any, Union, TypeGuard, TypedDict
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import tendrils.utils

# custom types
numeric = Union[float, int, u.Quantity]
tabular = tuple[SkyCoord, Table]
flows_args = TypedDict('args', {'rotate': float, 'target': Union[str, int], 'shifta': float, 'shiftd': float, 'plot': bool})

class MissingCoords(SkyCoord):
    def __init__(self):
        super().__init__(u.Quantity(0, u.deg), u.Quantity(0, u.deg), frame='icrs')
MISSING_COORDS = MissingCoords()

def api_token() -> str:
    """Try to get the API token from environment variable or from tendrils config.ini"""
    token = str(os.environ.get('FLOWS_API_TOKEN'))
    if token is not None and token != 'None':
        tendrils.utils.set_api_token(token=token, overwrite=True)
        return token

    cfg = tendrils.utils.load_config()
    token = cfg.get('api', 'token', fallback='None')
    if token.lower() in ['', 'none', 'bad', 'bad_token', 'badtoken']:
        token = None
    else:
        tendrils.utils.set_api_token(token=token, overwrite=True)
    if isinstance(token, str):
        return token
    raise TypeError(f'API token is of type {type(token)} instead of string')


def is_quantity(qt: list[Any]) -> TypeGuard[list[u.Quantity]]:
    for q in qt:
        if not isinstance(q, u.Quantity):
            return False
    return True


if __name__ == '__main__':
    token = api_token()
    tendrils.utils.set_api_token('None')