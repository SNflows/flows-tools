import os
import argparse
from typing import Union
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import tendrils

# custom types
numeric = Union[float, int, u.Quantity]
tabular = tuple[SkyCoord, Table]


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
    return token


if __name__ == '__main__':
    token = api_token()
    tendrils.utils.set_api_token('None')
