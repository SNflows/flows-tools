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


def parse():
    """Parse command line input to get target, position angle (rotate), alpha and delta offsets (shifta, shiftd)
    """
    parser = argparse.ArgumentParser(description='Calculate Brightest Star')
    parser.add_argument('-t', '--target', help="calculate for this targetname or targetid", type=str, default='None',
                        action='store')
    parser.add_argument('-r', '--rotate', help='rotation angle in degrees', type=float, default=0.0, action='store')
    parser.add_argument('-a', '--shifta', help='shift alpha in arcsec', type=float, default=0.0, action='store')
    parser.add_argument('-d', '--shiftd', help='shift delta in arcsec', type=float, default=0.0, action='store')
    parser.add_argument('-p', '--plot', help='whether to query images and plot', action='store_true')

    args = parser.parse_args()
    if args.target == 'None':
        parser.error('target id or name not provided, use -t <targetid> or <targetname>')
    elif args.target.isnumeric():
        args.target = int(args.target)
    return args.rotate, args.target, args.shifta, args.shiftd, args.plot


def api_token():
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
