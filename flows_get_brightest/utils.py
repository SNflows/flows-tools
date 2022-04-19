import argparse
from typing import Union
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

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

