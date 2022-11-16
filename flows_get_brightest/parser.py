import argparse
import astropy.units as u
from .instruments import Hawki, FixedSizeInstrument, Instrument


def parse() -> tuple[float, int | str, float, float, bool, type[Instrument]]:
    """Parse command line input to get target, position angle (rotate), alpha and delta offsets (shifta, shiftd)"""
    parser = argparse.ArgumentParser(description="Calculate Brightest Star")
    parser.add_argument(
        "-t", "--target", help="calculate for this targetname or targetid", type=str, default="None", action="store"
    )
    parser.add_argument("-r", "--rotate", help="rotation angle in degrees", type=float, default=0.0, action="store")
    parser.add_argument("-a", "--shifta", help="shift alpha in arcsec", type=float, default=0.0, action="store")
    parser.add_argument("-d", "--shiftd", help="shift delta in arcsec", type=float, default=0.0, action="store")
    parser.add_argument("-p", "--plot", help="whether to query images and plot", action="store_true")
    parser.add_argument(
        "-i",
        "--instrument",
        help="instrument name",
        choices=["Hawki", "FixedSize"],
        type=str,
        default="Hawki",
        action="store",
    )
    parser.add_argument(
        "--size",
        help="Instrument FoV in arcmin (if using FixedSize instrument), finder chart will be roughly twice the size.",
        type=float,
        default=7.5,
        action="store",
    )

    args = parser.parse_args()
    if args.target == "None":
        parser.error("target id or name not provided, use -t <targetid> or <targetname>")
    elif args.target.isnumeric():
        args.target = int(args.target)

    if args.instrument not in ["Hawki", "FixedSize"]:
        raise ValueError(f"Instrument {args.instrument} not supported, use Hawki or FixedSize")
    instrument: type[Instrument] = Hawki
    if args.instrument == "FixedSize":
        instrument = FixedSizeInstrument
        instrument.field_hw = args.size << u.arcmin

    return args.rotate, args.target, args.shifta, args.shiftd, args.plot, instrument
