import argparse
from typing import Any, TypeGuard

import astropy.units as u

from .catalogs import SkyViewSurveys
from .instruments import FixedSizeInstrument, Hawki, Instrument
from .utils import StrEnum


class Parsers(StrEnum):
    get_brightest = "get_brightest"
    make_fc = "make_fc"


def get_defaults(use_parser: Parsers) -> argparse.ArgumentParser:
    """Get default arguments for a given parser"""

    if use_parser == Parsers.get_brightest:
        parser = argparse.ArgumentParser(description="Calculate Brightest Star & Optionally plot Finder Chart")
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
            help=(
                "Instrument FoV in arcmin (if using FixedSize instrument), finder chart will be roughly twice the size."
            ),
            type=float,
            default=7.5,
            action="store",
        )
        return parser
    elif use_parser == Parsers.make_fc:
        parser = argparse.ArgumentParser(description="Make Finder chart")
        group = parser.add_mutually_exclusive_group(required=False)
        group.add_argument(
            "image",
            help="OPTIONAL: image to use for plot if not querying",
            nargs="?",
            type=str,
            action="store",
            default=None,
        )
        group.add_argument(
            "-s",
            "--survey",
            help=(
                f"survey to query for image {[e.value for e in SkyViewSurveys]}, else pass path to fits image with WCS."
            ),
            action="store",
            type=str,
            default="DSS",
        )
        parser.add_argument(
            "-t", "--target", help="FLOWS Targetname or targetid", type=str, default="None", action="store"
        )

        parser.add_argument(
            "-i",
            "--instrument",
            help="instrument name",
            choices=["Hawki", "FixedSize"],
            type=str,
            default="FixedSize",
            action="store",
        )
        parser.add_argument(
            "--size",
            "--fov",
            help="Instrument FoV in arcmin (if using FixedSize instrument), finder chart will be slightly bigger.",
            type=float,
            default=7.5,
            action="store",
        )
        parser.add_argument(
            "--cmap", help="matplotlib colormap to use for image", type=str, default="gist_yarg", action="store"
        )

        group1 = parser.add_argument_group("Plot Scaling:")
        group1.add_argument("--scale", help="How to scale image", type=str, default="linear", action="store")
        group1.add_argument("--sigma", help="Sigma for Z-scaling", type=float, default=2.5, action="store")
        group1.add_argument("--contrast", help="Contrast [0,1] for Z-scaling", type=float, default=0.25, action="store")

        group2 = parser.add_argument_group("Rotation & Shifts:")
        group2.add_argument("-r", "--rotate", help="rotation angle in degrees", type=float, default=0.0, action="store")
        group2.add_argument("-a", "--shifta", help="shift alpha in arcsec", type=float, default=0.0, action="store")
        group2.add_argument("-d", "--shiftd", help="shift delta in arcsec", type=float, default=0.0, action="store")
        return parser

    raise ValueError(f"Parser {use_parser} not supported")


def get_instrument(args: argparse.Namespace) -> type[Instrument]:
    if args.instrument not in ["Hawki", "FixedSize"]:
        raise ValueError(f"Instrument {args.instrument} not supported, use Hawki or FixedSize")
    instrument: type[Instrument] = Hawki
    if args.instrument == "FixedSize":
        instrument = FixedSizeInstrument
        instrument.field_hw = args.size << u.arcmin
    return instrument


def check_flows_target(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if args.target == "None":
        parser.error("target id or name not provided, use -t <targetid> or <targetname>")
    elif args.target.isnumeric():
        args.target = int(args.target)


def parse_brightest() -> tuple[float, int | str, float, float, bool, type[Instrument]]:
    parser = get_defaults(Parsers.get_brightest)
    args = parser.parse_args()

    check_flows_target(args, parser)

    instrument = get_instrument(args)

    return args.target, args.rotate, args.shifta, args.shiftd, args.plot, instrument


def parse_fc() -> (
    tuple[float, int | str, float, float, type[Instrument], str | None, str | None, str, str, float, float]
):
    parser = get_defaults(Parsers.get_brightest)
    args = parser.parse_args()

    check_flows_target(args, parser)
    instrument = get_instrument(args)

    if args.survey not in [e.value for e in SkyViewSurveys]:
        raise ValueError(f"Survey {args.survey} not supported, use {[e.value for e in SkyViewSurveys]}")

    if args.image is not None:
        args.survey = None

    if not (0 < args.contrast <= 1):
        raise ValueError(f"Contrast must be between 0 and 1, got {args.contrast}")

    return (
        args.rotate,
        args.target,
        args.shifta,
        args.shiftd,
        instrument,
        args.image,
        args.survey,
        args.cmap,
        args.scale,
        args.sigma,
        args.contrast,
    )


def parse(use_parser: Parsers = Parsers.get_brightest) -> Any:
    """Parse command line input to get target, position angle (rotate), alpha and delta offsets (shifta, shiftd)"""
    match use_parser:
        case Parsers.get_brightest:
            return parse_brightest()
        case Parsers.make_fc:
            return parse_fc()
