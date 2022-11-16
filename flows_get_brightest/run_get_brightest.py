from __future__ import annotations
from typing import cast
import warnings
from erfa import ErfaWarning
from .plots import Plotter
from .auth import test_connection
from .observer import get_flows_observer
from .parser import parse


# Most useless warnings ever spammed for every operation by this package!
warnings.filterwarnings("ignore", category=ErfaWarning, append=True)
warnings.filterwarnings("ignore", message="invalid value", category=RuntimeWarning, append=True)


def main():
    # Parse input
    rot, tid, shifta, shiftd, make_fc, inst = parse()

    # Test connection to flows:
    test_connection()

    # Print brightest star in (first) field
    obs = get_flows_observer(rot, tid, shifta, shiftd, inst)
    obs.check_bright_stars(region=obs.regions[0])

    # Make finding chart if requested
    radius = cast(float, inst.field_hw.value) * 2
    if make_fc:
        Plotter(obs).make_finding_chart(radius=radius)


if __name__ == "__main__":
    main()
