from __future__ import annotations

import warnings
from typing import cast

from erfa import ErfaWarning

from .argparser import parse_brightest
from .auth import test_connection
from .observer import get_flows_observer
from .plan import make_plan
from .plots import Plotter

# Most useless warnings ever spammed for every operation by this package!
warnings.filterwarnings("ignore", category=ErfaWarning, append=True)
warnings.filterwarnings("ignore", message="invalid value", category=RuntimeWarning, append=True)


def main():
    # Parse input
    rot, tid, shifta, shiftd, make_fc, inst = parse_brightest()

    # Test connection to flows:
    test_connection()

    # Print brightest star in (first) field
    plan = make_plan(rot, shifta, shiftd)
    obs = get_flows_observer(plan, tid, inst)
    obs.check_bright_stars(region=obs.regions[0])

    # Make finding chart if requested
    radius = cast(float, inst.field_hw.value) * 2.0
    if make_fc:
        Plotter(obs).make_finding_chart(radius=radius)


if __name__ == "__main__":
    main()
