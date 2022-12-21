from __future__ import annotations

import warnings
from typing import cast

from erfa import ErfaWarning

from .argparser import parse_fc
from .auth import test_connection
from .observer import get_flows_observer
from .plan import make_plan
from .plots import Plotter, get_zscaler

# Most useless warnings ever spammed for every operation by this package!
warnings.filterwarnings("ignore", category=ErfaWarning, append=True)
warnings.filterwarnings("ignore", message="invalid value", category=RuntimeWarning, append=True)


def main():
    # Parse input
    rot, tid, shifta, shiftd, inst, image, survey, cmap, scale, sigma, contrast = parse_fc()

    # Test connection to flows:
    test_connection()

    # Whether to query for image or use local image
    plan = make_plan(rot, shifta, shiftd, image=image, survey=survey, scale=scale)

    # Create observer
    obs = get_flows_observer(plan, tid, inst)

    # Make finding chart if requested
    radius = cast(float, inst.field_hw.value) * 1.2

    zscaler = get_zscaler(krej=sigma, contrast=contrast)
    Plotter(obs, cmap=cmap).make_finding_chart(radius=radius, zscaler=zscaler)


if __name__ == "__main__":
    main()
