from __future__ import annotations
import warnings
from erfa import ErfaWarning
from .plots import Plotter
from .auth import test_connection
from .observer import get_flows_observer
from .utils import parse

# Most useless warnings ever spammed for every operation by this package!
warnings.filterwarnings('ignore', category=ErfaWarning, append=True)


def main():
    # Parse input
    rot, tid, shifta, shiftd, make_fc = parse()

    # Test connection to flows:
    test_connection()

    # Print brightest star in field
    obs = get_flows_observer(rot, tid, shifta, shiftd)
    obs.check_bright_stars(region=obs.regions[0])

    # Make finding chart if requested
    if make_fc:
        Plotter(obs).make_finding_chart(radius=14)


if __name__ == '__main__':
    main()
