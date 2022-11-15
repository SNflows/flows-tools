from flows_get_brightest.utils import parse
from flows_get_brightest.plots import Plotter
from flows_get_brightest.observer import get_flows_observer, Observer
import pytest

ARGS1 = [0, 0, 0, 0, False]

def test_get_brightest():
    rot, tid, shifta, shiftd, make_fc = parse()

    # Print brightest star in field
    obs = get_flows_observer(rot, tid, shifta, shiftd)
    obs.check_bright_stars(region=obs.regions[0])

    # Make finding chart if requested
    if make_fc:
        Plotter(obs).make_finding_chart(radius=14)