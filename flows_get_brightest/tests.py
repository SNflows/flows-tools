import sys
from typing import Optional
from unittest import mock

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from .argparser import parse_brightest, parse_fc
from .auth import test_connection
from .catalogs import SkyViewSurveys
from .instruments import FixedSizeInstrument, Hawki, Instrument
from .observer import Observer, get_flows_observer
from .plan import Plan, make_plan
from .plots import Plotter
from .utils import api_token

rot, tid, shifta, shiftd = 30, 8, 10, 10
token = api_token()


@pytest.mark.skip(reason="Skip until destructive overwrite is fixed")
def test_auth(monkeypatch):
    """
    Test the auth module.
    """
    with pytest.warns(RuntimeWarning):
        monkeypatch.setattr("builtins.input", lambda _: "bad_token")
        test_connection()


@pytest.mark.parametrize("img", [None, "test.fits"])
@pytest.mark.parametrize("scale", ["linear", "log"])
@pytest.mark.parametrize("survey", [survey.value for survey in SkyViewSurveys])
def test_make_plan(img: Optional[str], scale: str, survey: str):
    """
    Test the plan module.
    """
    plan: Plan = make_plan(rot, shifta, shiftd, image=img)
    assert plan.rotation.value == rot
    assert plan.alpha.value == shifta
    assert plan.delta.value == shiftd
    assert plan.local_image == img
    assert plan.survey in [s for s in SkyViewSurveys]  # refactor to test output call of observer.
    assert plan.image_scale in ["linear", "log"]  # refactor to use enum and more meaningfull img test.


@pytest.fixture
def fixedinstrument_obs() -> Observer:
    plan = make_plan(rot, shifta, shiftd)
    plan.plan_finder_chart(image=None, survey=SkyViewSurveys.DSS.value, scale="linear")
    return get_flows_observer(plan, tid, instrument=FixedSizeInstrument)


@pytest.fixture
def Hawki_obs() -> Observer:
    plan = make_plan(rot, shifta, shiftd)
    return get_flows_observer(plan, tid, instrument=Hawki)


@pytest.fixture
def observer(request) -> Observer:
    return request.getfixturevalue(request.param)


@pytest.mark.parametrize("observer", ["fixedinstrument_obs", "Hawki_obs"], indirect=True)
def test_get_brightest(capsys, observer):
    observer.check_bright_stars(region=observer.regions[0])
    captured = capsys.readouterr()
    assert "Brightest star has" in captured.out
    sys.stdout.write(captured.out)
    sys.stderr.write(captured.err)


@pytest.mark.parametrize("observer", ["fixedinstrument_obs", "Hawki_obs"], indirect=True)
def test_plan(observer):
    assert observer.plan.rotation == rot * u.deg  # type: ignore
    assert observer.plan.alpha == shifta * u.arcsec  # type: ignore
    assert observer.plan.delta == shiftd * u.arcsec  # type: ignore
    assert observer.plan.shift is True
    assert observer.plan.rotate is True


@pytest.mark.parametrize("observer", ["fixedinstrument_obs", "Hawki_obs"], indirect=True)
def test_observer(observer):
    isinstance(observer, Observer)


@pytest.mark.slow
@pytest.mark.parametrize("observer", ["fixedinstrument_obs", "Hawki_obs"], indirect=True)
def test_make_finding_chart(observer, monkeypatch):
    import matplotlib.pyplot as plt

    monkeypatch.setattr(plt, "show", lambda: None)
    plotter = Plotter(observer)
    ax = plotter.make_finding_chart(observer, savefig=False)
    assert len(ax.get_images()) > 0
    title = ax.get_title()
    assert title.startswith(f"{observer.target.info['target_name']}")
    assert title.endswith("FC")


# End to end test with Hawki.
ARGS0 = (0, 8, 0, 0, False, 12.2)
ARGS1 = (30, 8, 0, 0, False, 11.5)
ARGS2 = (30, 8, -50, 100, False, 12.3)
argnames = ["-r", "-t", "-a", "-d", "-p"]


@pytest.mark.parametrize("args", [ARGS0, ARGS1, ARGS2])
def test_parse_brightest(args) -> None:
    cmd_args: list[str] = []
    for i, arg in enumerate(args[:-1]):
        if isinstance(arg, bool):
            if arg is False:
                continue
        cmd_args.append(f"{argnames[i]} {arg}")

    with mock.patch("sys.argv", ["flows_get_brightest"] + list(cmd_args)):
        rot, tid, shifta, shiftd, make_fc, inst = parse_brightest()
        assert rot == args[0]
        assert tid == args[1]
        assert shifta == args[2]
        assert shiftd == args[3]
        assert make_fc == args[4]
        assert inst == Hawki


@pytest.mark.slow
@pytest.mark.parametrize("args", [ARGS0, ARGS1, ARGS2])
def test_end_to_end(args: tuple[int, int | str, int, int, bool, float]) -> None:
    rot, tid, shifta, shiftd, make_fc, brightest = args

    # Print brightest star in field
    plan = make_plan(rot, shifta, shiftd)
    obs = get_flows_observer(plan, tid)
    stars = obs.check_bright_stars(region=obs.regions[0])

    assert pytest.approx(np.round(stars.min(), 1)) == brightest

    # Make finding chart if requested
    if make_fc:
        plotter = Plotter(obs)
        assert plotter.obs == obs


# End to end test with Hawki.
NARGS0 = (0, 8, 0, 0, "DSS", "Hawki", 7.5, "viridis", "linear")
NARGS1 = (30, 8, 0, 0, "DSS", "FixedSize", 8.5, "gist_yarg", "linear")
NARGS2 = (30, 8, -50, 100, "2MASS-H", "FixedSize", 9.5, "gist_yarg", "asinh")
nargnames = ["-r=", "-t=", "-a=", "-d=", "-s=", "-i=", "--size=", "--cmap=", "--scale="]


@pytest.mark.parametrize("args", [NARGS0, NARGS1, NARGS2])
def test_parse_fc(args):
    cmd_args: list[str] = []
    for i, arg in enumerate(args):
        cmd_args.append(f"{nargnames[i]}{arg}")

    with mock.patch("sys.argv", ["flows_get_brightest"] + list(cmd_args)):
        rot, tid, shifta, shiftd, inst, image, survey, cmap, scale, sigma, contrast = parse_fc()
        assert rot == args[0]
        assert tid == args[1]
        assert shifta == args[2]
        assert shiftd == args[3]
        assert survey == args[4]
        assert isinstance(inst(SkyCoord(0, 0, unit="deg")), Instrument)
        assert image is None
        assert scale == args[8]
        assert cmap == args[7]
        assert isinstance(sigma, float)
        assert isinstance(contrast, float)
        assert not isinstance(image, str)

        if args[5] == "FixedSize":
            assert inst == FixedSizeInstrument
            pytest.approx(inst.field_hw, u.Quantity(args[6], u.arcmin))
        elif args[5] == "Hawki":
            assert inst == Hawki
