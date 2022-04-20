import pytest
import astropy.units as u
import tendrils
from .observer import get_flows_observer, Observer
from .auth import test_connection
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
        monkeypatch.setattr('builtins.input', lambda _: 'bad_token')
        test_connection()


@pytest.fixture
def observer():
    return get_flows_observer(rot, tid, shifta, shiftd)


def test_get_brightest(capsys, observer):
    observer.check_bright_stars(region=observer.regions[0])
    captured = capsys.readouterr()
    assert "Brightest star has" in captured.out


def test_plan(observer):
    assert observer.plan.rotation == rot * u.deg
    assert observer.plan.alpha == shifta * u.arcsec
    assert observer.plan.delta == shiftd * u.arcsec
    assert observer.plan.shift is True
    assert observer.plan.rotate is True


def test_observer(observer):
    isinstance(observer, Observer)


@pytest.mark.slow
def test_make_finding_chart(observer, monkeypatch):
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, 'show', lambda: None)
    plotter = Plotter(observer)
    ax = plotter.make_finding_chart(observer, savefig=False)
    assert len(ax.get_images()) > 0
    title = ax.get_title()
    assert title.startswith(f"{observer.target.info['target_name']}")
    assert title.endswith("FC")
